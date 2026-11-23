// Implementation of the TRAMO automatic-fit assembly (item 23-F F.4): the
// TramoKernel.process composition of the landed units - the 23-A builder
// mirror and the 23-B transform mirror run LIVE and are cross-checked
// against the captured pre-stage facts, the F.3 round driver runs verbatim,
// the differencing/arma seams re-run their landed stage units live where
// the stage envelope qualifies (header "witness tiers"), and the final
// surface recomputes the D9 regression-block covariance
// V = unscaledCovariance * (ssq / (dim - nx - k)) - the
// ConcentratedLikelihood.covariance(k, true) chain, javap-pinned
// (RegSarimaModel.of passes getArimaSpec().freeParametersCount() and
// iconst_1 for the scaling flag; the chain computes s = ssq / (double)
// (dim - nx - k) once, then FastMatrix.mul(s) elementwise; nx == 0 answers
// FastMatrix.EMPTY) - against the captured coefficientsCovariance upper
// bits. Decoded tables: tmp_opencode/kalos-m5-pinning.md SECTION 32 with
// the header-dated amendments; corpus:
// validation/evidence/native/tramo-fit (Scratch23f, fit-capture-1/2/witness
// byte-identical, 34 identity all-1 rows + the pinned fit_const crash
// witness, which enters no test case).

#include <kalos/tramo_fit.hpp>

#include "kalos/detail/fdlibm.hpp"
#include "mlj_engine.hpp"

#include <kalos/tramo_arma_selection.hpp>
#include <kalos/tramo_differencing_decision.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace kalos {
namespace {

void require(bool condition, const std::string& what) {
    if (!condition)
        throw std::invalid_argument("tramo_fit: " + what);
}

bool bit_equal(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

std::string stage_word(TramoKernelStageResult r) {
    switch (r) {
        case TramoKernelStageResult::changed: return "Changed";
        case TramoKernelStageResult::unchanged: return "Unchanged";
        case TramoKernelStageResult::unprocessed: return "Unprocessed";
        case TramoKernelStageResult::failed: return "Failed";
    }
    return "?";
}

std::string orders_word(const TramoOrders& o) {
    return std::to_string(o.p) + "," + std::to_string(o.d) + "," +
           std::to_string(o.q) + "," + std::to_string(o.bp) + "," +
           std::to_string(o.bd) + "," + std::to_string(o.bq);
}

bool orders_equal(const TramoOrders& a, const TramoOrders& b) {
    return orders_word(a) == orders_word(b);
}

bool outliers_equal(const std::vector<std::string>& a,
                    const std::vector<std::string>& b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

// Seam kinds that mutate the description OUTSIDE the module seams (the
// round-0 restore, the pass-2 gate restore/lastSolution, the pass-1
// difference drop): a module ENTRY state may legitimately differ from the
// previous witnessed post-state across one of them.
bool kernel_transition_between(const std::vector<TramoKernelSeam>& seams,
                               std::size_t from, std::size_t to) {
    for (std::size_t k = from + 1; k <= to; ++k) {
        switch (seams[k].kind) {
            case TramoKernelSeamKind::round0:
            case TramoKernelSeamKind::pass2_gate:
            case TramoKernelSeamKind::differencing_drop:
                return true;
            default:
                break;
        }
    }
    return false;
}

// The module seams carry the post-run description state.
bool carries_state(TramoKernelSeamKind kind) {
    return kind == TramoKernelSeamKind::differencing ||
           kind == TramoKernelSeamKind::arma ||
           kind == TramoKernelSeamKind::scan;
}

const char* module_name(TramoKernelSeamKind kind) {
    switch (kind) {
        case TramoKernelSeamKind::differencing: return "tramo_differencing";
        case TramoKernelSeamKind::arma: return "tramo_arma_selection";
        default: return "tramo_outlier_rounds";
    }
}

}  // namespace

TramoFitResult tramo_fit(const TramoFitRequest& request) {
    TramoFitResult result;
    const TramoFitSpec& spec = request.spec;

    // ---- structural validation (throws) ---------------------------------
    require(!request.series.empty(), "empty series");
    require(request.start_month >= 1 && request.start_month <= 12,
            "start_month outside 1..12");

    // ---- fail-closed refusal matrix (before any state is built) ---------
    auto refuse = [&result](const char* reason, std::string note) {
        result.status = std::string("refused:") + reason;
        result.refusal_note = std::move(note);
        return result;
    };
    if (spec.annual_frequency != 12)
        return refuse("d11-frequency",
                      "annual frequency " +
                          std::to_string(spec.annual_frequency) +
                          " is outside the monthly-only envelope (D11)");
    for (const double v : request.series) {
        if (!std::isfinite(v))
            return refuse("missing-values",
                          "the fit envelope carries no interpolation plan "
                          "(qualified at unit level in 23-A/23-B)");
    }
    if (spec.td_automatic || spec.td_test || spec.easter_test ||
        spec.user_variables)
        return refuse("d10-calendar",
                      "the calendar family (trading days automatic/tested, "
                      "Easter tests, user variables) is the joint 21-I "
                      "surface (D10)");
    if (spec.auto_method != tramo_auto_method::unused)
        return refuse("i10-automethod",
                      "the automatic regression methods (FTEST/WALD/BIC/AIC) "
                      "are refused (calendar-roadmap I-10)");
    if (spec.outliers_correction)
        return refuse("oc-precorrection",
                      "the robust outliers pre-correction has never been "
                      "witnessed through the fit seam");
    if (spec.seats_decomposition_mode != -1)
        return refuse("d6-decomposition",
                      "the fit never decomposes (register D6); decomposition "
                      "policy belongs to the 24-D/24-E consumers");

    const TramoKernelRoundsPrestage& prestage = request.kernel.prestage;

    // ---- stage 1: the 23-A builder mirror, cross-checked ----------------
    tramo_build_request breq;
    breq.values = request.series;
    breq.annual_frequency = spec.annual_frequency;
    breq.start_year = request.start_year;
    breq.start_month = request.start_month;
    breq.using_auto_model = spec.using_auto_model;
    breq.transform_function = spec.transform_function;
    breq.transform_adjust = tramo_lop::none;
    breq.arima = spec.fixed_arima;  // read only when !using_auto_model
    breq.mean.used = spec.mean_used;
    breq.mean.trend_constant = spec.mean_used;
    breq.mean.test = spec.mean_is_test;
    breq.delta_tc = spec.delta_tc;
    breq.outliers = request.builder_outliers;
    const tramo_build_result built = build_tramo_model(breq);

    require(static_cast<std::size_t>(built.domain_length) ==
                request.series.size(),
            "builder domain length disagrees with the series");
    require(built.missing.empty(),
            "builder registered missings inside the fit envelope");
    require(built.log_transformation ==
                (spec.transform_function == tramo_transform_function::log),
            "builder log flag disagrees with the transform function");
    {
        // Jar-decoded pre-stage (javap TramoKernel.testSeasonality): the
        // description is rewritten to SarimaOrders.m011(freq) = (0,1,1)(0,0,0)
        // with ctx.seasonal := false only when autoModel is on and the
        // seasonality verdict toInt() < 2 (the fit_gs2 / fit_shortfail shape:
        // prestate orders (0,1,1)(0,0,0) with seasonal off). Every other row
        // keeps the registered orders, so the builder mirror is the authority
        // there; a mismatch is admitted only in the decoded downgrade shape.
        TramoOrders built_orders;
        built_orders.d = built.arima.d;
        built_orders.bd = built.arima.bd;
        built_orders.p = static_cast<int>(built.arima.phi.size());
        built_orders.q = static_cast<int>(built.arima.theta.size());
        built_orders.bp = static_cast<int>(built.arima.bphi.size());
        built_orders.bq = static_cast<int>(built.arima.btheta.size());
        if (!orders_equal(built_orders, prestage.orders)) {
            TramoOrders m011;
            m011.d = 1;
            m011.q = 1;
            require(spec.using_auto_model && !prestage.seasonal &&
                        orders_equal(prestage.orders, m011),
                    "captured prestage orders " + orders_word(prestage.orders) +
                        " match neither the builder registration " +
                        orders_word(built_orders) + " nor the seasonality-"
                        "downgraded (0,1,1)(0,0,0) shape");
        }
    }
    if (prestage.regression_test == "Unprocessed")
        require(built.mean == prestage.mean,
                "builder mean flag disagrees with the captured prestage "
                "mean");
    result.builder = built;

    // ---- stage 2: the 23-B transform mirror, cross-checked --------------
    // testTransformation runs its machine only under TransformationType.Auto
    // with outliersCorrection off (javap, TramoKernel); the +fs row leaves
    // the context defaults in force (seasonal true, log off).
    if (spec.transform_function == tramo_transform_function::automatic) {
        TramoTransformDecisionRequest treq;
        treq.values = built.interpolated;
        treq.period = spec.annual_frequency;
        treq.mean = built.mean;
        treq.detect_seasonality = spec.using_auto_model;
        treq.spec_seasonality = true;  // the ModellingContext default
        // LogSpec.getLogPreference(): the resolved fct through the pinned
        // log (the TramoSpec default is Math.log(0.95)); the resolved line
        // is the authority and the request carries the raw double.
        treq.log_preference = spec.log_preference;
        treq.precision = spec.intermediate_precision;
        result.transformation = decide_tramo_transform(treq);
    } else {
        result.transformation.seasonal = true;
        result.transformation.spec_pre = "(0,1,1)(0,1,1)";
    }
    result.transform_matches_capture =
        result.transformation.seasonal == prestage.seasonal &&
        result.transformation.log_transform == prestage.log;
    require(result.transformation.seasonal == prestage.seasonal,
            "seasonality flag disagrees with the captured prestage");
    require(result.transformation.log_transform == prestage.log,
            "log flag disagrees with the captured prestage");

    // ---- stage 3: the ami() round loop ------------------------------------
    // Jar-decoded ami() branch (javap TramoKernel.ami / isFullySpecified):
    // the round driver is skipped only when the transform is already decided
    // and nothing is automatic - function != Auto && !autoModel &&
    // !outliersDetection && !meanTest && !td.test && !td.automatic &&
    // !easter.test (the last three clauses are vacuous here: the D10 refusal
    // above gates the calendar family out). fit_nomean (auto transform,
    // templates otherwise fixed) therefore runs at least one iterate - the
    // seam stream is cross-checked against the decoded verdict.
    const bool fully_specified =
        spec.transform_function != tramo_transform_function::automatic &&
        !spec.using_auto_model && !spec.outliers_used && !spec.mean_is_test;
    if (fully_specified) {
        require(request.kernel.seams.empty(),
                "a fully-specified request carries round-loop seams");
    } else {
        require(!request.kernel.seams.empty(),
                "a round-driver request carries no round-loop seams");
        result.kernel = tramo_run_kernel_rounds(request.kernel);
    }

    // ---- stage 4: the per-seam live cross-checks --------------------------
    const std::vector<TramoKernelSeam>& seams = request.kernel.seams;
    // The evolving seasonality context (dated 2026-11-23 javap re-decode):
    // the TramoKernel module factories (differencingModule(), armaModule())
    // hand the builder seasonal(ctx.seasonal) - the context at module
    // creation, not the pre-stage value. The seasonality controller seams
    // carry the context mutations (TramoKernelSeam.seasonal_after, from the
    // per-step rsnap witness); the gs2 / shortfail pass-3 modules are created
    // after a "scontroller res Changed" flip back to seasonal.
    bool ctx_seasonal = prestage.seasonal;
    std::size_t next_input = 0;
    for (std::size_t s = 0; s < seams.size(); ++s) {
        const TramoKernelSeam& seam = seams[s];
        if (seam.kind == TramoKernelSeamKind::seasonality_controller) {
            ctx_seasonal = seam.seasonal_after;
            continue;
        }
        if (seam.kind != TramoKernelSeamKind::differencing &&
            seam.kind != TramoKernelSeamKind::arma &&
            seam.kind != TramoKernelSeamKind::scan)
            continue;
        while (next_input < request.seam_inputs.size() &&
               request.seam_inputs[next_input].seam < s)
            ++next_input;
        const bool has_input = next_input < request.seam_inputs.size() &&
                               request.seam_inputs[next_input].seam == s;

        TramoFitSeamCheck check;
        check.seam = s;
        check.unit = module_name(seam.kind);
        if (!has_input) {
            check.tier = "fact";
            check.note = "no fitseam entry record";
            result.seam_checks.push_back(check);
            ++result.seams_fact;
            continue;
        }
        const TramoFitSeamInputs& in = request.seam_inputs[next_input];

        // Zero-drift entry chain: the entry state must equal the previous
        // witnessed module post-state unless a kernel-owned transition sits
        // between them.
        for (std::size_t t = s; t-- > 0;) {
            if (!carries_state(seams[t].kind))
                continue;
            if (!kernel_transition_between(seams, t, s)) {
                require(orders_equal(seams[t].orders, in.entry_orders) &&
                            seams[t].mean == in.entry_mean &&
                            outliers_equal(seams[t].outliers,
                                           in.entry_outliers),
                        "seam " + std::to_string(s) +
                            " entry state disagrees with the previous "
                            "witnessed post-state");
            }
            break;
        }

        if (seam.kind == TramoKernelSeamKind::differencing &&
            in.entry_outliers.empty() && !spec.outliers_ml &&
            in.round == 1) {
            // Live 23-C re-run: the initial-round, outlier-free ladder over
            // the mirrored builder response. The differencingModule()
            // factory sets Builder.initial(round == 1) (javap offsets
            // 52-65), and every corpus diff seam is a round-1 call.
            TramoDifferencingRequest dreq;
            dreq.series = built.interpolated;
            dreq.log_input = prestage.log;
            dreq.period = spec.annual_frequency;
            dreq.d0 = in.entry_orders.d;
            dreq.bd0 = in.entry_orders.bd;
            dreq.p0 = in.entry_orders.p;
            dreq.q0 = in.entry_orders.q;
            dreq.bp0 = in.entry_orders.bp;
            dreq.bq0 = in.entry_orders.bq;
            dreq.mean0 = in.entry_mean;
            dreq.initial_round = true;
            dreq.seasonal = ctx_seasonal;
            dreq.ub1 = spec.ub1;
            dreq.ub2 = spec.ub2;
            dreq.cancel = spec.cancel;
            dreq.ml = false;
            dreq.nobs = in.entry_nobs > 0 ? in.entry_nobs : -1;
            const TramoDifferencingDecision d =
                decide_tramo_differencing(dreq);
            const TramoKernelStageResult mapped =
                d.result == TramoDifferencingDecision::Result::Changed
                    ? TramoKernelStageResult::changed
                    : d.result ==
                              TramoDifferencingDecision::Result::Unprocessed
                          ? TramoKernelStageResult::unprocessed
                          : TramoKernelStageResult::unchanged;
            TramoOrders post;
            post.p = d.spec_p;
            post.d = d.spec_d;
            post.q = d.spec_q;
            post.bp = d.spec_bp;
            post.bd = d.spec_bd;
            post.bq = d.spec_bq;
            require(mapped == seam.result,
                    "seam " + std::to_string(s) +
                        ": live differencing result " + stage_word(mapped) +
                        " disagrees with the seam fact " +
                        stage_word(seam.result));
            require(orders_equal(post, seam.orders),
                    "seam " + std::to_string(s) +
                        ": live differencing orders " + orders_word(post) +
                        " disagree with the seam fact " +
                        orders_word(seam.orders));
            require(d.mean == seam.mean,
                    "seam " + std::to_string(s) +
                        ": live differencing mean flag disagrees");
            check.tier = "live";
            check.note = "match " + stage_word(mapped) + " " +
                         orders_word(post);
        } else if (seam.kind == TramoKernelSeamKind::arma &&
                   !in.ols_residuals.empty()) {
            // Live 23-D re-run over the JVM-pinned olsResiduals hand-over.
            TramoArmaSelectionRequest areq;
            areq.residuals = in.ols_residuals;
            areq.nobs = in.nobs;
            areq.period = spec.annual_frequency;
            areq.p0 = in.entry_orders.p;
            areq.d0 = in.entry_orders.d;
            areq.q0 = in.entry_orders.q;
            areq.bp0 = in.entry_orders.bp;
            areq.bd0 = in.entry_orders.bd;
            areq.bq0 = in.entry_orders.bq;
            areq.seasonal = ctx_seasonal;
            areq.accept_white_noise = spec.accept_airline;
            const TramoArmaSelectionResult a = select_tramo_arma(areq);
            const TramoKernelStageResult mapped =
                a.result == TramoArmaSelectionResult::Result::Changed
                    ? TramoKernelStageResult::changed
                    : a.result == TramoArmaSelectionResult::Result::Failed
                          ? TramoKernelStageResult::failed
                          : a.result ==
                                    TramoArmaSelectionResult::Result::
                                        Unprocessed
                                ? TramoKernelStageResult::unprocessed
                                : TramoKernelStageResult::unchanged;
            TramoOrders post;
            post.p = a.spec_p;
            post.d = a.spec_d;
            post.q = a.spec_q;
            post.bp = a.spec_bp;
            post.bd = a.spec_bd;
            post.bq = a.spec_bq;
            require(mapped == seam.result,
                    "seam " + std::to_string(s) +
                        ": live arma result disagrees with the seam fact");
            require(orders_equal(post, seam.orders),
                    "seam " + std::to_string(s) +
                        ": live arma orders disagree with the seam fact");
            check.tier = "live";
            check.note = "match " + stage_word(mapped) + " " +
                         orders_word(post);
        } else {
            check.tier = "fact";
            if (seam.kind == TramoKernelSeamKind::scan)
                check.note =
                    "the outliers machine is qualified in the "
                    "tramo-outlier-rounds corpus; the fit seam rides as a "
                    "fact (curva witness)";
            else if (seam.kind == TramoKernelSeamKind::arma)
                check.note = "no olsResiduals hand-over captured";
            else
                check.note =
                    "outside the 23-C live envelope (outlier columns, "
                    "ml, or a re-round ladder)";
        }
        result.seam_checks.push_back(check);
        if (check.tier == "live")
            ++result.seams_live;
        else
            ++result.seams_fact;
    }

    // ---- stage 5: the final surface and the D9 covariance -----------------
    const TramoFitFinalRecord& fin = request.final;
    require(fin.dim > 0, "final record without a dim");
    require(fin.nx >= 0, "negative nx");
    require(fin.free_arma_parameters == fin.orders.p + fin.orders.q +
                                         fin.orders.bp + fin.orders.bq,
            "free_arma_parameters disagrees with the final orders");
    require(fin.parameters.size() ==
                static_cast<std::size_t>(fin.free_arma_parameters),
            "the final parameter vector length disagrees with k");
    require(fin.coefficients.size() == static_cast<std::size_t>(fin.nx),
            "the coefficient vector length disagrees with nx");
    require(fin.unscaled_covariance.size() ==
                static_cast<std::size_t>(fin.nx) * fin.nx,
            "unscaledCovariance is not nx x nx");
    require(fin.captured_covariance_upper.size() ==
                static_cast<std::size_t>(fin.nx) * (fin.nx + 1) / 2,
            "the captured covariance upper triangle has the wrong size");
    require(fin.scaling_factor, "the cl view is not the scaling-factor "
                                "branch covariance(k, true) consumes");
    require(fin.differenced_response.size() ==
                static_cast<std::size_t>(fin.dim),
            "the differenced response length disagrees with cl dim");
    require(fin.variables.size() == fin.outlier_names.size(),
            "the variable records disagree with the outlier names");
    {
        const std::size_t base = fin.mean ? 1 : 0;
        require(base + fin.variables.size() == fin.coefficients.size(),
                "the regression design (mean + outliers) disagrees with nx");
        for (std::size_t i = 0; i < fin.variables.size(); ++i) {
            require(fin.variables[i].name == fin.outlier_names[i],
                    "variable record " + std::to_string(i) +
                        " disagrees with the outlier name list");
            require(bit_equal(fin.variables[i].coefficient,
                              fin.coefficients[base + i]),
                    "the coefficient of variable " +
                        fin.variables[i].name +
                        " disagrees with the cl design slot");
        }
    }

    // The differenced response the landing machine consumed, recomputed
    // natively (pinned log, then the RegArimaUtility.differencingFilter /
    // BackFilter orientation of the 23-C pin: the JVM hands the BackFilter
    // the REVERSED polynomial, applied directly - the (-1)^(d + bd)
    // equivalence of the 23-C note holds entry-wise, but a cancelling sum
    // lands at JVM +0.0 where negating the ascending result would give
    // -0.0 (the fit_fct_air witness at index 99), so the reversed filter is
    // applied as the JVM builds it).
    std::vector<double> y = built.interpolated;
    if (fin.log) {
        for (double& v : y) {
            require(v > 0.0, "non-positive value under the log channel");
            v = detail::fdlibm_log(v);
        }
    }
    std::vector<double> filter = mlj::differencing_filter(
        spec.annual_frequency, fin.orders.d, fin.orders.bd);
    std::reverse(filter.begin(), filter.end());
    std::vector<double> w;
    if (filter.size() > 1) {
        require(y.size() > filter.size() - 1,
                "series too short for the differencing filter");
        w.assign(y.size() - (filter.size() - 1), 0.0);
        mlj::apply_backfilter(filter, y, 0, w);
    } else {
        w = y;
    }
    require(w.size() == fin.differenced_response.size(),
            "the recomputed differenced response length disagrees");
    for (std::size_t i = 0; i < w.size(); ++i) {
        require(bit_equal(w[i], fin.differenced_response[i]),
                "the recomputed differenced response differs from the "
                "capture at index " + std::to_string(i));
    }

    // The D9 covariance: V = uc * (ssq / (double)(dim - nx - k)), element
    // mul after the single scalar division (bytecode-conform sequence).
    result.ssq = fin.ssq;
    result.cl_dim = fin.dim;
    result.cl_nx = fin.nx;
    result.cl_nmissing = fin.nmissing;
    result.log_determinant = fin.log_determinant;
    result.covariance.captured_upper = fin.captured_covariance_upper;
    result.covariance.n = static_cast<std::size_t>(fin.nx);
    if (fin.nx > 0) {
        require(fin.dim - fin.nx - fin.free_arma_parameters > 0,
                "non-positive residual degrees of freedom");
        const double scale =
            fin.ssq / static_cast<double>(fin.dim - fin.nx -
                                          fin.free_arma_parameters);
        result.covariance_scale = scale;
        result.covariance.covariance.resize(
            static_cast<std::size_t>(fin.nx) * fin.nx);
        for (int i = 0; i < fin.nx; ++i) {
            for (int j = 0; j < fin.nx; ++j) {
                result.covariance
                    .covariance[static_cast<std::size_t>(i) * fin.nx + j] =
                    fin.unscaled_covariance[static_cast<std::size_t>(i) *
                                                fin.nx + j] *
                    scale;
            }
        }
        bool all = true;
        for (int i = 0; i < fin.nx && all; ++i) {
            for (int j = i; j < fin.nx; ++j) {
                const double v =
                    result.covariance
                        .covariance[static_cast<std::size_t>(i) * fin.nx + j];
                result.covariance.covariance_upper.push_back(v);
                // Row-major upper triangle: row r has nx - r entries, so
                // entry (i, j), i <= j, sits at i*nx - i*(i-1)/2 + (j - i).
                if (!bit_equal(v, fin.captured_covariance_upper[
                                      static_cast<std::size_t>(i) * fin.nx -
                                          static_cast<std::size_t>(i) *
                                              (i - 1) / 2 +
                                          (j - i)])) {
                    all = false;
                    break;
                }
            }
        }
        result.covariance.matches_capture = all;
        require(all, "the recomputed D9 covariance differs from the "
                     "captured coefficientsCovariance bits");
    } else {
        // The nx == 0 branch: ConcentratedLikelihood.covariance answers
        // FastMatrix.EMPTY (bytecode offsets 0-12).
        result.covariance.matches_capture =
            fin.captured_covariance_upper.empty();
        require(result.covariance.matches_capture,
                "an empty-regression record carries covariance bits");
    }

    result.log = fin.log;
    result.mean = fin.mean;
    result.orders = fin.orders;
    result.parameters = fin.parameters;
    result.variables = fin.variables;
    result.ls_obs = fin.ls_obs;
    result.ls_eff = fin.ls_eff;
    result.ls_npar = fin.ls_npar;
    result.ls_ll = fin.ls_ll;
    result.ls_tadj = fin.ls_tadj;
    result.ls_all = fin.ls_all;
    result.ls_ssq = fin.ls_ssq;
    result.ls_aic = fin.ls_aic;
    result.ls_aicc = fin.ls_aicc;
    result.ls_bic = fin.ls_bic;
    result.ls_bicc = fin.ls_bicc;
    result.ls_bic2 = fin.ls_bic2;
    result.ls_hq = fin.ls_hq;
    result.full_residuals = fin.full_residuals;
    result.status = "available";
    return result;
}

}  // namespace kalos
