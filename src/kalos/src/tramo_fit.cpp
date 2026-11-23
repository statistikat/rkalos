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
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/outlier_scoring.hpp"
#include "mlj_engine.hpp"

#include <kalos/tramo_arma_selection.hpp>
#include <kalos/tramo_differencing_decision.hpp>
#include <kalos/tramo_final_estimation.hpp>
#include <kalos/tramo_outlier_rounds.hpp>
#include <kalos/tramo_regression_tests.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <optional>
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

std::vector<double> differenced_response(std::vector<double> y, bool log, int period, int d, int bd) {
    if (log) {
        for (double& v : y) {
            require(v > 0.0, "non-positive value under the log channel");
            v = detail::fdlibm_log(v);
        }
    }
    auto filter = mlj::differencing_filter(period, d, bd);
    if (filter.size() <= 1) return y;
    std::reverse(filter.begin(), filter.end());
    require(y.size() >= filter.size(), "series too short for the differencing filter");
    std::vector<double> result(y.size() - filter.size() + 1);
    mlj::apply_backfilter(filter, y, 0, result);
    return result;
}

std::vector<double> arma_residuals(const tramo_build_result& built, bool log,
                                    const TramoFitSeamInputs& in) {
    const int period = built.annual_frequency;
    const int d = in.entry_orders.d, bd = in.entry_orders.bd;
    auto y = differenced_response(built.interpolated, log, period, d, bd);
    std::vector<std::vector<double>> columns;
    if (in.entry_mean) columns.emplace_back(y.size(), 1.0);
    require(in.entry_variables.size() == in.entry_outliers.size(), "incomplete ARMA entry design");
    for (std::size_t j = 0; j < in.entry_variables.size(); ++j) {
        const auto& var = in.entry_variables[j];
        require(var.name == in.entry_outliers[j] && var.position >= 0 &&
                    static_cast<std::size_t>(var.position) < built.interpolated.size(),
                "invalid ARMA entry outlier name or position");
        require(var.type == TramoOutlierKind::Ao || var.type == TramoOutlierKind::Ls,
                "ARMA entry design requires AO or LS");
        std::vector<double> column(built.interpolated.size());
        const auto type = var.type == TramoOutlierKind::Ao ? OutlierType::Additive : OutlierType::LevelShift;
        fill_outlier_column(type, var.position, period, 0., column);
        require(var.column.empty() || var.column == column, "ARMA outlier column disagrees with factory");
        columns.push_back(differenced_response(std::move(column), false, period, d, bd));
    }
    if (columns.empty()) return y;
    const auto beta = detail::pinned_ols_coefficients(y, columns);
    require(beta && beta->size() == columns.size(), "ARMA input regression failed");
    for (std::size_t j = 0; j < columns.size(); ++j)
        for (std::size_t i = 0; i < y.size(); ++i)
            y[i] += -(*beta)[j] * columns[j][i];
    return y;
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

struct DetectedOutlier {
    OutlierType type;
    int position;
    std::string name;
};

std::optional<std::vector<DetectedOutlier>> check_scan(const TramoFitSpec& spec, const tramo_build_result& built, bool log,
                const TramoFitSeamInputs& in, const TramoKernelSeam& expected) {
    TramoOutlierScanRequest request{};
    request.y = built.interpolated;
    request.frequency = built.annual_frequency;
    request.mean = in.entry_mean;
    request.log = log;
    const auto& o = in.entry_orders;
    request.p = o.p;
    request.d = o.d;
    request.q = o.q;
    request.bp = o.bp;
    request.bd = o.bd;
    request.bq = o.bq;
    request.ao = spec.ao;
    request.ls = spec.ls;
    request.tc = spec.tc;
    request.so = spec.so;
    request.tc_rate = spec.delta_tc;
    request.precision = spec.intermediate_precision;
    request.maximum_likelihood = spec.outliers_ml;
    // The pinned builder stores its two default caps in crossed fields.
    request.max_round = 30;
    request.max_outliers = 50;
    request.cv = expected.critical_value;
    request.upper_bound = built.estimation_length;
    if (o.d == 0 && o.bd == 0 && !tramo_stationary_scan_supported(request)) return std::nullopt;
    const auto scan = tramo_outlier_rounds(request);
    const auto outcome = scan.status == TramoOutlierRoundsResult::Status::Failed ? TramoKernelStageResult::failed :
        scan.outliers.empty() ? TramoKernelStageResult::unchanged : TramoKernelStageResult::changed;
    require(outcome == expected.result, "live scan outcome disagrees with the module witness");
    std::vector<OutlierType> factories;
    if (spec.ao) factories.push_back(OutlierType::Additive);
    if (spec.ls) factories.push_back(OutlierType::LevelShift);
    if (spec.tc) factories.push_back(OutlierType::TransitoryChange);
    if (spec.so) factories.push_back(OutlierType::SeasonalOutlier);
    std::vector<std::string> names;
    std::vector<DetectedOutlier> variables;
    for (const auto& selected : scan.outliers) {
        const auto type = factories.at(selected[1]);
        names.push_back(monthly_outlier_name(type, selected[0], built.domain_year, built.domain_position + 1));
        variables.push_back({type, selected[0], names.back()});
    }
    require(names == expected.outliers, "live scan selection disagrees with the ordered outlier witness");
    require(orders_equal(o, expected.orders) && expected.mean == in.entry_mean,
            "live scan description disagrees with the module witness");
    // ModelDescription.sortVariables stably orders automatic outliers by date
    // before estimating the regression. The scan witness uses detection order.
    std::stable_sort(variables.begin(), variables.end(),
                     [](const auto& a, const auto& b) { return a.position < b.position; });
    return variables;
}

struct EstimationEntry {
    TramoOrders orders;
    bool mean;
    std::vector<DetectedOutlier> variables;
};

struct FittedOutlierEffects {
    TramoOrders orders;
    bool mean;
    std::vector<double> parameters;
    std::vector<std::string> names;
    std::vector<double> corrected_response;
    TramoFinalEstimatorRunResult fit;
};

TramoTestStatistics fitted_test_statistics(const FittedOutlierEffects& effects) {
    const auto& fit = effects.fit;
    TramoTestStatistics test;
    test.coefficients = fit.cl_coefs;
    const auto nx = test.coefficients.size();
    require(fit.cl_unscaled_covariance.size() == nx * nx, "invalid fitted regression covariance");
    for (std::size_t i = 0; i < nx; ++i)
        test.cov_diagonal.push_back(fit.cl_unscaled_covariance[i * nx + i]);
    test.ssq = fit.cl_ssq;
    test.dim = fit.cl_dim;
    test.free_parameters = static_cast<int>(effects.parameters.size());
    return test;
}

void check_round_statistics(const FittedOutlierEffects& effects, const TramoKernelStatsRow& row,
                            int frequency, int n, std::size_t seam) {
    auto test = fitted_test_statistics(effects);
    TramoModelStatisticsRequest request;
    request.residuals = effects.fit.cl_residuals;
    request.coefficients = std::move(test.coefficients);
    request.cov_diagonal = std::move(test.cov_diagonal);
    request.ssq = test.ssq;
    request.dim = test.dim;
    request.sigma2 = test.ssq / test.dim;
    request.degrees_of_freedom = test.dim - static_cast<int>(request.coefficients.size());
    request.annual_frequency = frequency;
    request.estimation_domain = n;
    request.ami_outliers = static_cast<int>(effects.names.size());
    request.free_parameters = test.free_parameters;
    const auto stats = tramo_model_statistics(request);
    const std::string stage = " at kernel seam " + std::to_string(seam);
    require(std::equal(row.residuals.begin(), row.residuals.end(), request.residuals.begin(),
                       request.residuals.end(), bit_equal), "fitted likelihood residuals disagree" + stage);
    require(stats.built && row.dim == request.dim && row.free_parameters == request.free_parameters &&
            row.observations == stats.observations_count && row.outliers_count == stats.outliers_count &&
            row.effective_observations == stats.effective_observations_count, "fitted statistics counts disagree" + stage);
    require(bit_equal(row.ljung_box, stats.ljung_box) && bit_equal(row.ljung_box_pvalue, stats.ljung_box_pvalue),
            "fitted Ljung-Box statistics disagree" + stage);
    require(bit_equal(row.bic, stats.bic), "fitted BIC disagrees" + stage);
    require(bit_equal(row.se, stats.se), "fitted standard error disagrees" + stage);
}

bool check_regression_test(const FittedOutlierEffects& effects, bool test_mean,
                           const TramoKernelSeam& seam) {
    TramoRegressionTestRequest request;
    request.mean = effects.mean;
    request.test_mean = test_mean;
    request.mean_threshold = seam.tsig;
    request.statistics = fitted_test_statistics(effects);
    const auto decision = run_tramo_regression_tests(request);
    const auto stage = decision.estimation_cleared ? TramoKernelStageResult::changed
                                                   : TramoKernelStageResult::unchanged;
    require(seam.result == stage && seam.mean == decision.mean, "fitted regression test decision disagrees");
    // The probe asks for tstats after the test: clearing the model also clears its likelihood.
    const auto t = decision.estimation_cleared ? std::vector<double>{}
                                               : tramo_test_tstats(request.statistics, true);
    require(std::equal(t.begin(), t.end(), seam.tstats.begin(), seam.tstats.end(), bit_equal),
            "fitted regression t-statistics disagree");
    return decision.estimation_cleared;
}

FittedOutlierEffects check_final_estimation(const TramoFitSpec& spec, const std::vector<double>& y,
                            EstimationEntry& entry, const TramoKernelSeam& expected) {
    TramoFinalEstimatorRunRequest options;
    const auto& o = entry.orders;
    options.orders = {o.p, o.d, o.q, o.bp, o.bd, o.bq};
    options.mean = entry.mean;
    options.eps = spec.precision;
    options.ur = spec.final_unit_root_threshold;
    options.pass = expected.pass_at_call;
    options.ami = spec.using_auto_model;
    options.outliers = spec.outliers_used;
    std::vector<std::vector<double>> columns;
    std::vector<std::string> names;
    for (const auto& var : entry.variables) {
        names.push_back(var.name);
        columns.emplace_back(y.size());
        fill_outlier_column(var.type, var.position, spec.annual_frequency, spec.delta_tc, columns.back());
    }
    auto estimate = tramo_final_estimator_run(options, y, columns);
    const auto& t = estimate.orders_after;
    require(estimate.returned == (expected.returned != 0) && estimate.nnsig == expected.nnsig,
            "live final estimation decision disagrees with the module witness");
    entry.orders = {t.p, t.d, t.q, t.bp, t.bd, t.bq};
    entry.mean = estimate.mean_out;
    require(orders_equal(entry.orders, expected.orders) &&
            names == expected.outliers, "live final estimation specification disagrees with the witness");
    require(std::equal(estimate.parameters_after.begin(), estimate.parameters_after.end(),
                       expected.parameters.begin(), expected.parameters.end(), bit_equal),
            "live final estimation parameters disagree with the module witness");
    FittedOutlierEffects effects{entry.orders, entry.mean, estimate.parameters_after, std::move(names), {}, {}};
    if (!columns.empty() && estimate.cl_coefs.size() == columns.size() + static_cast<std::size_t>(entry.mean)) {
        // RegArimaUtility.regressionEffect accumulates the effects before
        // the differencing module subtracts their sum from the model response.
        std::vector<double> total(y.size(), 0.0);
        for (std::size_t j = 0; j < columns.size(); ++j)
            for (std::size_t i = 0; i < y.size(); ++i)
                total[i] += estimate.cl_coefs[j + entry.mean] * columns[j][i];
        effects.corrected_response = y;
        for (std::size_t i = 0; i < y.size(); ++i) effects.corrected_response[i] -= total[i];
    }
    effects.fit = std::move(estimate);
    return effects;
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
    const bool has_fixed_regressors = std::any_of(built.variables.begin(), built.variables.end(),
        [](const auto& var) { return var.kind != tramo_variable_kind::trend_constant; });
    bool ctx_seasonal = prestage.seasonal;
    std::size_t next_input = 0;
    std::optional<EstimationEntry> final_entry;
    std::optional<FittedOutlierEffects> effects, reference_effects;
    for (std::size_t s = 0; s < seams.size(); ++s) {
        const TramoKernelSeam& seam = seams[s];
        if (effects && effects->fit.cl_state != 0 &&
            (seam.kind == TramoKernelSeamKind::round0 || seam.kind == TramoKernelSeamKind::pass2_gate)) {
            check_round_statistics(*effects, seam.stats, built.annual_frequency,
                                   static_cast<int>(built.interpolated.size()), s);
            ++result.statistics_checked;
        }
        if (seam.kind == TramoKernelSeamKind::round0)
            reference_effects = effects;
        if (seam.kind == TramoKernelSeamKind::pass2_gate && seam.reverted)
            effects = reference_effects;
        if (seam.kind == TramoKernelSeamKind::estimation && seam.estimation_taken)
            effects.reset();
        if (effects && effects->fit.cl_state != 0 && seam.kind == TramoKernelSeamKind::regression_test) {
            if (check_regression_test(*effects, result.kernel.ami_check_mu, seam))
                effects.reset();
            ++result.regression_tests_checked;
        }
        if (effects && (seam.kind == TramoKernelSeamKind::seasonality_controller ||
                        seam.kind == TramoKernelSeamKind::model_controller) &&
            (seam.parameters != effects->parameters || seam.outliers != effects->names ||
             !orders_equal(seam.orders, effects->orders) || seam.mean != effects->mean))
            effects.reset();
        if (s == 0 && seam.kind == TramoKernelSeamKind::estimation && prestage.outliers.empty()) {
            final_entry = EstimationEntry{prestage.orders, prestage.mean, {}};
            continue;
        }
        if (final_entry && seam.kind == TramoKernelSeamKind::regression_test && seam.retry_site &&
            s > 0 && seams[s - 1].kind == TramoKernelSeamKind::estimation_model) {
            final_entry->mean = seam.mean;
            continue;
        }
        if (seam.kind == TramoKernelSeamKind::seasonality_controller) {
            ctx_seasonal = seam.seasonal_after;
            final_entry.reset();
            continue;
        }
        if (seam.kind == TramoKernelSeamKind::estimation_model) {
            TramoFitSeamCheck check{s, "final_estimator_run", "fact", "entry regression design unavailable"};
            if (!has_fixed_regressors && final_entry) {
                const auto y = differenced_response(built.interpolated, prestage.log, built.annual_frequency, 0, 0);
                effects = check_final_estimation(spec, y, *final_entry, seam);
                check.tier = "live";
                check.note = "match final estimation decision, orders and parameters";
                ++result.seams_live;
            } else {
                effects.reset();
                ++result.seams_fact;
            }
            result.seam_checks.push_back(std::move(check));
            continue;
        }
        // An intervening stage can replace the description. Carry an entry
        // only across its immediately verified predecessor or mean-test retry.
        final_entry.reset();
        if (seam.kind == TramoKernelSeamKind::scan && seam.outliers.empty())
            final_entry = EstimationEntry{seam.orders, seam.mean, {}};
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

        const bool can_derive_arma = !has_fixed_regressors &&
            (in.entry_outliers.empty() || !in.entry_variables.empty());
        const bool fitted_effects = effects && !effects->corrected_response.empty() &&
            effects->names == in.entry_outliers && orders_equal(effects->orders, in.entry_orders) &&
            effects->mean == in.entry_mean;
        if (seam.kind == TramoKernelSeamKind::differencing && !has_fixed_regressors &&
            (in.entry_outliers.empty() || fitted_effects)) {
            TramoDifferencingRequest dreq;
            if (in.entry_outliers.empty()) {
                dreq.series = built.interpolated;
                dreq.log_input = prestage.log;
            } else {
                dreq.series = effects->corrected_response;
            }
            dreq.period = spec.annual_frequency;
            dreq.d0 = in.entry_orders.d;
            dreq.bd0 = in.entry_orders.bd;
            dreq.p0 = in.entry_orders.p;
            dreq.q0 = in.entry_orders.q;
            dreq.bp0 = in.entry_orders.bp;
            dreq.bq0 = in.entry_orders.bq;
            dreq.mean0 = in.entry_mean;
            dreq.initial_round = in.round == 1;
            dreq.eps = spec.intermediate_precision;
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
        } else if (seam.kind == TramoKernelSeamKind::scan && !has_fixed_regressors && in.entry_outliers.empty()) {
            auto selected = check_scan(spec, built, prestage.log, in, seam);
            if (selected) {
                final_entry = EstimationEntry{in.entry_orders, in.entry_mean, std::move(*selected)};
                check.tier = "live";
                check.note = "match " + stage_word(seam.result) + ", ordered outlier selection";
            } else {
                check.tier = "fact";
                check.note = "outside the qualified stationary scan design";
            }
        } else if (seam.kind == TramoKernelSeamKind::arma &&
                   (can_derive_arma || !in.ols_residuals.empty())) {
            TramoArmaSelectionRequest areq;
            if (can_derive_arma) {
                areq.residuals = arma_residuals(built, prestage.log, in);
                areq.nobs = built.estimation_length;
                if (!in.ols_residuals.empty()) {
                    require(areq.nobs == in.nobs && areq.residuals.size() == in.ols_residuals.size() &&
                            std::equal(areq.residuals.begin(), areq.residuals.end(),
                                       in.ols_residuals.begin(), bit_equal),
                            "seam " + std::to_string(s) + ": native ARMA residuals differ from the witness");
                }
            } else {
                areq.residuals = in.ols_residuals;
                areq.nobs = in.nobs;
            }
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
            if (in.entry_outliers.empty())
                final_entry = EstimationEntry{post, in.entry_mean, {}};
            check.note = "match " + stage_word(mapped) + " " +
                         orders_word(post);
        } else {
            check.tier = "fact";
            if (seam.kind == TramoKernelSeamKind::scan)
                check.note =
                    "outside the live scan design or differencing envelope";
            else if (seam.kind == TramoKernelSeamKind::arma)
                check.note = "no olsResiduals hand-over captured";
            else
                check.note =
                    "no current fitted outlier effects for the differencing input";
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
    const auto w = differenced_response(built.interpolated, fin.log, spec.annual_frequency,
                                        fin.orders.d, fin.orders.bd);
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
