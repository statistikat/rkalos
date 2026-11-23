// Item 23-F F.1: transcription of the private decision methods of the
// pinned jdplus.tramoseats.base.core.tramo.FinalEstimator (bytecode
// offsets in tmp_opencode/kalos-m5-pinning.md SECTION 26.2/26.6; probe
// adjudication out/fe1e against the 3.9.0 jars).

#include "kalos/tramo_final_estimation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/detail/ucarima_engine.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/final_information.hpp"
#include "kalos/hannan_rissanen.hpp"
#include "kalos/regression_likelihood.hpp"

namespace kalos {
namespace {

using detail::cplx;
using detail::Poly;
using detail::p_roots;
using detail::gt_l; // dcmpl + ifgt: strict, NaN excluded
using detail::gt_g; // dcmpg + ifgt: NaN counts as above
using detail::ge_g; // dcmpg + ifge: NaN counts as at-or-above
using detail::lt_g; // dcmpg + iflt: strict, NaN excluded

[[noreturn]] void refuse(const std::string& why) {
    throw std::invalid_argument("tramo_final_estimation: " + why);
}

// Complex.abs: the scaled modulus m * sqrt(1 + (n / m) ^ 2) with the
// (0, 0) and n == 0 shortcuts, NOT the correctly rounded std::hypot
// (23-C pin of the same JVM primitive).
double jvm_complex_abs(double re, double im) {
    if (re == 0.0 && im == 0.0) {
        return 0.0;
    }
    const double a = std::fabs(re), b = std::fabs(im);
    const double m = std::max(a, b), n = std::min(a, b);
    if (n == 0.0) {
        return m;
    }
    const double t = n / m;
    return m * std::sqrt(1.0 + t * t);
}

// getRegularAR()/getSeasonalAR(): Polynomial.ofInternal([1 | block]) - the
// trimming variant (leading zeros of the block shrink the effective
// degree; mr06 of the corpus pins it) - optionally followed by mirror()
// (plain reversal of the trimmed array, itself untrimmed) and roots().
std::vector<cplx> ar_value_roots(
    const std::vector<double>& block, const char* what, bool mirror) {
    std::vector<double> coefficients;
    coefficients.reserve(block.size() + 1);
    coefficients.push_back(1.0);
    coefficients.insert(coefficients.end(), block.begin(), block.end());
    Poly poly = detail::p_of(std::move(coefficients));
    if (mirror) {
        std::reverse(poly.c.begin(), poly.c.end());
    }
    std::vector<cplx> roots;
    try {
        roots = p_roots(std::move(poly), what);
    } catch (const detail::jvm_decompose_abort& abort) {
        throw std::runtime_error(
            std::string("tramo_final_estimation: ") + what + ": " + abort.reason);
    }
    return roots;
}

std::vector<cplx> ar_mirror_roots(const std::vector<double>& block, const char* what) {
    return ar_value_roots(block, what, true);
}

} // namespace

std::string tramo_final_orders_to_string(const TramoFinalOrders& orders) {
    return "(" + std::to_string(orders.p) + "," + std::to_string(orders.d) + ","
        + std::to_string(orders.q) + ")(" + std::to_string(orders.bp) + ","
        + std::to_string(orders.bd) + "," + std::to_string(orders.bq) + ")";
}

std::vector<TramoFinalComplex> tramo_mirror_roots(
    const std::vector<double>& ascending_coefficients) {
    Poly poly = detail::p_of(std::vector<double>(ascending_coefficients));
    std::reverse(poly.c.begin(), poly.c.end());
    std::vector<cplx> roots;
    try {
        roots = p_roots(std::move(poly), "mirror_roots");
    } catch (const detail::jvm_decompose_abort& abort) {
        throw std::runtime_error(
            "tramo_final_estimation: mirror_roots: " + abort.reason);
    }
    std::vector<TramoFinalComplex> out;
    out.reserve(roots.size());
    for (const cplx& z : roots) {
        out.push_back({z.re, z.im});
    }
    return out;
}

int tramo_searchur_roots(const std::vector<TramoFinalComplex>& roots, double ur) {
    int count = 0;
    for (const TramoFinalComplex& z : roots) {
        const double im = std::fabs(z.im);
        const double modulus = jvm_complex_abs(z.re, z.im);
        if (gt_l(modulus, ur) && !gt_g(im, 0.05) && gt_l(z.re, 0.0)) {
            ++count;
        }
    }
    return count;
}

int tramo_searchur(const std::vector<double>& ar_block, double ur) {
    const std::vector<cplx> roots = ar_mirror_roots(ar_block, "searchur");
    std::vector<TramoFinalComplex> converted;
    converted.reserve(roots.size());
    for (const cplx& z : roots) {
        converted.push_back({z.re, z.im});
    }
    return tramo_searchur_roots(converted, ur);
}

TramoUnitRootsResult tramo_check_unit_roots(const TramoUnitRootsRequest& request) {
    TramoUnitRootsResult result;
    result.orders_after = request.orders;
    if (request.orders.p < 0 || request.orders.d < 0 || request.orders.q < 0
        || request.orders.bp < 0 || request.orders.bd < 0 || request.orders.bq < 0) {
        refuse("negative orders");
    }
    // offsets 16-70: if (p > 0 && d < 2 && searchur(AR mirror roots) != 0)
    // { setP(p - 1); setD(d + 1); ans = false; } The length guards sit
    // inside the JVM read sites: a block is required only when its scan
    // actually runs (a gated-off scan never touches the polynomial).
    if (request.orders.p > 0 && request.orders.d < 2) {
        if (request.phi.size() < static_cast<std::size_t>(request.orders.p)) {
            refuse("phi block shorter than p");
        }
        result.regular_checked = true;
        result.regular_count = tramo_searchur(request.phi, request.ur);
        if (result.regular_count != 0) {
            result.regular_fired = true;
            --result.orders_after.p;
            ++result.orders_after.d;
            result.unchanged = false;
        }
    }
    // offsets 72-126: the seasonal twin with the bd < 1 gate.
    if (request.orders.bp > 0 && request.orders.bd < 1) {
        if (request.bphi.size() < static_cast<std::size_t>(request.orders.bp)) {
            refuse("bphi block shorter than bp");
        }
        result.seasonal_checked = true;
        result.seasonal_count = tramo_searchur(request.bphi, request.ur);
        if (result.seasonal_count != 0) {
            result.seasonal_fired = true;
            --result.orders_after.bp;
            ++result.orders_after.bd;
            result.unchanged = false;
        }
    }
    return result;
}

TramoCommonRootsResult tramo_check_common_roots(const TramoCommonRootsRequest& request) {
    TramoCommonRootsResult result;
    result.orders_after = request.orders;
    if (request.orders.p < 0 || request.orders.d < 0 || request.orders.q < 0
        || request.orders.bp < 0 || request.orders.bd < 0 || request.orders.bq < 0) {
        refuse("negative orders");
    }
    // offsets 16-115: p != 0 && q != 0 computes both root arrays (roots()
    // on getRegularAR()/getRegularMA(), NO mirror) and inverts the entries
    // in place; nothing reads the arrays afterwards - but the roots() calls
    // can still raise, so they execute faithfully and only the results are
    // dropped. Complex.inv on a zero root yields NaN without raising, so
    // the inversion itself is unobservable.
    if (request.orders.p != 0 && request.orders.q != 0) {
        if (request.phi.size() < static_cast<std::size_t>(request.orders.p)
            || request.theta.size() < static_cast<std::size_t>(request.orders.q)) {
            refuse("regular blocks shorter than p/q");
        }
        (void)ar_value_roots(request.phi, "checkCommonRoots regular AR", false);
        (void)ar_value_roots(request.theta, "checkCommonRoots regular MA", false);
    }
    // offsets 116-165: bp == 1 && bq == 1 && |bphi(1) - btheta(1)| < cancel.
    const bool rule_eligible = request.orders.bp == 1 && request.orders.bq == 1;
    if (rule_eligible) {
        if (request.bphi.empty() || request.btheta.empty()) {
            refuse("bphi/btheta block missing for the bp == bq == 1 rule");
        }
        const double delta = std::fabs(request.bphi[0] - request.btheta[0]);
        if (lt_g(delta, request.cancel)) { // dcmpg, ifge: NaN never cancels
            result.cancelled = true;
            result.unchanged = false;
            result.orders_after.bp = 0;
            result.orders_after.bq = 0;
        }
    }
    return result;
}

TramoFinalEstimatorTestResult tramo_final_estimator_test(
    const TramoFinalEstimatorTestRequest& request) {
    const TramoFinalOrders& orders = request.orders;
    if (orders.p < 0 || orders.d < 0 || orders.q < 0 || orders.bp < 0 || orders.bd < 0
        || orders.bq < 0) {
        refuse("negative orders");
    }
    const std::size_t dim = static_cast<std::size_t>(orders.p + orders.bp + orders.q
        + orders.bq);
    if (request.parameters.size() < dim) {
        refuse("parameters shorter than p + bp + q + bq");
    }
    if (dim > 0 && request.asy_covariance_diagonal.size() < dim) {
        refuse("asymptotic covariance diagonal shorter than the parameter count");
    }
    const double tsig = request.tsig;
    const double sig2 = request.sample_length <= 150 ? 0.15 : 0.1; // offset 17-34

    TramoFinalEstimatorTestResult result;
    result.orders_after = orders;
    double ratio = 99999.0; // offset 76; the x13 fork starts at 999
    bool ratio_set = false;
    int index = -1; // offset 81; running SarimaModel parameter index

    // One block: offsets 112-183 (P), 193-275 (Bp), 285-370 (Q), 380-464 (Bq).
    struct Block {
        int order;
        bool* flag;
        bool gated;      // P has no ratio gate
        bool clears_p;
        bool clears_bp;
        bool clears_q;
        bool updates_ratio;
    };
    Block blocks[4] = {{orders.p, &result.fired_p, false, false, false, false, true},
        {orders.bp, &result.fired_bp, true, true, false, false, true},
        {orders.q, &result.fired_q, true, true, true, false, true},
        {orders.bq, &result.fired_bq, true, true, true, true, false}};
    for (const Block& block : blocks) {
        index += block.order;
        if (block.order <= 0) {
            continue;
        }
        const double a = std::fabs(request.parameters[static_cast<std::size_t>(index)]);
        const double v = request.asy_covariance_diagonal[static_cast<std::size_t>(index)];
        if (!gt_l(v, 0.0)) { // dcmpl, ifle: NaN variance skips the block
            continue;
        }
        const double t = a / std::sqrt(v);
        if (ge_g(t, tsig) || ge_g(a, sig2)) { // dcmpg, ifge on both
            continue;
        }
        if (block.gated && !gt_l(ratio, t)) { // dcmpl, ifle: strict, NaN t blocks
            continue;
        }
        *block.flag = true;
        if (block.clears_p) {
            result.fired_p = false;
        }
        if (block.clears_bp) {
            result.fired_bp = false;
        }
        if (block.clears_q) {
            result.fired_q = false;
        }
        if (block.updates_ratio) {
            ratio = t;
            ratio_set = true;
        }
    }
    result.ratio_final = ratio_set ? ratio : 0.0;
    result.total = (result.fired_p ? 1 : 0) + (result.fired_bp ? 1 : 0)
        + (result.fired_q ? 1 : 0) + (result.fired_bq ? 1 : 0);
    if (result.total == 0) { // offsets 477-483: no clone, no event
        return result;
    }
    if (result.fired_p) {
        --result.orders_after.p;
    }
    if (result.fired_bp) {
        --result.orders_after.bp;
    }
    if (result.fired_q) {
        --result.orders_after.q;
    }
    if (result.fired_bq) {
        --result.orders_after.bq;
    }
    result.spec_changed = true;
    result.info_detail = tramo_final_orders_to_string(result.orders_after);
    return result;
}

// ---------------------------------------------------------------------------
// Item 23-F F.2: the estimate() round-loop driver. SECTION 29 decodes; the
// Scratch23n corpus (n23-capture-*.out) confirms the decoded driver against
// the live FinalEstimator bit for bit on all ten fe_est cases.
// ---------------------------------------------------------------------------

namespace {

// The fe_est family rebuild (SECTION 29.6): the probe's StrictMath series,
// bit-identical through the fdlibm mirror; the operation grouping follows
// the Java source exactly (left-associative chain, 2 * PI * t / 12 inside
// the trig call).
std::vector<double> est_family_series(TramoEstFamily family, int n) {
    constexpr double kPi = 3.141592653589793; // StrictMath.PI
    std::vector<double> y(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) {
        const double t = static_cast<double>(k + 1);
        y[static_cast<std::size_t>(k)] = 50.0 + 0.05 * t
            + 5.0 * detail::fdlibm_sin(2.0 * kPi * t / 12.0)
            + 3.0 * detail::fdlibm_cos(4.0 * kPi * t / 12.0)
            + 0.5 * detail::fdlibm_sin(14.0 * kPi * t / 12.0);
    }
    if (family == TramoEstFamily::LevelShift) {
        for (int k = 48; k < n; ++k) {
            y[static_cast<std::size_t>(k)] += 20.0;
        }
    } else if (family == TramoEstFamily::TwoOutliers) {
        y[10] += 12.0;
        y[60] -= 18.0;
    }
    return y;
}

// SarimaMapping2-style HannanRissanen start in JVM coordinates
// ([phi | bphi | theta | btheta], AR blocks negated from the Kalos
// convention); the mirror of the scan machine's hr_to_jvm.
std::vector<double> hr_start_to_jvm(const std::array<std::vector<double>, 4>& c) {
    std::vector<double> jvm;
    for (const double v : c[0]) {
        jvm.push_back(-v);
    }
    for (const double v : c[1]) {
        jvm.push_back(-v);
    }
    for (const double v : c[2]) {
        jvm.push_back(v);
    }
    for (const double v : c[3]) {
        jvm.push_back(v);
    }
    return jvm;
}

TramoFinalEstimatorRunResult run_final_estimator(
    const TramoFinalEstimatorRunRequest& request, const std::vector<double>* observations) {
    if (observations && (observations->size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        (request.n != 0 && static_cast<std::size_t>(request.n) != observations->size()))) {
        refuse("response length does not match the supplied observations");
    }
    const int n = observations ? static_cast<int>(observations->size()) : request.n;
    if (n <= 0) {
        refuse("the response length must be positive");
    }
    if (!observations && request.family == TramoEstFamily::TwoOutliers && n <= 60) {
        refuse("the TwoOutliers family reaches past the response");
    }
    if (!std::isfinite(request.eps) || request.eps <= 0.0
        || !std::isfinite(request.ur)) {
        refuse("eps must be positive and finite and ur finite");
    }
    TramoFinalEstimatorRunResult result;
    TramoFinalOrders orders = request.airline
        ? TramoFinalOrders{0, 1, 1, 0, 1, 1}
        : request.orders;
    if (!request.airline) {
        const int dims[6] = {orders.p, orders.d, orders.q, orders.bp, orders.bd, orders.bq};
        for (const int v : dims) {
            if (v < 0) {
                refuse("negative specification order");
            }
        }
    }
    if (orders.d + 12LL * orders.bd >= n) {
        refuse("the differencing consumes the response");
    }
    if (request.mean && n - orders.d - 12 * orders.bd <= 1) {
        refuse("the differenced design needs a non-trivial response");
    }

    // RegSarimaModelling state: the response, the optional raw constant
    // column, the specification (orders + arima parameters, empty while the
    // arima_ field is nulled behind a setSpecification) and the estimation
    // slot (0 null, 1 full with point, 2 zero-parameter).
    std::vector<double> family_series;
    if (observations) {
        for (double v : *observations)
            if (!std::isfinite(v))
                refuse("the model-scale observations must be finite");
    } else {
        family_series = est_family_series(request.family, n);
    }
    const auto& y = observations ? *observations : family_series;
    const int period = 12;
    const auto differencing = [&] {
        return detail::pinned_differencing_polynomial(period, orders.d, orders.bd);
    };
    std::vector<double> pars; // JVM [phi | bphi | theta | btheta]

    int nnsig = 0;
    int estimation_state = 0;
    int cl_dim = 0;
    double cl_ssq = 0.0;
    std::vector<double> cl_coefs, covariance_diagonal;
    auto null_estimation = [&] {
        estimation_state = 0;
        cl_dim = 0;
        cl_ssq = 0.0;
        cl_coefs.clear();
        covariance_diagonal.clear();
        result.landing_information.clear();
    };
    auto set_specification = [&] {
        pars.clear(); // arima_ := null; parameters() reads empty
        null_estimation();
    };

    // ProcessingLog mirror: push emits before the depth rises, pop after it
    // falls (the probe Observer), info/remark at the current depth.
    int depth = 0;
    result.events.emplace_back(1, depth); // push("full estimation")
    ++depth;
    auto finish = [&](bool returned, int remark) {
        if (remark != 0) {
            result.events.emplace_back(remark, depth); // 3 remark, 4 info
        }
        --depth;
        result.events.emplace_back(2, depth); // pop
        result.returned = returned;
        result.nnsig = nnsig;
        result.mean_out = request.mean;
        result.orders_after = orders;
        result.parameters_after = pars;
        result.cl_state = estimation_state;
        result.cl_dim = cl_dim;
        result.cl_ssq = cl_ssq;
        result.cl_coefs = cl_coefs;
        result.covariance_diagonal = covariance_diagonal;
        return result;
    };
    // The pinned catch (Exception): remark("estimation faied") [jar typo
    // verbatim], pop, return false; the estimation slot keeps whatever the
    // last completed round left.
    auto estimation_failed = [&] { return finish(false, 3); };

    int continues = 0; // "while (i++ < 5)" post-test: at most five bodies

    for (int round = 0;; ++round) {
        TramoFinalEstimatorRoundRecord rec;
        rec.round = round;
        const int dim = orders.p + orders.bp + orders.q + orders.bq;
        rec.dim = dim;

        // m.estimate(computer): ModelDescription.estimate over
        // desc.regarima() with RegSarimaComputer
        // (minimizer LM, precision eps, computeExactFinalDerivatives true).
        const int dly_size = n - orders.d - 12 * orders.bd;
        const auto dly = detail::pinned_backfilter_apply(differencing(), y);
        std::vector<std::vector<double>> design_cols; // machine (differenced-domain) view
        if (request.mean) {
            // The pinned ConcentratedLikelihoodComputer's meanCorrection
            // constant: 1.0 over the differenced domain (PinnedSsqProvider's
            // qualified view; the differenced ones column is identically 0).
            design_cols.push_back(
                std::vector<double>(static_cast<std::size_t>(dly_size), 1.0));
        }

        if (dim == 0) {
            // The zero-parameter shortcut: DEFAULT_COMPUTER.compute(m4).
            auto cl = request.mean
                ? detail::pinned_deviances_evaluate_multi(dly, 0, 0, 0, 0, period,
                    design_cols, {}, detail::PinnedComputerShape::Fixed)
                : detail::pinned_deviances_evaluate(dly, 0, 0, 0, 0, period, {}, {});
            if (!cl) {
                result.rounds.push_back(rec);
                return estimation_failed(); // "Gls failed" et al.
            }
            estimation_state = 2; // no landing point
            cl_dim = static_cast<int>(cl->e.size());
            cl_ssq = cl->ssq_err;
            cl_coefs = request.mean ? cl->b : std::vector<double>{};
            rec.shortcircuit = true;
            rec.nnsig_after = nnsig;
            rec.outcome = "return_true";
            result.rounds.push_back(rec);
            return finish(true, 0); // the dim == 0 short-circuit
        }

        // HannanRissanenInitializer (stabilize, useDefaultIfFailed) over the
        // differenced model, then the RAPC chain warm-started from it.
        RegressionDesign hr_design;
        hr_design.columns = design_cols.empty() ? 0 : design_cols.size();
        if (!design_cols.empty()) {
            hr_design.values = design_cols[0];
            hr_design.coefficients.assign(design_cols.size(), 0.0);
            hr_design.fixed.assign(design_cols.size(), 0);
        }
        const HannanRissanenStart hr = hannan_rissanen_start(dly, hr_design,
            HannanRissanenOrders{orders.p, orders.bp, orders.q, orders.bq, period});
        rec.hr_source = hr.source == HannanRissanenStart::Source::Fitted
            ? "fitted"
            : hr.source == HannanRissanenStart::Source::ZeroVariance ? "zerovariance"
                                                                     : "defaulted";
        const std::vector<double> warm = hr_start_to_jvm(hr.coefficients);

        SarimaEstimationRequest est_request;
        est_request.series.values = y;
        est_request.series.missing.assign(static_cast<std::size_t>(n), 0);
        if (request.mean) {
            est_request.design.columns = 1;
            est_request.design.values.assign(static_cast<std::size_t>(n), 1.0);
            est_request.design.coefficients.assign(1, 0.0);
            est_request.design.fixed.assign(1, 0);
        } else {
            est_request.design.columns = 0;
        }
        const std::size_t block_sizes[4] = {static_cast<std::size_t>(orders.p),
            static_cast<std::size_t>(orders.q), static_cast<std::size_t>(orders.bp),
            static_cast<std::size_t>(orders.bq)};
        for (std::size_t b = 0; b < 4; ++b) {
            est_request.blocks[b].initial.assign(block_sizes[b], 0.0);
            est_request.blocks[b].fixed.assign(block_sizes[b], 0);
        }
        est_request.period = period;
        est_request.d = orders.d;
        est_request.seasonal_d = orders.bd;
        est_request.control = EstimationControl{request.eps, 1e-4, true, false, true,
            true, EstimationStartingPoint::HannanRissanen,
            EstimationSolver::PinnedLevenbergMarquardt};
        est_request.request_landing_point = true;
        std::array<std::vector<double>, 4> start;
        {
            std::size_t at = 0;
            for (int i = 0; i < orders.p; ++i) {
                start[0].push_back(-warm[at++]);
            }
            for (int i = 0; i < orders.bp; ++i) {
                start[2].push_back(-warm[at++]);
            }
            for (int i = 0; i < orders.q; ++i) {
                start[1].push_back(warm[at++]);
            }
            for (int i = 0; i < orders.bq; ++i) {
                start[3].push_back(warm[at++]);
            }
        }
        const SarimaEstimationResult est = optimize_sarima_from(est_request, start);
        if (!est.has_landing_point) {
            result.rounds.push_back(rec);
            return estimation_failed(); // the pinned Evaluation threw
        }
        pars.clear();
        for (const double v : est.coefficients[0]) {
            pars.push_back(-v);
        }
        for (const double v : est.coefficients[2]) {
            pars.push_back(-v);
        }
        for (const double v : est.coefficients[1]) {
            pars.push_back(v);
        }
        for (const double v : est.coefficients[3]) {
            pars.push_back(v);
        }
        rec.landing = pars;

        auto cl = request.mean
            ? detail::pinned_deviances_evaluate_multi(dly,
                static_cast<std::size_t>(orders.p), static_cast<std::size_t>(orders.bp),
                static_cast<std::size_t>(orders.q), static_cast<std::size_t>(orders.bq),
                period, design_cols, pars, detail::PinnedComputerShape::Full)
            : detail::pinned_deviances_evaluate(dly,
                static_cast<std::size_t>(orders.p), static_cast<std::size_t>(orders.bp),
                static_cast<std::size_t>(orders.q), static_cast<std::size_t>(orders.bq),
                period, {}, pars);
        if (!cl) {
            result.rounds.push_back(rec);
            return estimation_failed();
        }
        estimation_state = 1;
        cl_dim = static_cast<int>(cl->e.size());
        cl_ssq = cl->ssq_err;
        cl_coefs = request.mean ? cl->b : std::vector<double>{};

        // The Point tail, the RegArmaProcessor.compute exact-derivatives
        // branch (computeExactFinalDerivatives, the FE's computer flag):
        // LinearModel.calcResiduals regresses the landing allCoefficients off
        // the differenced response (the mean column included; the pinned b is
        // the ConcentratedLikelihood allCoefficients the JVM feeds in), the
        // deviance machine re-runs DESIGN-FREE on those residuals with the
        // SarimaMapping2 stencil (NumericalDerivatives over
        // DefaultLikelihoodEvaluation.deviance = Likelihood.ssq * factor, the
        // pinned ssq_e), and score/information scale with 0.5 * nobs / fx on
        // the cdf nobs y_d.length - X.getColumnsCount() - dim. The
        // concentrated deviance can be indefinite at the landing (es01
        // vd1 < 0), which the pinned SymmetricMatrix.inverse Cholesky-then-
        // Doolittle fallback absorbs; a singular inverse mirrors the
        // MatrixException into the round catch.
        std::vector<double> tail_y = dly;
        if (cl->b.size() != design_cols.size()) {
            result.rounds.push_back(rec);
            return estimation_failed(); // calcResiduals would throw
        }
        for (std::size_t k = 0; k < design_cols.size(); ++k) {
            for (std::size_t i = 0; i < tail_y.size(); ++i) {
                tail_y[i] -= design_cols[k][i] * cl->b[k];
            }
        }
        final_information_request fi;
        fi.orders = SarimaOrders{orders.p, orders.bp, orders.q, orders.bq};
        fi.point = pars;
        fi.nobs = dly_size - static_cast<int>(design_cols.size()) - dim;
        {
            const auto tail_c = tail_y;
            const SarimaOrders o{orders.p, orders.bp, orders.q, orders.bq};
            fi.deviance =
                [tail_c, o, period](
                    const std::vector<double>& x) -> std::optional<double> {
                    const auto pt = detail::pinned_deviances_evaluate(tail_c,
                        static_cast<std::size_t>(o.p),
                        static_cast<std::size_t>(o.bp),
                        static_cast<std::size_t>(o.q),
                        static_cast<std::size_t>(o.bq), period, {}, x);
                    if (!pt) {
                        return std::nullopt;
                    }
                    return pt->ssq_e;
                };
        }
        const final_information_result tail = final_information(fi);
        rec.tail_ssq_e = tail.fx;
        rec.tail_score = tail.score;
        rec.tail_information = tail.information;
        if (tail.status != final_information_status::available) {
            result.rounds.push_back(rec);
            return estimation_failed(); // the Point build throws; the round catch
        }
        result.landing_information = tail.information;
        covariance_diagonal = tail.asymptotic_covariance_diagonal;

        // if (!ami) { pop; return true; }: the auto-modeling gate, after the
        // estimation is assigned. nnsig keeps its previous value; es03 shows
        // it is indistinguishable from the nnsig == 0 exit for a first round.
        if (!request.ami) {
            rec.nnsig_after = nnsig;
            rec.outcome = "return_true";
            result.rounds.push_back(rec);
            return finish(true, 0);
        }

        // The decision cascade on the fresh landing.
        std::vector<double> phi, bphi, theta, btheta;
        {
            std::size_t at = 0;
            phi.assign(pars.begin() + at, pars.begin() + at + orders.p);
            at += orders.p;
            bphi.assign(pars.begin() + at, pars.begin() + at + orders.bp);
            at += orders.bp;
            theta.assign(pars.begin() + at, pars.begin() + at + orders.q);
            at += orders.q;
            btheta.assign(pars.begin() + at, pars.begin() + at + orders.bq);
        }
        const TramoUnitRootsResult cur = tramo_check_unit_roots(
            TramoUnitRootsRequest{orders, phi, bphi, request.ur});
        rec.checkur = cur.unchanged ? 1 : 0;
        if (!cur.unchanged) {
            orders = cur.orders_after;
            set_specification();
            rec.nnsig_after = nnsig;
            rec.outcome = "continue";
            result.rounds.push_back(rec);
            if (++continues >= 5) {
                return finish(false, 0); // the loop test failed
            }
            continue;
        }
        const TramoCommonRootsResult cc = tramo_check_common_roots(
            TramoCommonRootsRequest{orders, phi, bphi, theta, btheta, 0.044});
        rec.checkcr = cc.unchanged ? 1 : 0;
        if (cc.unchanged) {
            TramoFinalEstimatorTestRequest treq;
            treq.orders = orders;
            treq.parameters = pars;
            treq.asy_covariance_diagonal = covariance_diagonal;
            treq.sample_length = dly_size;
            treq.tsig = 1.0;
            const TramoFinalEstimatorTestResult tst = tramo_final_estimator_test(treq);
            rec.test_total = tst.total;
            nnsig = tst.total;
            if (tst.spec_changed) { // setSpecification + the info event
                orders = tst.orders_after;
                set_specification();
                result.events.emplace_back(4, depth);
            }
        } else {
            nnsig = 2; // the seasonal AR/MA cancellation
            orders = cc.orders_after;
            set_specification();
        }
        rec.nnsig_after = nnsig;
        if (nnsig == 0) {
            rec.outcome = "return_true";
            result.rounds.push_back(rec);
            return finish(true, 0);
        }
        if (nnsig >= 2 && request.outliers && request.pass <= 1) {
            rec.outcome = "return_false";
            result.rounds.push_back(rec);
            return finish(false, 0);
        }
        rec.outcome = "continue";
        result.rounds.push_back(rec);
        if (++continues >= 5) {
            return finish(false, 0);
        }
    }
}

} // namespace

TramoFinalEstimatorRunResult tramo_final_estimator_run(
    const TramoFinalEstimatorRunRequest& request) {
    return run_final_estimator(request, nullptr);
}

TramoFinalEstimatorRunResult tramo_final_estimator_run(
    const TramoFinalEstimatorRunRequest& request, const std::vector<double>& model_series) {
    return run_final_estimator(request, &model_series);
}

} // namespace kalos
