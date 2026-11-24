// Item 23-B: the round-0 TRAMO outlier scan of the pinned jdplus 3.9.0
// stack (transcription map: docs/contracts/tramo-outlier-scan-contract.md).
//
// Machine order of a legal request: the OutliersDetectionModule.make size
// gate (comatip with TramoUtility.autlar), the ModelDescription pre-stage
// (AverageInterpolator.cleanMissings with the retained missing positions,
// the log transformation), the differenced RegArmaModel (response through
// the pinned differencing backfilter, design = interpolated point columns,
// the constant column when the mean is estimated, the pre-registered
// indicator columns), FastOutliersDetector.estimateModel at round 0 (Ols
//Computer residuals, HannanRissanen start with stabilization and the round-0
// installation, the warm maximum-likelihood seam when the tail demands it),
// the FastOutlierDetector initmodel (DEFAULT_COMPUTER generalized least
// squares on the whitened design, the stationary-model Kalman innovations,
// Mad2(50)), the processOutlier T/coefficient tables over the pi * filter
// weight columns, the searchMax sweep and the |maxT| >= cv gate.
#include "kalos/tramo_outlier_scan.hpp"

#include "tramo_scan_machine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {
namespace {

namespace machine = kalos::detail;

constexpr int kPeriod = machine::kScanPeriod;

using machine::clean_missings;
using machine::ScanOrders;
using machine::tramo_autlar;
using machine::tramo_comatip;
using machine::mad2_default;
using machine::calc_residuals;
using machine::raw_regular;
using machine::raw_seasonal;
using machine::hr_to_jvm;
using machine::ScanFactory;
using machine::register_factories;

double qnan() { return machine::scan_qnan(); }

} // namespace

TramoOutlierScanDecision scan_tramo_outliers(const TramoOutlierScanRequest& request) {
    if (request.frequency != kPeriod) {
        throw std::invalid_argument(
            "scan_tramo_outliers: monthly data only (register D11)");
    }
    const int p = request.p, d = request.d, q = request.q;
    const int bp = request.bp, bd = request.bd, bq = request.bq;
    if (p < 0 || d < 0 || q < 0 || bp < 0 || bd < 0 || bq < 0) {
        throw std::invalid_argument("scan_tramo_outliers: negative orders are rejected");
    }
    if (d + bd == 0) {
        throw std::invalid_argument(
            "scan_tramo_outliers: a positive differencing order is required");
    }
    if (request.y.empty()) {
        throw std::invalid_argument("scan_tramo_outliers: an empty series is rejected");
    }
    if (request.tc_rate <= 0.0 || !(request.tc_rate < 1.0) || !std::isfinite(request.tc_rate)) {
        throw std::invalid_argument("scan_tramo_outliers: the tc rate must lie in (0, 1)");
    }
    if (request.lower_bound > request.upper_bound) {
        throw std::invalid_argument(
            "scan_tramo_outliers: resolved span bounds must be ordered");
    }
    detail::validate_scan_regressors(request);
    for (const auto& pre : request.prespec) {
        if (pre.second == TramoOutlierKind::So) {
            throw std::invalid_argument("scan_tramo_outliers: a pre-registered SO is "
                                        "outside the transcribed pinned surface");
        }
        const int kind = static_cast<int>(pre.second);
        if (kind < 0 || kind > 3 || pre.first < 0) {
            throw std::invalid_argument("scan_tramo_outliers: unknown prespec variable");
        }
    }
    if (!request.interpolate) {
        for (const double v : request.y) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument("scan_tramo_outliers: missing observations "
                                            "require interpolation");
            }
        }
    }

    TramoOutlierScanDecision decision;
    decision.status = TramoOutlierScanDecision::Status::Unchanged;
    decision.position = -1;
    decision.type = -1;
    decision.max_t = qnan();
    decision.coefficient = qnan();
    decision.factories = 0;
    decision.lower_bound = request.lower_bound;
    decision.upper_bound = request.upper_bound;
    decision.critical_value = request.cv;

    // Registered factories in the pinned order (echoed even for a rejected
    // series: the module would build the same mix).
    const std::vector<ScanFactory> factories = register_factories(
        request.ao, request.ls, request.tc, request.tc_rate, request.so);
    decision.factories = static_cast<int>(factories.size());

    const std::size_t n = request.y.size();

    // make(): the comatip gate over the finite count of the raw series; a
    // positive verdict also forces the maximum-likelihood flag.
    int n_finite = 0;
    for (const double v : request.y) {
        if (std::isfinite(v)) {
            ++n_finite;
        }
    }
    const ScanOrders orders{p, d, q, bp, bd, bq};
    const int matip = tramo_comatip(n_finite, orders, request.mean,
        request.maximum_likelihood);
    if (matip < 0) {
        decision.status = TramoOutlierScanDecision::Status::Failed;
        return decision;
    }
    const bool mvx = matip > 0;

    // calc() switches the scan off on an empty mix or collapsed bounds (the
    // module then reports the empty selection: Unchanged).
    if (factories.empty() || request.upper_bound <= request.lower_bound) {
        return decision; // Unchanged
    }

    // The pre-stage: AverageInterpolator.cleanMissings (the positions are
    // retained as the model missings) and the log transformation.
    std::vector<double> z = request.y;
    std::vector<int> missing_at;
    {
        bool any_missing = false;
        for (const double v : z) {
            if (!std::isfinite(v)) {
                any_missing = true;
                break;
            }
        }
        if (any_missing) {
            if (!clean_missings(z, missing_at)) {
                decision.status = TramoOutlierScanDecision::Status::Failed;
                return decision;
            }
            std::sort(missing_at.begin(), missing_at.end()); // JVM IntList: descending
        }
    }
    if (request.log) {
        for (const double v : z) {
            if (!(v > 0.0)) {
                throw std::invalid_argument(
                    "scan_tramo_outliers: the log requires strictly positive observations");
            }
        }
        for (double& v : z) {
            v = detail::fdlibm_log(v);
        }
    }

    // The differenced model: response and design through the pinned
    // differencing; design = missing-point columns, the constant when the
    // mean is estimated, the pre-registered variable columns (in-domain
    // only; an out-of-domain variable is tagged excluded and never enters).
    const std::vector<double> kernel =
        detail::pinned_differencing_polynomial(kPeriod, d, bd);
    const std::size_t deg = kernel.size() - 1;
    const std::vector<double> dly = detail::pinned_backfilter_apply(kernel, z);
    if (dly.empty()) {
        decision.status = TramoOutlierScanDecision::Status::Failed;
        return decision;
    }
    std::vector<std::vector<double>> xcols;
    for (const int m : missing_at) {
        std::vector<double> col(z.size(), 0.0);
        col[static_cast<std::size_t>(m)] = 1.0;
        xcols.push_back(detail::pinned_backfilter_apply(kernel, col));
    }
    if (request.mean) {
        xcols.push_back(std::vector<double>(dly.size(), 1.0));
    }
    auto user_columns = request.calendar_regressors;
    std::vector<std::pair<std::size_t, int>> prespec_cells; // (row, factory column)
    for (const auto& pre : request.prespec) {
        const int kind = static_cast<int>(pre.second);
        if (static_cast<std::size_t>(pre.first) >= n) {
            continue; // excluded by ModelDescription; the module exclude no-ops
        }
        user_columns.push_back(machine::factory_raw(kind, request.tc_rate, pre.first, n));
        const int row = pre.first;
        int col = -1;
        for (std::size_t o = 0; o < factories.size(); ++o) {
            if (factories[o].kind == kind) {
                col = static_cast<int>(o);
                break;
            }
        }
        if (col >= 0) { // outlierType < 0 (factory missing) leaves exclude no-op
            prespec_cells.emplace_back(static_cast<std::size_t>(row), col);
        }
    }

    for (const auto& column : user_columns)
        xcols.push_back(detail::pinned_backfilter_apply(kernel, column));

    // estimateModel at round 0: the OlsComputer residuals (a null result is
    // "Ols failed"), then the HannanRissanen start with stabilization; a
    // fitted or zero-variance start is installed (round 0 always installs)
    // and a stabilizer rewrite (or the maximum-likelihood flag) sends the
    // tail through the warm re-estimation seam.
    std::vector<double> res = dly;
    if (!xcols.empty()) {
        const auto ols = detail::pinned_ols_coefficients(dly, xcols);
        if (!ols) {
            decision.status = TramoOutlierScanDecision::Status::Failed;
            return decision;
        }
        res = calc_residuals(dly, xcols, *ols);
    }

    // The SarimaSpec builder fill at JVM coordinates (phi -0.1, theta and
    // btheta -0.2), the model estimateModel keeps when HannanRissanen fails.
    std::vector<double> jvm;
    for (int i = 0; i < p + bp; ++i) jvm.push_back(-0.1);
    for (int i = 0; i < q + bq; ++i) jvm.push_back(-0.2);

    bool need_optimize = mvx;
    if (p + bp + q + bq > 0) {
        const HannanRissanenStart hr = hannan_rissanen_start(
            res, RegressionDesign{0, {}, {}, {}},
            HannanRissanenOrders{p, bp, q, bq, kPeriod});
        if (hr.source != HannanRissanenStart::Source::Defaulted) {
            // hr.coefficients is the stabilized landing (the JVM installs
            // SarimaMapping.stabilize(hr.getModel())); the landing itself may
            // sit outside the stability region, in which case the stabilize
            // was not the identity (flag6 = 0) and the detector tail runs the
            // warm Levenberg-Marquardt optimizeModel from the installed model.
            jvm = hr_to_jvm(hr.coefficients);
            if (hr.stabilized) {
                need_optimize = true;
            }
        }
    }
    if (need_optimize) {
        SarimaEstimationRequest warm{};
        machine::set_scan_warm_design(warm, z, missing_at, request.mean, user_columns);
        warm.blocks[0].initial.assign(static_cast<std::size_t>(p), 0.0);
        warm.blocks[0].fixed.assign(static_cast<std::size_t>(p), 0);
        warm.blocks[1].initial.assign(static_cast<std::size_t>(q), 0.0);
        warm.blocks[1].fixed.assign(static_cast<std::size_t>(q), 0);
        warm.blocks[2].initial.assign(static_cast<std::size_t>(bp), 0.0);
        warm.blocks[2].fixed.assign(static_cast<std::size_t>(bp), 0);
        warm.blocks[3].initial.assign(static_cast<std::size_t>(bq), 0.0);
        warm.blocks[3].fixed.assign(static_cast<std::size_t>(bq), 0);
        warm.period = kPeriod;
        warm.d = d;
        warm.seasonal_d = bd;
        warm.control = EstimationControl{request.precision, 1e-4, true, false, true, true,
            EstimationStartingPoint::HannanRissanen,
            EstimationSolver::PinnedLevenbergMarquardt};
        std::array<std::vector<double>, 4> start;
        {
            const auto at = [&](std::size_t o) {
                return jvm[static_cast<std::size_t>(o)];
            };
            for (int i = 0; i < p; ++i) start[0].push_back(-at(i));
            for (int i = 0; i < bp; ++i) start[2].push_back(-at(p + i));
            for (int i = 0; i < q; ++i) start[1].push_back(at(p + bp + i));
            for (int i = 0; i < bq; ++i) start[3].push_back(at(p + bp + q + i));
        }
        const SarimaEstimationResult est = optimize_sarima_from(warm, start);
        const auto& c = est.coefficients; // [phi, theta, bphi, btheta], Kalos
        for (const auto& block : c) {
            for (const double value : block) {
                if (!std::isfinite(value)) {
                    decision.status = TramoOutlierScanDecision::Status::Failed;
                    return decision;
                }
            }
        }
        jvm.clear();
        for (const double v : c[0]) jvm.push_back(-v);
        for (const double v : c[2]) jvm.push_back(-v);
        for (const double v : c[1]) jvm.push_back(v);
        for (const double v : c[3]) jvm.push_back(v);
    }

    // FastOutlierDetector.initmodel: the DEFAULT_COMPUTER regression over
    // the whitened design (a rank failure is "Gls failed"), the regression
    // residuals filtered through the stationary model (the design-free
    // Kalman machine when allCoefficients is empty) and the Mad2 scale.
    const std::size_t p_sz = static_cast<std::size_t>(p);
    const std::size_t bp_sz = static_cast<std::size_t>(bp);
    const std::size_t q_sz = static_cast<std::size_t>(q);
    const std::size_t bq_sz = static_cast<std::size_t>(bq);
    std::vector<double> el;
    if (!xcols.empty()) {
        const auto gls = detail::pinned_gls_coefficients(
            dly, p_sz, bp_sz, q_sz, bq_sz, kPeriod, xcols, jvm);
        if (!gls) {
            decision.status = TramoOutlierScanDecision::Status::Failed;
            return decision;
        }
        const std::vector<double> wr = calc_residuals(dly, xcols, *gls);
        const auto innovations = detail::pinned_fast_kalman_residuals(
            wr, p_sz, bp_sz, q_sz, bq_sz, kPeriod, jvm);
        if (!innovations) {
            decision.status = TramoOutlierScanDecision::Status::Failed;
            return decision;
        }
        el = *innovations;
    } else {
        const auto innovations = detail::pinned_kalman_innovations(
            dly, p_sz, bp_sz, q_sz, bq_sz, kPeriod, jvm);
        if (!innovations) {
            decision.status = TramoOutlierScanDecision::Status::Failed;
            return decision;
        }
        el = *innovations;
    }
    const double mad = mad2_default(el);

    // The scan tables: SingleOutlierDetector.prepare rows (the estimation
    // domain length, matching el + differencing degree), zones per factory
    // and the exclusions (missing cells for every factory column, the
    // in-domain pre-registered cells at their own column).
    const std::size_t m = el.size() + deg;
    if (m != n) {
        decision.status = TramoOutlierScanDecision::Status::Failed;
        return decision; // defensive: the pinned prepare uses the domain length
    }
    const std::size_t k = factories.size();
    std::vector<char> allowed(m * k, 0);
    for (std::size_t o = 0; o < k; ++o) {
        const long long zs = factories[o].zone_start;
        const long long ze = factories[o].zone_end;
        for (std::size_t row = 0; row < m; ++row) {
            allowed[row + o * m] = static_cast<long long>(row) >= zs
                && static_cast<long long>(row) < static_cast<long long>(m) - ze;
        }
    }
    for (const int mp : missing_at) {
        for (std::size_t o = 0; o < k; ++o) {
            allowed[static_cast<std::size_t>(mp) + o * m] = 0; // exclude guard holds
        }
    }
    for (const auto& cell : prespec_cells) {
        allowed[cell.first + static_cast<std::size_t>(cell.second) * m] = 0;
    }

    // The table fill and the winner search run on the shared machine
    // (src/tramo_scan_machine.hpp), the same kernels the item 23-F F.2 rounds
    // driver calls every round; live folds the zone/exclusion mask with the
    // domain bounds.
    std::vector<char> live(m * k, 0);
    for (std::size_t o = 0; o < k; ++o) {
        for (std::size_t row = 0; row < m; ++row) {
            live[row + o * m] = static_cast<char>(allowed[row + o * m]
                && static_cast<long long>(row) >= request.lower_bound
                && static_cast<long long>(row) < request.upper_bound);
        }
    }
    machine::ScanTables tables;
    machine::build_scan_tables(factories, jvm, p_sz, bp_sz, q_sz, bq_sz, kernel,
        el, mad, live, true, tables);
    int idx = -1;
    double dmax = 0.0;
    machine::search_max(tables, m, idx, dmax);
    if (idx == -1) {
        decision.status = TramoOutlierScanDecision::Status::Failed;
        return decision;
    }
    const std::vector<double>& t_table = tables.t_table;

    const int pos_max = idx % static_cast<int>(m);
    const int o_max = idx / static_cast<int>(m);
    const double max_t = t_table[static_cast<std::size_t>(idx)];
    // FastOutliersDetector.process gate: Math.abs(maxT) < cv breaks (the
    // dcmpg comparison lets a NaN pass through).
    if (std::isless(std::abs(max_t), request.cv)) {
        return decision; // Unchanged
    }
    decision.status = TramoOutlierScanDecision::Status::Selected;
    decision.position = pos_max;
    decision.type = o_max;
    decision.max_t = max_t;
    decision.coefficient = tables.coef_table[static_cast<std::size_t>(idx)];
    return decision;
}

} // namespace kalos
