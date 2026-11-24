#include "kalos/tramo_outlier_rounds.hpp"

#include <algorithm>
#include <cmath>

#include "tramo_scan_machine.hpp"

#include <stdexcept>

namespace kalos {

bool tramo_stationary_scan_supported(const TramoOutlierScanRequest& request) {
    const bool regular_ar = (request.p == 1 || request.p == 2) && request.q == 0 &&
                            request.bp == 0 && request.bq == 0;
    const bool seasonal_arma = ((request.p == 2 && request.q == 0) || (request.p == 1 && request.q == 1)) &&
                               request.bp == 1 && request.bq == 1;
    const bool short_ma = request.p == 0 && request.bp == 0 && request.bq == 0 && (request.q == 0 || request.q == 1);
    const bool seasonal_ar = request.p == 0 && request.q == 0 && request.bp == 1 && request.bq == 0;
    const bool last_chance = request.p == 3 && request.q == 1 && request.bp == request.bq &&
                             (request.bp == 0 || request.bp == 1);
    return request.frequency == 12 && request.d == 0 && request.bd == 0 &&
           (regular_ar || seasonal_arma || short_ma || seasonal_ar || last_chance) &&
           request.mean && !request.maximum_likelihood && !request.so && request.prespec.empty() &&
           request.calendar_regressors.empty() &&
           std::all_of(request.y.begin(), request.y.end(), [](double y) { return std::isfinite(y); });
}

TramoOutlierRoundsResult tramo_outlier_rounds(
    const TramoOutlierScanRequest& request) {
    // The 23-B static gate verbatim (same checks, same order, unit prefix).
    if (request.frequency != detail::kScanPeriod) {
        throw std::invalid_argument(
            "tramo_outlier_rounds: monthly data only (register D11)");
    }
    const int p = request.p, d = request.d, q = request.q;
    const int bp = request.bp, bd = request.bd, bq = request.bq;
    if (p < 0 || d < 0 || q < 0 || bp < 0 || bd < 0 || bq < 0) {
        throw std::invalid_argument(
            "tramo_outlier_rounds: negative orders are rejected");
    }
    if (d == 0 && bd == 0 && !tramo_stationary_scan_supported(request)) {
        throw std::invalid_argument(
            "tramo_outlier_rounds: unsupported stationary scan design");
    }
    if (request.y.empty()) {
        throw std::invalid_argument(
            "tramo_outlier_rounds: an empty series is rejected");
    }
    if (request.tc_rate <= 0.0 || !(request.tc_rate < 1.0)
        || !std::isfinite(request.tc_rate)) {
        throw std::invalid_argument(
            "tramo_outlier_rounds: the tc rate must lie in (0, 1)");
    }
    if (request.lower_bound > request.upper_bound) {
        throw std::invalid_argument(
            "tramo_outlier_rounds: resolved span bounds must be ordered");
    }
    detail::validate_scan_regressors(request);
    for (const auto& pre : request.prespec) {
        if (pre.second == TramoOutlierKind::So) {
            throw std::invalid_argument("tramo_outlier_rounds: a pre-registered "
                                        "SO is outside the transcribed pinned "
                                        "surface");
        }
        const int kind = static_cast<int>(pre.second);
        if (kind < 0 || kind > 3 || pre.first < 0) {
            throw std::invalid_argument(
                "tramo_outlier_rounds: unknown prespec variable");
        }
    }
    if (!request.interpolate) {
        for (const double v : request.y) {
            if (!std::isfinite(v)) {
                throw std::invalid_argument(
                    "tramo_outlier_rounds: missing observations require "
                    "interpolation");
            }
        }
    }

    TramoOutlierRoundsResult result;
    result.status = TramoOutlierRoundsResult::Status::Failed;
    result.note = "make";
    result.round = 0;

    // make(): the comatip size gate over the finite count of the raw series;
    // a positive verdict also forces the maximum-likelihood flag.
    int n_finite = 0;
    for (const double v : request.y) {
        if (std::isfinite(v)) {
            ++n_finite;
        }
    }
    const detail::ScanOrders orders{p, d, q, bp, bd, bq};
    const int matip = detail::tramo_comatip(
        n_finite, orders, request.mean, request.maximum_likelihood);
    if (matip < 0) {
        return result; // "make null": the module is Failed without a driver
    }

    detail::ScanMachine sm;
    sm.orders = orders;
    sm.mean = request.mean;
    sm.log = request.log;
    sm.mvx = request.maximum_likelihood || matip > 0;
    sm.precision = request.precision;
    sm.tc_rate = request.tc_rate;
    sm.cv = request.cv;
    sm.max_round = request.max_round;
    sm.max_outliers = request.max_outliers;
    sm.lower_bound = request.lower_bound;
    sm.upper_bound = request.upper_bound;
    sm.factories = detail::register_factories(
        request.ao, request.ls, request.tc, request.tc_rate, request.so);
    sm.user_raw = request.calendar_regressors;
    detail::prestage(sm, request.y, request.prespec,
        "tramo_outlier_rounds: the log requires strictly positive observations");
    if (sm.failed) {
        return result;
    }

    // The driver. Unlike the round-0 unit, an empty factory mix or collapsed
    // bounds does NOT short-circuit here: the witness pins estimateModel and
    // the drain running before SingleOutlierDetector.process reports false
    // (sc_spannone, note sodfalse).
    detail::process_rounds(sm, result.trace);
    for (const auto& ev : result.trace) {
        if (ev.tag == "manret") {
            result.status = ev.manret_returned
                ? TramoOutlierRoundsResult::Status::Ok
                : TramoOutlierRoundsResult::Status::Failed;
            result.note = ev.note;
            result.round = ev.manret_round;
        }
    }
    result.outliers = sm.outliers;
    return result;
}

} // namespace kalos
