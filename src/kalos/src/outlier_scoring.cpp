#include "kalos/outlier_scoring.hpp"

#include "outlier_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {

namespace {

using outlier_detail::dcmpl;
using outlier_detail::dcmpg;
using outlier_detail::java_max;
using outlier_detail::robust_lu_solve;

// Pinned jdplus Constants.getEpsilon() / 1/sqrt(2*pi).
constexpr double kEpsilon = 1.1102230246251565e-16;
constexpr double kSqrpi = 0.3989422804014327;

// ---- Utility.intProbability (lower tail), Cody regions ----------------------

double int_probability_lower(double x) {
    static constexpr double A[5] = {2.2352520354606837, 161.02823106855587,
        1067.6894854603709, 18154.98125334356, 0.06568233791820745};
    static constexpr double B[4] = {47.202581904688245, 976.0985517377767,
        10260.932208618979, 45507.78933502673};
    static constexpr double C[9] = {0.39894151208813466, 8.883149794388377,
        93.50665613217785, 597.2702763948002, 2494.5375852903726,
        6848.190450536283, 11602.65143764735, 9842.714838383978,
        1.0765576773720192e-8};
    static constexpr double D[8] = {22.266688044328117, 235.387901782625,
        1519.3775994075547, 6485.558298266761, 18615.571640885097,
        34900.95272114598, 38912.00328609327, 19685.429676859992};
    // Pinned Utility.P/Q (as filled in <clinit>): the asymptotic numerator
    // P carries 0.023073441764940174 as its leading (index 5) coefficient.
    static constexpr double P[6] = {0.215898534057957, 0.12740116116024736,
        0.022235277870649807, 0.0014216191932278934, 2.9112874951168793e-5,
        0.023073441764940174};
    static constexpr double Q[5] = {1.284260096144911, 0.4682382124808651,
        0.06598813786892856, 0.0037823963320275824, 7.297515550839662e-5};

    const double ax = std::fabs(x);
    double res;
    if (ax <= 0.66291) { // THRSH
        const double z2 = x == 0.0 ? 0.0 : x * x;
        double num = A[4] * z2;
        double den = z2;
        for (int k = 0; k < 3; ++k) {
            num = (num + A[k]) * z2;
            den = (den + B[k]) * z2;
        }
        res = x * (num + A[3]) / (den + B[3]) + 0.5;
    } else if (ax <= std::sqrt(32.0)) { // ROOT32
        double num = C[8] * ax;
        double den = ax;
        for (int k = 0; k < 7; ++k) {
            num = (num + C[k]) * ax;
            den = (den + D[k]) * ax;
        }
        double r = (num + C[7]) / (den + D[7]);
        const double i = std::floor(1.6 * ax) / 1.6;
        const double f = (ax - i) * (ax + i);
        r = std::exp(-i * i * 0.5) * std::exp(-f * 0.5) * r;
        res = x > 0 ? 1.0 - r : r;
    } else {
        const double t = 1.0 / (x * x);
        double num = P[5] * t;
        double den = t;
        for (int k = 0; k < 4; ++k) {
            num = (num + P[k]) * t;
            den = (den + Q[k]) * t;
        }
        double r = t * (num + P[4]) / (den + Q[4]);
        r = (kSqrpi - r) / ax;
        const double i = std::floor(1.6 * x) / 1.6; // Signed x in this region.
        const double f = (x - i) * (x + i);
        r = std::exp(-i * i * 0.5) * std::exp(-f * 0.5) * r;
        res = x > 0 ? 1.0 - r : r;
    }
    return res;
}

// ---- Utility.intProbabilityInverse (AS 241 seed + Newton) --------------------

double calc_poly(const double* a, std::size_t len, double x) {
    double acc = a[len - 1];
    for (std::size_t i = len - 1; i-- > 0;) acc = acc * x + a[i];
    return acc;
}

double stn_val(double t) {
    static constexpr double XNUM[5] = {-0.322232431088, -1.0, -0.342242088547,
        -0.0204231210245, -4.53642210148e-5};
    static constexpr double XDENOM[5] = {0.099348462606, 0.588581570495,
        0.531103462366, 0.10353775285, 0.0038560700634};
    double sgn;
    double u;
    if (dcmpg(t, 0.5) > 0) {
        sgn = 1.0;
        u = 1.0 - t;
    } else {
        sgn = -1.0;
        u = t;
    }
    const double w = std::sqrt(-2.0 * std::log(u));
    return (w + calc_poly(XNUM, 5, w) / calc_poly(XDENOM, 5, w)) * sgn;
}

double int_probability_inverse(double p) {
    constexpr double TOL = 1e-13;
    const double t = dcmpg(p, 0.5) > 0 ? 1.0 - p : p;
    double z = stn_val(t);
    for (int it = 0; it < 200; ++it) {
        const double err = int_probability_lower(z) - t;
        const double step = err / (kSqrpi * std::exp(-0.5 * z * z));
        z = z - step;
        if (!(dcmpg(std::fabs(step / z), TOL) > 0)) return p > 0.5 ? -z : z;
    }
    throw std::overflow_error("Too many iterations");
}

// ---- CriticalValueComputer ---------------------------------------------------

double calc_val(int n, double alpha) {
    if (n == 1) return 1.96;
    const double dn = static_cast<double>(n);
    const double k = 2.0 - std::sqrt(1.0 + alpha);
    const double d7 = std::sqrt(2.0 * std::log(dn));
    // Pinned stack order: the divisor pushed last (2 * v7) divides the earlier
    // S, and v9 = v7 - S / (2 * v7).
    const double d9 = (std::log(std::log(dn)) + std::log(12.566370614359172)) / (2.0 * d7);
    const double d11 = -std::log(-0.5 * std::log(k));
    return d11 / d7 + (d7 - d9);
}

double advanced_critical_value(int n, double alpha) {
    if (n == 1) {
        // Normal.getProbabilityInverse(alpha / 2, Upper): p <- 1 - p, then the
        // lower-tail inverse.
        return int_probability_inverse(1.0 - alpha / 2.0);
    }
    const double dn = static_cast<double>(n);
    static constexpr int nodes[3] = {2, 100, 200};
    double v[3] = {0.0, 0.0, 0.0};
    // M is 3x3 column-major storage [col * 3 + row].
    std::vector<double> M(9, 0.0);
    for (int r = 0; r < 3; ++r) {
        M[0 * 3 + static_cast<std::size_t>(r)] = 1.0;
        M[2 * 3 + static_cast<std::size_t>(r)] =
            std::sqrt(2.0 * std::log(static_cast<double>(nodes[r])));
        M[1 * 3 + static_cast<std::size_t>(r)] =
            (std::log(std::log(static_cast<double>(nodes[r]))) + std::log(12.566370614359172)) /
            (2.0 * M[2 * 3 + static_cast<std::size_t>(r)]);
    }
    v[0] = int_probability_inverse((1.0 + std::sqrt(1.0 - alpha)) / 2.0);
    for (int r = 1; r < 3; ++r) v[r] = calc_val(nodes[r], alpha);
    std::vector<double> rhs(v, v + 3);
    robust_lu_solve(M, 3, rhs);
    const double d9 = std::sqrt(2.0 * std::log(dn));
    const double d11 = (std::log(std::log(dn)) + std::log(12.566370614359172)) / (2.0 * d9);
    return rhs[0] + rhs[1] * d11 + rhs[2] * d9;
}

} // namespace

double standard_normal_lower_cdf(double x) {
    return int_probability_lower(x);
}

double standard_normal_lower_quantile(double p) {
    // Pinned getProbabilityInverse guard: dcmpg on both bounds, so a NaN p
    // passes the guard and only escapes through the 200-iteration loop.
    if (!(dcmpg(p, 1e-15) >= 0 && dcmpg(1.0 - p, 1e-15) >= 0)) {
        throw std::invalid_argument("probability out of range for the pinned inverse");
    }
    return int_probability_inverse(p);
}

double default_critical_value(std::size_t observations) {
    if (observations == 0) {
        throw std::invalid_argument("no critical value is defined for zero observations");
    }
    if (observations > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("observation count exceeds the pinned int domain");
    }
    const double advanced =
        advanced_critical_value(static_cast<int>(observations), 0.05); // X13Utility alpha.
    return java_max(2.8, advanced);                                    // X13Utility calcCv.
}

// ---- Factories ----------------------------------------------------------------

double transitory_change_rate(int frequency, double monthly_rate) {
    if (frequency == 12) return monthly_rate;
    if (frequency <= 0 || 12 % frequency != 0) {
        throw std::invalid_argument("Invalid frequnecy"); // Pinned message, sic.
    }
    return std::pow(monthly_rate, 12 / frequency);
}

OutlierZone outlier_excluding_zone(OutlierType type, int period) {
    switch (type) {
    case OutlierType::Additive: return OutlierZone{0, 0};
    case OutlierType::LevelShift: return OutlierZone{1, 1};
    case OutlierType::TransitoryChange: return OutlierZone{0, 1};
    case OutlierType::SeasonalOutlier:
        if (period <= 1) {
            throw std::invalid_argument("the seasonal outlier factory requires period > 1");
        }
        return OutlierZone{static_cast<std::size_t>(2 * period),
            static_cast<std::size_t>(2 * period)};
    }
    throw std::invalid_argument("unknown outlier type");
}

void fill_outlier_column(OutlierType type, std::size_t position, int period,
    double transitory_rate, std::vector<double>& column) {
    const std::size_t n = column.size();
    if (position >= n) throw std::invalid_argument("outlier position outside the column");
    switch (type) {
    case OutlierType::Additive: // AdditiveOutlierFactory.fill: block.set(pos, 1).
        column[position] = 1.0;
        return;
    case OutlierType::LevelShift: // LevelShiftFactory (zeroEnded): -1 on [0, pos), 0 on [pos, n).
        for (std::size_t i = 0; i < position; ++i) column[i] = -1.0;
        for (std::size_t i = position; i < n; ++i) column[i] = 0.0;
        return;
    case OutlierType::SeasonalOutlier: {
        if (period <= 1) throw std::invalid_argument("the seasonal factory requires period > 1");
        const double neg = -1.0 / static_cast<double>(period - 1);
        // PeriodicOutlierFactory (zeroEnded): backward from position - 1 in
        // runs of period - 1 entries of neg then one 1.0 entry.
        std::size_t pos = position;
        int k = 1;
        while (true) {
            if (k < period && pos != 0) {
                --pos;
                if (pos < n) column[pos] = neg;
                ++k;
                continue;
            }
            if (pos == 0) return;
            --pos;
            if (pos < n) column[pos] = 1.0;
            k = 1;
        }
    }
    case OutlierType::TransitoryChange: {
        // TransitoryChangeFactory.fill: 1 at the position, then v *= coef;
        // a next value below ZERO = 1e-15 in magnitude is never written.
        constexpr double ZERO = 1.0e-15;
        double v = 1.0;
        for (std::size_t i = position; i < n; ++i) {
            column[i] = v;
            v *= transitory_rate;
            if (dcmpg(std::fabs(v), ZERO) < 0) return;
        }
        return;
    }
    }
    throw std::invalid_argument("unknown outlier type");
}

// ---- t statistics ---------------------------------------------------------------

std::vector<double> regression_t_statistics(
    const RegressionLikelihood& likelihood, std::size_t concentrated_dimension) {
    if (likelihood.status != "available") {
        throw std::invalid_argument(
            "t statistics require an available regression likelihood (status \"" +
            likelihood.status + "\")");
    }
    // Pinned ConcentratedLikelihood.tstats(0, false) vector path: stdev[i] =
    // sqrt(covdiag[i] * ssq / dim) for every coefficient (no per-entry zero
    // or NaN special cases; NaN and zeros propagate through the division).
    std::vector<double> tstats(likelihood.coefficients.size(), 0.0);
    const double df = static_cast<double>(concentrated_dimension);
    for (std::size_t i = 0; i < likelihood.coefficients.size(); ++i) {
        const double c = likelihood.unscaled_covariance.empty()
            ? std::nan("")
            : likelihood.unscaled_covariance[i * likelihood.coefficients.size() + i];
        const double stdev = std::sqrt(c * likelihood.sum_squared_standardized / df);
        tstats[i] = likelihood.coefficients[i] / stdev;
    }
    return tstats;
}

// ---- Driver (ExactOutliersDetector) ----------------------------------------------

namespace {

struct ActiveOutlier {
    std::size_t position;
    OutlierType type;
};

} // namespace

OutlierScoringResult detect_outliers(const OutlierScoringRequest& request,
    const OutlierEstimateFunction& estimate, const OutlierOptimizeFunction& optimize,
    const Interrupt& interrupt) {
    if (request.types.empty()) {
        throw std::invalid_argument("outlier detection requires at least one candidate type");
    }
    for (std::size_t i = 0; i < request.types.size(); ++i) {
        for (std::size_t j = i + 1; j < request.types.size(); ++j) {
            if (request.types[i] == request.types[j]) {
                throw std::invalid_argument("duplicate outlier candidate type");
            }
        }
    }
    const std::size_t rows = request.base.series.values.size();
    if (std::any_of(request.base.series.missing.begin(), request.base.series.missing.end(),
            [](std::uint8_t flag) { return flag != 0; })) {
        throw std::invalid_argument(
            "the qualified outlier driver is defined for complete responses only");
    }
    if (std::any_of(request.base.design.fixed.begin(), request.base.design.fixed.end(),
            [](std::uint8_t flag) { return flag != 0; })) {
        throw std::invalid_argument(
            "the qualified outlier driver requires a fully free regression design");
    }
    if (request.upper_bound < request.lower_bound || request.upper_bound > rows) {
        throw std::invalid_argument("outlier scan bounds must satisfy l <= r <= length");
    }
    if (request.max_round == 0 || request.max_outliers == 0) {
        throw std::invalid_argument("the pinned round and outlier limits are positive");
    }

    OutlierScoringResult result;
    RegressionDesign design = request.base.design;
    std::array<std::vector<double>, 4> coefficients = request.initial_coefficients;
    std::vector<double> tstats; // Latest available estimation, free-column order.
    std::vector<ActiveOutlier> outliers;
    std::vector<std::uint8_t> excluded(rows * request.types.size(), 0);

    const auto call_estimate = [&](bool warm, const std::array<std::vector<double>, 4>& start) {
        SarimaEstimationRequest candidate = request.base;
        candidate.design = design;
        OutlierEstimation est =
            warm ? optimize(candidate, start, interrupt) : estimate(candidate, interrupt);
        if (est.available && est.t_statistics.size() != design.columns) {
            throw std::logic_error(
                "the estimation callback returned a t statistic vector that does not cover "
                "every design column");
        }
        return est;
    };
    const auto adopt = [&](const OutlierEstimation& est) {
        coefficients = est.coefficients;
        tstats = est.t_statistics;
    };
    const auto remove_outlier = [&](std::size_t index) {
        const std::size_t column = design.columns - outliers.size() + index;
        design.values.erase(design.values.begin() +
                static_cast<std::ptrdiff_t>(column) * static_cast<std::ptrdiff_t>(rows),
            design.values.begin() +
                static_cast<std::ptrdiff_t>(column + 1) * static_cast<std::ptrdiff_t>(rows));
        design.coefficients.erase(design.coefficients.begin() +
            static_cast<std::ptrdiff_t>(column));
        design.fixed.erase(design.fixed.begin() + static_cast<std::ptrdiff_t>(column));
        --design.columns;
        outliers.erase(outliers.begin() + static_cast<std::ptrdiff_t>(index));
        result.changed = true;
    };

    // process(): the opening estimateModel(map, true) whose result is
    // discarded; a failure keeps the incoming model, exactly pinned.
    {
        const auto est = call_estimate(false, coefficients);
        if (est.available) adopt(est);
    }

    // execute(): the forward-add loop.
    std::size_t rounds = 0;
    while (true) {
        OutlierScanRequest scan;
        scan.series = request.base.series;
        scan.design = design;
        scan.arma = coefficients;
        scan.period = request.base.period;
        scan.d = request.base.d;
        scan.seasonal_d = request.base.seasonal_d;
        scan.lower_bound = request.lower_bound;
        scan.upper_bound = request.upper_bound;
        scan.types = request.types;
        scan.transitory_rate = request.transitory_rate;
        scan.excluded = excluded;
        const auto found = scan_outlier_candidates(scan, interrupt);
        if (!found.available) {
            result.stop_reason = "scan_failed";
            break;
        }
        excluded = found.excluded_after;
        if (!(dcmpl(std::fabs(found.max_t_stat), request.critical_value) > 0)) {
            result.stop_reason = "below_cv";
            break;
        }
        ++rounds;
        // addNewOutlier: append the factory column last, remember and exclude.
        std::vector<double> column(rows, 0.0);
        fill_outlier_column(request.types[found.max_type], found.max_position,
            request.base.period, request.transitory_rate, column);
        design.values.insert(design.values.end(), column.begin(), column.end());
        design.coefficients.push_back(0.0);
        design.fixed.push_back(0);
        ++design.columns;
        outliers.push_back(ActiveOutlier{found.max_position, request.types[found.max_type]});
        excluded[found.max_position + found.max_type * rows] = 1;
        result.changed = true;

        OutlierDecision decision;
        decision.position = found.max_position;
        decision.type = request.types[found.max_type];
        decision.t_statistic = found.max_t_stat;
        decision.coefficient = found.coefficient_table[found.max_position +
            found.max_type * rows];

        const auto est = call_estimate(false, coefficients);
        if (!est.available) {
            // Pinned rollback: remove the column (the add's exclude() stays),
            // warm-reoptimize (discarded), then straight to the verifymodel
            // drain.
            remove_outlier(outliers.size() - 1);
            result.rejections.push_back(OutlierRejection{decision.position, decision.type,
                std::fabs(decision.t_statistic), rounds, OutlierRejection::Origin::AddRollback});
            const auto rollback = call_estimate(true, coefficients);
            if (rollback.available) adopt(rollback);
            result.stop_reason = "estimate_failed";
            break;
        }
        adopt(est);
        result.outliers.push_back(decision);
        if (!(rounds < request.max_round && outliers.size() < request.max_outliers)) {
            result.stop_reason = rounds >= request.max_round ? "max_rounds" : "max_outliers";
            break;
        }
    }

    // END: the verifymodel drain, warm re-estimating after every removal.
    while (!outliers.empty()) {
        const std::size_t off = design.columns - outliers.size();
        if (tstats.size() < design.columns) {
            throw std::logic_error("the driver lost track of the estimation t statistics");
        }
        std::size_t best = 0;
        for (std::size_t i = 1; i < outliers.size(); ++i) {
            if (dcmpg(std::fabs(tstats[off + i]), std::fabs(tstats[off + best])) < 0) {
                best = i; // Pinned dcmpg/ifge: a NaN never replaces the best; ties
                          // keep the earlier entry. A NaN incumbent stays best.
            }
        }
        const double abest = std::fabs(tstats[off + best]);
        if (dcmpl(abest, request.critical_value) >= 0) break; // NaN drops (dcmpl -1).
        OutlierRejection rejection;
        rejection.position = outliers[best].position;
        rejection.type = outliers[best].type;
        rejection.abs_t_statistic = abest;
        rejection.round = rounds;
        const std::size_t type_index = static_cast<std::size_t>(
            [&] {
                for (std::size_t j = 0; j < request.types.size(); ++j) {
                    if (request.types[j] == outliers[best].type) return j;
                }
                throw std::logic_error("unreachable");
            }());
        excluded[rejection.position + type_index * rows] = 0; // allow()
        rejection.origin = OutlierRejection::Origin::VerifyModel;
        result.rejections.push_back(rejection);
        // The kept-decision trace mirrors the internal list (adds appended in
        // column order, rollbacks never recorded), so the same index drops it.
        result.outliers.erase(result.outliers.begin() + static_cast<std::ptrdiff_t>(best));
        remove_outlier(best);
        const auto est = call_estimate(true, coefficients);
        if (!est.available) break;
        adopt(est);
    }

    result.rounds = rounds;
    result.coefficients = coefficients;
    result.t_statistics = tstats;
    return result;
}

} // namespace kalos
