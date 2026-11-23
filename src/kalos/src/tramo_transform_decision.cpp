// Item 23-B: the TRAMO pre-stage transform decision. Bytecode-pinned
// transcription of TramoUtility.calcCv, TramoSeasonalityDetector (the
// SeasonalityTests test(seq, freq, 1, true) machine: DifferencingResult.of,
// the Qs LjungBox variant at the seasonal lags, the Friedman month-rank test
// and the toInt ladder) and the TramoKernel pre-stage wiring (the m011
// downgrade remark and the LogLevelModule log-versus-level decision on the
// airline or m011 fixed-order model). The contract with the pinned-semantics
// notes and the fail-closed matrix is docs/contracts/tramo-transform-decision-
// contract.md; the JVM corpus is validation/evidence/native/tramo-transform.
#include <kalos/tramo_transform_decision.hpp>

#include <kalos/arma_diagnostics.hpp>
#include <kalos/detail/fdlibm.hpp>
#include <kalos/detail/special_functions.hpp>
#include <kalos/detail/ssq_deviances.hpp>
#include <kalos/estimation_control.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {

double tramo_critical_value(std::size_t observations) {
    if (observations > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("tramo_critical_value: observation count exceeds int");
    }
    // TramoUtility.calcCv delegates to CriticalValueComputer.simpleComputer
    // and clamps with Math.max(2.0, v); every branch below already exceeds
    // the floor, so the max is folded away (documented in the contract).
    const double n = static_cast<double>(static_cast<int>(observations));
    if (observations <= 50) return 3.3;
    if (observations < 450) return 3.3 + 0.0025 * (n - 50.0);
    return 4.3;
}

namespace {

constexpr int kPeriod = 12; // D11: the monthly-only envelope.

// DoubleSeq.delta(1, 1) followed by removeMean(): the first difference in
// original order, then every value less the running mean (the pinned
// DoubleSeq.average sum order).
std::vector<double> differenced_mean_removed(const std::vector<double>& y) {
    std::vector<double> d;
    d.reserve(y.size() > 0 ? y.size() - 1 : 0);
    for (std::size_t i = 1; i < y.size(); ++i) d.push_back(y[i] - y[i - 1]);
    double sum = 0.0;
    for (const double value : d) sum += value;
    const double mean = d.empty() ? 0.0 : sum / static_cast<double>(d.size());
    for (double& value : d) value -= mean;
    return d;
}

// LjungBox(seq).lag(12).autoCorrelationsCount(2).usePositiveAutoCorrelations
// ().value(): the sum over i = 1, 2 of ac(12i)^2 / (n - 12i) over the
// strictly positive autocorrelations, with the pinned early exit to exactly
// 0.0 when the FIRST seasonal autocorrelation is not strictly positive
// (dcmpl: a NaN autocorrelation takes this exit too).
double qs_statistic(const std::vector<double>& x) {
    const double n = static_cast<double>(x.size());
    double sum = 0.0;
    for (int i = 1; i <= 2; ++i) {
        const double ac = auto_correlation(x, 0.0, i * kPeriod);
        if (!(ac > 0.0)) { // Pinned dcmpl/ifgt; NaN skips.
            if (i == 1) return 0.0;
            continue;
        }
        sum += (ac * ac) / (n - static_cast<double>(i * kPeriod));
    }
    return sum * n * (n + 2.0); // LjungBox.build: value * n * (n + 2)
}

// Ranking.sort: the pinned midrank assignment. Items sort by Double.compare
// ascending (NaN last, -0.0 before +0.0) and runs of tied values (the pinned
// dcmpl==0 equality, where -0.0 == +0.0 and NaN never ties) share the average
// rank; a run of one takes its sorted position plus one. OUT is the pinned
// column-major matrix row: entry (r, c) lives at r + c*stride.
void ranking_sort(const double* row, int length, double* out, int stride) {
    std::vector<std::pair<double, int>> items(static_cast<std::size_t>(length));
    for (int i = 0; i < length; ++i) items[static_cast<std::size_t>(i)] = {row[i], i};
    std::stable_sort(items.begin(), items.end(),
        [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
            if (std::isnan(a.first)) return false;
            if (std::isnan(b.first)) return !std::isnan(a.first) ? true : false;
            if (a.first != b.first) return a.first < b.first;
            return std::signbit(a.first) && !std::signbit(b.first); // -0.0 < +0.0
        });
    int i = 0;
    while (i < length) {
        const int start = i;
        ++i;
        while (i < length && items[static_cast<std::size_t>(i)].first
                   == items[static_cast<std::size_t>(i - 1)].first) {
            ++i;
        }
        const int count = i - start;
        // Pinned: rank = (i0 as double) + (count + 1) * 0.5; a run of one
        // takes its sorted position plus one (identical value, pinned split).
        const double rank = count == 1 ? static_cast<double>(start + 1)
                                       : static_cast<double>(start)
                                             + static_cast<double>(count + 1) * 0.5;
        for (int k = start; k < i; ++k) {
            out[static_cast<std::size_t>(items[static_cast<std::size_t>(k)].second)
                * static_cast<std::size_t>(stride)] = rank;
        }
    }
}

// Friedman(seq, 12).build() with the default useChi2: the rank matrix over
// whole blocks (a NaN t yields a NaN p and an invalid test exactly like the
// pinned StatisticalTest.isFinite gate).
TramoSeasonalityTest friedman_test(const std::vector<double>& x) {
    const int m = static_cast<int>(x.size());
    const int n = m / kPeriod; // drop(m - n*period, 0): leading partial block
    const int skip = m - n * kPeriod;
    // The n x period FastMatrix in column-major storage: entry (r, c) sits at
    // r + c*n; rows are the ranking blocks, columns the months.
    std::vector<double> ranks(static_cast<std::size_t>(n) * kPeriod, 0.0);
    for (int r = 0; r < n; ++r) {
        ranking_sort(x.data() + skip + r * kPeriod, kPeriod, ranks.data() + r, n);
    }
    // FastMatrix.sum() and the column sums accumulate column-major.
    double total = 0.0;
    std::vector<double> colsum(static_cast<std::size_t>(kPeriod), 0.0);
    for (int c = 0; c < kPeriod; ++c) {
        for (int r = 0; r < n; ++r) {
            const double v = ranks[static_cast<std::size_t>(r)
                + static_cast<std::size_t>(c) * static_cast<std::size_t>(n)];
            total += v;
            colsum[static_cast<std::size_t>(c)] += v;
        }
    }
    const double nd = static_cast<double>(n);
    const double grand = total / (nd * static_cast<double>(kPeriod));
    double sst = 0.0;
    for (int c = 0; c < kPeriod; ++c) {
        const double d = colsum[static_cast<std::size_t>(c)] / nd - grand;
        sst += d * d;
    }
    sst *= nd;
    double ssq = 0.0;
    for (int c = 0; c < kPeriod; ++c) {
        for (int r = 0; r < n; ++r) {
            const double d = ranks[static_cast<std::size_t>(r)
                + static_cast<std::size_t>(c) * static_cast<std::size_t>(n)] - grand;
            ssq += d * d;
        }
    }
    const double sse = ssq / (nd * static_cast<double>(kPeriod - 1));
    TramoSeasonalityTest t;
    t.value = sst / sse;
    t.p_value = detail::chi2_upper_p(static_cast<double>(kPeriod - 1), t.value);
    t.valid = std::isfinite(t.value) && std::isfinite(t.p_value);
    return t;
}

// One log/level branch: the fixed-order (0,1,1)(0,1,1) or (0,1,1)(0,0,0)
// model with the meanCorrection constant column, estimated through the
// qualified pinned machine (TramoUtility.processor: HannanRissanen start with
// stabilization, sum-of-squares Levenberg-Marquardt at the module precision,
// maximum likelihood), and then re-evaluated at the landing coordinates
// through the pinned deviances machine to obtain ln(ssq * factor) - the
// ConcentratedLikelihoodWithMissing pair the module reads.
double transform_branch_ll(
    const std::vector<double>& y, bool seasonal, double precision) {
    SarimaEstimationRequest request;
    request.series.values = y;
    request.series.missing.assign(y.size(), 0);
    request.design.columns = 1;
    request.design.values.assign(y.size(), 1.0);
    request.design.coefficients = {0.0};
    request.design.fixed = {0};
    request.blocks[1].initial = {0.0}; // theta free
    request.blocks[1].fixed = {0};
    if (seasonal) {
        request.blocks[3].initial = {0.0}; // btheta free
        request.blocks[3].fixed = {0};
    }
    request.period = kPeriod;
    request.d = 1;
    request.seasonal_d = seasonal ? 1 : 0;
    request.control = EstimationControl{precision, 1e-4, true, false, true, true,
        EstimationStartingPoint::HannanRissanen, EstimationSolver::PinnedLevenbergMarquardt};

    const SarimaEstimationResult est = estimate_sarima_fixed(request, Interrupt{});
    const auto& c = est.coefficients; // [phi, theta, bphi, btheta], Kalos form
    for (const auto& block : c) {
        for (const double value : block) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("tramo transform decision: non-finite landing coefficients");
            }
        }
    }
    // JVM-convention [phi, bphi, theta, btheta] coordinates of the landing
    // point (the AR blocks cross the sign flip; here both are empty).
    std::vector<double> x;
    for (const double value : c[1]) x.push_back(value);
    for (const double value : c[3]) x.push_back(value);
    const std::vector<double> machine_y = detail::pinned_backfilter_apply(
        detail::pinned_differencing_polynomial(kPeriod, 1, seasonal ? 1 : 0), y);
    std::vector<double> design(machine_y.size(), 1.0);
    const auto point = detail::pinned_deviances_evaluate(
        machine_y, 0, 0, 1, seasonal ? 1 : 0, kPeriod, design, x);
    if (!point || !std::isfinite(point->ssq_e) || point->ssq_e <= 0.0) {
        throw std::runtime_error("tramo transform decision: deviances machine rejected the point");
    }
    return detail::fdlibm_log(point->ssq_e);
}

} // namespace

TramoSeasonalityReport tramo_seasonality(const std::vector<double>& values) {
    if (values.size() < 2) {
        throw std::invalid_argument("tramo_seasonality: at least two observations are required");
    }
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("tramo_seasonality: non-finite observations are rejected");
        }
    }
    TramoSeasonalityReport report;
    report.differenced = differenced_mean_removed(values);

    report.qs.value = qs_statistic(report.differenced);
    report.qs.p_value = detail::chi2_upper_p(2.0, report.qs.value);
    report.qs.valid = std::isfinite(report.qs.value) && std::isfinite(report.qs.p_value);
    report.friedman = friedman_test(report.differenced);

    // TramoSeasonalityDetector.hasSeasonality: the dcmpg ladder on the two
    // p-values (a NaN p takes neither branch).
    int qscore = 0, nscore = 0, score = 0, ost = 0;
    if (report.qs.p_value < 0.01) {
        qscore = 2;
        ++ost;
        ++score;
    } else if (report.qs.p_value < 0.05) {
        qscore = 1;
        ++score;
    }
    if (report.friedman.p_value < 0.01) {
        nscore = 2;
        ++ost;
        ++score;
    } else if (report.friedman.p_value < 0.05) {
        nscore = 1;
        ++score;
    }
    if (qscore == 2) {
        report.verdict = "Strong";
        report.to_int = 3;
    } else if (score == 2) {
        report.verdict = "Moderate";
        report.to_int = 2;
    } else if (ost == 1 || score == 1) {
        report.verdict = "Weak";
        report.to_int = 1;
    } else {
        report.verdict = "NotObservable";
        report.to_int = 0;
    }
    (void)nscore; // The non-parametric score point only feeds the ladder.
    return report;
}

TramoTransformDecisionResult decide_tramo_transform(
    const TramoTransformDecisionRequest& request) {
    if (request.period != kPeriod) {
        throw std::invalid_argument("decide_tramo_transform: only period 12 is qualified (D11)");
    }
    if (request.values.size() < 2) {
        throw std::invalid_argument("decide_tramo_transform: at least two observations are required");
    }
    for (const double value : request.values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "decide_tramo_transform: non-finite observations are rejected");
        }
    }
    if (!std::isfinite(request.precision) || request.precision <= 0.0) {
        throw std::invalid_argument("decide_tramo_transform: precision must be finite and positive");
    }
    if (!std::isfinite(request.log_preference)) {
        throw std::invalid_argument("decide_tramo_transform: log preference must be finite");
    }

    TramoTransformDecisionResult result;
    result.log_preference = request.log_preference;

    bool downgraded = false;
    if (request.detect_seasonality) {
        result.seasonality = tramo_seasonality(request.values);
        if (result.seasonality.to_int < 2) { // The TramoKernel downgrade path.
            result.remarks.push_back("No identifiable seasonality");
            result.seasonal = false;
            downgraded = true;
        } else {
            result.seasonal = request.spec_seasonality;
        }
    } else {
        result.seasonal = request.spec_seasonality;
    }
    // spec_pre reports the ModelDescription: only the detection downgrade
    // rewrites it (setSpecification(m011)); a fixed non-seasonal spec leaves
    // the description airline while the LogLevelModule flag is false.
    result.spec_pre = downgraded ? "(0,1,1)(0,0,0)" : "(0,1,1)(0,1,1)";

    // The pinned anyMatch gate: DoubleSeq values are rejected by v <= 0.0,
    // where NaN never triggers (dcmpg), and no estimation runs.
    const bool negative = std::any_of(request.values.begin(), request.values.end(),
        [](double v) { return v <= 0.0; });
    if (negative) {
        result.nonpositive = true;
        result.loglevel_status = "Unchanged";
        result.remarks.push_back("negative values, levels are chosen");
        return result;
    }

    // LogLevelModule.process: the LEVEL branch first (fields keep the
    // clear() defaults 1.0 / 0.0 / 0.0 until a branch writes them), then the
    // log transform (Math.log per value, the running mean as slog), then the
    // LOG branch with log = ln(ssq * factor) + 2 * slog. Any estimation
    // exception surfaces as Failed exactly like the pinned
    // RegSarimaModelling overload's RuntimeException catch.
    bool log_estimated = false;
    try {
        result.level_ll = transform_branch_ll(request.values, result.seasonal, request.precision);
        double slog = 0.0;
        std::vector<double> logs;
        logs.reserve(request.values.size());
        for (const double value : request.values) {
            const double lv = detail::fdlibm_log(value);
            logs.push_back(lv);
            slog += lv;
        }
        slog /= static_cast<double>(logs.size());
        result.log_correction = slog;
        const double ln = transform_branch_ll(logs, result.seasonal, request.precision);
        result.log_ll = ln + 2.0 * slog;
        log_estimated = true;
    } catch (const std::exception&) {
        result.loglevel_status = "Failed";
        result.log_transform = false;
        return result;
    }

    result.choosing_log = log_estimated
        && (result.log_ll + request.log_preference) < result.level_ll; // dcmpg-ifge
    result.loglevel_status = result.choosing_log ? "Changed" : "Unchanged";
    result.log_transform = result.choosing_log;
    return result;
}

} // namespace kalos
