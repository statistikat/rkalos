// Monthly TRAMO regression-test layer, including the Fast calendar tests,
// of the pinned jdplus 3.9.0 stack (item 23-E). The bytecode pin lives in
// docs/contracts/tramo-regression-tests-contract.md (decoded facts distilled
// in tmp_opencode/kalos-m5-pinning.md SECTION 25); the double-captured JVM
// corpus under validation/evidence/native/tramo-regression-tests/ qualifies
// every decision and every echoed double bit-for-bit.
//
// Reuse map (no new kernels): the Ljung-Box autocorrelation machine is the
// landed arma_diagnostics::auto_correlation (AutoCovariances
// .autoCorrelationFunction), its length rule arma_diagnostics::calc_lb_length
// (LjungBox.defaultAutoCorrelationsCount), the Chi-square upper tail is
// detail::chi2_upper_p (23-C.2), Math.log routes to detail::fdlibm_log
// (capture JVM runs with the dexp/dlog intrinsics disabled) and Math.sqrt is
// the exact IEEE-754 root. The t-statistic pair (ser/tstats) is transcribed
// here because the 21-G RegressionTestStatistics::tstats() computes the
// PLAIN arm while every TRAMO accept site reads the ADDMU arm (25.3).

#include "kalos/tramo_regression_tests.hpp"

#include <cmath>
#include <limits>

#include "kalos/arma_diagnostics.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/special_functions.hpp"

namespace kalos {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// JVM double comparison instructions: dcmpg pushes 1 and dcmpl pushes -1
// when either operand is NaN.
int dcmpg(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return 1;
    return a > b ? 1 : (a < b ? -1 : 0);
}

int dcmpl(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return -1;
    return a > b ? 1 : (a < b ? -1 : 0);
}

// java.lang.Math.max propagates a NaN argument (std::max would return the
// other operand).
double jvm_max(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return kNaN;
    return a > b ? a : b;
}

void require_valid_hand_over(const TramoTestStatistics& statistics) {
    if (statistics.cov_diagonal.size() != statistics.coefficients.size())
        throw std::invalid_argument(
            "the unscaledCovariance diagonal must match the coefficient count");
    if (statistics.dim < 0 || statistics.free_parameters < 0)
        throw std::invalid_argument(
            "negative likelihood counts are outside the tramo regression-tests envelope");
}

// TRegressionTest.accept over the given t-vector with the bytecode bucket
// rules: dcmpl(|t|, thigh) >= 0 counts high, dcmpl(|t|, tlow) >= 0 counts
// low; a NaN |t| compares below both thresholds and counts in NEITHER.
TramoTAcceptance t_accept(const std::vector<double>& t, int first, int dim,
                          double tlow, double thigh, int nsig) {
    if (first < 0 || dim < 0 || static_cast<std::size_t>(dim) > t.size() ||
        static_cast<std::size_t>(first) > t.size() - static_cast<std::size_t>(dim))
        throw std::invalid_argument(
            "the tested t-window lies outside the likelihood hand-over");
    TramoTAcceptance out;
    for (int j = 0; j < dim; ++j) {
        const double a = std::fabs(t[static_cast<std::size_t>(first + j)]);
        if (dcmpl(a, thigh) >= 0) {
            ++out.high;
        } else if (dcmpl(a, tlow) >= 0) {
            ++out.low;
        }
    }
    const int nsig_eff = nsig < dim ? nsig : dim;  // java.lang.Math.min(II)I
    out.accepted = out.high > 0 || out.low >= nsig_eff;
    return out;
}

}  // namespace

bool tramo_ami_check_mu(bool using_auto_model, bool mean_is_test) {
    return using_auto_model || mean_is_test;
}

bool tramo_round0_accepts(double ljung_box_pvalue, double limit,
                          bool accept_airline) {
    return accept_airline &&
           dcmpg(1.0 - ljung_box_pvalue, limit) < 0;
}

TramoReduceVa tramo_reduce_va(double curva, double reduce_va) {
    if (dcmpg(curva, 2.0) == 0) return {false, curva};
    return {true, jvm_max(2.0, curva * (1.0 - reduce_va))};
}

TramoOrders tramo_last_solution_orders(const TramoOrders& orders,
                                       bool seasonal) {
    TramoOrders next = orders;
    next.p = 3;
    if (orders.bd > 0) next.bp = 0;
    next.q = 1;
    if (seasonal) next.bq = 1;
    return next;
}

TramoEstimateModelRetry tramo_estimate_model_retry(
    const TramoEstimateModelRetryRequest& request) {
    if (request.pass < 0)
        throw std::invalid_argument(
            "the estimate-model retry pass must be non-negative");
    TramoEstimateModelRetry out;
    int retry = 0;  // the for-counter; increments on every loop-head pass
    for (;;) {
        if (out.estimations_consumed >= request.estimation_results.size())
            throw std::invalid_argument(
                "the estimate-model retry request ran out of estimation results");
        const bool estimated =
            request.estimation_results[out.estimations_consumed++];
        if (!estimated) {
            if (request.pass == 1 && request.need_outliers_desc &&
                request.need_auto_desc)
                out.reduce = tramo_reduce_va(request.curva, request.reduce_va);
            return out;  // completed == false
        }
        if (request.pass3 || request.pass != 0) {
            out.completed = true;
            return out;
        }
        if (retry >= 5) {
            out.completed = true;
            return out;
        }
        ++retry;
        if (out.tests_consumed >= request.test_results.size())
            throw std::invalid_argument(
                "the estimate-model retry request ran out of regression-test results");
        if (request.test_results[out.tests_consumed++]) {
            out.completed = true;
            return out;
        }
        // ifeq: a failed mean test re-enters the estimation (goto L0).
    }
}

TramoModelStatistics tramo_model_statistics(
    const TramoModelStatisticsRequest& request) {
    if (request.annual_frequency != 12)
        throw std::invalid_argument(
            "the tramo regression-tests envelope is monthly frequency only (D11)");
    if (request.residuals.empty())
        throw std::invalid_argument(
            "tramo model statistics needs a non-empty residual block");
    for (const double v : request.residuals) {
        if (!std::isfinite(v))
            throw std::invalid_argument(
                "non-finite residuals are outside the tramo model statistics envelope");
    }
    require_valid_hand_over({request.coefficients, request.cov_diagonal,
                             request.ssq, request.dim,
                             request.free_parameters});
    if (request.degrees_of_freedom < 0 || request.estimation_domain < 0 ||
        request.ami_outliers < 0)
        throw std::invalid_argument(
            "negative likelihood or description counts are outside the tramo "
            "regression-tests envelope");

    TramoModelStatistics s;
    const int n = static_cast<int>(request.residuals.size());
    const int k = request.free_parameters;
    s.lb_length = calc_lb_length(request.annual_frequency, n);
    s.outliers_count = request.ami_outliers;
    s.observations_count = request.estimation_domain;
    s.effective_observations_count = request.dim;

    const int df = s.lb_length - k;
    if (df <= 0) {
        // new LjungBox(...).build() feeds Chi2(df <= 0), whose constructor
        // throws with no catch on the main path of ModelStatistics.of: the
        // whole row set fails (the JVM witness is
        // `ms ms_craft_dfneg built 0 threw ...`). No row is invented.
        s.built = false;
        s.ljung_box = s.ljung_box_pvalue = s.se = s.bic = kNaN;
        s.seasonal_ljung_box = s.seasonal_ljung_box_pvalue = kNaN;
        return s;
    }
    s.built = true;

    // LjungBox(acf, n).autoCorrelationsCount(m).hyperParametersCount(k)
    // .build(): value = sum_{i=1..m} r_i^2 / (double)(n - i) (sign filter 0),
    // statistic = value * (double)n * (double)(n + 2), df = m - k.
    double value = 0.0;
    for (int i = 1; i <= s.lb_length; ++i) {
        const double r = auto_correlation(request.residuals, 0.0, i);
        value += (r * r) / static_cast<double>(n - i);
    }
    s.ljung_box =
        value * static_cast<double>(n) * static_cast<double>(n + 2);
    s.ljung_box_pvalue = detail::chi2_upper_p(static_cast<double>(df),
                                              s.ljung_box);

    // bic = Math.log(sigma2) + ((double)(k + nx) * Math.log(dim)) / dim
    s.bic = detail::fdlibm_log(request.sigma2) +
            (static_cast<double>(k + static_cast<int>(request.coefficients.size())) *
             detail::fdlibm_log(static_cast<double>(request.dim))) /
                static_cast<double>(request.dim);

    // se = Math.sqrt(ssq / (double)(degreesOfFreedom() - k))
    s.se = std::sqrt(request.ssq /
                     static_cast<double>(request.degrees_of_freedom - k));

    // Seasonal block (annual_frequency > 1, always 12 in this envelope):
    // new LjungBox(acf, n).autoCorrelationsCount(2).lag(freq).build(); the
    // df is 2 (lag != 1 skips the hyper-parameter subtraction), so the build
    // never throws and the JVM catch surface is unreachable for finite
    // residuals. An out-of-range lag contributes a zero autocorrelation
    // (witnessed live by slb = 0, slbp = 1 on n = 7, freq = 12).
    double seasonal_value = 0.0;
    for (int i = 1; i <= 2; ++i) {
        const double r = auto_correlation(
            request.residuals, 0.0, i * request.annual_frequency);
        seasonal_value +=
            (r * r) / static_cast<double>(n - i * request.annual_frequency);
    }
    s.seasonal_ljung_box = seasonal_value * static_cast<double>(n) *
                           static_cast<double>(n + 2);
    s.seasonal_ljung_box_pvalue =
        detail::chi2_upper_p(2.0, s.seasonal_ljung_box);
    return s;
}

TramoRefStatsGateResult tramo_ref_stats_gate(
    const TramoRefStatsGateRequest& request) {
    if (request.pass < 0 || request.round < 0)
        throw std::invalid_argument(
            "the refStats gate pass and round must be non-negative");
    if (request.reference_outliers < -1 || request.current_outliers < 0)
        throw std::invalid_argument(
            "the outlier tallies nb/nr must be non-negative (nr = -1 marks a null reference)");

    TramoRefStatsGateResult r;
    double lbp = request.ljung_box_complement;
    const double se = request.standard_error;

    // The current-model state (restore and the last-chance rewrite below).
    TramoOrders orders = request.orders;
    bool mean = request.mean;
    int outliers = request.current_outliers;

    if (request.reference_outliers >= 0) {  // refAuto != null
        const double rlbp = request.reference_ljung_box_complement;
        const double rse = request.reference_standard_error;
        bool revert = false;
        if (request.reference_outliers <= request.current_outliers) {
            const TramoOrders& o = request.orders;
            const bool rule_a =
                dcmpg(lbp, 0.95) < 0 && dcmpg(rlbp, 0.75) < 0 &&
                dcmpg(rse, se) < 0;
            const bool rule_b = request.pass == 1 &&
                                dcmpl(lbp, 0.95) >= 0 &&
                                dcmpg(rlbp, 0.95) < 0;
            const bool rule_c =  // dead subset of rule_a in the pinned build
                dcmpg(lbp, 0.95) < 0 && dcmpg(rlbp, 0.75) < 0 &&
                dcmpg(rlbp, lbp) < 0 && dcmpg(rse, se) < 0;
            const bool rule_d = dcmpl(lbp, 0.95) >= 0 &&
                                dcmpg(rlbp, 0.95) < 0 &&
                                dcmpg(rse, se) < 0;
            const bool shape_a = o.d == 0 && o.bd == 1 && o.p == 1 &&
                                 dcmpg(request.phi1, -0.82) <= 0 &&
                                 o.q <= 1 && o.bp == 0 && o.bq == 1;
            const bool shape_b = o.d == 1 && o.bd == 0 && o.p == 0 &&
                                 o.q == 1 && o.bp == 1 &&
                                 dcmpg(request.bphi1, -0.65) <= 0 &&
                                 o.bq <= 1;
            revert = rule_a || rule_b || rule_c || rule_d || shape_a ||
                     shape_b;
        }
        if (revert) {
            // restore(m): m.set(copyOf(refAuto.desc), refAuto.est)
            r.reverted = true;
            orders = request.reference_orders;
            mean = request.reference_mean;
            outliers = request.reference_outliers;
            lbp = rlbp;
        }
        // !revert: refAuto := copyOf(m), refStats := stats - the state
        // simply persists; nothing to model.
    }

    r.effective_ljung_box_complement = lbp;
    r.pcr = request.pcr;
    if (request.pass == 1)
        r.pcr += 0.025;
    else if (request.pass >= 2)
        r.pcr += 0.015;

    if (dcmpg(lbp, r.pcr) <= 0) {
        // Accepted: pass/round/curva and the need flags survive unchanged.
        r.accepted = true;
        r.pass = request.pass;
        r.round = request.round;
        r.curva = request.curva;
        r.orders_after = orders;
        r.mean_after = mean;
        r.outliers_after = outliers;
        return r;
    }
    if (request.pass == 1 && request.outliers_detection) {
        r.reduce = tramo_reduce_va(request.curva, request.reduce_va);
        r.remark = true;  // the bytecode calls log.remark unconditionally
                          // after the (possibly no-op) reduceVa()
        r.curva = r.reduce.curva;
    } else {
        r.curva = request.curva;
    }
    r.round = request.round + 1;
    r.pass = request.pass + 1;
    if (r.pass <= 2) {
        // Bytecode 534-544: needAutoModelling = !changed (iload_1; ifne ->
        // push 0, else push 1), needOutliers = isOutliersDetection().
        r.need_auto_modelling = !request.changed;
        r.need_outliers = request.outliers_detection;
    } else {
        // lastSolution(m): orders := {p 3, bp 0 if bd > 0, q 1, bq 1 if
        // seasonal}; removeVariable(isOutlier(v, true)); warning "last
        // chance model"; round := 1; needOutliers := isOutliersDetection();
        // needAutoModelling := false.
        r.last_solution = true;
        orders = tramo_last_solution_orders(orders, request.seasonal);
        outliers = 0;
        r.round = 1;
        r.need_outliers = request.outliers_detection;
        r.need_auto_modelling = false;
    }
    r.orders_after = orders;
    r.mean_after = mean;
    r.outliers_after = outliers;
    return r;
}

std::vector<double> tramo_test_tstats(const TramoTestStatistics& statistics,
                                      bool addmu) {
    require_valid_hand_over(statistics);
    const int nx = static_cast<int>(statistics.coefficients.size());
    std::vector<double> t(static_cast<std::size_t>(nx));
    if (nx == 0) return t;  // ConcentratedLikelihood.EMPTYARRAY
    // ser(nobs, addmu): n = addmu ? dim - nx - nobs : dim;
    // ser[j] = sqrt((vdiag[j] * ssq) / (double) n); t[j] = coef[j] / ser[j].
    const int nobs = statistics.free_parameters;
    const int n = addmu ? statistics.dim - nx - nobs : statistics.dim;
    for (int j = 0; j < nx; ++j) {
        const double ser = std::sqrt((statistics.cov_diagonal[static_cast<std::size_t>(j)] *
                                      statistics.ssq) /
                                     static_cast<double>(n));
        t[static_cast<std::size_t>(j)] =
            statistics.coefficients[static_cast<std::size_t>(j)] / ser;
    }
    return t;
}

TramoTAcceptance tramo_t_accepts(const TramoTestStatistics& statistics,
                                 int first, int dim, double tlow, double thigh,
                                 int nsig) {
    return t_accept(tramo_test_tstats(statistics, true), first, dim, tlow,
                    thigh, nsig);
}

TramoRegressionTestDecision run_tramo_regression_tests(
    const TramoRegressionTestRequest& request) {
    if (request.kind == TramoRegressionTestKind::fast) {
        if ((request.auto_trading_days && (request.trading_days_position < 0 || request.trading_days_dimension < 1))
            || (request.auto_length_of_period && request.length_of_period_position < 0)
            || (request.auto_easter && request.easter_position < 0))
            throw std::invalid_argument("fast calendar tests require the fitted column positions");
        if (request.auto_moving_holiday)
            throw std::invalid_argument(
                "a moving-holiday variable has no filter in the pinned fast regression "
                "test machine and no corpus witness; refused fail closed");
    } else {
        if (request.has_trading_days)
            throw std::invalid_argument(
                "a trading-days builder field in the default regression test machine "
                "belongs to item 21-I.3 (calendar branch)");
        if (request.has_length_of_period)
            throw std::invalid_argument(
                "a length-of-period builder field in the default regression test machine "
                "belongs to item 21-I.3 (calendar branch)");
        if (request.has_easter)
            throw std::invalid_argument(
                "an Easter builder field in the default regression test machine belongs "
                "to item 21-I.3 (calendar branch)");
    }
    if (request.period != 12)
        throw std::invalid_argument(
            "the tramo regression-tests envelope is monthly frequency only (D11)");
    require_valid_hand_over(request.statistics);

    TramoRegressionTestDecision decision;
    decision.mean = request.mean;
    decision.inert_trading_days = request.plain_trading_days;
    decision.inert_length_of_period = request.plain_length_of_period;
    decision.inert_easter = request.plain_easter;

    if (request.kind == TramoRegressionTestKind::def && !request.test_mean) {
        // DefaultRegressionTest: no calendar field (refused above otherwise)
        // and meanTest == null: the machine returns Unprocessed without
        // touching the model.
        decision.result = TramoRegressionTestResult::unprocessed;
        return decision;
    }

    const bool calendar = request.kind == TramoRegressionTestKind::fast &&
        (request.auto_trading_days || request.auto_length_of_period || request.auto_easter);
    const bool test_mean = request.test_mean && request.mean;
    if (!calendar && !test_mean) return decision;
    const auto t = tramo_test_tstats(request.statistics, true);
    if (calendar) {
        auto passes = [&](int first, int width) {
            return t_accept(t, first, width, request.variable_threshold, request.variable_threshold, 1).accepted;
        };
        bool drop_days = request.auto_trading_days && !passes(request.trading_days_position,
                                                               request.trading_days_dimension);
        // A significant LP retains the associated TD variable even when its
        // own t test failed. LP is retested only after an unsuccessful TD test.
        if (drop_days && request.auto_length_of_period) {
            drop_days = !passes(request.length_of_period_position, 1);
            decision.removed_length_of_period = drop_days;
        }
        decision.removed_trading_days = drop_days;
        decision.removed_easter = request.auto_easter && !passes(request.easter_position, 1);
        if (drop_days || decision.removed_easter) {
            decision.result = TramoRegressionTestResult::changed;
            decision.estimation_cleared = true;
        }
    }
    if (!test_mean) return decision;
    if (t.empty()) throw std::invalid_argument("the mean test requires a fitted mean coefficient");
    decision.mean_tstat = t[0];
    decision.mean_test_ran = true;
    const TramoTAcceptance accepted =
        t_accept(t, 0, 1, request.mean_threshold, request.mean_threshold, 1);
    decision.mean_test_accepted = accepted.accepted;

    if (!accepted.accepted) {
        decision.result = TramoRegressionTestResult::changed;
        decision.mean = false;       // newDesc.setMean(false)
        decision.estimation_cleared = true;  // m.set(newDesc, null) /
                                             // m.clearEstimation()
    }
    return decision;
}

}  // namespace kalos
