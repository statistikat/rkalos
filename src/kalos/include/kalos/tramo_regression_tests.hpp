#ifndef KALOS_TRAMO_REGRESSION_TESTS_HPP
#define KALOS_TRAMO_REGRESSION_TESTS_HPP

// Transcription of the NO-CALENDAR branch of the tramo regression-test layer
// of the pinned jdplus 3.9.0 stack (item 23-E of
// docs/tramo-seats-roadmap.md): the rows of ModelStatistics.of that the refStats
// gate consumes, the Ljung-Box/se/BIC machines behind them, the mean-rule
// layer (AmiOptions.checkMu, the round-0 acceptAirline rule, the mean-test
// thresholds the two regression-test modules read) and the non-automatic
// members of FastRegressionTest / DefaultRegressionTest that a request can
// reach without a trading-day, length-of-period or Easter column.
//
// Boundary (docs/contracts/tramo-regression-tests-contract.md): every request
// that would drive a calendar-variable kernel - the Fast machine's
// auto-identified td/lp/easter tests and the Default machine's calendar
// builder fields and calendar accept paths - belongs to item 21-I.3
// (maintainer ruling of 2026-11-22 in docs/calendar-roadmap.md, I-4: no
// second trading-day implementation here). Those requests are refused fail
// closed with the pinned reasons below; nothing is approximated. Calendar-
// free requests run the transcription in full.
//
// The likelihood hand-overs are REQUIRED REQUEST DATA (21-G precedent): the
// decisions consume the coefficients, the unscaledCovariance diagonal, ssq,
// dim and the free-parameter count of an already-completed estimation, never
// an estimation this unit runs. The pass2 gate takes both ModelStatistics
// rows as they stand (the current-side row is what ModelStatistics.of
// produces over the hand-over, so callers build it with
// tramo_model_statistics; the reference-side row arrives from the previous
// pass or from a stored reference model).
//
// Numeric notes (pinned in tmp_opencode/kalos-m5-pinning.md SECTION 25):
// - TRegressionTest.accept(cl, nobs, first, dim) calls
//   cl.tstats(nobs >= 0 ? nobs : 0, nobs >= 0), so every TRAMO site
//   (nobs = k = freeParametersCount >= 0) reads the ADDMU arm of
//   ConcentratedLikelihood.ser/tstats (denominator dim - nx - k), NOT the
//   plain arm the 21-G stage-R sites read (nobs = -1, denominator dim). The
//   unit carries its own operation-for-operation pair for both arms; the
//   21-G RegressionTestStatistics::tstats() (plain arm) is NOT reused.
// - The bucket comparisons of accept() are dcmpl-based: dcmpl(a, thigh) >= 0
//   counts high, dcmpl(a, tlow) >= 0 counts low, and a NaN |t| compares
//   below both thresholds under dcmpl and so counts in NEITHER bucket (the
//   same C++ `a >= t` semantics the landed 21-G t_test_accepts uses; the
//   live 21-G call sites never carry NaN t-values).
// - Math.log routes through the fdlibm transcription (capture JVM runs with
//   the dexp/dlog intrinsics disabled); Math.sqrt is the exact IEEE-754
//   sqrt, matching std::sqrt bit-for-bit; p-values reuse detail::chi2_upper_p
//   (no new distribution kernel).

#include <stdexcept>
#include <vector>

namespace kalos {

// ---------------------------------------------------------------------------
// AmiOptions layer (TramoKernel.readAmiOptions / iterate / reduceVa /
// lastSolution)
// ---------------------------------------------------------------------------

// AmiOptions.from(spec).checkMu = usingAutoModel || meanIsTest: the mean test
// is wired whenever automatic modelling is on or the user asked for the mean
// to be tested. Live AMI witnesses: TR0-TR2 give false, TR3 and TRfull true.
bool tramo_ami_check_mu(bool using_auto_model, bool mean_is_test);

// The round-0 acceptance rule of TramoKernel.iterate:
// acceptAirline && dcmpg(1 - refStats.getLjungBoxPvalue(), ljungBoxLimit) < 0.
// `ljung_box_pvalue` is the stored reference row's p-value, `limit` is
// options.ljungBoxLimit (the AutoModelSpec.pcr start value); the dcmpg
// ordering means a NaN p-value does NOT accept.
bool tramo_round0_accepts(double ljung_box_pvalue, double limit,
                          bool accept_airline);

// TramoKernel.reduceVa(): if (dcmpg(curva, 2.0) == 0) return false; else
// curva = Math.max(2.0, curva * (1 - options.reduceVa)) and return true.
// Math.max propagates NaN (unlike std::max), so a NaN curva stays NaN.
struct TramoReduceVa {
    bool changed = false;
    double curva = 0.0; // value after the update (unchanged when !changed)
};

TramoReduceVa tramo_reduce_va(double curva, double reduce_va);

// TramoKernel.lastSolution() orders rewrite: p := 3; if (bd > 0) bp := 0;
// q := 1; if (context.seasonal) bq := 1; d and bd pass through. The caller
// owns the companion effects (removeVariable(isOutlier), round := 1).
struct TramoOrders {
    int p = 0, d = 0, q = 0, bp = 0, bd = 0, bq = 0;
};

TramoOrders tramo_last_solution_orders(const TramoOrders& orders,
                                       bool seasonal);

// TramoKernel.estimateModel() retry ladder as a pure decision over request
// data: pass and pass3 are kernel fields, need_outliers_desc /
// need_auto_desc are the two desc predicates evaluated at the failure site,
// estimation_results / test_results are the finite sequences of
// FinalEstimator.estimate and testRegression(m, 1.96) outcomes the loop
// consumes (true = success). The loop: an estimation failure returns
// !completed (after the pass==1 reduceVa arm); with pass3 || pass != 0 the
// first success completes; otherwise the Fast mean test may fail at most
// five times (the for-counter increments on every pass through the loop
// head, so five failures then completion), each failure re-entering the
// estimation. Request data exhaustion throws std::invalid_argument.
struct TramoEstimateModelRetryRequest {
    int pass = 0;
    bool pass3 = false;
    bool need_outliers_desc = false;
    bool need_auto_desc = false;
    double curva = 2.0;
    double reduce_va = 0.12;
    std::vector<bool> estimation_results;
    std::vector<bool> test_results;
};

struct TramoEstimateModelRetry {
    bool completed = false; // return value of estimateModel
    std::size_t estimations_consumed = 0;
    std::size_t tests_consumed = 0;
    TramoReduceVa reduce; // the failure-arm reduceVa (identity when unused)
};

TramoEstimateModelRetry tramo_estimate_model_retry(
    const TramoEstimateModelRetryRequest& request);

// ---------------------------------------------------------------------------
// ModelStatistics.of
// ---------------------------------------------------------------------------

// The ConcentratedLikelihood hand-over ModelStatistics.of reads, plus the
// description fields it reads (annual frequency, the arimaSpec free-parameter
// count, the estimation-domain length and the ami-outlier tally).
// cov_diagonal is the unscaledCovariance DIAGONAL only: the consumed rows
// need nothing else (bic uses nx = coefficients.size(), se uses
// ssq/(dof - k)). The Ljung-Box n is residuals.size(), not dim.
struct TramoModelStatisticsRequest {
    std::vector<double> residuals;     // cl.e(), length n >= 1, finite
    std::vector<double> coefficients;  // cl.coefficients() (nx)
    std::vector<double> cov_diagonal;  // cl.unscaledCovariance().diagonal()
    double ssq = 0.0;                  // cl.ssq()
    int dim = 0;                       // cl.dim()
    int degrees_of_freedom = 0;        // cl.degreesOfFreedom()
    double sigma2 = 0.0;               // cl.sigma2()
    int annual_frequency = 12;
    int free_parameters = 0;           // desc.getArimaSpec().freeParametersCount()
    int estimation_domain = 0;         // desc.getEstimationDomain().getLength()
    // ModelStatistics.outliersCount: the NUMBER OF ami-flagged outlier
    // VARIABLES (variables().filter(isOutlier(v,true)).count()), distinct
    // from the pass2 nb = countRegressors(isOutlier(v,true)) that sums
    // core.dim() (the gate request carries nb/nr separately).
    int ami_outliers = 0;
};

// The rows the tramo decision layer consumes, plus the record-only block it
// never reads (the skewness/stableMean/stableVariance surfaces stay at the
// builder defaults here; the JVM kernels are not landed and the tramo gate
// never reads them - decision-inert, documented, not transcribed).
// `built == false` is the JVM throw surface (the main Ljung-Box build feeds
// Chi2(df = m - k <= 0), which throws out of of() with no catch): no row is
// invented, every double is NaN and the caller must not consume the row.
// The seasonal pair is computed whenever annual_frequency > 1 (catch surface
// unreachable for finite residuals; an out-of-range seasonal lag contributes
// a zero autocorrelation, witnessed live by slb = 0, slbp = 1 on n = 7).
struct TramoModelStatistics {
    bool built = false;
    double ljung_box = 0.0;
    double ljung_box_pvalue = 0.0;
    double se = 0.0;
    double bic = 0.0;
    double seasonal_ljung_box = 0.0;
    double seasonal_ljung_box_pvalue = 0.0;
    int lb_length = 0;          // LjungBox.defaultAutoCorrelationsCount(freq, n)
    int outliers_count = 0;
    int observations_count = 0; // desc.getEstimationDomain().getLength()
    int effective_observations_count = 0; // cl.dim()
};

TramoModelStatistics tramo_model_statistics(
    const TramoModelStatisticsRequest& request);

// ---------------------------------------------------------------------------
// The refStats gate (TramoKernel.pass2)
// ---------------------------------------------------------------------------

// Gate inputs as the bytecode reads them: `ljung_box_complement` and
// `reference_ljung_box_complement` are the 1 - pvalue forms the two rows are
// converted into at the top of pass2 (the current side from
// tramo_model_statistics, the reference side from the stored refStats),
// `reference_outliers` is the nr countRegressors tally, -1 = refAuto null.
// `phi1`/`bphi1` are arima().phi(1) / arima().bphi(1); the shape guards read
// them only when the matching order is 1.
//
// The revert rules translate the bytecode chain exactly (rules A..D then the
// two shape guards, evaluated in order, first hit reverts):
//   A: lbp < .95 && rlbp < .75 && rse < se                      (dcmpg)
//   B: pass == 1 && lbp >= .95 && rlbp < .95                    (NaN lbp out, dcmpl)
//   C: A's guards && rlbp < lbp && rse < se                     (dead subset of A)
//   D: lbp >= .95 && rlbp < .95 && rse < se                     (NaN lbp out, dcmpl)
//   shapeA: d==0 bd==1 p==1 phi(1) <= -.82 q <= 1 bp==0 bq==1   (dcmpg)
//   shapeB: d==1 bd==0 p==0 q==1 bp==1 bphi(1) <= -.65 bq <= 1  (dcmpg)
// nr > nb skips ALL rules and lands on keep.
struct TramoRefStatsGateRequest {
    int pass = 0;                    // kernel pass before the call
    int round = 0;
    double pcr = 0.975;              // evolved AutoModelSpec.pcr
    bool changed = false;            // iterate's changed flag
    bool seasonal = true;            // TramoContext.seasonal
    double curva = 2.0;              // outliers-detection critical value
    double reduce_va = 0.12;         // AmiOptions.reduceVa
    bool outliers_detection = true;  // isOutliersDetection()

    // Current model state (the pass2 `m` operand).
    TramoOrders orders;
    bool mean = false;
    double phi1 = 0.0;
    double bphi1 = 0.0;
    double ljung_box_complement = 0.0;
    double standard_error = 0.0;
    int current_outliers = 0;        // nb = countRegressors(isOutlier(v,true))

    // Reference state: -1 in reference_outliers means refAuto == null.
    TramoOrders reference_orders;    // refAuto.desc.specification() on revert
    bool reference_mean = false;     // refAuto.desc.isMean() on revert
    int reference_outliers = -1;     // nr
    double reference_ljung_box_complement = 0.0;
    double reference_standard_error = 0.0;
};

struct TramoRefStatsGateResult {
    bool accepted = false;      // return value true: the ladder stops here
    bool reverted = false;      // restore(m) ran, lbp := reference lbp
    bool remark = false;        // "reduction of the critical value ..." remark
    bool last_solution = false; // the pass > 2 last-chance rewrite ran
    TramoReduceVa reduce;       // the pass-1 continue arm (identity when not run)
    int pass = 0;               // pass after the call
    int round = 0;              // round after the call (lastSolution pins 1)
    double pcr = 0.0;           // pcr after the pass-dependent increment
    double curva = 0.0;         // curva after the pass-1 reduceVa arm
    bool need_auto_modelling = false;
    bool need_outliers = false;
    // Current-model state after the call (restore / last-chance rewrite
    // applied; the pgm witness block).
    TramoOrders orders_after;
    bool mean_after = false;
    int outliers_after = 0;
    double effective_ljung_box_complement = 0.0; // lbp at the return test
};

TramoRefStatsGateResult tramo_ref_stats_gate(
    const TramoRefStatsGateRequest& request);

// ---------------------------------------------------------------------------
// The regression-test modules
// ---------------------------------------------------------------------------

// Which pinned machine a request replays.
enum class TramoRegressionTestKind {
    fast,  // tramo.FastRegressionTest.test
    def,   // tramo.DefaultRegressionTest.test
};

enum class TramoRegressionTestResult {
    changed,
    unchanged,
    unprocessed,  // Default only: no calendar field and no mean test
};

// One likelihood hand-over with the arimaSpec free-parameter count the accept
// sites pass as nobs (the addmu arm of cl.tstats, denominator
// dim - nx - nobs). coefficients are in JVM column order: the mean
// coefficient, if any, is coefficient 0 (meanTest.accept(cl, k, 0, 1)).
struct TramoTestStatistics {
    std::vector<double> coefficients;
    std::vector<double> cov_diagonal;
    double ssq = 0.0;
    int dim = 0;
    int free_parameters = 0;
};

// Bytecode-exact ConcentratedLikelihood.tstats(nobs, addmu) over the
// hand-over with nobs = free_parameters: addmu selects the denominator
// dim - nx - free_parameters, plain uses dim.
// t_j = coef_j / sqrt((vdiag_j * ssq) / (double) n). The JVM gate on
// isScalingFactor() is assumed satisfied by the hand-over (every captured
// estimation satisfies it).
std::vector<double> tramo_test_tstats(const TramoTestStatistics& statistics,
                                      bool addmu);

// Bytecode-exact TRegressionTest(tlow, thigh, nsig).accept over
// statistics via the ADDMU arm: high counts dcmpl(|t[first+j]|, thigh) >= 0,
// low counts the remaining dcmpl(|t[first+j]|, tlow) >= 0 (a NaN |t| counts
// in NEITHER bucket), accept = high > 0 || low >= min(nsig, dim).
struct TramoTAcceptance {
    bool accepted = false;
    int high = 0;
    int low = 0;
};

TramoTAcceptance tramo_t_accepts(const TramoTestStatistics& statistics,
                                 int first, int dim, double tlow,
                                 double thigh, int nsig);

struct TramoRegressionTestRequest {
    TramoRegressionTestKind kind = TramoRegressionTestKind::fast;

    // Description state the machine reads and mutates. For Default this is
    // the TEST-MODEL flag pair (copyOf(desc) without calendar fields, which
    // is what the no-calendar request builds).
    bool mean = true;            // regarima().isMean()
    int period = 12;

    // The likelihood the accept sites consume: the live estimation's cl for
    // Fast, the test-model estimation's cl for Default (the machine runs
    // RegArimaUtility.processor(true, precision) over its own copy of the
    // description; that estimation is request data, 21-G precedent).
    TramoTestStatistics statistics;

    // Mean-test layer: testMean off leaves the mean branch absent (the
    // machine has no mean test at all). The Default kernel never moves the
    // mean threshold off its 1.96 builder default; the corpus pins both.
    bool test_mean = true;
    double mean_threshold = 1.96;

    // Default-machine builder fields. A non-null calendar field is a 21-I.3
    // refusal; adjust with no calendar field cannot run its arm and
    // joinTest with no calendar field never builds its F kernel (both inert,
    // corpus-witnessed, accepted here).
    bool has_trading_days = false;
    bool has_length_of_period = false;
    bool has_easter = false;
    bool adjust = false;
    bool join_test = false;

    // Variables the Fast machine filters: a calendar-kind variable that
    // carries the ami attribute is the 21-I.3 seam and is refused; the
    // plain (non-ami) ones are inert in the pinned filters and ride along.
    // auto_moving_holiday has no filter in the pinned machine at all (no
    // corpus witness) and is refused as outside the qualified envelope.
    bool plain_trading_days = false;
    bool plain_length_of_period = false;
    bool plain_easter = false;
    bool auto_trading_days = false;
    bool auto_length_of_period = false;
    bool auto_easter = false;
    bool auto_moving_holiday = false;
};

struct TramoRegressionTestDecision {
    TramoRegressionTestResult result = TramoRegressionTestResult::unchanged;
    bool mean = true;                // desc.isMean() after the run
    bool estimation_cleared = false; // changed: m.set(newDesc, null) /
                                     // clearEstimation()
    // The tested t-statistic and the acceptance of the mean test, echoed for
    // the corpus (default values when the machine never reads the mean).
    bool mean_test_ran = false;
    bool mean_test_accepted = false;
    double mean_tstat = 0.0;
    // Inert echo: plain calendar variables seen on the description.
    bool inert_trading_days = false;
    bool inert_length_of_period = false;
    bool inert_easter = false;
};

TramoRegressionTestDecision run_tramo_regression_tests(
    const TramoRegressionTestRequest& request);

}  // namespace kalos

#endif
