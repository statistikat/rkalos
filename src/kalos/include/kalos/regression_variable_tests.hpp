#ifndef KALOS_REGRESSION_VARIABLE_TESTS_HPP
#define KALOS_REGRESSION_VARIABLE_TESTS_HPP
// Regression variable tests and mean rules (item 21-G, stage R): the
// bytecode transcription of the pinned jdplus 3.9.0 stage-R decision
// machines - RegressionVariablesTest.process (the decoder-built
// TradingDays/LengthOfPeriod/MovingHolidays passes with the
// DerivedRegressionTest second life and the all-fail-or-none gate, plus the
// meanTest removal), MeanController.test (the checkMu call sites at 1.96 and
// 1.6) and the RegArimaKernel final mean rule (rtval > 2.5). The three
// decoder scope facts are transcription constants with ledger evidence: the
// t-statistic cutoffs 1.96 (initial test, round-0 checkMu), 1.0 (final test
// meanTest) and 1.6 (redoEstimation checkMu), and CVALFINAL 2.5; the pinned
// pipeline always builds the tests through RegressionVariablesTest.Builder
// with tdTest (t, derived), never tdJointTest, so the FRegressionTest branch
// of the constructor is dead for every pinned call site and the unit
// rejects a joint cutoff fail-closed. The concentrated-likelihood
// statistics, the per-variable layout (findPosition, dim) and - where the
// pinned machine (re-)estimates for its decision - the estimation
// statistics are required request data (like the 21-F regression
// coefficients); no estimation runs inside this unit. The outlier-fix
// integration and user-variable elimination are NOT stage-R code in the
// pinned 3.9.0 bytecode (docs/contracts/regression-tests-contract.md).
// There is no Rcpp boundary; the unit is C++-only.
#include <cmath>
#include <cstddef>
#include <vector>

namespace kalos {

// MeanController.CVALFINAL, the strict cutoff of the final mean rule
// (MeanController.class constant pool; javap -constants).
inline constexpr double kMeanFinalCutoff = 2.5;

// The ConcentratedLikelihoodWithMissing statistics the stage-R machines
// read: coefficients(), unscaledCovariance() (row-major nx x nx), ssq() and
// dim(). Every stage-R call site uses the nobs = 0, ml = false forms, so
// only the (0, 0, false) t-statistics and tstats(0, false) vector are
// transcribed here; the ml forms and the covariance scaling belong to the
// reporting paths, not to these decisions.
struct RegressionTestStatistics {
    std::vector<double> coefficients;        // cl.coefficients(), nx >= 0
    std::vector<double> unscaled_covariance; // cl.unscaledCovariance(), row-major
    double ssq;                              // cl.ssq()
    int dim;                                 // cl.dim(), effective observations
    void validate() const;

    // cl.tstats(0, false): ser_j = sqrt((uvarDiag_j * ssq) / dim) then
    // coef_j / ser_j, no zero guards (the bytecode vector form has none).
    std::vector<double> tstats() const;
    // cl.tstat(0, 0, false) as consumed by MeanController: the ser(i, n, ml)
    // scalar form with its pinned guards (zero covariance entry -> NaN,
    // zero coefficient -> 0/0 -> NaN).
    double mean_tstat() const;
};

// DerivedRegressionTest(tsig, ml = true).accept(cl, -1, first, dim) as the
// RegressionVariablesTest derived branch invokes it: tcur = vcur /
// sqrt((blockSum * ssq) / dim) with vcur the negated left-fold of
// coefficients[first, first + dim) and blockSum the row-major left-fold of
// the unscaledCovariance block, accept = |tcur| > tsig.
double derived_test_tstat(const RegressionTestStatistics& statistics,
                          int first, int dim);

// TRegressionTest(tlow, thigh, nsig).accept(cl, -1, first, dim) on the
// precomputed tstats(0, false) vector: high counts |t| >= thigh, low counts
// the remaining |t| >= tlow (the bytecode else-if chain), accept =
// high > 0 || low >= min(nsig, dim). All stage-R tests are the single-
// argument form (tlow = thigh = t, nsig = 1).
bool t_test_accepts(const std::vector<double>& tstats, double tlow,
                    double thigh, int nsig, int first, int dim);

// Variable.Core classification as the three process() passes filter it
// (ModellingUtility.isTradingDays / isLengthOfPeriod / isMovingHoliday);
// user and other variables are never tested by the live machine.
enum class RegressionVariableKind {
    trading_days,
    length_of_period,
    moving_holiday,
    user,
    other,
};

struct RegressionVariableLayout {
    RegressionVariableKind kind;
    int position; // desc.findPosition(core) at process() time
    int dim;      // core.dim() (1 for LoP and moving holidays)
};

// One RegressionVariablesTest.process(m) run: the builder scalars exactly
// as X13SpecDecoder.readRegression builds them (meanTest(tmu),
// tdTest(ttd, derived), movingHolidaysTest(tmh), ftd always 0), the tested
// variable layout and the description state the machine reads and mutates.
struct RegressionVariableTestRequest {
    RegressionTestStatistics statistics; // the estimation cl reads
    std::vector<RegressionVariableLayout> variables; // desc.variables()
    double mean_test;  // builder meanTest cutoff (1.96 / 1.0; 0 = no test)
    double td_test;    // builder tdTest cutoff (1.96; 0 = no test)
    double mh_test;    // builder movingHolidaysTest cutoff (1.96; 0 = none)
    bool derived_test; // builder derived flag (tdTest(t, derived))
    double td_joint;   // builder ftd; must be 0 - dead decoder path
    bool mean;         // desc.isMean()
    bool adjusted;     // desc.isAdjusted() (log && lp != None)
    int preadjustment; // desc.getPreadjustment(): 0 None, 1 LeapYear,
                       // 2 GoodFriday
    void validate() const;
};

struct RegressionVariableTestResult {
    bool changed;                 // ProcessingResult.Changed (clearEstimation)
    bool mean;                    // desc.isMean() after the run
    int preadjustment;            // desc.getPreadjustment() after the run
    std::vector<bool> removed;    // per request variable
    // Audit mirrors of the machine internals (for bit-exact sweeps):
    std::vector<double> tstats;            // cl.tstats(0, false) as used
    std::vector<bool> tested;              // variable entered a test pass
    std::vector<bool> accepted;            // its t-test accepted
    std::vector<bool> derived_tested;      // the derived second life ran
    std::vector<double> derived_tstats;    // tcur (quiet_NaN when untested)
};

RegressionVariableTestResult run_regression_variable_tests(
    const RegressionVariableTestRequest& request);

// One MeanController(cval).test(m) run (the checkMu machine). The
// statistics are those of the estimation the machine uses: the current
// estimation when the mean flag is on, the mean-added estimation through
// RegArimaUtility.processor(true, 1.0E-5) when it is off (required request
// data; the estimation itself is 21-H wiring).
struct MeanControllerRequest {
    RegressionTestStatistics statistics;
    bool mean; // desc.isMean()
    double cval; // MeanController ctor: 1.96 (round 0), 1.6 (redoEstimation)
    void validate() const;
};

struct MeanControllerResult {
    bool changed; // Changed implies clearEstimation and the flag flip
    bool mean;    // desc.isMean() after the run
    double tstat; // the consumed tstat(0, 0, false)
};

MeanControllerResult run_mean_controller(const MeanControllerRequest& request);

// The RegArimaKernel final mean rule (calc, immediately before the final
// estimate): isAutoModelling() && !desc.isMean() && options.checkMu &&
// rtval0 > 2.5 -> setMean(true) + clearEstimation(). rtval0 is the
// ModelController rTval residual-mean z-statistic, a required scalar (its
// computation is 21-H round-structure machinery); the strict > with the
// pinned NaN ordering (dcmpl) is what this unit qualifies.
bool final_mean_rule_adds_mean(bool auto_modelling, bool mean, bool check_mu,
                               double rtval);

} // namespace kalos
#endif
