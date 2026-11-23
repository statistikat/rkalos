#ifndef KALOS_TRAMO_KERNEL_ROUNDS_HPP
#define KALOS_TRAMO_KERNEL_ROUNDS_HPP

// Transcription of the round-driver seam of TramoKernel from the pinned
// jdplus 3.9.0 stack (item 23-F.3 of docs/tramo-seats-roadmap.md): the
// TramoKernel.ami() round loop (the ten iterate() attempts), the
// iterate() body (estimate site, differencing / arma / outliers guards,
// the change accumulator, the estimateModel() retry ladder, the round-0
// reference snapshot with its acceptAirline rule, the refStats gate call,
// the seasonality-controller pass3 reset and the controller sweep), the
// pass2() continue arms (reduceVa, the pcr ladder, lastSolution) and the
// ProcessingLog event stream the loop produces.
//
// This unit owns the kernel DECISIONS and the LOG STRUCTURE only. Every
// heavy stage the loop consults is REQUIRED REQUEST DATA carried as an
// ordered seam fact (the 21-G precedent; the F.2 driver precedent of
// tramo_final_estimator_run): the differencing / arma / outliers-module
// outcomes and the model states they leave behind, the FinalEstimator
// attempts, the FastRegressionTest sites, the ModelStatistics.of rows,
// the seasonality- and model-controller results, and the module-owned log
// events inside each module block. The kernel-level machines the driver
// runs itself are the landed ones: the refStats gate, the round-0 accept,
// reduceVa, the lastSolution rewrite and the estimateModel retry ladder
// through kalos/tramo_regression_tests.hpp, and the initProcessing
// critical value through tramo_critical_value (23-B). Nothing is
// approximated: a fact that arrives where the state machine does not
// expect one is a desync and throws std::invalid_argument, and facts the
// machine consumed are cross-checked against the JVM's own echoed
// decisions (estimate taken, round0 accept, gate verdicts, gate return,
// iterate end) so a divergence between this transcription and the pinned
// build cannot pass silently.
//
// Boundary (docs/contracts/tramo-kernel-rounds-contract.md): the unit is
// monthly-frequency only (D11) - a non-12 annual frequency refuses with
// the pinned reason, witnessed by kr_q_air in the replay corpus. The
// pre-stage (testSeasonality / testTransformation / the regression
// module) and the post-loop ModelEstimator build are outside this unit:
// the pre-stage enters as the prestage facts (including its log block)
// and the final model block stays in the probe corpus identity witness.
// Calendar-variable requests never reach a calendar kernel here because
// the pre-stage machine is not invoked (its result rides as a fact); the
// 21-I.3 routing of the landed 23-E seam is untouched.
//
// Semantics pinned in tmp_opencode/kalos-m5-pinning.md SECTION 25.1
// (iterate / estimateModel / pass2 / lastSolution / reduceVa decode),
// SECTION 26.4 (the ami() outer decode) and SECTION 30.8, adjudicated
// byte-exact by the double-captured JVM corpus in
// validation/evidence/native/tramo-kernel-rounds/ (Scratch23r manual
// driver identical to the jar kernel on all 35 cases):
// - the step text is "round " + pass (NOT round); kr_q_air's pass-3
//   entry with round pinned to 1 by lastSolution witnesses the slot,
// - needDifferencing = needAutoModelling && !(round == 2 and no
//   ami-outlier variable on the description) && (d < 2 || bd >= 1);
//   needAutoModelling(desc) drops the first conjunct's orders test;
//   needOutliers(desc) = isOutliersDetection() && needOutliersField,
// - isAutoModelling() is spec.isUsingAutoModel() (kr_tr1_air, resolved
//   autoModel 0, skips differencing, arma, the gate AND the seasonality
//   controller while the outliers selection keeps running),
// - the change accumulator assigns: differencing result, then the
//   post-arma specification/mean comparison, then |= the outliers
//   result (bytecode offsets 63-79, 129-166, 194-218),
// - the gate is entered with the NEGATED accumulator: pass2(!changed,
//   m); its continue arm then sets needAutoModelling = !operand = the
//   accumulator value read at the call (kr_q_air's post-pass-1
//   snapshot keeps needAuto = 1 with a changed arma module),
// - the outliers block clears the estimation whenever an ami-outlier
//   variable is on the description BEFORE scanning, so each kernel
//   round feeds a clean design,
// - lastSolution rewrites from the CURRENT specification (p := 3, bp := 0
//   when bd > 0, q := 1, bq := 1 when the context is seasonal), clears
//   the outliers, warns "last chance model", pins round to 1 and
//   invalidates the estimation (the kr_q_air pass-3 estimate site
//   re-estimates),
// - reduceVa fires in the estimateModel failure arm (pass == 1 with both
//   desc predicates) and in the pass-1 gate continue arm, whose remark
//   fires unconditionally after it.
//
// Numeric note: every double in a fact is the raw IEEE-754 value the
// pinned build produced (the probe printed %.17g and the replay restores
// the exact bits); the driver only adds 1.0 - x complements, the pcr and
// curva ladders and the calcCv init, all bit-identical to the JVM.

#include <kalos/tramo_regression_tests.hpp>

#include <string>
#include <vector>

namespace kalos {

// ---------------------------------------------------------------------------
// The ProcessingLog surface
// ---------------------------------------------------------------------------

// One ProcessingLog line. The kernel-emitted events (the "tramo" and
// "automatic model identification" pushes/pops, the "round N" steps, the
// gate-continue remark and the last-chance warning) are produced by this
// unit; everything inside a module block is carried on the seam fact that
// the block belongs to (module-owned text, masked payloads included - the
// probe masked locale-dependent info payloads as "-", and the replay
// compares the masked text as data).
enum class TramoKernelLogKind { push, pop, info, remark, warning, step };

struct TramoKernelLogEvent {
    TramoKernelLogKind kind = TramoKernelLogKind::info;
    int level = 0;
    std::string text;
};

// ---------------------------------------------------------------------------
// Seam facts (request data)
// ---------------------------------------------------------------------------

// ProcessingResult as the pinned modules return it (the outliers machine
// also answers Failed on too-short designs, kr_shortfail's pass-1 scan).
enum class TramoKernelStageResult { changed, unchanged, unprocessed, failed };

// The seam kinds, in the probe's capture vocabulary. The driver expects a
// strict ordering derived from its own state; an unexpected kind is a
// desync.
enum class TramoKernelSeamKind {
    estimation,             // m.estimate() at the iterate() head
    differencing,           // differencing module block
    differencing_drop,      // pass-1 difference-order change: outliers cleared
    arma,                   // arma module block
    arma_q1,                // pass >= 3 zero-parameter q := 1 rewrite
    scan,                   // outliers module block (carries the curva)
    estimation_model,       // FinalEstimator.estimate attempt
    estimation_reduce,      // estimateModel failure-arm reduceVa echo
    regression_test,        // FastRegressionTest site (retry ladder / iterate)
    round0,                 // round-0 reference snapshot and accept
    pass2_gate,             // the refStats gate call (carries the fresh row)
    pass2_return,           // the gate's return value and evolved pcr
    seasonality_controller, // scontroller (carries the context mutations)
    model_controller,       // one controller sweep entry
    iterate_end,            // the iterate() return site
};

enum class TramoKernelIterateEnd {
    done,          // fell through the controller sweep, true
    pass3_reset,   // scontroller Changed with pass3 off, false
    estimation_fail, // estimateModel returned false, false
    regression_fail, // the raw-tsig mean test failed, false
};

// The ModelStatistics.of row as the probe dumped it (freq, k, dim, obs,
// lbp, se, lb, bic, oc, eff plus the e-array). The Ljung-Box pair is
// reproducible bit-for-bit from the e-array through the landed 23-E
// primitives (the replay witnesses that); se and bic need the likelihood
// sigma2/ssq/degreesOfFreedom fields the capture does NOT record, so
// the standalone round replay carries both. The integrated fit regenerates
// all row fields from its live final estimator.
struct TramoKernelStatsRow {
    int free_parameters = 0;   // desc.getArimaSpec().freeParametersCount()
    int dim = 0;               // cl.dim()
    int observations = 0;      // desc.getEstimationDomain().getLength()
    double ljung_box_pvalue = 0.0; // the row's lbp slot (the p-value)
    double se = 0.0;
    double ljung_box = 0.0;
    double bic = 0.0;
    int outliers_count = 0;    // ModelStatistics.outliersCount()
    int effective_observations = 0; // cl.dim() echoed as eff
    std::vector<double> residuals; // cl.e()
};

// One ordered seam fact. Only the fields the kind documents are
// meaningful (the replay's generator fills exactly those; the driver
// reads exactly those).
struct TramoKernelSeam {
    TramoKernelSeamKind kind = TramoKernelSeamKind::estimation;
    std::vector<TramoKernelLogEvent> events; // module-owned block, push..pop

    // differencing / arma / scan / regression_test / seasonality_ and
    // model_controller: the ProcessingResult the module answered.
    TramoKernelStageResult result = TramoKernelStageResult::unchanged;
    // The description state the module left behind (post-run
    // desc.specification() / isMean() / outlier names; "empty" prints as
    // "-"). differencing, scan and both controllers carry parameters
    // (desc.arima().parameters()); arma carries neither outliers nor
    // parameters.
    TramoOrders orders;
    bool mean = false;
    std::vector<std::string> outliers;
    std::vector<double> parameters;
    double critical_value = 0.0; // scan: the curva the block ran with

    // estimation: the kernel's needEstimation() answer and the precision
    // m.estimate() consumed (0 taken => no precision recorded).
    bool estimation_taken = false;
    double intermediate_precision = 0.0;

    // estimation_model: FinalEstimator.estimate return, the kernel pass
    // at construction, getChangedParametersCount(), post state.
    int returned = 1;
    int pass_at_call = 0;
    int nnsig = 0;

    // estimation_reduce: the reduceVa curva echo (never fires in the
    // corpus).
    // regression_test: the probe site ("retry<i>" inside the ladder,
    // "iterate" for the raw-tsig test), the threshold, post mean, the
    // cl.tstats(k, true) dump.
    bool retry_site = false;
    int retry_index = 0;
    double tsig = 0.0;
    std::vector<double> tstats;

    // round0: the refstats row; the JVM accept echo.
    // pass2_gate: the fresh row, the JVM-echoed gate operands and verdict
    // (cross-check), and the ARIMA orders the shape guards read
    // (desc.arima().orders(), from the probe's p2dbg line).
    TramoKernelStatsRow stats;
    bool accept = false;      // round0 echo
    int outlier_regressors = 0;    // nb = countRegressors(isOutlier)
    int reference_regressors = -1; // nr (-1: refAuto null)
    bool reverted = false;    // JVM verdict echo
    std::string rule = "-";   // JVM rule echo (A..F, "skip-nrnb", "-")
    double gate_ljung_box_complement = 0.0;       // gate line lbp (post-restore)
    double gate_standard_error = 0.0;             // gate line se
    double gate_reference_ljung_box_complement = 0.0; // gate line rlbp
    double gate_reference_standard_error = 0.0;       // gate line rse
    TramoOrders arima_orders;

    // pass2_return: the gate return value and the evolved pcr.
    double returned_pcr = 0.0;

    // seasonality_controller: the TramoContext mutations the controller
    // left (context.seasonal and originalSeasonalityTest - the kr_gs2 and
    // kr_q_air pass-3 snapshots read both).
    bool seasonal_after = false;
    int original_seasonality_after = 0;

    // iterate_end.
    TramoKernelIterateEnd end = TramoKernelIterateEnd::done;
};

// ---------------------------------------------------------------------------
// Request
// ---------------------------------------------------------------------------

// The resolved TramoSpec facts the kernel reads (the probe's
// "KER <id> resolved" line).
struct TramoKernelRoundsSpec {
    int annual_frequency = 12;
    int observations = 0;        // regarima.getActualObservationsCount()
                                 // (the calcCv init input; below the series
                                 // length when the series carries missings)
    bool using_auto_model = true;
    bool outliers_detection = true;
    bool accept_airline = false;
    double pcr = 0.95;           // options.ljungBoxLimit (AutoModelSpec.pcr)
    double reduce_va = 0.12;     // options.reduceVa (AutoModelSpec.pc)
    double tsig = 1.0;           // AutoModelSpec.getTsig()
    double critical_value = 0.0; // options.va (OutlierSpec.criticalValue);
                                 // 0.0 triggers the calcCv init (23-B)
    double intermediate_precision = 1e-7; // AmiOptions.getIntermediatePrecision
    bool mean_is_test = false;   // RegressionSpec.getMean().isTest()
    int controller_count = 0;    // the ModelController list size
};

// The ami() pre-stage: its log block, the context mutations
// (testSeasonality sets seasonal and originalSeasonalityTest; the
// not-run case keeps the pinned context defaults seasonal = true,
// originalSeasonalityTest = 0 - kr_tr1_air), the regression module result
// and the resulting description state (the probe's prestate lines).
struct TramoKernelRoundsPrestage {
    std::vector<TramoKernelLogEvent> events;
    bool seasonal = true;
    int original_seasonality_test = 0;
    std::string regression_test; // ProcessingResult name, echoed
    TramoOrders orders;
    bool mean = false;
    bool log = false;            // the prestate log flag (echo-only)
    std::vector<std::string> outliers;
};

struct TramoKernelRoundsRequest {
    TramoKernelRoundsSpec spec;
    TramoKernelRoundsPrestage prestage;
    std::vector<TramoKernelSeam> seams;
};

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

// The per-iterate-entry kernel state (the probe's rsnap line), taken at
// the "round N" step event, exactly where the probe reflected it.
struct TramoKernelRoundSnapshot {
    int pass = 0;
    int round = 0;
    bool pass3 = false;
    bool need_outliers = false;
    bool need_auto_modelling = false;
    double curva = 0.0;
    double pcr = 0.0;
    bool ref_auto = false;
    bool ref_stats = false;
    bool seasonal = false;
    int original_seasonality_test = 0;
};

// The decision the unit took at each seam fact (index-aligned with
// request.seams). Only the kind's own fields are meaningful; they are
// RECOMPUTED here (the fact's echoed verdicts are cross-checks, and the
// replay compares these echoes against the capture).
struct TramoKernelSeamDecision {
    TramoKernelSeamKind kind = TramoKernelSeamKind::estimation;

    int flag = 0;        // estimation: taken; round0: accept; pass2_return:
                         // return value; iterate_end: end code
    double value = 0.0;  // pass2_return: the evolved pcr

    // pass2_gate echo (the landed gate's operands and verdict):
    double effective_ljung_box_complement = 0.0; // lbp (post-restore)
    double current_standard_error = 0.0;
    double reference_ljung_box_complement = 0.0;
    double reference_standard_error = 0.0;
    int current_regressors = 0;
    int reference_regressors = -1;
    bool reverted = false;
    bool skipped_rules = false; // nr > nb (the JVM "skip-nrnb" echo)
    bool accepted = false;      // the gate returned true
    bool remark = false;        // the reduction-of-curva remark fired
    bool last_solution = false; // the last-chance rewrite fired
    int pass_after = 0;
    int round_after = 0;
    double pcr_after = 0.0;
    double curva_after = 0.0;
    bool need_auto_modelling_after = false;
    bool need_outliers_after = false;
    TramoOrders orders_after;
    bool mean_after = false;
    int outliers_after = 0;
};

struct TramoKernelRoundsRun {
    // initProcessing (curva after the calcCv trigger, the pcr start).
    double initial_curva = 0.0;
    double initial_pcr = 0.0;
    bool ami_check_mu = false; // tramo_ami_check_mu(using_auto_model,
                               // mean_is_test) - the wiring echo
    int iterations = 0;        // iterate() attempts used (of the ten)
    bool converged = false;    // an iterate() returned true
    std::vector<TramoKernelRoundSnapshot> snapshots;
    std::vector<TramoKernelLogEvent> events;
    std::vector<TramoKernelSeamDecision> decisions; // aligned with seams
    // Final state (for the contract's witness tables).
    TramoKernelRoundSnapshot final_state;
};

// Replay the ami() round loop over the request. Throws
// std::invalid_argument on an out-of-envelope request (non-12 annual
// frequency - D11 - non-finite facts, negative counts) and on any desync
// between the seam facts and the state machine's expectations.
TramoKernelRoundsRun tramo_run_kernel_rounds(
    const TramoKernelRoundsRequest& request);

}  // namespace kalos

#endif
