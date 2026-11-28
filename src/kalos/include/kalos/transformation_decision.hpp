#ifndef KALOS_TRANSFORMATION_DECISION_HPP
#define KALOS_TRANSFORMATION_DECISION_HPP

#include "estimation_control.hpp"

#include <cstddef>
#include <string>

namespace kalos {

// Literal transcription of the pinned jdplus LogLevelModule comparison
// (isChoosingLog: dcmpl against aiccLog + correction, NaN on either side
// selects levels). With the pinned call-site correction -2.0 this reads
// "levels are kept only when the log AICC is at least 2.0 worse"; the
// boundary is strict in favour of levels (a difference of exactly 2.0
// selects levels, a difference just below 2.0 selects logs). The correction
// is request data, never a default: no caller of decide_transformation
// receives an implicit threshold.
bool log_level_choose_log(double aicc_level, double aicc_log, double correction);

// One fixed-order estimation request for both branches of the decision.
// branch.series carries the LEVEL data; the log branch is the pinned
// copyOf + setLogTransformation(true) twin built by the unit. The control
// must carry the pinned decision processor
// (DemetraUtility.processor(true, 1e-7)): maximum_likelihood true and
// starting_point Multiple, with the builder-default feps 1e-4 and the
// pinned call-site eps 1e-7 kept as required request data. Any other
// control is a rejected request, not a silent fallback, because the pinned
// log/level module builds no other computer.
struct TransformationDecisionRequest {
    SarimaEstimationRequest branch;
    // Pinned LogLevelModule Builder aiccLogCorrection (TransformSpec
    // DEF_AICDIFF, -2.0). Not validated: the pinned comparison has defined
    // results for every double, including NaN, which selects levels.
    double aicc_log_correction;
    // 21-I stage 1 (m17-td1): the calendar-carrying variant of the same
    // pinned module. The X13 decoder puts the spec trading-days variables in
    // the initial description and hands LogLevelModule$Builder.preadjust
    // (autoAdjust ? lp : None, readTransformation offsets 43-60), so the
    // level branch estimates with [td, lp], the log branch with the "lp"
    // variable removed (LogLevelModule offsets 146-153) and, under a
    // preadjustment, the length-of-period factors folded into the data
    // (LengthOfPeriodTransformation before the log; the Jacobian gains
    // sum log f over the estimation window, the same structure the CEDM
    // arms use). Empty vectors keep the pinned shape byte-identical.
    // Columns carry the RAW domain (the estimation unit's unit-root filter
    // maps them, 23-F F.2 multi-column machine); every column must have
    // branch.series size, factors likewise with strictly positive entries.
    std::vector<std::vector<double>> level_columns;
    std::vector<std::vector<double>> log_columns;
    std::vector<double> log_preadjust_factors;
};

// One branch as the pinned module leaves it after ModelDescription.estimate
// with the Jacobian-corrected statistics (llAdjustment, LikelihoodStatistics
// AICC form). The branch is unavailable when the pinned estimate() would
// have thrown (the module keeps an Estimation there as null): every field
// then keeps the pinned post-clear() value and estimation stays default.
struct TransformationBranchReport {
    bool available;
    SarimaEstimationResult estimation;
    // Concentrated maximum likelihood at the final candidate; NaN when the
    // retained candidate is unavailable (the pinned likelihood_only guard).
    double log_likelihood;
    // Pinned LogJacobian value: exact 0.0 for the level branch, and
    // -sum log(y_t) over the differenced estimation sample [d + D*period, L)
    // for the log branch. NaN is sanitized to 0.0 by the pinned statistics
    // builder before the adjustment.
    double jacobian;
    // log_likelihood after the adjustment (pinned adjustedLogLikelihood:
    // the raw likelihood whenever the sanitized adjustment is zero).
    double adjusted_log_likelihood;
    // Pinned effectiveObservationsCount L - missing - (d + D*period) and
    // estimatedParametersCount free-arma + design-columns + 1.
    std::size_t observations;
    std::size_t parameters;
    // Pinned LikelihoodStatistics AICC: -2 * (adjusted - n*k/(n-k-1)), in
    // that arithmetic order. NaN propagates; it is never guarded.
    double aicc;
};

struct TransformationDecisionResult {
    // "nonpositive" (a strictly non-positive value, the pinned levels
    // short-circuit; both branches stay unrun), "failed" (both branch
    // estimates are unavailable, the pinned ProcessingResult.Failed keeps
    // levels), or "decided".
    std::string outcome;
    // Pinned isChoosingLog: a missing log branch selects levels, a missing
    // level branch selects logs, otherwise log_level_choose_log decides.
    bool log_chosen;
    // Pinned module fields aiccLevel/aiccLog: 0.0 when the branch did not
    // run or was unavailable, the computed AICC otherwise.
    double aicc_level;
    double aicc_log;
    TransformationBranchReport level;
    TransformationBranchReport log;
};

// Stage-T log/level decision transcribed from the pinned sa-core
// LogLevelModule.process over the qualified 21-A control layer: the short-
// circuit scan (a strictly non-positive value; the pinned dcmpg predicate
// lets NaN pass), the finite-response rejection (a positive NaN or infinity
// takes the pinned missing-data staging, which the qualified replay does
// not reproduce), the level estimate, the copyOf + log-estimate twin with
// the LogTransformation data rule (log only on finite values, non-finite
// left untouched), the pinned Jacobian window and statistics adjustment,
// the AICC forms and the comparison above. The branch is otherwise a plain
// fixed-order 21-A request: masked responses, fixed regression coefficients
// and fixed ARMA entries exceed the qualified scope and throw instead of
// approximating the pinned behaviour. A branch exception is recorded as an
// unavailable branch, mirroring the pinned null handling and its
// exception-cleanup branches; it never escapes as a failure of the decision.
// This unit decides; it does not mutate the modelling and feeds no fit
// seam yet (docs/contracts/transformation-decision-contract.md).
TransformationDecisionResult decide_transformation(
    const TransformationDecisionRequest& request, const Interrupt& interrupt = {});

} // namespace kalos
#endif
