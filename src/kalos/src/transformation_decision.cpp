#include "kalos/transformation_decision.hpp"

#include "kalos/detail/fdlibm.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace kalos {
namespace {

double unavailable() { return std::numeric_limits<double>::quiet_NaN(); }

// Pinned post-clear() branch: the module fields keep their defaults and no
// estimate exists. The aicc field keeps the pinned 0.0, not NaN.
TransformationBranchReport unavailable_branch() {
    TransformationBranchReport branch;
    branch.available = false;
    branch.log_likelihood = unavailable();
    branch.jacobian = 0.0;
    branch.adjusted_log_likelihood = unavailable();
    branch.observations = 0;
    branch.parameters = 0;
    branch.aicc = 0.0;
    return branch;
}

// Pinned LikelihoodStatistics.calc AICC assignment, in the exact bytecode
// arithmetic order: -2 * (adjusted - n*k/((n - k) - 1)) with the int counts
// widened to double. Divisions by zero keep IEEE results; nothing guards.
double pinned_aicc(double adjusted_log_likelihood, std::size_t observations,
    std::size_t parameters) {
    const double n = static_cast<double>(observations);
    const double k = static_cast<double>(parameters);
    return -2.0 * (adjusted_log_likelihood - (n * k) / ((n - k) - 1.0));
}

// Pinned LikelihoodStatistics build(): a NaN adjustment is sanitized to
// exact 0.0 first, and adjustedLogLikelihood copies the raw likelihood
// whenever the sanitized adjustment is zero (Inf is not sanitized).
double pinned_adjusted(double log_likelihood, double& jacobian) {
    if (std::isnan(jacobian)) jacobian = 0.0;
    return jacobian == 0.0 ? log_likelihood : log_likelihood + jacobian;
}

// Runs one branch. The request passed in has already cleared the decision
// unit's own preflight (pinned control, fully free blocks, complete series,
// enough room to difference), so an escaping exception here is the pinned
// estimate() throw path and makes the branch unavailable.
TransformationBranchReport run_branch(SarimaEstimationRequest request, bool log_branch,
    std::size_t warmup, std::size_t free_count, const Interrupt& interrupt) {
    auto branch = unavailable_branch();
    if (log_branch) {
        // Pinned LogTransformation: log on finite values, others untouched.
        for (auto& value : request.series.values) {
            if (std::isfinite(value)) value = detail::fdlibm_log(value);
        }
        // Pinned LogJacobian on the estimation window with no missing
        // positions: sum the transformed values left-to-right over
        // [warmup, L), skipping non-finite entries, then negate.
        double sum = 0.0;
        for (std::size_t i = warmup; i < request.series.values.size(); ++i) {
            const double value = request.series.values[i];
            if (std::isfinite(value)) sum += value;
        }
        branch.jacobian = -sum;
    }
    try {
        branch.estimation = estimate_sarima_fixed(request, interrupt);
    } catch (const std::exception&) {
        return unavailable_branch(); // Pinned null branch.
    } catch (...) {
        return unavailable_branch();
    }
    branch.available = true;
    branch.log_likelihood = branch.estimation.log_likelihood;
    branch.adjusted_log_likelihood = pinned_adjusted(branch.log_likelihood, branch.jacobian);
    branch.observations = request.series.values.size() - warmup;
    branch.parameters = free_count + request.design.columns + 1;
    branch.aicc = pinned_aicc(branch.adjusted_log_likelihood, branch.observations,
        branch.parameters);
    return branch;
}

} // namespace

bool log_level_choose_log(double aicc_level, double aicc_log, double correction) {
    // Pinned dcmpl chain: a NaN on either side is unordered and selects
    // levels; equality (the strict > fails) selects levels.
    return std::isgreater(aicc_level, aicc_log + correction);
}

TransformationDecisionResult decide_transformation(
    const TransformationDecisionRequest& request, const Interrupt& interrupt) {
    TransformationDecisionResult decision;
    const auto& series = request.branch.series;

    // Pinned step 1 (after clear(), which leaves the 0.0 field defaults):
    // any value with dcmpg(v, 0) <= 0 - a strictly non-positive double; NaN
    // compares +1 under dcmpg and passes the scan - short-circuits to levels
    // before a computer is even built, so this runs ahead of all request
    // validation and both branches stay unrun.
    for (const double value : series.values) {
        if (value <= 0.0) {
            decision.outcome = "nonpositive";
            decision.log_chosen = false;
            decision.aicc_level = 0.0;
            decision.aicc_log = 0.0;
            decision.level = unavailable_branch();
            decision.log = unavailable_branch();
            return decision;
        }
    }
    // A positive non-finite value passes the pinned scan and then enters the
    // pinned interpolation/missing-data staging inside buildTransformation,
    // which the qualified replay does not reproduce. Reject it explicitly,
    // after the scan (a series carrying both a zero and a NaN still takes
    // the pinned short-circuit) and ahead of the control checks.
    for (const double value : series.values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "the log/level decision is qualified for finite responses; a positive "
                "NaN or infinity takes the pinned missing-data staging");
        }
    }

    // Pinned call-site scope: the module builds exactly
    // DemetraUtility.processor(true, precision), and the qualified replay
    // covers the estimate()-to-process() path (fully free blocks) over the
    // complete-window, fully free-design requests the 21-A initializer is
    // qualified for. Each rejection is explicit; nothing is approximated.
    request.branch.control.validate();
    if (request.branch.control.starting_point != EstimationStartingPoint::Multiple) {
        throw std::invalid_argument(
            "the pinned log/level processor uses the Multiple starting point");
    }
    for (const auto& block : request.branch.blocks) {
        if (std::any_of(block.fixed.begin(), block.fixed.end(),
                [](std::uint8_t flag) { return flag != 0; })) {
            throw std::invalid_argument(
                "the log/level decision is qualified for fully free SARIMA blocks; the "
                "pinned estimate() dispatches fixed specifications to the unqualified "
                "optimize path");
        }
    }
    if (std::any_of(series.missing.begin(), series.missing.end(),
            [](std::uint8_t flag) { return flag != 0; })) {
        throw std::invalid_argument(
            "the log/level decision is qualified for complete series");
    }
    if (std::any_of(request.branch.design.fixed.begin(),
            request.branch.design.fixed.end(),
            [](std::uint8_t flag) { return flag != 0; })) {
        throw std::invalid_argument(
            "the log/level decision is qualified for fully free regression designs");
    }
    const std::size_t warmup = static_cast<std::size_t>(request.branch.d) +
        static_cast<std::size_t>(request.branch.seasonal_d) *
            static_cast<std::size_t>(request.branch.period);
    if (series.values.size() <= warmup) {
        throw std::invalid_argument(
            "the response is too short to difference for the log/level decision");
    }

    std::size_t free_count = 0;
    for (const auto& block : request.branch.blocks) {
        for (const auto flag : block.fixed) free_count += flag ? 0 : 1;
    }

    // Pinned steps: estimate the level description, copyOf + setLog-
    // Transformation(true) of the same description, estimate the copy. The
    // level estimate's setFreeParameters write-back is invisible here: the
    // pinned Multiple start overwrites every coefficient in both branches.
    decision.level = run_branch(request.branch, false, warmup, free_count, interrupt);
    decision.log = run_branch(request.branch, true, warmup, free_count, interrupt);

    if (!decision.level.available && !decision.log.available) {
        decision.outcome = "failed";
    } else {
        decision.outcome = "decided";
    }
    decision.aicc_level = decision.level.available ? decision.level.aicc : 0.0;
    decision.aicc_log = decision.log.available ? decision.log.aicc : 0.0;

    // Pinned isChoosingLog, null rules first: a missing log selects levels,
    // a missing level selects logs, otherwise the comparison decides. With
    // both missing the pinned ProcessingResult.Failed keeps levels.
    if (!decision.log.available) {
        decision.log_chosen = false;
    } else if (!decision.level.available) {
        decision.log_chosen = true;
    } else {
        decision.log_chosen = log_level_choose_log(
            decision.aicc_level, decision.aicc_log, request.aicc_log_correction);
    }
    return decision;
}

} // namespace kalos
