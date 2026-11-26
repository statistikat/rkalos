#ifndef KALOS_ESTIMATION_CONTROL_HPP
#define KALOS_ESTIMATION_CONTROL_HPP

#include "arma_objective.hpp"
#include "arma_parameters.hpp"
#include "optimizer.hpp"
#include "regression_likelihood.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

// The starting-point policy of the pinned jdplus RegSarimaComputer. All four
// starts are qualified: Zero and Default are the pinned fills, HannanRissanen
// runs the transcribed HannanRissanenInitializer and Multiple adds the two
// internal-precision passes with the pinned selection. The two initializer
// starts require a complete response and a fully free regression design (see
// estimate_sarima_fixed); selecting them otherwise is a rejected request, not
// a silent fallback. Promotion of any pinned value into configuration is held
// by the maintainer decision list of docs/item21-roadmap.md.
enum class EstimationStartingPoint { Zero, Default, HannanRissanen, Multiple };

// Solver seam (item 21-D.2): the qualified Kalos BFGS over the mapped −ll
// objective (the 21-A/B/C corpus), or the transcribed pinned sum-of-squares
// Levenberg-Marquardt machine of RegArmaSsqFunction (detail/ssq_levmar.hpp)
// over the physical JVM-convention coordinates of the request surface. Fully
// free coefficient blocks run the machine directly; blocks with fixed entries
// run the same machine over the reduced free-coordinate seam (the pinned
// SarimaFixedMapping, detail::SsqFixedSeam) — the eligibility the M9-H
// qualification measured (tmp_opencode/m9-perf/lmext, 2026-11-25).
enum class EstimationSolver { Bfgs, PinnedLevenbergMarquardt };

// Pinned RegSarimaComputer control fields, kept as required request data:
// there is no default-constructed control and the core applies no implicit
// value. Numerically active in this unit: eps, starting_point. The remaining
// fields record the pinned call-site fidelity and are validated below.
struct EstimationControl {
    // Final optimization precision (pinned DEF_EPS 1e-7; the pinned
    // PROCESSOR singleton uses 1e-9). Forwarded to the solver as its
    // function tolerance; the solver keeps its qualified gradient criterion.
    double eps;
    // Pinned internal precision (DEF_INTERNAL_EPS 1e-4) of the Multiple-start
    // passes. Positive and finite; read only by the Multiple starting policy.
    double feps;
    // Pinned ml: objective in maximum-likelihood form. Only true is qualified
    // by the existing objective; false is rejected rather than approximated.
    bool maximum_likelihood;
    // Pinned mt (builder useParallelProcessing): a stored, never-read field
    // of the pinned processor. Accepted at either value; no Kalos effect.
    bool unit_root_search;
    // Pinned cdf (builder useCorrectedDegreesOfFreedom): reduces the
    // score/information sample by the mapping dimension. This unit exposes no
    // score or information, so either value is accepted with no effect.
    bool corrected_degrees_of_freedom;
    // Pinned fast (builder computeExactFinalDerivatives sets its negation):
    // derivative source of the pinned processor. The qualified Kalos solver
    // owns its derivative policy, so either value is accepted.
    bool fast_derivatives;
    EstimationStartingPoint starting_point;
    EstimationSolver solver;
    void validate() const;
};

// One fixed-order Sarima estimation request for the control layer. The
// blocks use the shared Kalos coefficient convention (AR 1 - phi(B), MA
// 1 + theta(B); [regular AR, regular MA, seasonal AR, seasonal MA] block
// order) and the mean/drift enter through explicit regression columns exactly
// as in the qualified regression-likelihood contract; no implicit mean exists
// here.
struct SarimaEstimationRequest {
    MaskedSeries series;
    RegressionDesign design;
    // Entries flagged fixed keep their initial value exactly. Free entries
    // must be admissible as map templates; their values never act as a
    // starting point because the pinned starting-point policy overwrites
    // every coefficient before optimization.
    std::array<ArmaParameterBlock, 4> blocks;
    int period;
    int d;
    int seasonal_d;
    // When the response has already been unit-root-filtered by the caller
    // (d = seasonal_d = 0), the spec SarimaOrders the pinned machine still
    // reads its unit-root-cancellation gates from (orders.getDifferenceOrder
    // () and orders.getBd()) are carried here instead; negative inherits
    // from d/seasonal_d.
    int spec_d = -1;
    int spec_bd = -1;
    EstimationControl control;
    // Qualified 20b solver controls; function_tolerance is replaced by eps.
    OptimizerControl optimizer;
    // Item 23-F F.2: ask the pinned sum-of-squares path to attach the landing
    // point derivatives (SsqLmResult tail) to the result, the machine behind
    // the pinned gradientAtMinimum()/curvatureAtMinimum() used by the
    // FinalEstimator round driver. Ignored by every other entry point.
    bool request_landing_point = false;
    // GlsSarimaComputer selects the unit-circle SarimaMapping for its refit.
    bool gls_sarima_mapping = false;
};

// One decision record of the control layer: the starting-point provenance,
// each pinned estimate() chain (internal-precision Multiple-start chains when
// the Multiple policy runs, then the final eps chain) and the final chain's
// two unit-root-cancellation restart attempts.
struct EstimationPassRecord {
    // estimate | ur_regular | ur_seasonal, prefixed by "multiple_first." or
    // "multiple_second." for the internal-precision chains, plus "initialize"
    // for the HannanRissanen/Multiple starting point.
    std::string stage;
    bool attempted;      // The pinned gates fired and an optimization ran.
    bool gates_passed;   // Both declared polynomials fell within UR_LIMIT 0.1.
    bool accepted;       // Retained (estimate/initialize) or strictly improved (ur_*).
    std::string optimizer_status; // Empty when no optimization ran.
    std::size_t evaluations; // Objective evaluations of this pass.
    double log_likelihood; // Adjusted concentrated ML; NaN when unavailable.
    // Provenance of the initialize record: the HannanRissanen start source
    // ("fitted", "zerovariance" or "defaulted") plus ",stabilized" when the
    // pinned stabilization rewrote a block. Empty for every other stage.
    std::string detail;
};

struct SarimaEstimationResult {
    // converged | likelihood_only | the producing pass's solver status.
    std::string status;
    bool converged;
    // Final coefficients, Kalos convention; fixed entries preserved exactly.
    std::array<std::vector<double>, 4> coefficients;
    double log_likelihood; // NaN when the retained candidate is unavailable.
    // Landing seam of the pinned sum-of-squares path (item 23-F F.2), set
    // only when request_landing_point was asked and the pinned solver ran;
    // it doubles as the machine-failure seam for the FinalEstimator round
    // driver. The Point tail itself lands in the final_information producer
    // (the residuals-based deviance, not here).
    bool has_landing_point = false;
    // Chronological decision records: initialize (HannanRissanen starts
    // only), the internal-precision multiple_first./multiple_second. chains
    // (Multiple only), then the final eps chain estimate, ur_regular,
    // ur_seasonal.
    std::vector<EstimationPassRecord> passes;
};

// Fixed-order estimation control over the qualified objective and optimizer:
// the pinned Zero/Default starting fills, the transcribed HannanRissanen
// initializer (item 21-A.2) for the HannanRissanen and Multiple starts, the
// pinned Multiple selection between its two internal-precision chains and the
// eps-precision final pass, and the pinned unit-root-cancellation restart
// passes (regular pass gated on p, q, d and |1 + phi(B)|, |1 + theta(B)| at
// B = 1; seasonal pass gated on p, q, seasonal_d and the same polynomials at
// B = -1; a deflated restart is kept only when its likelihood is strictly
// better, and any candidate failure keeps the incumbent, mirroring the pinned
// exception branches). The HannanRissanen-dependent starts require a complete
// response and a fully free regression design, and fail explicitly otherwise;
// the pinned Point score/information are not produced by this unit. Invalid
// static inputs throw; unqualified controls fail explicitly.
SarimaEstimationResult estimate_sarima_fixed(
    const SarimaEstimationRequest& request, const Interrupt& interrupt = {});

// Warm-start counterpart of estimate_sarima_fixed, transcribed from the
// pinned SarimaEstimator.optimize path used by the outlier driver's per-
// iteration re-estimation: the same static validation and objective seam as
// estimate_sarima_fixed, but NO starting-point policy and NO unit-root-
// cancellation restarts - one eps-precision optimization run from the given
// physical coefficients (Kalos convention; fixed entries are re-imposed from
// the request blocks). The start must cover the request exactly and be
// admissible for its map; otherwise the request is rejected. The record
// trace is the single "estimate" pass. A fully fixed model reports the
// likelihood-only result exactly as estimate_sarima_fixed does.
SarimaEstimationResult optimize_sarima_from(const SarimaEstimationRequest& request,
    const std::array<std::vector<double>, 4>& start, const Interrupt& interrupt = {});

} // namespace kalos
#endif
