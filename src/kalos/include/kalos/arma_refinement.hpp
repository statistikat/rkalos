#ifndef KALOS_ARMA_REFINEMENT_HPP
#define KALOS_ARMA_REFINEMENT_HPP

// ARMA refinement and comparison (item 21-F, stage M part 2): the bytecode
// transcription of the pinned jdplus.x13.base.core.x13.regarima.
// AutoModellingModule.process - the refinement + comparison stage that
// follows the initial identification in the pinned round-0 pipeline. The
// unit linearizes the data of a given estimation (RegArimaUtility.
// linearizedData: user regressors removed with the concentrated-likelihood
// coefficients, the mean kept), re-enters the differencing decision on the
// linearized response (decide_differencing, 21-D), re-differences it with
// the pinned UnitRoots.D(d) * UnitRoots.D(period, bd) filter, refines the
// stationary orders through the ArmaModule.select machine on that sequence
// (select_arma_orders, 21-E) and compares the resulting specification and
// mean flag against the current ones (the setSpecification/setMean
// Changed/Unchanged decision). Any failure inside that guarded region takes
// the pinned airline fallback (SarimaOrders.airline, mean off) instead of
// throwing. It is a decision unit feeding (not yet wired into) the fit
// seams; the round-structure wiring is stage Z (21-H). The roadmap's
// SarimaComponent hill-climb and AICc-comparator vocabulary are dead or
// stage-R code in the pinned 3.9.0 bytecode (see the contract). There is no
// Rcpp boundary; the unit is C++-only and every pinned value is a
// transcription constant or required request data, never configuration
// (docs/contracts/arma-refinement-contract.md).
#include "differencing_decision.hpp"
#include "initial_arma_identification.hpp"
#include "numerics.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace kalos {

// The upstream model state the pinned module consumes: the raw response and
// regressors, the mean flag, the concentrated-likelihood regression
// coefficients of the CURRENT estimation (cl.coefficients(): the mean
// coefficient first when the mean flag is on, then one per user column) and
// the current SarimaOrders. Together these are exactly what
// AutoModellingModule.process reads from RegSarimaModelling (est + desc);
// the ARMA landing itself never enters the machine (linearizedData reads the
// regression coefficients, not the ARMA parameters). The two module states
// (DifferencingModule, ArmaModule) arrive as their scalar fields, split
// below exactly as X13SpecDecoder.readAutoModel builds them. Every field is
// required request data; the core applies no implicit value.
struct ArmaRefinementRequest {
    std::vector<double> series;                // response y, complete and finite
    std::vector<std::vector<double>> regressors; // user columns, full span
    bool mean;                                 // ModelDescription.isMean()
    // ConcentratedLikelihoodWithMissing.coefficients() of the current
    // estimation: (mean ? 1 : 0) + regressors.size() finite values.
    std::vector<double> regression_coefficients;
    int p, q, bp, bq;  // current specification, stationary part
    int d, seasonal_d; // current specification differencing orders
    int period;        // qualified domain {4, 12}

    // DifferencingModule state (the decoder builds exactly these from the
    // AutoModelSpec): clamps, root-scan boundaries, cancellation start and
    // the module estimation precision.
    int max_d;
    int max_seasonal_d;
    double ub1;
    double ub2;
    double cancel;
    double precision;

    // ArmaModule state: panel size, grid bounds, the balanced/mixed policy,
    // the stored-never-read acceptWhiteNoise field and the candidate
    // estimation precision (independent of the differencing module's).
    int models_count;
    int max_p;
    int max_q;
    int max_bp;
    int max_bq;
    bool balanced;
    bool mixed;
    bool accept_white_noise;
    double eps;
    void validate() const;
};

// The full module outcome. On the pinned normal path the airline fallback
// flags are false, differencing holds the IDifferencingModule$Info record
// and the select records mirror ArmaModule.getPreferedModels plus the
// candidate evaluations and refinement trace; on the pinned RuntimeException
// path fallback is set, the specification is airline(period), mean is off,
// and the records stop at the last completed stage.
struct ArmaRefinement {
    int p = 0, q = 0, bp = 0, bq = 0;  // resulting specification, stationary part
    int d = 0, seasonal_d = 0;         // resulting differencing orders
    bool mean = false;                 // resulting mean flag
    bool specification_changed = false; // the pinned !spec0.equals(spec) branch
    bool mean_changed = false;          // the pinned mean != mean0 branch
    bool changed = false;               // clearEstimation + ProcessingResult.Changed
    bool fallback = false;              // pinned "default model selected" remark
    // Which guarded stage raised into the pinned catch ("differencing" or
    // "select"); empty without a fallback.
    std::string fallback_stage;
    // RegArimaUtility.linearizedData output (full span; set whenever the
    // pure arithmetic ran).
    std::vector<double> linearized;
    // The sequence handed to ArmaModule.select: the linearized data filtered
    // by the pinned (1 - B^d)(1 - B^period)^bd BackFilter, or the linearized
    // data itself when both orders are 0. Empty when the region stopped
    // before the filter.
    std::vector<double> re_differenced;
    std::optional<DifferencingDecision> differencing; // 21-D record (d, bd,
                                                      // mean, tmean, trace)
    std::vector<ArmaCandidateEvaluation> evaluations; // select grids, all stages
    std::vector<ArmaPreferredModel> preferred_models; // merged panel (nmod)
    std::vector<ArmaSelectionStep> refinement;        // eleven-rule scan trace
};

// Runs the pinned AutoModellingModule.process on the given model and
// estimation data. The guarded region maps its RuntimeException catch onto
// std::runtime_error from the inner qualified units (differencing decision,
// candidate machine, unfilled panel): such a failure produces the airline
// fallback result, it never throws. Static input violations throw
// std::invalid_argument (request validation composes the 21-D and 21-E
// scalar validators; a linearized response the differencing unit rejects
// statically is a request-data error, not a pinned path). Interrupt
// exceptions always propagate: the first firing of the callback arms the
// guarded region to rethrow instead of converting to the fallback.
ArmaRefinement refine_arma_model(
    const ArmaRefinementRequest& request, const Interrupt& interrupt = {});

} // namespace kalos
#endif
