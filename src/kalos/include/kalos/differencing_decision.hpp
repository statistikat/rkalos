#ifndef KALOS_DIFFERENCING_DECISION_HPP
#define KALOS_DIFFERENCING_DECISION_HPP

// Transcription of the pinned jdplus x13 DifferencingModule (stage F, item
// 21-D): the d/D selection machine step0/initstep/maincondition/cond1/
// allcond/finalcond plus the mean-correction test calcMean. The unit maps a
// (linearized) response and frequency onto differencing orders and the mean
// flag, exactly as IDifferencingModule$Info records them. Every threshold
// and constant is pinned by bytecode (docs/contracts/differencing-decision-
// contract.md); promotion into configuration is held by the maintainer
// decision list of docs/item21-roadmap.md, so the core applies no defaults.
#include "regression_likelihood.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

// One differencing decision run. All fields are required request data; the
// pinned x13 decoder call site supplies: max_d = 2, max_seasonal_d = 1,
// ub1 = 1/1.0309278350515465, ub2 = 1/1.1363636363636365 (the module stores
// the reciprocals of the x13 spec boundaries), cancel = 0.1, and
// precision = 1e-7. The pinned module ml flag is never set by the decoder
// and is not exposed: the ML branch condition transcribes as
// (failure || useml).
struct DifferencingDecisionRequest {
    // The response entering DifferencingModule.process (in the pinned
    // pipeline the linearized data of the current estimation). Complete and
    // finite; missing data belongs to the earlier preparation stage.
    std::vector<double> series;
    int period;          // qualified period: 4 or 12 (estimation-seam domain)
    int max_d;           // finalcond clamp for d
    int max_seasonal_d;  // finalcond clamp for bd and the hasSeas gate
    double ub1;          // searchur root-modulus boundary
    double ub2;          // maincondition unit-root boundary companion
    double cancel;       // initial c_ (cancellation threshold, less 0.002 per round)
    double precision;    // module eps: trial-fit and mean-test estimation precision
    void validate() const;
};

// One iteration record of the machine for root-cause inspection. The JVM
// observer can only see the final (d, bd, mean, tmean) quadruple, so these
// records qualify structurally (unit-level tier), not by corpus values.
// Coefficients are reported in the pinned JVM convention (1 + phi(B)).
struct DifferencingIterationRecord {
    // 0 = step0 (initstep(true) + root scan), 1.. = nextstep rounds.
    int iteration;
    int p, q, bp, bq;                          // orders set by initstep
    std::string hr_source;                     // "fitted" or "defaultfill"
    bool hr_stabilized;                        // unstable-gate stabilization rewrote a block
    bool failure;                              // pinned failure flag after the unstable gate
    bool trial_fit;                            // ML branch executed an optimize()
    std::string trial_status;                  // warm-reentry status, empty without trial fit
    double phi, theta, bphi, btheta;           // lastModel coefficients consumed
    double rmax, rsmax;                        // searchur maxima (step0) or allcond fills
    double cancel_threshold;                   // c_ as seen by maincondition (step0: c_ before)
    int maincondition;                         // maincondition() return (nextstep rounds; 0 at step0)
    int carried;                               // cond1 return (round 0), else maincondition result
    bool useml;                                // useml leaving the iteration
    bool mlused;                               // mlused leaving the iteration
    int d, seasonal_d;                         // spec d/bd leaving the iteration
};

// The decision; mirrors IDifferencingModule$Info (d, bd, mean, tmean) plus
// the transcribed internals needed to band-check the estimation couplings.
struct DifferencingDecision {
    int d;
    int seasonal_d;
    bool mean_correction;
    double t_mean;
    bool ml_used;      // mlused at the end of the run (not JVM-observable)
    int nextstep_calls;
    std::vector<DifferencingIterationRecord> trace;
};

// The pinned calc() sequence: step0, the nextstep loop (at most six calls;
// the loop breaks when finalcond settles or the iteration guard trips), and
// calcMean. Static input violations throw std::invalid_argument. The pinned
// "Non convergence in IDDIF" condition of the trial fit throws
// std::runtime_error; pipeline numerical failures propagate (the pinned
// raw pipeline has no catch-all at this call site). The mean test follows
// the pinned branches: pure algebra when the spec is stationary, one
// Multiple-start estimation otherwise; a rank-deficient mean regressor
// after whitening reproduces the pinned out-of-bounds access as
// std::runtime_error. Interrupt exceptions always propagate.
DifferencingDecision decide_differencing(const DifferencingDecisionRequest& request,
    const Interrupt& interrupt = {});

} // namespace kalos
#endif
