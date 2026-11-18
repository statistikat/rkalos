#ifndef KALOS_INITIAL_ARMA_IDENTIFICATION_HPP
#define KALOS_INITIAL_ARMA_IDENTIFICATION_HPP

// Initial ARMA identification (item 21-E, stage M part 1): the bytecode
// transcription of the pinned jdplus.x13.base.core.x13.regarima.ArmaModule
// initial-identification machine composed onto the qualified estimation
// seams of items 21-A/21-D. The unit reproduces the OLS-residual cursor of
// RegArimaUtility.olsResiduals, the three-stage candidate enumeration over
// (bp, bq) / (p, q) / (bp, bq) grids with RegArmaBic sorting and merging,
// and the BIC-margin refinement select(cursor, d, bd) with its eleven
// ordered acceptance rules. It is a decision unit feeding (not yet wired
// into) the fit seams; the AutoModellingModule refinement of stage M part 2
// lives in kalos/arma_refinement.hpp and consumes the select_arma_orders
// seam below (the pinned AutoModellingModule calls ArmaModule.select
// directly on the re-differenced linearized data, without an OLS cursor);
// the round-structure wiring belongs to later slices. There is no Rcpp
// boundary; the unit is C++-only and every pinned threshold is a
// transcription constant, not
// configuration (docs/contracts/initial-arma-identification-contract.md).
#include "numerics.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

// The whole pinned module state plus the upstream model the residual cursor
// is computed from. Every field is required request data: the struct has no
// default-constructed form and the core applies no implicit value. The
// decoder-pinned call-site values (models_count 5, max_p/max_q 2,
// max_bp/max_bq 1, balanced false, mixed true, accept_white_noise false,
// eps 1.0E-7) arrive from the call site; see the contract.
struct InitialArmaIdentificationRequest {
    std::vector<double> series;            // response, complete and finite
    std::vector<std::vector<double>> regressors; // pre-differencing user columns
    bool mean_correction;                  // pinned RegArimaModel.isMean()
    int period;                            // qualified domain {4, 12}
    int d;                                 // spec regular differencing order
    int seasonal_d;                        // spec seasonal differencing order
    int models_count;                      // pinned nmod (estimations slots)
    int max_p;
    int max_q;
    int max_bp;
    int max_bq;
    bool balanced;                         // rules 1-7 gate on it
    bool mixed;                            // false drops p>0&&q>0 candidates
    // Pinned acceptWhiteNoise: a stored, never-read field of the pinned
    // 3.9.0 bytecode (constructor putfield, no getfield). Accepted at either
    // value with no effect; the stage-3 (bp) clamp keys on the seasonal
    // differencing order, as pinned.
    bool accept_white_noise;
    double eps;                            // candidate estimation precision
    void validate() const;
};

// One RegArmaBic candidate evaluation (HR initialization + one pinned LM
// chain over the deviances machine + the BIC2 statistics assembly).
struct ArmaCandidateEvaluation {
    int stage;   // 1 seasonal (bp, bq) grid, 2 regular (p, q) grid,
                 // 3 seasonal refinement grid
    int slot;    // 0-based evaluation order inside the stage (grid order)
    int p, q, bp, bq;
    // HannanRissanen start source "fitted" | "zerovariance" | "defaulted",
    // suffixed ",stabilized"; empty for a zero-parameter candidate, where
    // the qualified seam short-circuits before initialization and no start
    // exists.
    std::string initial_source;
    std::string fit_status;      // estimate-pass solver status
    double log_likelihood;       // machine logLikelihood at the landing
    double bic;                  // BIC2 = (-2 ll + k log n) / n
    std::array<std::vector<double>, 4> coefficients; // landing, Kalos blocks
};

// One merged estimations[] entry (the pinned getPreferedModels panel), in
// ascending BIC order.
struct ArmaPreferredModel {
    int p, q, bp, bq;
    double bic;
    double log_likelihood;
    std::array<std::vector<double>, 4> coefficients;
};

// One refinement-scan iteration over the merged list (every index, not only
// the accepted ones).
struct ArmaSelectionStep {
    int candidate; // loop index i (1-based over the merged list)
    int rule;      // accepted rule 1..11, 0 when no rule accepted
    int incumbent; // resulting incumbent index
    int moves;     // resulting accepted-move count
};

struct InitialArmaIdentification {
    int p, q, bp, bq;                      // selected stationary orders
    std::vector<double> residuals;         // the OLS-residual cursor
    std::vector<ArmaCandidateEvaluation> evaluations; // all stages, grid order
    std::vector<ArmaPreferredModel> preferred_models; // merged top list
    std::vector<ArmaSelectionStep> refinement; // empty when no refinement ran
};

// The pinned ArmaModule.select(cursor, period, d, bd) on a caller-owned
// cursor (the module entry point the 21-F AutoModellingModule transcription
// drives directly: the pinned module never builds an OLS cursor there). The
// grid bounds, panel size and refinement flags are the same ArmaModule
// fields as in InitialArmaIdentificationRequest and stay required request
// data; the cursor replaces the whole upstream model (no regressors, no
// mean column: the pinned candidate fits are design-free on the sequence
// handed to select, exactly as in the initial-identification grid).
struct ArmaSelectRequest {
    std::vector<double> cursor;            // response, complete and finite
    int period;                            // qualified domain {4, 12}
    int d;                                 // refinement polarity: spec d
    int seasonal_d;                        // refinement polarity: spec bd
    int models_count;                      // pinned nmod (estimations slots)
    int max_p;
    int max_q;
    int max_bp;
    int max_bq;
    bool balanced;                         // rules 1-7 gate on it
    bool mixed;                            // false drops p>0&&q>0 candidates
    bool accept_white_noise;               // stored, never-read (21-E)
    double eps;                            // candidate estimation precision
    void validate() const;
};

// The select outcome: the selected stationary orders, every candidate
// evaluation of the three stages, the merged estimations[] panel and the
// refinement trace (the same records the full initial identification
// produces; the cursor itself stays caller-owned).
struct ArmaSelectResult {
    int p, q, bp, bq;
    std::vector<ArmaCandidateEvaluation> evaluations;
    std::vector<ArmaPreferredModel> preferred_models;
    std::vector<ArmaSelectionStep> refinement;
};

// Runs the pinned select(cursor, period, d, bd) machine on the given cursor.
// Semantics, tiers and the pinned-abort equivalences are those of the
// select section of the contract: an empty grid, a candidate the pinned
// machine rejects (start or landing) or a panel that cannot fill all
// models_count slots throw std::runtime_error (the module-level callers
// decide between Failed and the airline fallback); invalid static inputs
// throw std::invalid_argument. Interrupt exceptions always propagate.
ArmaSelectResult select_arma_orders(
    const ArmaSelectRequest& request, const Interrupt& interrupt = {});

// Runs the pinned initial identification. The cursor is the OLS residual of
// the differenced upstream model (the filtered response regressed on fresh
// ones for mean_correction plus the BackFilter-filtered user columns, via
// the pinned calcResiduals column-major accumulation); each candidate is
// estimated with the GlsSarimaComputer chain (HannanRissanen initialization,
// one eps-precision pinned-Levenberg-Marquardt chain over the qualified
// deviances machine) on the cursor with no regression design (the pinned
// candidate model carries no mean; the mean correction lives in the cursor),
// and scored by BIC2 = ((-2 logLikelihood + (p+q+bp+bq+1) log n) / n with
// fdlibm log; the selection then follows the pinned sort/merge/refinement
// rules exactly.
//
// Fail-closed scope rejections (std::invalid_argument): non-finite or empty
// series, regressor length/finite violations, negative orders, period
// outside {4, 12} (the qualified estimation-seam domain; the pinned JVM
// accepts any period), models_count < 1, negative grid bounds, non-finite or
// non-positive eps. Pinned-abort equivalences (std::runtime_error, the
// pinned ProcessingResult.Failed paths, never a silent fallback): an OLS
// cursor with no observations after differencing ("Not enough observations"
// - the pinned EcoException), a degenerate OLS cursor response below the
// pinned norm2 epsilon (the pinned null-results NPE), a candidate whose
// start or landing the pinned machine rejects (the pinned pre-minimize or
// final-evaluation throw), and a merged list that cannot fill all
// models_count slots (the pinned refinement NPE on the empty slot).
// Interrupt exceptions always propagate.
InitialArmaIdentification identify_initial_arma(
    const InitialArmaIdentificationRequest& request, const Interrupt& interrupt = {});

} // namespace kalos
#endif
