#ifndef KALOS_TRAMO_DIFFERENCING_DECISION_HPP
#define KALOS_TRAMO_DIFFERENCING_DECISION_HPP

// Transcription of the live round-1 differencing decision of the pinned
// jdplus 3.9.0 tramo driver: tramo.internal.DifferencingModule.process over
// a fully-estimated RegSarimaModelling (item 23-C of
// docs/tramo-seats-roadmap.md). The unit owns the whole round-1 ladder -
// the comespd size gate with its airline fallback, the step0 / nextstep /
// initstep / maincondition / cond1 / finalcond loop with its per-step
// HannanRissanen start and warm maximum-likelihood re-estimation, the
// searchur root scans feeding rmax/rsmax, the FastKalmanFilter tmean and the
// isMeanCorrection length table - and decides the Changed/Unchanged result
// of the module stage (specification and mean updates of the description).
//
// The request carries the module inputs as they stand after the driver's
// entry estimation: `series` is the interpolated data minus the outlier
// regression effect (RegArimaUtility.interpolatedData and regressionEffect
// in the pinned process(RegSarimaModelling)); on the complete,
// variable-free corpus these reduce to the observed series itself, which the
// golden facts confirm case by case through the echoed ladder panel.
//
// The qualification boundary (docs/contracts/
// tramo-differencing-decision-contract.md): the warm maximum-likelihood seam
// of initstep covers the design-free/constant-column machine qualified by
// item 23-B; if its landing reaches a pinned unit-root-cancellation restart
// (tryUrpCancelling / tryUrmCancelling gates of the private
// RegSarimaComputer.optimize tail) the unit fails closed instead of
// approximating - the corpus never reaches that state. Battery-stage replay
// (SeasonalityController passes, RegularUnder/SeasonalUnder(+2), the
// SeasonalOverDifferencingTest comparison and the residual-ladder spectral
// peaks) is item 23-C.2; its landing facts are pinned inert in the golden.

#include <string>
#include <vector>

namespace kalos {

struct TramoDifferencingRequest {
    // Interpolated, outlier-effect-free response (module field y before the
    // log channel). Complete and finite; the corpus enters it as the observed
    // series together with the description's log flag: with log_input the
    // unit applies the pinned natural-log channel (fdlibm log, the same
    // application as transform_series) so the ladder runs on
    // regarima().y exactly as the module does.
    std::vector<double> series;
    bool log_input = false;
    int period = 12;

    // Entry specification and mean flag of the description (the ladder
    // starts from 0,0 when initial_round, else from d0,bd0). The ARMA
    // orders p0,q0,bp0,bq0 only feed the airline isAirline guard.
    int d0 = 0, bd0 = 0;
    int p0 = 0, q0 = 0, bp0 = 0, bq0 = 0;
    bool mean0 = true;
    bool initial_round = true;

    // DifferencingModule builder pin (TR defaults). maxbd = seasonal ? 1 : 0
    // with maxd = 2 (MAXD/MAXBD constants); the thresholds are raw
    // AutoModelSpec values, not reciprocals.
    bool seasonal = true;
    double ub1 = 0.97, ub2 = 0.91, cancel = 0.05, eps = 1e-7;
    bool ml = false;

    // regarima().getObservationsCount() behind the comespd gate; -1 uses the
    // series length (equal on the corpus).
    int nobs = -1;
};

struct TramoDifferencingDecision {
    enum class Result {
        Changed,      // specification or mean updated (or airline applied)
        Unchanged,    // ladder kept the entry specification and mean
        Unprocessed,  // airline guard: description already airline
    };
    Result result = Result::Unchanged;

    // Which fallback fired, if any: default_model is the !comespd remark
    // "default model selected (not enough obs.)", selection_failed the
    // RuntimeException catch remark "differencing selection failed".
    bool default_model = false;
    bool selection_failed = false;
    bool airline_applied = false;

    // Description after the module (what desc.setSpecification/setMean left).
    int spec_p = 0, spec_d = 0, spec_q = 0;
    int spec_bp = 0, spec_bd = 0, spec_bq = 0;
    bool mean = true;

    // IDifferencingModule$Info quadruple.
    bool mean_correction = false;
    double tmean = 0.0;

    // Private field panel (the probe's pinning surface).
    double rmax = 0.0, rsmax = 0.0, c = 0.0, din = 0.0;
    int iter = 0;
    bool ml = false, useml = false, mlused = false;
    bool has_last_model = false;
    // Last model parameters in JVM block order [phi | bphi | theta | btheta].
    std::vector<double> last_model;
    // Module spec (mspec) orders of the ladder's SarimaOrders field.
    int mspec_p = 0, mspec_d = 0, mspec_q = 0;
    int mspec_bp = 0, mspec_bd = 0, mspec_bq = 0;

    // ProcessingLog sequence of the module section: "push differencing
    // selection", optional "remark <text>", optional "info diffInfo", "pop".
    std::vector<std::string> events;
};

// Runs the round-1 differencing ladder. Throws std::invalid_argument on
// degenerate requests (empty or non-finite series, period < 2). Internal
// pinned failures (non-converging ESPDIF estimation, missing identification)
// follow the JVM exception path into the airline fallback and never throw;
// the unit-root-cancellation seam throws std::runtime_error.
TramoDifferencingDecision decide_tramo_differencing(const TramoDifferencingRequest& request);

} // namespace kalos

#endif // KALOS_TRAMO_DIFFERENCING_DECISION_HPP
