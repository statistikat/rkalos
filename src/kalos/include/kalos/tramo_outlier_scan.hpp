#ifndef KALOS_TRAMO_OUTLIER_SCAN_HPP
#define KALOS_TRAMO_OUTLIER_SCAN_HPP

#include <cstddef>
#include <utility>
#include <vector>

namespace kalos {

// Item 23-B decision unit: the ROUND-0 TRAMO outlier scan of the pinned
// jdplus 3.9.0 stack, observed through the live
// tramo.internal.OutliersDetectionModule entry point. The transcription runs
// the OutliersDetectionModule.make size gate (TramoUtility.autlar / comatip),
// the ModelDescription pre-stage (AverageInterpolator.cleanMissings, the log
// transformation, the differenced RegArmaModel design with the interpolated
// point columns, the constant column, calendar regressors and pre-registered outlier
// columns), the FastOutliersDetector round-0 estimateModel (OlsComputer
// residuals, HannanRissanen start with stabilization and installation, the
// maximum-likelihood warm seam when requested), the FastOutlierDetector
// initmodel (ConcentratedLikelihoodComputer.DEFAULT_COMPUTER regression and
// FastKalmanFilter innovations over the stationary model), the Mad2(50)
// robust scale, the per-factory T and coefficient tables of
// FastOutlierDetector.processOutlier (pi / (pi * filter) weight columns with
// the incremental ssq accumulator), the SingleOutlierDetector.searchMax
// sweep and the |maxT| >= cv gate of FastOutliersDetector.process. The
// bytecode pin lives in docs/contracts/tramo-outlier-scan-contract.md.
//
// Round-0 scope envelope (stop gate, register D3): this unit observes the
// FIRST scan round only - exactly what the committed corpus round-0 facts
// pin. The multi-round driver (further estimateModel/sod rounds, the backward
// removal step, the maxRound/maxOutliers exits and the robust AO+LS
// correction pre-stage) stays open under D3; its corpus traces serve as the
// decision-brief evidence, and two corpus cases (so, spanout) already show
// the final multi-round set diverging from the round-0 selection.
//
// Register D11: monthly data only. Static request violations (frequency
// other than 12, negative orders, an empty series, a differencing degree of
// zero, missing observations without interpolation, a pre-registered SO)
// throw std::invalid_argument before any state; the pinned SO prespec fill
// loop is outside the transcribed surface, so a pre-registered SO is
// rejected rather than approximated. Pinned runtime failures inside the
// round-0 chain (a rejected size gate, "Ols failed", a failed whitening, a
// table search without any admissible cell) map onto Status::Failed, the
// module's ProcessingResult.Failed, never onto an approximation.

// The pinned factory registration order: AO (filter 1/1, no exclusion
// zones), LS (1/(1-B), one-point zones), TC (1/(1-cB), no start zone and a
// one-point end zone) and SO (1/(1-B^12), two-period zones); the scan table
// column of a factory is its ordinal here, matching the corpus outlier type
// codes and the module's Info types [AO, LS, TC, SO].
enum class TramoOutlierKind { Ao = 0, Ls = 1, Tc = 2, So = 3 };

struct TramoOutlierScanRequest {
    // Raw series observations; a NaN observation is a missing value and is
    // interpolated and carried as a free differenced-model column exactly
    // like the pinned pre-stage. Interpolation is required in the pinned
    // corpus recipe; missings without it are rejected.
    std::vector<double> y;

    int frequency;   // 12 is the only accepted value (register D11)
    bool mean;       // ModelDescription.setMean
    bool log;        // the log transformation of the decision pre-stage

    // Fixed SarimaOrders of the specification; all stationary coefficients
    // are free and carry the pinned SarimaSpec builder fill
    // (phi = -0.1, theta = btheta = -0.2) until estimateModel installs the
    // HannanRissanen start. d + bd must be positive.
    int p, d, q, bp, bd, bq;

    // OutliersDetectionModule factory mix, tc decay, estimation precision,
    // the maximum-likelihood flag (comatip's own verdict is added to it
    // inside make) and the caps. Only round 0 is ever observed, so max_round
    // and max_outliers only echo through make.
    bool ao, ls, tc, so;
    double tc_rate;
    double precision;
    bool maximum_likelihood;
    int max_round;
    int max_outliers;

    // The scan critical value (the TramoUtility.calcCv grid of the module
    // recipe feeds it; it enters make and the round-0 gate unchanged).
    double cv;

    // AverageInterpolator.cleanMissings on the raw observations, the pinned
    // corpus recipe; missing observations without it are rejected.
    bool interpolate;

    // Selection span already resolved to estimation-domain coordinates:
    // all -> [0, length), TimeSelector.between(l, h) -> [l, h + 1], an
    // out-of-range selector -> the collapsed bounds that switch the scan
    // off. make takes the domain-relative start and end of the selected
    // sub-domain; the calendar arithmetic belongs to the caller.
    int lower_bound, upper_bound;

    // Pre-registered outlier variables in insertion order (position in
    // estimation-domain coordinates, kind). An in-domain variable enters
    // the differenced design (indicator convolved with the differencing,
    // the TC column by the pinned repeated-multiply decay) and excludes its
    // own scan cell; an out-of-domain variable is tagged excluded by
    // ModelDescription and the module exclusion is a guarded no-op, so it
    // neither enters the design nor excludes anything.
    std::vector<std::pair<int, TramoOutlierKind>> prespec;

    // Undifferenced calendar design in the description's fitted column order.
    // Each column spans the observations. These precede outlier indicators.
    std::vector<std::vector<double>> calendar_regressors;
};

struct TramoOutlierScanDecision {
    // Failed: make rejected the series (comatip < 0) or the round-0 chain
    // failed the way the module maps onto ProcessingResult.Failed.
    // Unchanged: the scan ran but selected nothing (bounds collapsed, the
    // scan skipped, no admissible cell or the gate failed). Selected: the
    // round-0 winner passed the |maxT| >= cv gate.
    enum class Status { Failed, Unchanged, Selected };
    Status status;

    // Round-0 winner (Selected only): estimation-domain row and factory
    // column, the raw T table value and the coefficient table cell.
    int position;
    int type;
    double max_t;
    double coefficient;

    // Echoes of the machine state the corpus events pin.
    int factories;            // registered scan-table columns
    int lower_bound, upper_bound; // effective scan bounds
    double critical_value;    // the request cv, echoed
};

// Runs the round-0 scan. Static request violations throw
// std::invalid_argument; a request outside the transcribed pinned surface
// (a pre-registered SO, missing observations without interpolation) throws
// std::invalid_argument as the documented fail-closed choice; pinned
// numerical dead ends inside a legal request report Status::Failed.
TramoOutlierScanDecision scan_tramo_outliers(const TramoOutlierScanRequest& request);

} // namespace kalos
#endif
