#ifndef KALOS_OUTLIER_SCORING_HPP
#define KALOS_OUTLIER_SCORING_HPP

#include "estimation_control.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace kalos {

// Pinned x13 outlier-factory order (OutliersDetectionModule.make): AO, LS, TC
// and SO. The ordinal is the pinned factory index j used by the scan tables
// and by the detected (position, type) pairs.
enum class OutlierType : int {
    Additive = 0,
    LevelShift = 1,
    TransitoryChange = 2,
    SeasonalOutlier = 3,
};

// ---- Critical values -------------------------------------------------------
// standard_normal_lower_cdf is the pinned jdplus Utility.intProbability
// Cody/branch transcription for the lower tail (THRSH = 0.66291, the 0/1/2
// region series and the two-factor exp split, transcribed operation for
// operation including the signed floor(1.6*x) split in the third region).
// standard_normal_lower_quantile is the pinned getProbabilityInverse path:
// stnVal (AS 241-style seed) plus the Newton refinement at TOL 1e-13 over at
// most 200 iterations, with the pinned guard on p and 1 - p. Both reproduce
// the pinned doubles bit for bit on the qualified range; the quantile throws
// std::invalid_argument where the pinned code throws DStatException (p
// strictly below 1e-15 or 1 - p strictly below it; a NaN p passes the pinned
// guard), and std::overflow_error on the pinned "Too many iterations" escape
// (which is where a NaN p ends up).
double standard_normal_lower_cdf(double x);
double standard_normal_lower_quantile(double p);

// Pinned X13Utility.calcCv(n) = max(2.8, advancedComputer(0.05).apply(n))
// including the 3x3 basis interpolation solved with the pinned robust solver
// (LULinearSystemSolver2 at 1e-13). n == 0 throws std::invalid_argument; the
// pinned chain has no defined value there (the module never asks).
double default_critical_value(std::size_t observations);

// ---- Outlier factories -----------------------------------------------------
// Pinned TransitoryChangeFactory.rate(freq, monthly): the monthly rate itself
// for freq == 12, pow(monthly, 12 / freq) when 12 % freq == 0, and a rejected
// request otherwise (the pinned "Invalid frequnecy" TsException).
double transitory_change_rate(int frequency, double monthly_rate);

// Pinned excluding zones (start, end) of each factory: AO (0, 0), LS (1, 1),
// SO (2P, 2P) and TC (0, 1). Position p is in the initial allowed zone when
// start <= p and p < length - end.
struct OutlierZone {
    std::size_t start;
    std::size_t end;
};
OutlierZone outlier_excluding_zone(OutlierType type, int period);

// Pinned factory fills over a fresh zero column of the response length:
// AO sets column[position] = 1; the zero-ended LS writes -1 on [0, position];
// the zero-ended SO writes the periodic backward pattern (-1/(P-1) repeated
// P-1 times then +1, repeating) ending at position; TC writes the ramp
// coefficient^(j - position) from position upward, stopping the write loop at
// the first magnitude below the pinned ZERO cutoff 1e-15 (the cutoff value
// itself is still written). The trailing zero entries stay zero. period must
// exceed 1 for the SO type; position must index the column.
void fill_outlier_column(OutlierType type, std::size_t position, int period,
                         double transitory_rate, std::vector<double>& column);

// ---- t-stats (ConcentratedLikelihood.tstats(0, false) recipe) -------------
// One t-statistic per free design column, in free-column order:
// coefficient[i] / sqrt(unscaled_covariance[i,i] * ssq / dim) in that
// arithmetic order, with dim the concentrated dimension and ssq the
// concentrated sum of squares (the rc=false vector path: no per-coefficient
// special cases, NaN and zero propagate through the division). Fixed columns
// contribute no entry (the pinned nx covers the estimated coefficients).
// likelihood must be an available fit (its status "available"); other
// statuses throw std::invalid_argument.
std::vector<double> regression_t_statistics(const RegressionLikelihood& likelihood,
                                            std::size_t concentrated_dimension);

// ---- Scan engine (ExactSingleOutlierDetector.calc) --------------------------
// The scan request carries the current state exactly once: the response, the
// current design (base columns first, appended outlier columns last), the
// fitted stationary ARMA blocks in the shared Kalos convention (unit
// innovation variance, the pinned SarimaModel constant) and the exclusion
// table. The differenced model, Ansley whitening, pivoted-QR initialization,
// Mad(50) robust scale, MLJ residuals and the T/coefficient tables follow the
// pinned bytecode; see docs/contracts/outlier-scoring-contract.md.
struct OutlierScanRequest {
    // Complete response over the estimation domain; masked values throw (the
    // qualified scan reproduces the no-missing differenced model only).
    MaskedSeries series;
    // All coefficients must be free and all values finite; the differenced
    // design is built from these columns in order.
    RegressionDesign design;
    // Stationary blocks [ar, ma, sar, sma]; the differencing filter is built
    // from period, d and seasonal_d.
    std::array<std::vector<double>, 4> arma;
    int period;
    int d;
    int seasonal_d;
    // Pinned setBounds(l, r) scan window over the undifferenced positions.
    std::size_t lower_bound;
    std::size_t upper_bound;
    // Factory order for the T-table columns; duplicates are rejected.
    std::vector<OutlierType> types;
    double transitory_rate;
    // Pinned exclusion table (rows x types, column-major, 1 = excluded).
    // Excluded cells keep NaN in T and are never recomputed.
    std::vector<std::uint8_t> excluded;
};

struct OutlierScanResult {
    // false mirrors a pinned calc() return of false; the tables keep the
    // cleared state (all zero) then, like the pinned T after clear(false).
    bool available = false;
    // "ok" | "empty_bounds" | "too_short" | "rank_deficient" |
    // "whitening_failed". Invalid inputs throw instead.
    std::string status;
    std::size_t rows = 0; // Undifferenced response length (T-table rows).
    std::vector<double> t_table;            // rows x types, column-major.
    std::vector<double> coefficient_table;  // rows x types, column-major.
    // Pinned searchMax winner (posMax, oMax); max_type == max_type_sentinel
    // when no cell beat zero (the pinned oMax == -1 escape).
    static constexpr std::size_t max_type_sentinel = static_cast<std::size_t>(-1);
    std::size_t max_position = 0;
    std::size_t max_type = 0;
    double max_t_stat = 0.0; // Signed T(max_position, max_type); 0 with no winner.
    // Updated exclusion table (rows x types, column-major): the request's
    // exclusions plus every cell the pinned processOutlier excluded through
    // the non-positive residual-variance rule (the pinned exclude(p, j) is
    // permanent for later scans until an explicit allow). The driver feeds
    // this back into the next scan.
    std::vector<std::uint8_t> excluded_after;
};

OutlierScanResult scan_outlier_candidates(const OutlierScanRequest& request,
                                          const Interrupt& interrupt = {});

// ---- Driver (ExactOutliersDetector) -----------------------------------------
// One estimateModel() outcome: the retained candidate of one pinned
// processor.process()/processor.optimize() call and its coefficient
// t-statistics. available == false mirrors a null pinned estimation and
// triggers the pinned rollback paths: pinned estimateModel returns false
// ONLY for a null RegArimaEstimation, so a stalled-but-finite fit is
// available and its tstats drive the decisions (convergence itself is not
// the gate). The pinned "Gls failed" module failure has no available=false
// analogue and stays an exception from the callback.
struct OutlierEstimation {
    bool available;
    std::array<std::vector<double>, 4> coefficients;
    std::vector<double> t_statistics; // One entry per design column, order.
    std::string status;               // Diagnostic; never read by the driver.
};

// The driver passes the request with the current design (base columns first,
// appended outlier columns last) and the current ARMA order template; the
// callbacks own the estimation policy translation (21-A control layer). The
// optimize callback additionally receives the current fitted coefficients as
// the pinned warm start.
using OutlierEstimateFunction = std::function<OutlierEstimation(
    const SarimaEstimationRequest&, const Interrupt&)>;
using OutlierOptimizeFunction = std::function<OutlierEstimation(
    const SarimaEstimationRequest&, const std::array<std::vector<double>, 4>& start,
    const Interrupt&)>;

struct OutlierScoringRequest {
    // Fixed-order base request without outlier columns. All design
    // coefficients must be free (the pinned t-statistics cover every design
    // column); masked responses are rejected.
    SarimaEstimationRequest base;
    // Fitted coefficients of the incoming model (the pinned regarima state):
    // the first scan whitens with them and the first optimize warm start uses
    // them when the opening process-chain estimation fails.
    std::array<std::vector<double>, 4> initial_coefficients;
    std::vector<OutlierType> types; // Non-empty; duplicates rejected.
    // Pinned det.criticalValue, passed through unchanged (the module never
    // recomputes it; X13Utility.calcCv is exposed as default_critical_value).
    double critical_value;
    double transitory_rate;
    // Pinned setBounds(l, r) scan window over estimation-domain positions.
    std::size_t lower_bound;
    std::size_t upper_bound;
    // Pinned maxRound / maxOutliers (x13 module decoder defaults 50 / 30).
    std::size_t max_round;
    std::size_t max_outliers;
};

struct OutlierDecision {
    std::size_t position; // Estimation-domain index of the kept column.
    OutlierType type;
    double t_statistic;   // Signed selection T.
    double coefficient;   // LS coefficient at selection.
};

struct OutlierRejection {
    // Where the removal happened: a selected add whose re-estimation failed
    // (the pinned removeOutlier + discarded warm optimize + stop) or a
    // verifymodel drop. The pinned allow() (verifymodel) re-enables the
    // (position, type) cell for later scans; the add rollback leaves the
    // pinned exclude() in place.
    enum class Origin : int { AddRollback = 0, VerifyModel = 1 };
    std::size_t position;
    OutlierType type;
    double abs_t_statistic; // |t| observed by verifymodel at removal.
    std::size_t round;      // Selection round the removal followed.
    Origin origin;
};

struct OutlierScoringResult {
    // Pinned changed flag: set by every addNewOutlier/removeOutlier and
    // never reset; pinned process() initializes it to false.
    bool changed = false;
    // Exit of the selection loop: "scan_failed" | "below_cv" | "max_rounds" |
    // "max_outliers" | "estimate_failed". The verifymodel drain afterwards can
    // only remove outliers (rejections below).
    std::string stop_reason;
    std::size_t rounds = 0; // Selection rounds completed (pinned round).
    std::vector<OutlierDecision> outliers;          // Kept, in column order.
    std::vector<OutlierRejection> rejections;       // verifymodel trace.
    std::array<std::vector<double>, 4> coefficients; // Last available state.
    std::vector<double> t_statistics;                // Same state, column order.
    // One entry per completed (available) scan, in driver order. The argmax
    // triple reproduces the pinned searchMax outcome; beyond_cv is the pinned
    // |maxT| > cv gate. The first and last entries qualify against the JVM
    // odmscan measurements (entry-state replay and stop scan).
    struct ScanRecord {
        std::size_t position = 0;    // searchMax row
        std::size_t type = 0;        // factory index (or max_type_sentinel)
        double t_stat = 0.0;         // Signed argmax T
        bool beyond_cv = false;      // Pinned |maxT| > cv decision
    };
    std::vector<ScanRecord> scans;
    // Signed max T of the last completed scan; NaN when no scan completed
    // (additive 21-H.3 field: the oracle stop-scan comparison quantity).
    double last_scan_max_t = std::numeric_limits<double>::quiet_NaN();
};

// Stage-O forward selection transcribed from the pinned ExactOutliersDetector
// over the qualified 21-A control layer: the opening estimateModel(true), the
// scan / |maxT| > cv (the pinned dcmpl/iflt stops on a NaN) / add / re-estimate
// loop with the pinned remove-and-warm-reoptimize rollback, the round and
// outlier limits, and the closing verifymodel drain (min-|t| drop while the
// pinned dcmpg scan keeps the earlier entry on ties and NaN, and a NaN minimum
// is dropped) with its warm re-estimations. Callbacks must keep the estimation contract above;
// any exception they throw propagates, mirroring the pinned module-failure
// paths. This unit decides which outliers would be selected and records the
// decision trace; turning them into regression specifications belongs to the
// stage-Z integration (docs/item21-roadmap.md).
OutlierScoringResult detect_outliers(const OutlierScoringRequest& request,
                                     const OutlierEstimateFunction& estimate,
                                     const OutlierOptimizeFunction& optimize,
                                     const Interrupt& interrupt = {});

} // namespace kalos
#endif
