#ifndef KALOS_TRAMO_DIFFERENCING_BATTERY_HPP
#define KALOS_TRAMO_DIFFERENCING_BATTERY_HPP

// Transcription of the round-1 differencing BATTERY of the pinned jdplus
// 3.9.0 tramo driver (item 23-C.2 of docs/tramo-seats-roadmap.md): the
// SeasonalityController first pass (computeReferenceModel with the
// SeasonalityTests ladder, the SeasonalFTest gate and the
// ctx.originalSeasonalityTest score), the RegularUnderDifferencingTest /
// SeasonalUnderDifferencingTest(+2) loop and the second
// (compareReferenceModels) pass, plus the per-state ports the battery
// consumes: the SeasonalOverDifferencingTest code, the residual-seasonality
// ladder score, the SpectralPeaks seasonal-peak gate and the
// incomplete-beta F p-value (detail::f_upper_p, validated bit-exactly in
// tests/test_special_functions.cpp against the double-captured
// validation/reference/fcdf_golden.txt).
//
// Qualification boundary (docs/contracts/tramo-differencing-decision-
// contract.md, battery section): estimation, the maximum-likelihood /
// JointTest algebra and the periodogram stay pinned at the JVM boundary -
// every candidate-model panel and every ladder/peak/F-test value is a
// captured getter of the live run (validation/evidence/native/
// tramo-differencing-battery). The unit replays the DECISIONS over those
// pinned inputs: candidate construction and comparator direction are
// asserted against the captured candidate stream and processing results.
// Anything outside the ported domain - outlier-module-on estimation (the
// 23-F multi-round driver), calendar regressors in the
// SeasonalUnderDifferencingTest2 branch, unobservable candidate states -
// fails closed with std::runtime_error; nothing is approximated.

#include <cstdint>
#include <vector>

namespace kalos {

// One SpectralPeaks record (public final fields ar/tu of the pinned class).
// Encoded from the echoed enum names: AR in {A, a, none, undef}, Tukey in
// {T, t, none, undef}.
enum class PeakAr : unsigned char { A, a, none, undef };
enum class PeakTukey : unsigned char { T, t, none, undef };

struct PeakRecord {
    PeakAr ar;
    PeakTukey tukey;
};

// SeasonalityTests stage facts of one state (qs / nonParametricTest /
// differenced length / spectral-peaks records; *_null mirrors a null
// StatisticalTest).
struct BatteryTests {
    bool qs_null = false, np_null = false;
    double qs_v = 0.0, qs_p = 0.0, np_v = 0.0, np_p = 0.0;
    int differenced_length = 0;
    bool tukey_null = false, ar_null = false;
    // getSpectralPeaks(): empty vector pins null; otherwise the records in
    // array order (length 1, 2 or 6 - the dispatch domain of
    // SpectralPeaks.hasSeasonalPeaks).
    std::vector<PeakRecord> peaks;
};

// The SeasonalFTest triple (StatisticalTest value, JVM-computed p-value and
// the two degrees of freedom read from the F description string).
struct BatteryFTest {
    bool present = false;
    double value = 0.0, pvalue = 0.0, df1 = 0.0, df2 = 0.0;
};

// The ten ModelStatistics getters the ModelComparator / controllers read
// (getLjungBox, getLjungBoxPvalue, getSeasonalLjungBox(+Pvalue), getBic,
// getObservationsCount, getOutliersCount, getSkewnessPvalue,
// getSkewnessAbsvalue, getStableMeanPvalue, getStableVariancePvalue).
struct BatteryPanel {
    double ljung_box = 0.0, ljung_box_p = 0.0, season_box = 0.0, season_box_p = 0.0;
    double bic = 0.0, skew_abs = 0.0, skew_p = 0.0, stable_mean_p = 0.0, stable_var_p = 0.0;
    int nobs = 0, nout = 0;
};

// One observable state of the driver segment (a ":e"-style tag of the
// capture): specification, mean, the coefficient front values the gates
// read, the residual vector behind the ru2 needProcessing autocovariance
// count, and the observation panels.
struct BatteryState {
    int p = 0, d = 0, q = 0, bp = 0, bd = 0, bq = 0;
    bool mean = false;
    int nvar = 0;
    std::vector<double> phi;  // phi(1..p), the regular AR coefficients
    double bphi1 = 0.0;       // bphi(1) (meaningful when bp >= 1)
    double btheta1 = 0.0;     // btheta(1) (meaningful when bq >= 1)
    std::vector<double> resid;
    BatteryFTest ftest;
    BatteryPanel mstat;
    BatteryTests ctlst;  // SeasonalityTests.test(linearized, freq, 1, true)
    BatteryTests resl;   // SeasonalityTests.residualSeasonalityTest(e, freq)
    int resl_score = 0;  // the captured ladder score
    int over_code = 0;   // SeasonalOverDifferencingTest.test(m) code
};

// One ObservingEstimator candidate snapshot (specification, mean flag,
// estimation result and the comparator panel of the estimated candidate).
struct BatteryCandidate {
    int p = 0, d = 0, q = 0, bp = 0, bd = 0, bq = 0;
    bool mean = false;
    bool ok = true;
    BatteryPanel mstat;
};

// The full case: TramoContext entry flags and the pinned state stream in
// capture order (:e, :s1, :ru1, :su1, :ru2, :su2) plus the candidate
// snapshots in estimation order (candidates[0] is the initial estimation
// behind ":e", not a battery candidate).
struct BatteryCase {
    int freq = 12;
    bool seasonal = true;
    int ost = 0;
    BatteryState after_entry;      // :e
    BatteryState after_pass1;      // :s1
    BatteryState after_ru1;        // :ru1
    BatteryState after_su1;        // :su1
    BatteryState after_ru2;        // :ru2
    BatteryState after_su2;        // :su2 (the final state)
    std::vector<BatteryCandidate> candidates;
};

struct BatteryDecision {
    enum class Result { Unchanged, Changed, Unprocessed, Failed };
    Result seas1 = Result::Unchanged;
    Result ru1 = Result::Unchanged, su1 = Result::Unchanged;
    Result ru2 = Result::Unchanged, su2 = Result::Unchanged;
    Result seas2 = Result::Unchanged;
    // Description and context after the whole segment (the ":su2" state).
    int spec_p = 0, spec_d = 0, spec_q = 0, spec_bp = 0, spec_bd = 0, spec_bq = 0;
    bool mean = false;
    bool seasonal = false;
    int ost = 0;
    // Candidate ordinals consumed by the pass-1 re-estimation and by the
    // battery steps, in consumption order (the test asserts this stream
    // equals the captured candidate sequence). The per-step ordinals (-1
    // when the step never reached its candidate) locate each step's snapshot.
    std::vector<int> consumed;
    int cand_pass1 = -1, cand_ru1 = -1, cand_su1 = -1, cand_ru2 = -1, cand_su2 = -1;
    // (spec, mean, ctx) snapshot after seas1 and after each battery step, in
    // execution order (5 entries: after the controller pass, ru1, su1, ru2,
    // su2); the replay compares them against the captured step/seas lines.
    struct Snap {
        int p, d, q, bp, bd, bq;
        bool mean;
        bool seasonal;
        int ost;
    };
    std::vector<Snap> snaps;
};

// Replays the SeasonalityController pass, the four battery steps and the
// second pass over the pinned inputs. The pass-1 re-estimation and every
// step candidate consume the next captured candidate snapshot and assert
// (by throwing std::runtime_error) that the constructed specification,
// mean flag and estimation success match the JVM snapshot; an unconsumed
// snapshot at the end is reported by the replay test.
BatteryDecision decide_tramo_differencing_battery(const BatteryCase& c);

// --- the per-state ports (each asserts against captured facts) -------------

// SpectralPeaks.hasSeasonalPeaks over the record array (null pins the empty
// vector): the length dispatch (6 monthly / 2 quarterly / 1 half-yearly)
// with the exact count/ALL/NONE branches; any other length is false, as in
// the JVM default arm.
bool battery_has_seasonal_peaks(const std::vector<PeakRecord>& peaks);

// The SeasonalityTests score built by residualSeasonalityTest (and by the
// 5-argument seasonalityTest wrapper the over-side uses) over stage facts:
// qs and nonParametricTest significance at 0.01 plus the spectral-peaks
// contribution behind the differenced-length gate (>= 80, or >= 60 with a
// non-12 frequency).
int battery_ladder_score(const BatteryTests& t, int freq);

// SeasonalOverDifferencingTest.test over one state: the specification gate
// (bd >= 1, bq >= 1, btheta(1) < -0.7, no bphi guard on a bp = 1 model),
// the 5-argument ladder score on the controller-side stage facts plus the
// F significance, returning the 0/1/2 preference code of the JVM.
int battery_over_differencing_code(const BatteryState& s, int freq);

} // namespace kalos

#endif // KALOS_TRAMO_DIFFERENCING_BATTERY_HPP
