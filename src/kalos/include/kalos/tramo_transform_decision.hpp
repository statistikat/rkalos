#ifndef KALOS_TRAMO_TRANSFORM_DECISION_HPP
#define KALOS_TRAMO_TRANSFORM_DECISION_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

// Item 23-B decision unit for the TRAMO pre-stage: the TramoUtility.calcCv
// critical-value grid, the TramoSeasonalityDetector machine (the (1-B)
// removeMean differencing fed to the Qs Ljung-Box variant at the seasonal
// lags and the Friedman month-rank test, the p-value verdict ladder and the
// m011 downgrade), and the LogLevelModule log-versus-level decision run on
// the airline (or m011) fixed-order model through the pinned GlsSarimaComputer
// (HR start, sum-of-squares Levenberg-Marquardt, maximum likelihood) with the
// ln(ssq * factor) comparison and the log-preference offset. The bytecode pin
// lives in docs/contracts/tramo-transform-decision-contract.md.
//
// Scope envelope (register D11): monthly data only; the series must be
// complete and finite. The robust AO+LS correction that TramoKernel can run
// ahead of the decision belongs to the multi-round outlier driver and is held
// by open register D3: a request that would need it is rejected, never
// approximated. Violations raise std::invalid_argument before any state.

// TramoUtility.calcCv(n) = max(2.0, CriticalValueComputer.simpleComputer()):
// n <= 50 gives 3.3, 50 < n < 450 gives the exact 3.3 + 0.0025 * (n - 50)
// ramp, n >= 450 gives 4.3 (all branches exceed the 2.0 floor). A negative
// observation count throws; the pinned driver never asks for one.
double tramo_critical_value(std::size_t observations);

// One raw test readout: value, upper-tail p, and the validity flag the pinned
// StatisticalTest carries (TestsUtility.testOf marks a test valid exactly
// when the distribution was constructed with a positive degrees of freedom).
struct TramoSeasonalityTest {
    double value = 0.0;
    double p_value = 0.0;
    bool valid = false;
};

// TramoSeasonalityDetector.hasSeasonality(values, period) with the raw
// getters: the (1 - B) differenced, mean-removed sequence (the single
// differenced block the pinned DifferencingResult.of(seq, period, 1, true)
// yields), the Qs statistic (LjungBox with zero-mean autocorrelations, lag
// equal to the period, two seasonal-lag autocorrelations, positive
// autocorrelations only - the sum early-exits to exactly 0.0 when the first
// seasonal autocorrelation is non-positive), the Friedman statistic on the
// same differenced block (midranks over the trailing whole months), and the
// verdict ladder: Qs p < 0.01 is Strong; both tests fired is Moderate; one
// fired (or one score point) is Weak; otherwise NotObservable. The ladder
// comparisons use the pinned dcmpg thresholds 0.01 and 0.05.
struct TramoSeasonalityReport {
    std::string verdict;   // Strong | Moderate | Weak | NotObservable
    int to_int = 0;        // the pinned Seasonality.toInt()
    TramoSeasonalityTest qs;       // Chi2 with 2.0 degrees of freedom
    TramoSeasonalityTest friedman; // Chi2 with period - 1 degrees of freedom
    std::vector<double> differenced;
};

// The detector machine over one complete finite monthly series.
TramoSeasonalityReport tramo_seasonality(const std::vector<double>& values);

// One TRAMO pre-stage transform decision request. The robust-correction
// input is intentionally absent: with no pre-correction the decision runs on
// the description as built.
struct TramoTransformDecisionRequest {
    std::vector<double> values; // complete, finite, size >= 2 * period
    int period;                 // 12 (D11)
    // ModelDescription flags as TramoKernel sets them: setAirline(true) and
    // setMean(mean). The mean flag never reaches the log/level model - the
    // pinned LogLevelModule always estimates with the meanCorrection constant
    // column (bytecoded-verified) - but it selects the regression design of
    // the later stages and is carried as request data.
    bool mean = true;
    bool detect_seasonality = true; // the auto-detection path
    bool spec_seasonality = true;   // the seasonal flag when detection is off
    // TramoSpec.log().logPreference(); the JVM default is Math.log(0.95).
    double log_preference = 0.0;
    // TramoSpec.log().estimationPrecision(); the processor precision the
    // module builds its GlsSarimaComputer with. Positive and finite.
    double precision = 0.0;
};

struct TramoTransformDecisionResult {
    // Seasonality pre-stage (verdict/qs/np) and the resulting flags.
    TramoSeasonalityReport seasonality;
    std::string spec_pre;       // "(0,1,1)(0,1,1)" | "(0,1,1)(0,0,0)"
    bool seasonal = false;
    bool nonpositive = false;   // the negative-values short-circuit fired

    // ProcessingResult of LogLevelModule.process: Unchanged | Changed |
    // Failed (an estimation exception, mirroring the pinned RuntimeException
    // catch). The nonpositive path is Unchanged with nonpositive set, like
    // the pinned remark-and-return.
    std::string loglevel_status;
    bool choosing_log = false;
    // Module field defaults after clear() are log = 1.0, level = 0.0,
    // slog = 0.0; a failed branch keeps its default exactly like the pinned
    // module fields (the ln(ssq * factor) of the successful branch survives).
    double log_ll = 1.0;
    double level_ll = 0.0;
    double log_correction = 0.0; // mean of the log-transformed response
    double log_preference = 0.0;
    bool log_transform = false;  // desc.setLogTransformation(choosingLog)

    // Chronological ProcessingLog remarks of the pre-stage: "No identifiable
    // seasonality" (downgrade) and "negative values, levels are chosen".
    std::vector<std::string> remarks;
};

TramoTransformDecisionResult decide_tramo_transform(
    const TramoTransformDecisionRequest& request);

} // namespace kalos
#endif
