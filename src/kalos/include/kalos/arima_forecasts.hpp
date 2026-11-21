#ifndef KALOS_ARIMA_FORECASTS_HPP
#define KALOS_ARIMA_FORECASTS_HPP

// Exact-diffuse ARIMA point forecasts/backcasts for the automatic pipeline's
// nonzero-horizon route (item 21-H.4). Faithful port of the pinned JVM
// `jdplus.toolkit.base.core.ssf.arima.ExactArimaForecasts` consumer stack
// (CkmsFilter steady-state filtering with the diffuse fast initializer) as
// `X13Kernel.initialStep` drives it: the model is the landing SarimaModel
// point (unit innovation variance, unit-root polynomials kept in the AR) and
// the series is the preadjusted linearized data on the estimation calendar.
//
// Coefficients follow the JVM polynomial convention throughout (the AR and
// MA polynomials are 1 + c1 B + ...; the driver-side export sign convention
// from the pipeline contract applies at the call site). The mean-correction
// flag mirrors `isMeanCorrection` - a mean/TrendConstant variable in the
// landing model; it folds one regular difference out of the filtered model
// (the AR unit roots and the MA polynomial are both multiplied by 1-B) and
// reports the bar-polynomial mean of the final state.
//
// Acceptance: measured-tier against the JVM reference captures (the pinning
// notes kalos-21h4-pinning.md section E); the filter trajectory is algorithm-
// faithful, not bit-locked. Statuses follow the pipeline vocabulary:
// "available", or "failed:arima_forecasts" when the steady-state initializer
// rejects the model (the JVM "Invalid fast filter" SsfException).

#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

struct ArimaForecastModel
{
    int period = 12;
    int p = 0, d = 0, q = 0;
    int bp = 0, bd = 0, bq = 0;
    // JVM-surface coefficients, sizes p, bp, q, bq.
    std::vector<double> phi, bphi, theta, btheta;
    bool mean_correction = false;
};

struct ArimaForecastsResult
{
    std::string status = "available";
    std::string message;
    // Chronological: backcasts ascending to the first observation, then
    // forecasts after the last one. Empty when the status is not "available".
    std::vector<double> backcasts;
    std::vector<double> forecasts;
    double mean = 0.0;
};

// `series` is the estimation-domain response (the JVM r = linv - v27
// alignment); it must be finite and non-empty. Negative or zero horizon
// requests are caller errors: the JVM never constructs the filter for them
// and this function simply returns empty vectors.
ArimaForecastsResult arima_forecasts(const ArimaForecastModel& model,
    const std::vector<double>& series, std::size_t backcast_horizon,
    std::size_t forecast_horizon);

} // namespace kalos

#endif
