#pragma once

#include "kalos/tramo_seats_model_validator.hpp"

#include <cstdint>

namespace kalos {

// Public preprocessing residuals (residuals.tsres and residuals.ser).
struct tramo_tr0_residuals {
    std::int64_t start_year = 0;
    int start_month = 0;
    std::vector<double> values;
    double standard_error = 0;
};

// Fitted projection consumed by the monthly RSA0 decomposition. The model
// carries the two estimated MA coefficients, date, observations and variance.
struct tramo_tr0_result {
    tramo_seats_model_request model;
    double mean_coefficient = 0;
    double log_likelihood = 0;
    double residual_ssq = 0;
    int effective_observations = 0;
    int estimated_parameters = 0;
    tramo_tr0_residuals residuals;
};

// The fully specified TR0 airline model estimates both MA coefficients and
// the mean correction from a complete monthly series. It skips the automatic
// transformation, outlier and model-identification stages.
tramo_tr0_result estimate_tramo_tr0(const std::vector<double>& series, int start_year, int start_month);

struct tramo_tr0_forecast_errors_result {
    std::string status;
    std::vector<double> standard_errors;
};

// Original-series forecast standard errors reported by RegSarimaModel (y_ef).
// Months are a nonnegative count on the future domain. A zero horizon returns
// no_requested_predictions. Numerical failures retain the covariance status.
tramo_tr0_forecast_errors_result tramo_tr0_forecast_errors(const tramo_tr0_result& fit, int months);

struct TramoFitResult;
// Public y_ef on a complete monthly stock TR1 fit with an estimated drift,
// including AO, LS and stock TC effects and the lognormal back-transformation.
tramo_tr0_forecast_errors_result tramo_tr1_forecast_errors(const TramoFitResult& fit, int months);

// Public original-series errors for monthly, automatically identified TR3 models.
tramo_tr0_forecast_errors_result tramo_tr3_forecast_errors(const TramoFitResult& fit, int months);

}  // namespace kalos
