#pragma once

#include "kalos/tramo_seats_model_validator.hpp"

namespace kalos {

// Fitted projection consumed by the monthly RSA0 decomposition. The model
// carries the two estimated MA coefficients, date, observations and variance.
struct tramo_tr0_result {
    tramo_seats_model_request model;
    double mean_coefficient = 0;
    double log_likelihood = 0;
    double residual_ssq = 0;
    int effective_observations = 0;
    int estimated_parameters = 0;
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

}  // namespace kalos
