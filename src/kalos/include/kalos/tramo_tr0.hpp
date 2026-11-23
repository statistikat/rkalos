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

}  // namespace kalos
