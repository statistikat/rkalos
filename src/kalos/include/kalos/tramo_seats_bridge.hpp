#pragma once

#include "kalos/tramo_fit.hpp"
#include "kalos/tramo_seats_model_validator.hpp"

namespace kalos {

// 24-D: convert the available 23-F fitted model to a fixed-parameter SEATS
// request on the original monthly calendar. Observed data must be complete.
// AO, LS, and TC with the stock monthly decay of 0.7 are qualified on
// both log and level fits.
// A refused fit or an unsupported design throws invalid_argument.
tramo_seats_model_request make_tramo_seats_model(const TramoFitResult& fit);

}  // namespace kalos
