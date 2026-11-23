#pragma once

#include "kalos/seats_extract_wk.hpp"
#include "kalos/tramo_fit.hpp"

namespace kalos {

// 24-D component merge at the fixed fitted-model boundary. Expects all five
// value components and complete component rows for each supplied forecast or
// backcast interval. Conditional WK standard deviations remain in the WK result.
std::vector<seats_component_series> merge_tramo_seats_components(
    const TramoFitResult& fit, const std::vector<seats_component_series>& wk_final);

struct tramo_seats_pipeline_result {
    seats_extract_wk_result extraction;
    std::vector<seats_component_series> finals;
};

// Runs the bridge, the WK kernel, and the deterministic-effects merge for
// an available monthly TRAMO fit. Counts follow the SEATS convention (a
// negative count is in years). A WK error has its existing status and no
// finals. Unsupported fits and D6 approximation-required models refuse.
tramo_seats_pipeline_result run_tramo_seats_pipeline(
    const TramoFitResult& fit, int backcast_count = 0, int forecast_count = -2);

}  // namespace kalos
