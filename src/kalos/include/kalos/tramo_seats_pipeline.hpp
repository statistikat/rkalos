#pragma once

#include "kalos/seats_extract_wk.hpp"
#include "kalos/tramo_fit.hpp"
#include "kalos/tramo_tr0.hpp"

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
// finals. Unsupported fits and unqualified approximations raise an error.
tramo_seats_pipeline_result run_tramo_seats_pipeline(
    const TramoFitResult& fit, int backcast_count = 0, int forecast_count = -2);

struct tramo_seats_rsa0_result {
    tramo_tr0_result preprocessing;
    tramo_seats_pipeline_result decomposition;
};

// Public monthly RSA0 template, from original observations to product finals.
// Defaults to the directly observed JVM horizon of two years of forecasts.
tramo_seats_rsa0_result run_tramo_seats_rsa0(const std::vector<double>& series,
    int start_year, int start_month, int backcast_count = 0, int forecast_count = -2);

struct tramo_seats_automatic_result {
    TramoFitResult preprocessing;
    tramo_seats_pipeline_result decomposition;
};

// Monthly RSA1 from observations. SEATS preserves its existing status or
// refusal for unsupported model and deterministic-effect shapes.
tramo_seats_automatic_result run_tramo_seats_rsa1(const std::vector<double>& series,
    int start_year, int start_month, int backcast_count = 0, int forecast_count = -2);

// Monthly raw RSA3 for the identification/restart envelope of estimate_tramo_tr3.
tramo_seats_automatic_result run_tramo_seats_rsa3(const std::vector<double>& series,
    int start_year, int start_month, int backcast_count = 0, int forecast_count = -2);

// Monthly stock RSAfull with automatic TRAMO calendar selection.
tramo_seats_automatic_result run_tramo_seats_rsafull(const std::vector<double>& series,
    int start_year, int start_month, int backcast_count = 0, int forecast_count = -2);

}  // namespace kalos
