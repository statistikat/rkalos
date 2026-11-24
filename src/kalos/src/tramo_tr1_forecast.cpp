#include "kalos/tramo_tr0.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/outlier_scoring.hpp"
#include "kalos/tramo_seats_bridge.hpp"
#include "tramo_airline_forecast.hpp"

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace kalos {
namespace {

void require(bool valid, const char* reason) {
    if (!valid) throw std::invalid_argument(std::string("tramo_tr1_forecast_errors: ") + reason);
}

// (1-B)(1-B^12) on the observed and forecast domains.
double difference(const std::vector<double>& series, std::size_t t) {
    return (series[t] - series[t - 1]) - (series[t - 12] - series[t - 13]);
}

double projection(const detail::AirlineForecastBoundary::State& weights,
                  const detail::AirlineForecastBoundary::State& state) {
    return std::inner_product(weights.begin(), weights.end(), state.begin(), 0.0);
}

}  // namespace

tramo_tr0_forecast_errors_result tramo_tr1_forecast_errors(const TramoFitResult& fit, int months) {
    require(months >= 0, "negative forecast count");
    // The fitted bridge validates the complete monthly domain, coefficients,
    // registered effects, and joint positive-variance likelihood.
    const auto model = make_tramo_seats_model(fit);
    require(model.mean_correction && model.phi.empty() && model.bphi.empty() && model.theta.size() == 1 &&
            model.btheta.size() == 1 && model.d == 1 && model.bd == 1, "requires the fitted TR1 airline model");
    const auto n = model.series.size();
    const auto p = fit.regression_coefficients.size();
    require(p == fit.variables.size() + 1 && p < n && static_cast<std::size_t>(fit.ls_npar) == p + 3 &&
            n == static_cast<std::size_t>(fit.ls_eff) + 13, "invalid fitted likelihood or design shape");
    require(fit.covariance.n == p && fit.covariance.covariance.size() / p == p &&
            fit.covariance.covariance.size() % p == 0, "invalid fitted coefficient covariance");
    if (months == 0) return {"no_requested_predictions", {}};
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max() - months))
        throw std::length_error("tramo_tr1_forecast_errors: extended domain exceeds the index range");
    const double theta = model.theta.front(), seasonal = model.btheta.front();
    // The fit has already selected the MA coefficients. Finite-MA conditioning
    // does not invert their polynomial, including at the near-unit fitted AO.
    auto response = fit.builder.interpolated;
    if (fit.log) for (double& x : response) x = detail::fdlibm_log(x);
    std::vector<std::vector<double>> columns;
    for (const auto& effect : fit.variables) {
        auto& column = columns.emplace_back(n + months);
        fill_outlier_column(static_cast<OutlierType>(effect.type), effect.position, 12, effect.tc_rate, column);
    }
    auto errors = detail::airline_forecast_uncertainty(std::move(response), columns, true, fit.log,
        theta, seasonal, fit.regression_coefficients, fit.covariance.covariance, model.innovation_variance, months);
    if (!errors) return {"covariance_range", {}};
    return {"available", std::move(*errors)};
}

std::optional<std::vector<double>> detail::airline_forecast_uncertainty(
    std::vector<double> response, const std::vector<std::vector<double>>& columns,
    bool mean, bool log, double theta, double seasonal, const std::vector<double>& coefficients,
    const std::vector<double>& covariance, double innovation_variance, int months) {
    const auto n = response.size(), p = coefficients.size();
    response.reserve(n + months);
    auto design = [&](std::size_t j, std::size_t t) { return mean && j == 0 ? 1.0 : difference(columns[j - mean], t); };
    detail::AirlineForecastBoundary boundary(theta, seasonal);
    using State = detail::AirlineForecastBoundary::State;
    State response_state{};
    std::vector<State> design_states(p);
    for (std::size_t t = 13; t < n; ++t) {
        const double variance = boundary.advance();
        if (!(variance > 0.0) || !std::isfinite(variance)) return std::nullopt;
        if (log) boundary.observe(response_state, difference(response, t), variance);
        for (std::size_t j = 0; j < p; ++j) boundary.observe(design_states[j], design(j, t), variance);
    }
    // The regression sensitivity is future X minus its integrated conditional
    // prediction. Propagate the ordinary and seasonal integrals for each column.
    std::vector<double> running(p), seasonal_sums(12 * p), sensitivity(p), errors(months);
    double future_variance = 0.0;
    for (int h = 1; h <= months; ++h) {
        State weights{};
        const double process = boundary.forecast_variance(h, weights, future_variance);
        const auto t = n + h - 1;
        for (std::size_t j = 0; j < p; ++j) {
            running[j] += design(j, t);
            auto& integrated = seasonal_sums[((h - 1) % 12) * p + j];
            integrated += running[j];
            sensitivity[j] = integrated - projection(weights, design_states[j]);
        }
        double variance = innovation_variance * process;
        for (std::size_t i = 0; i < p; ++i)
            for (std::size_t j = 0; j < p; ++j)
                variance += sensitivity[i] * covariance[i * p + j] * sensitivity[j];
        double error = std::sqrt(variance);
        if (log) {
            response.push_back(response[t - 1] + response[t - 12] - response[t - 13]);
            const double raw = response.back() + projection(weights, response_state) +
                std::inner_product(sensitivity.begin(), sensitivity.end(), coefficients.begin(), 0.0);
            // JD+ LogNormal.stdev uses the raw forecast and its total variance.
            const double log_variance = error * error;
            error = detail::fdlibm_exp(raw + 0.5 * log_variance) * std::sqrt(detail::fdlibm_exp(log_variance) - 1);
        }
        if (!std::isfinite(error)) return std::nullopt;
        errors[h - 1] = error;
    }
    return errors;
}

}  // namespace kalos
