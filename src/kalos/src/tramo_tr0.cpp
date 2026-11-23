#include "kalos/tramo_tr0.hpp"

#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/tramo_model_builder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {

tramo_tr0_result estimate_tramo_tr0(const std::vector<double>& series, int start_year, int start_month) {
    const auto& spec = *tramo_template("TR0");
    const int period = 12;
    const int difference_order = spec.arima.d + period * spec.arima.bd;
    // Two MA parameters, one regression coefficient and the variance.
    const int estimated = static_cast<int>(spec.arima.theta.size() + spec.arima.btheta.size()) + 2;
    if (series.size() <= static_cast<std::size_t>(difference_order + estimated) ||
        series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("estimate_tramo_tr0: positive residual degrees of freedom required");
    if (!std::all_of(series.begin(), series.end(), [](double x) { return std::isfinite(x); }))
        throw std::invalid_argument("estimate_tramo_tr0: a complete finite monthly series is required");

    tramo_build_request build;
    build.values = series;
    build.start_year = start_year;
    build.start_month = start_month;
    build.arima = spec.arima;
    build.mean = spec.mean;
    build.transform_function = spec.transform_function;
    auto description = build_tramo_model(build);
    if (!description.mean || description.log_transformation || !description.arima.is_airline())
        throw std::logic_error("estimate_tramo_tr0: inconsistent TR0 model registration");

    SarimaEstimationRequest request{};
    request.series.values = description.interpolated;
    request.series.missing.assign(series.size(), 0);
    request.design.columns = 1;
    request.design.values.assign(series.size(), 1.0);
    request.design.coefficients = {0.0};
    request.design.fixed = {0};
    for (int block : {1, 3}) {
        request.blocks[block].initial = {0.0};
        request.blocks[block].fixed = {0};
    }
    request.period = period;
    request.d = description.arima.d;
    request.seasonal_d = description.arima.bd;
    request.control = {1e-7, 1e-4, true, false, true, true,
        EstimationStartingPoint::HannanRissanen, EstimationSolver::PinnedLevenbergMarquardt};
    const auto fit = estimate_sarima_fixed(request);
    if (!std::isfinite(fit.log_likelihood) || fit.coefficients[1].size() != 1 || fit.coefficients[3].size() != 1)
        throw std::runtime_error("estimate_tramo_tr0: estimation produced no likelihood point");

    const auto dy = detail::pinned_backfilter_apply(
        detail::pinned_differencing_polynomial(period, request.d, request.seasonal_d), series);
    const auto& theta = fit.coefficients[1];
    const auto& btheta = fit.coefficients[3];
    const auto likelihood = detail::pinned_deviances_evaluate_multi(dy, 0, 0, 1, 1, period,
        {std::vector<double>(dy.size(), 1.0)}, {theta.front(), btheta.front()},
        detail::PinnedComputerShape::Full);
    if (!likelihood || likelihood->b.size() != 1 || !std::isfinite(likelihood->log_likelihood) ||
        !std::isfinite(likelihood->ssq_err) || likelihood->ssq_err <= 0)
        throw std::runtime_error("estimate_tramo_tr0: no finite positive-variance fitted likelihood");

    tramo_tr0_result result;
    result.residual_ssq = likelihood->ssq_err;
    result.mean_coefficient = likelihood->b.front();
    result.log_likelihood = likelihood->log_likelihood;
    result.effective_observations = static_cast<int>(dy.size());
    result.estimated_parameters = estimated;
    auto& model = result.model;
    model.series = std::move(description.interpolated);
    model.start_year = start_year;
    model.start_month = start_month;
    model.d = request.d;
    model.bd = request.seasonal_d;
    model.mean_correction = description.mean;
    model.theta = theta;
    model.btheta = btheta;
    model.innovation_variance = result.residual_ssq / (result.effective_observations - estimated);
    if (model.innovation_variance == 0)
        throw std::runtime_error("estimate_tramo_tr0: fitted innovation variance underflow");
    return result;
}

}  // namespace kalos
