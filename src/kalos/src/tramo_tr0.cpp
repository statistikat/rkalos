#include "kalos/tramo_tr0.hpp"

#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/forecast_information.hpp"
#include "kalos/preparation.hpp"
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

    auto dy = detail::pinned_backfilter_apply(
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

    // RegArimaUtility.fullResiduals: remove drift from the stationary
    // response, then run FastKalmanFilter.process on the fitted MA model.
    for (double& v : dy) v -= result.mean_coefficient;
    auto residuals = detail::pinned_fast_kalman_residuals(dy, 0, 0, 1, 1, period,
                                                          {theta.front(), btheta.front()});
    if (!residuals || residuals->size() != dy.size())
        throw std::runtime_error("estimate_tramo_tr0: residual filter failed");
    const auto start = TimeDomain(start_year, start_month, period, series.size()).at(difference_order);
    result.residuals.start_year = start.first;
    result.residuals.start_month = start.second;
    result.residuals.values = std::move(*residuals);
    // TsResiduals excludes the regression and free ARMA parameters from
    // its degrees of freedom. The fitted variance counts separately.
    result.residuals.standard_error = std::sqrt(result.residual_ssq /
        (result.effective_observations - (estimated - 1)));
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

tramo_tr0_forecast_errors_result tramo_tr0_forecast_errors(const tramo_tr0_result& fit, int months) {
    const auto& model = fit.model;
    if (months < 0) throw std::invalid_argument("tramo_tr0_forecast_errors: negative forecast count");
    if (model.frequency != 12 || model.period != 12 || model.d != 1 || model.bd != 1 || model.log ||
        !model.mean_correction || !model.phi.empty() || !model.bphi.empty() ||
        model.theta.size() != 1 || model.btheta.size() != 1)
        throw std::invalid_argument("tramo_tr0_forecast_errors: requires the fitted TR0 airline model");
    if (fit.estimated_parameters != 4 || fit.effective_observations <= fit.estimated_parameters ||
        model.series.size() != static_cast<std::size_t>(fit.effective_observations) + 13 ||
        !std::isfinite(fit.residual_ssq) || fit.residual_ssq <= 0.0)
        throw std::invalid_argument("tramo_tr0_forecast_errors: inconsistent TR0 likelihood");
    if (months == 0) return {"no_requested_predictions", {}};
    if (model.series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() - months))
        throw std::length_error("tramo_tr0_forecast_errors: extended domain exceeds the index range");

    std::array<ArmaParameterBlock, 4> blocks{};
    blocks[1] = {model.theta, {1}};
    blocks[3] = {model.btheta, {1}};
    for (const double value : {model.theta.front(), model.btheta.front()}) {
        const auto check = arma_admissibility({value}, true);
        if (check.status != "admissible") return {check.status, {}};
    }
    const ArmaParameterMap map(blocks);
    const TimeDomain observed(model.start_year, model.start_month, 12, model.series.size());
    auto drift = prepare_regression(align_regressors(observed, {}, 0, months), {}, 1, 1, 0, false, true);
    auto context = std::move(drift.design);
    context.values = std::move(drift.regressors.values);
    const MaskedSeries series{model.series, std::vector<std::uint8_t>(model.series.size(), 0)};
    const auto projected = forecast_information(series, context, 0, months, map, {}, 12, 1, 1);
    if (projected.status != "available") return {projected.status, {}};

    // RegSarimaModel.internalForecasts scales the process innovation by the
    // estimated-parameter count (including variance); the independently
    // fitted mean block uses ConcentratedLikelihood.covariance(k,true).
    const int df = fit.effective_observations - fit.estimated_parameters;
    const double process_scale = fit.residual_ssq / df;
    const double mean_scale = fit.residual_ssq / (df + 1);
    std::vector<double> errors(months);
    for (int i = 0; i < months; ++i) {
        const double process = (*projected.known_regression_covariance)(i, i) / projected.variance;
        const double mean = (*projected.regression_covariance)(i, i) / projected.variance;
        errors[i] = std::sqrt(process * process_scale + mean * mean_scale);
        if (!std::isfinite(errors[i])) return {"covariance_range", {}};
    }
    return {"available", std::move(errors)};
}

}  // namespace kalos
