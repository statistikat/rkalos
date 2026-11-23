#include "kalos/tramo_tr0.hpp"

#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/arma_parameters.hpp"
#include "kalos/tramo_model_builder.hpp"

#include <algorithm>
#include <array>
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

    TimeDomain(model.start_year, model.start_month, model.frequency, model.series.size());
    if (!std::all_of(model.series.begin(), model.series.end(), [](double v) { return std::isfinite(v); }))
        throw std::invalid_argument("tramo_tr0_forecast_errors: complete finite observations required");

    const double theta = model.theta.front(), seasonal_theta = model.btheta.front();
    const double interaction = theta * seasonal_theta;
    for (const double value : {theta, seasonal_theta}) {
        const auto check = arma_admissibility({value}, true);
        if (check.status != "admissible") return {check.status, {}};
    }

    // The differenced airline series is an MA polynomial with nonzero lags
    // 0, 1, 12 and 13. Condition the last 13 innovations on the observations.
    // Filtering a unit response in parallel gives the sensitivity to the
    // estimated drift and its unit-variance information.
    constexpr int lags = 13;
    using State = std::array<double, lags>;
    std::array<State, lags> covariance{};
    State drift_state{};
    for (int i = 0; i < lags; ++i) covariance[i][i] = 1.0;
    double drift_information = 0.0;
    for (int t = 0; t < fit.effective_observations; ++t) {
        State cross{};
        for (int i = 0; i < lags; ++i)
            cross[i] = theta * covariance[i][0] + seasonal_theta * covariance[i][11]
                     + interaction * covariance[i][12];
        const double variance = 1.0 + theta * cross[0] + seasonal_theta * cross[11]
                              + interaction * cross[12];
        if (!(variance > 0.0) || !std::isfinite(variance)) return {"covariance_range", {}};
        const double innovation = 1.0 - theta * drift_state[0] - seasonal_theta * drift_state[11]
                                - interaction * drift_state[12];
        drift_information += innovation * innovation / variance;
        // Shift in a new independent unit innovation, then condition the
        // shifted state on the observed value by a rank-one covariance update.
        State gain{};
        gain[0] = 1.0;
        std::copy_n(cross.begin(), lags - 1, gain.begin() + 1);
        for (int i = lags - 1; i >= 0; --i) {
            drift_state[i] = (i ? drift_state[i - 1] : 0.0) + gain[i] * innovation / variance;
            for (int j = lags - 1; j >= 0; --j) {
                const double prior = i && j ? covariance[i - 1][j - 1] : double(i == 0 && j == 0);
                covariance[i][j] = prior - gain[i] * gain[j] / variance;
            }
        }
    }
    if (!(drift_information > 0.0) || !std::isfinite(drift_information)) return {"covariance_range", {}};

    // g is the impulse response of 1 / ((1-B)(1-B^12)). Future
    // innovations are independent of the conditioned boundary state.
    const auto g = [](int k) { return k < 0 ? 0.0 : double(1 + k / 12); };
    const int df = fit.effective_observations - fit.estimated_parameters;
    const double process_scale = fit.residual_ssq / df;
    const double drift_scale = (fit.residual_ssq / (df + 1)) / drift_information;
    double future_variance = 0.0, drift_sum = 0.0;
    std::vector<double> errors(months);
    for (int h = 1; h <= months; ++h) {
        const double psi = g(h - 1) + theta * g(h - 2) + seasonal_theta * g(h - 13)
                         + interaction * g(h - 14);
        future_variance += psi * psi;
        drift_sum += g(h - 1);
        State weights{};
        for (int k = 0; k < lags; ++k) {
            if (k == 0) weights[k] += theta * g(h - 1);
            if (k < 12) weights[k] += seasonal_theta * g(h - 12 + k);
            weights[k] += interaction * g(h - 13 + k);
        }
        double process = future_variance, drift = drift_sum;
        for (int i = 0; i < lags; ++i) {
            drift -= weights[i] * drift_state[i];
            for (int j = 0; j < lags; ++j) process += weights[i] * covariance[i][j] * weights[j];
        }
        errors[h - 1] = std::sqrt(process * process_scale + drift * drift * drift_scale);
        if (!std::isfinite(errors[h - 1])) return {"covariance_range", {}};
    }
    return {"available", std::move(errors)};
}

}  // namespace kalos
