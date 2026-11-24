#include "kalos/tramo_tr0.hpp"
#include "kalos/calendar_generation.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/seats_exact_forecasts.hpp"
#include "kalos/outlier_scoring.hpp"
#include "kalos/tramo_seats_bridge.hpp"

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace kalos {
namespace {

namespace xf = detail::seats_xf;
using Vector = std::vector<double>;

// Propagate the first d observations with the unit-root polynomial to
// initialize the diffuse mean. The stationary contribution starts at zero.
Vector initial_state(const Vector& values, const Vector& difference, int dimension) {
    const auto degree = difference.size() - 1;
    Vector history(values.begin(), values.begin() + degree);
    for (int i = 0; i < dimension; ++i) {
        const auto t = history.size();
        double value = 0;
        for (std::size_t k = 1; k <= degree; ++k) value -= difference[k] * history[t - k];
        history.push_back(value);
    }
    return Vector(history.begin() + degree, history.end());
}

}  // namespace

tramo_tr0_forecast_errors_result tramo_selected_forecast_errors(const TramoFitResult& fit, int months) {
    if (months < 0) throw std::invalid_argument("tramo_selected_forecast_errors: negative forecast count");
    const auto model = make_tramo_seats_model(fit);
    const auto n = model.series.size(), p = fit.regression_coefficients.size();
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max() - months))
        throw std::length_error("tramo_selected_forecast_errors: extended domain exceeds the index range");
    const auto differencing_order = static_cast<long long>(model.d) + 12LL * model.bd;
    if (differencing_order != static_cast<long long>(n) - fit.ls_eff)
        throw std::invalid_argument("tramo_selected_forecast_errors: invalid effective sample size");
    const auto covariance_size = fit.covariance.covariance.size();
    const bool square = p ? covariance_size / p == p && covariance_size % p == 0 : covariance_size == 0;
    if (p != fit.variables.size() + fit.calendar.coefficients.size() + static_cast<std::size_t>(fit.mean) || p >= n ||
        fit.covariance.n != p || !square ||
        static_cast<std::size_t>(fit.ls_npar) != fit.parameters.size() + p + 1)
        throw std::invalid_argument("tramo_selected_forecast_errors: invalid fitted coefficients or covariance");

    if (months == 0) return {"no_requested_predictions", {}};

    detail::Sarima arima;
    arima.period = 12;
    arima.d = model.d;
    arima.bd = model.bd;
    arima.phi = model.phi;
    arima.bphi = model.bphi;
    arima.theta = model.theta;
    arima.btheta = model.btheta;
    const auto difference = arima.non_stationary_ar();
    const auto degree = static_cast<std::size_t>(difference.degree());


    const auto stationary = arima.stationary_ar(), ma = arima.ma();
    const auto full = detail::p_times(stationary, difference);
    const int dim = std::max(full.degree(), ma.degree() + 1);
    const xf::ArimaDyn transition{full.c, dim};
    auto acf = xf::Acf::of_polys(ma, stationary, 1);
    const auto stationary_psi = xf::RatFunc::of(ma, stationary).coefficients(dim);
    Vector autocovariances;
    try {
        autocovariances = acf.values(dim);
    } catch (const detail::seats_acf::arima_exception&) {
        return {"covariance_range", {}};
    }
    auto covariance = xf::p0_matrix(1, autocovariances, stationary_psi);
    if (degree) covariance = xf::xsxt(covariance, xf::sigma_matrix(difference.c, dim));
    const auto shock = xf::RatFunc::of(ma, full).coefficients(dim);
    const auto disturbance = xf::v_matrix(shock, 1);

    // Response followed by the original-domain drift, calendar and outlier columns.
    // A unit drift integrates a column of ones on the differenced domain.
    const auto length = n + months;
    std::vector<Vector> data{fit.builder.interpolated};
    if (fit.log) for (auto& value : data.front()) value = detail::fdlibm_log(value);
    if (fit.mean) {
        auto& drift = data.emplace_back(length);
        for (std::size_t t = degree; t < length; ++t) {
            drift[t] = 1;
            for (std::size_t k = 1; k <= degree; ++k) drift[t] -= difference.c[k] * drift[t - k];
        }
    }
    auto calendar = tramo_stock_calendar_columns(fit.calendar.selection,
        {model.start_year, model.start_month, 12, static_cast<int>(length)});
    if (calendar.size() != fit.calendar.coefficients.size())
        throw std::invalid_argument("tramo_selected_forecast_errors: inconsistent calendar design");
    for (auto& column : calendar) data.push_back(std::move(column));
    for (const auto& effect : fit.variables) {
        auto& column = data.emplace_back(length);
        fill_outlier_column(static_cast<OutlierType>(effect.type), effect.position, 12, effect.tc_rate, column);
    }
    std::vector<Vector> states;
    for (const auto& column : data) states.push_back(initial_state(column, difference.c, dim));

    Vector errors(months), gain(dim), sensitivity(p);
    for (std::size_t t = degree; t < length; ++t) {
        const double variance = covariance.get(0, 0);
        if (!(variance > 0) || !std::isfinite(variance)) return {"covariance_range", {}};
        if (t < n) {
            for (int i = 0; i < dim; ++i) gain[i] = covariance.get(i, 0);
            for (std::size_t j = 0; j < data.size(); ++j) {
                const double innovation = data[j][t] - states[j][0];
                xf::db_adday(states[j], innovation / variance, gain.data(), 1);
            }
            xf::add_xaxt(covariance, -1 / variance, gain);
        } else {
            double total = model.innovation_variance * variance;
            for (std::size_t j = 0; j < p; ++j) sensitivity[j] = data[j + 1][t] - states[j + 1][0];
            for (std::size_t i = 0; i < p; ++i)
                for (std::size_t j = 0; j < p; ++j)
                    total += sensitivity[i] * fit.covariance.covariance[i * p + j] * sensitivity[j];
            double error = std::sqrt(total);
            if (fit.log) {
                const double location = states[0][0] + std::inner_product(sensitivity.begin(), sensitivity.end(),
                                                                   fit.regression_coefficients.begin(), 0.0);
                const double log_variance = error * error;
                error = detail::fdlibm_exp(location + 0.5 * log_variance) * std::sqrt(detail::fdlibm_exp(log_variance) - 1);
            }
            if (!std::isfinite(error)) return {"covariance_range", {}};
            errors[t - n] = error;
        }
        for (auto& state : states) transition.tx(state.data());
        transition.tvt(covariance);
        for (std::size_t i = 0; i < covariance.a.size(); ++i) covariance.a[i] += disturbance.a[i];
    }
    return {"available", std::move(errors)};
}

tramo_tr0_forecast_errors_result tramo_tr3_forecast_errors(const TramoFitResult& fit, int months) {
    return tramo_selected_forecast_errors(fit, months);
}

}  // namespace kalos
