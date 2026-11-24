#include "kalos/tramo_differencing_battery.hpp"
#include "kalos/tramo_regression_tests.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "mlj_engine.hpp"
#include "pinned_cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {

BatteryFTest tramo_seasonal_f_test(const std::vector<double>& series, int start_year, int start_month,
    const TramoOrders& orders, const std::vector<std::vector<double>>& regressors) {
    constexpr int period = 12, contrasts = period - 1;
    if (start_year < 1 || start_year > 9999 || start_month < 1 || start_month > period ||
        orders.p < 0 || orders.d < 0 || orders.d > 2 || orders.q < 0 ||
        series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("invalid monthly domain or regular orders for seasonal F test");
    const auto finite = [](const std::vector<double>& v) {
        return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
    };
    if (!finite(series)) throw std::invalid_argument("seasonal F test requires complete observations");
    for (const auto& col : regressors) {
        if (col.size() != series.size() || !finite(col))
            throw std::invalid_argument("seasonal F test regressor does not match the observation domain");
    }
    const int n = static_cast<int>(series.size());
    const auto df = static_cast<long long>(n) - orders.d - orders.p - orders.q
        - static_cast<long long>(regressors.size()) - period;
    if (df <= 0) return {};
    // JDemetra drops the seasonal ARIMA orders and forces a mean after differencing.
    const auto filter = detail::pinned_differencing_polynomial(period, orders.d, 0);
    const auto y = detail::pinned_backfilter_apply(filter, series);
    std::vector<std::vector<double>> design{std::vector<double>(y.size(), 1.0)};
    SarimaEstimationRequest request{};
    request.series.values = series;
    request.series.missing.resize(series.size());
    request.period = period;
    request.gls_sarima_mapping = true;
    request.d = orders.d;
    request.design.values.assign(series.size(), 1.0);
    auto append = [&](const std::vector<double>& col) {
        request.design.values.insert(request.design.values.end(), col.begin(), col.end());
        design.push_back(detail::pinned_backfilter_apply(filter, col));
    };
    for (const auto& col : regressors) append(col);
    for (int month = 1; month < period; ++month) {
        std::vector<double> col(series.size(), 0.0);
        for (int t = 0; t < n; ++t) {
            const int phase = (start_month - 1 + t % period) % period + 1;
            if (phase == month) col[t] = 1.0;
            else if (phase == period) col[t] = -1.0;
        }
        append(col);
    }
    request.design.columns = design.size();
    request.design.coefficients.resize(design.size());
    request.design.fixed.resize(design.size());
    request.control = {1e-7, 1e-4, true, false, true, true,
        EstimationStartingPoint::Default, EstimationSolver::PinnedLevenbergMarquardt};
    std::array<std::vector<double>, 4> start;
    start[0].assign(orders.p, 0.1);
    start[1].assign(orders.q, -0.2);
    for (int b : {0, 1}) {
        request.blocks[b].initial = start[b];
        request.blocks[b].fixed.resize(start[b].size());
    }
    std::vector<double> jvm;
    if (orders.p || orders.q) {
        const auto fitted = optimize_sarima_from(request, start);
        if (!std::isfinite(fitted.log_likelihood)) return {};
        for (double ar : fitted.coefficients[0]) jvm.push_back(-ar);
        for (double ma : fitted.coefficients[1]) jvm.push_back(ma);
    }
    detail::PinnedRegressionDiagnostics diagnostics;
    const auto fitted = detail::pinned_deviances_evaluate_multi(y, orders.p, 0, orders.q, 0, period,
        design, jvm, detail::PinnedComputerShape::Fixed, &diagnostics);
    if (!fitted) return {};
    // JointTest takes the trailing block of unscaled coefficient covariance.
    const auto first = design.size() - contrasts;
    std::vector<double> beta(fitted->b.begin() + first, fitted->b.end());
    std::vector<double> covariance(contrasts * contrasts);
    for (int col = 0; col < contrasts; ++col)
        for (int row = 0; row < contrasts; ++row)
            covariance[col * contrasts + row] = diagnostics.unscaled_covariance[
                (first + row) * design.size() + first + col];
    if (!detail::pinned_lcholesky(covariance, contrasts)) return {};
    try {
        mlj::solve_lx(covariance, contrasts, beta);
    } catch (const std::runtime_error&) {
        return {};
    }
    const double statistic = (mlj::ssq(beta.data(), beta.size()) / contrasts) /
        (fitted->ssq_err / static_cast<double>(df));
    return {true, statistic, detail::f_upper_p(statistic, contrasts, df), contrasts, static_cast<double>(df)};
}

} // namespace kalos
