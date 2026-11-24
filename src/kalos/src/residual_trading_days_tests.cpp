// M7 (item 27 cutover gap): the generic-SA trading-day residual diagnostics
// (diagnostics.td-sa-all / td-sa-last / td-i-all / td-i-last), transcribed
// from the pinned jdplus.sa.base.core.diagnostics.ResidualTradingDaysTests:
// TradingDaysTest.maTest(series, false) over the six GenericTradingDays TD7
// contrast columns (Regression.matrix, no mean), the monthly ARIMA
// (0,1,1)(0,0,0)[12] — one regular difference, one free MA(1) coefficient —
// fitted through RegSarimaComputer.PROCESSOR (precision 1e-9, the same
// machine detail::fit_calendar_airline exposes under the reg_sarima
// processor with seasonal=false), and the pinned JointTest with
// hyperParametersCount(1) over the concentrated likelihood.
#include "kalos/residual_trading_days_tests.hpp"
#include "kalos/calendar_generation.hpp"
#include "kalos/detail/calendar_airline_likelihood.hpp"
#include "kalos/detail/special_functions.hpp"
#include "mlj_engine.hpp"
#include "pinned_cholesky.hpp"

#include <cmath>
#include <stdexcept>

namespace kalos {
namespace {
// GenericTradingDays.contrasts(DayClustering.TD7) yields six columns; the
// pinned JointTest counts the one hyper-parameter (the residual variance)
// outside the six contrast coefficients.
constexpr int contrast_count = 6;
constexpr int hyper_parameters = 1;
// The PROCESSOR singleton precision (not the RegArimaUtility 1e-7 default).
constexpr double processor_precision = 1e-9;

TradingDaysFTest ma_test(const std::vector<double>& series, int start_year, int start_month) {
    if (start_year < 1 || start_year > 9999 || start_month < 1 || start_month > 12)
        throw std::invalid_argument("residual trading-day test requires a valid monthly domain");
    // One differencing loss, six contrasts and the one hyper-parameter.
    if (series.size() <= static_cast<std::size_t>(1 + contrast_count + hyper_parameters))
        throw std::invalid_argument("residual trading-day test requires positive degrees of freedom");
    for (const double value : series)
        if (!std::isfinite(value))
            throw std::invalid_argument("residual trading-day test requires finite observations");
    const calendar_domain domain{start_year, start_month, 12, static_cast<int>(series.size())};
    const auto columns = trading_days_contrast(day_clustering::td7, domain);
    const auto fit = detail::fit_calendar_airline(series, columns, processor_precision,
                                                  detail::CalendarAirlineProcessor::reg_sarima, true, false);
    if (!std::isfinite(fit.ssq) || !(fit.ssq > 0.0))
        throw std::runtime_error("residual trading-day fit has no usable residual variance");
    if (fit.coefficients.size() != static_cast<std::size_t>(contrast_count) ||
        fit.unscaled_covariance.size() != static_cast<std::size_t>(contrast_count) * contrast_count)
        throw std::runtime_error("residual trading-day fit has no full contrast covariance");
    // JointTest.build(): the covariance is the symmetric
    // cl.unscaledCovariance(); lcholesky and solveLx consume its lower
    // triangle, identical in either storage order for a symmetric matrix.
    const double degrees_of_freedom = static_cast<double>(fit.dimension) - contrast_count - hyper_parameters;
    if (!(degrees_of_freedom > 0.0))
        throw std::invalid_argument("residual trading-day test requires positive degrees of freedom");
    auto covariance = fit.unscaled_covariance;
    auto beta = fit.coefficients;
    if (!detail::pinned_lcholesky(covariance, contrast_count))
        throw std::runtime_error("residual trading-day covariance is not factorable");
    try {
        mlj::solve_lx(covariance, contrast_count, beta);
    } catch (const std::runtime_error&) {
        throw std::runtime_error("residual trading-day covariance solve failed");
    }
    const double statistic = (mlj::ssq(beta.data(), beta.size()) / contrast_count) /
        (fit.ssq / degrees_of_freedom);
    return {statistic, detail::f_upper_p(statistic, contrast_count, degrees_of_freedom)};
}
} // namespace

ResidualTradingDaysTests residual_trading_days_tests(const std::vector<double>& sa,
    int sa_start_year, int sa_start_month, const std::vector<double>& irregular,
    int irregular_start_year, int irregular_start_month) {
    const auto sa_test = ma_test(sa, sa_start_year, sa_start_month);
    const auto irregular_test = ma_test(irregular, irregular_start_year, irregular_start_month);
    // The recorded option is flast = 0, where td2 tests the whole component:
    // ResidualTradingDaysTests computes one statistic per component and the
    // dashboard publishes it under both names.
    return {sa_test, sa_test, irregular_test, irregular_test};
}

} // namespace kalos
