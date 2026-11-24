#include "kalos/x13_estimated_calendar.hpp"
#include "kalos/detail/calendar_airline_likelihood.hpp"
#include "kalos/detail/fdlibm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {
namespace {
constexpr std::array<int, 4> durations{0, 1, 8, 15};

struct Input {
    std::vector<double> series;
    double adjustment = 0;
};

Input transform(const std::vector<double>& observations, const calendar_domain& domain,
                bool log, bool leap_adjustment) {
    Input input{observations};
    if (!log) return input;
    constexpr int difference_degree = 13;
    if (leap_adjustment) {
        const auto factors = x13_rsa2c_leap_factors(domain, false);
        for (std::size_t i = 0; i < factors.size(); ++i) {
            const double factor = factors[i];
            if (factor == 1.0) continue;
            input.series[i] *= factor;
            if (i >= difference_degree) input.adjustment += detail::fdlibm_log(factor);
        }
    }
    double sum = 0;
    for (std::size_t i = 0; i < input.series.size(); ++i) {
        input.series[i] = detail::fdlibm_log(input.series[i]);
        if (i >= difference_degree) sum += input.series[i];
    }
    input.adjustment -= sum;
    return input;
}

X13Rsa2cCalendarCandidate fit_candidate(const Input& input, const X13Rsa2cCalendar& calendar,
                                        bool log, const calendar_domain& domain) {
    const auto columns = x13_rsa2c_calendar_columns(calendar, log, domain);
    const auto fit = detail::fit_calendar_airline(input.series, columns, 1e-5);
    const int n = domain.length - 13;
    const int k = 3 + static_cast<int>(columns.size()); // two moving averages and the variance
    const double aicc = -2 * (fit.log_likelihood + input.adjustment -
                             static_cast<double>(n) * k / (n - k - 1));
    return {n, k, input.adjustment, fit.log_likelihood, fit.ssq, aicc};
}
} // namespace

std::vector<double> x13_rsa2c_leap_factors(const calendar_domain& domain, bool reverse) {
    if (domain.frequency != 12) throw std::invalid_argument("RSA2c requires monthly leap-year factors");
    auto days = length_of_period_column(calendar_lop::leap_year, domain);
    for (double& d : days) {
        const double actual = 28.25 + d;
        d = reverse ? actual / 28.25 : 28.25 / actual;
    }
    return days;
}

std::vector<std::vector<double>> x13_rsa2c_calendar_columns(const X13Rsa2cCalendar& calendar,
        bool log, const calendar_domain& domain) {
    if (domain.frequency != 12 || std::find(durations.begin(), durations.end(), calendar.easter_duration) == durations.end())
        throw std::invalid_argument("RSA2c requires monthly dates and a stock Easter duration");
    std::vector<std::vector<double>> columns;
    if (calendar.working_days) {
        columns = trading_days_contrast(day_clustering::td2, domain);
        if (!log) columns.push_back(length_of_period_column(calendar_lop::leap_year, domain));
    }
    if (calendar.easter_duration) {
        columns.push_back(easter_column(calendar.easter_duration, easter_correction::pre_computed, -1, domain));
    }
    return columns;
}

X13Rsa2cCalendarTest x13_rsa2c_calendar_test(const std::vector<double>& observations,
        int start_year, int start_month, bool log) {
    if (start_month < 1 || start_month > 12 || start_year < 1 || start_year > 9999 ||
        observations.size() <= 20 || observations.size() > std::numeric_limits<int>::max() ||
        std::any_of(observations.begin(), observations.end(), [log](double v) {
            return !std::isfinite(v) || (log && v <= 0);
        })) throw std::invalid_argument("RSA2c requires a complete monthly series on its selected scale");
    const calendar_domain domain{start_year, start_month, 12, static_cast<int>(observations.size())};
    const auto unadjusted = transform(observations, domain, log, false);
    const auto adjusted = transform(observations, domain, log, true);
    X13Rsa2cCalendarTest result;
    result.trading_days[0] = fit_candidate(unadjusted, {}, log, domain);
    result.trading_days[1] = fit_candidate(adjusted, {true, 0}, log, domain);
    result.selected.working_days = result.trading_days[1].aicc < result.trading_days[0].aicc;
    result.easter[0] = result.trading_days[result.selected.working_days];
    const auto& chosen = result.selected.working_days ? adjusted : unadjusted;
    std::size_t best = 0;
    for (std::size_t i = 1; i < durations.size(); ++i) {
        result.easter[i] = fit_candidate(chosen, {result.selected.working_days, durations[i]}, log, domain);
        if (result.easter[i].aicc < result.easter[best].aicc) best = i;
    }
    result.transformed = chosen.series;
    result.selected.easter_duration = durations[best];
    return result;
}

} // namespace kalos
