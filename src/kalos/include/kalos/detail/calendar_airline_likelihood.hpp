#ifndef KALOS_DETAIL_CALENDAR_AIRLINE_LIKELIHOOD_HPP
#define KALOS_DETAIL_CALENDAR_AIRLINE_LIKELIHOOD_HPP

#include <vector>

namespace kalos::detail {

enum class CalendarAirlineProcessor { gls, reg_sarima };

struct CalendarAirlineLikelihood {
    double log_likelihood;
    double ssq;
    std::vector<double> ma;
};

// RegArimaUtility.processor(true, precision), monthly airline without a mean.
CalendarAirlineLikelihood fit_calendar_airline(const std::vector<double>& series,
        const std::vector<std::vector<double>>& columns, double precision,
        CalendarAirlineProcessor processor = CalendarAirlineProcessor::gls);

} // namespace kalos::detail
#endif
