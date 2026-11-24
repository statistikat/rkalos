#ifndef KALOS_DETAIL_CALENDAR_AIRLINE_LIKELIHOOD_HPP
#define KALOS_DETAIL_CALENDAR_AIRLINE_LIKELIHOOD_HPP

#include <vector>

namespace kalos::detail {

struct CalendarAirlineLikelihood {
    double log_likelihood;
    double ssq;
};

// RegArimaUtility.processor(true, precision), monthly airline without a mean.
CalendarAirlineLikelihood fit_calendar_airline(const std::vector<double>& series,
        const std::vector<std::vector<double>>& columns, double precision);

} // namespace kalos::detail
#endif
