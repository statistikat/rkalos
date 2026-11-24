#ifndef KALOS_DETAIL_CALENDAR_AIRLINE_LIKELIHOOD_HPP
#define KALOS_DETAIL_CALENDAR_AIRLINE_LIKELIHOOD_HPP

#include <vector>

namespace kalos::detail {

enum class CalendarAirlineProcessor { gls, reg_sarima };

struct CalendarAirlineLikelihood {
    double log_likelihood;
    double ssq;
    std::vector<double> ma;
    // ConcentratedLikelihoodWithMissing surface (cl.dim()); always filled.
    int dimension = 0;
    // cl.coefficients() and cl.unscaledCovariance() (row major, square),
    // filled only when the call asked for with_diagnostics.
    std::vector<double> coefficients;
    std::vector<double> unscaled_covariance;
};

// RegArimaUtility.processor(true, precision), monthly model without a mean:
// seasonal=true is the airline (0,1,0)(0,1,1)[12], seasonal=false the regular
// (0,1,1)(0,0,0)[12] (one difference, one free MA(1) coefficient).
// with_diagnostics additionally surfaces the concentrated-likelihood
// coefficients and covariance (a collection pointer; it never changes the
// machine or any reported number).
CalendarAirlineLikelihood fit_calendar_airline(const std::vector<double>& series,
        const std::vector<std::vector<double>>& columns, double precision,
        CalendarAirlineProcessor processor = CalendarAirlineProcessor::gls,
        bool with_diagnostics = false, bool seasonal = true);

} // namespace kalos::detail
#endif
