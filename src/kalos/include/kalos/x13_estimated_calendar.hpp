#ifndef KALOS_X13_ESTIMATED_CALENDAR_HPP
#define KALOS_X13_ESTIMATED_CALENDAR_HPP

#include "calendar_generation.hpp"
#include <array>
#include <vector>

namespace kalos {

struct X13Rsa2cCalendar {
    bool working_days = false;
    int easter_duration = 0; // 0, 1, 8 or 15
};

struct X13Rsa2cCalendarCandidate {
    int nobs = 0;
    int parameters = 0;
    double transformation_adjustment = 0;
    double log_likelihood = 0;
    double ssq = 0;
    double aicc = 0;
};

struct X13Rsa2cCalendarTest {
    X13Rsa2cCalendar selected;
    std::vector<double> transformed; // selected estimation-scale observations
    std::array<X13Rsa2cCalendarCandidate, 2> trading_days; // absent, present
    std::array<X13Rsa2cCalendarCandidate, 4> easter; // absent, 1, 8 and 15 days
};

// February factors used by the stock log working-day model. Reverse factors
// reproduce the separate calendar back-transformation on a forecast domain.
std::vector<double> x13_rsa2c_leap_factors(const calendar_domain& domain, bool reverse);

std::vector<std::vector<double>> x13_rsa2c_calendar_columns(const X13Rsa2cCalendar& calendar,
        bool log, const calendar_domain& domain);

// Stock RSA2c calendar tests following transformation selection, before outlier detection.
X13Rsa2cCalendarTest x13_rsa2c_calendar_test(const std::vector<double>& observations,
        int start_year, int start_month, bool log);

} // namespace kalos
#endif
