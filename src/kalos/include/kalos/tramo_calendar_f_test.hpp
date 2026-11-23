#ifndef KALOS_TRAMO_CALENDAR_F_TEST_HPP
#define KALOS_TRAMO_CALENDAR_F_TEST_HPP

#include <vector>

namespace kalos {

// TRAMO's AutomaticFRegressionTest candidate comparison. The observations
// are on the model's estimation scale. This entry fits a monthly airline
// model with optional mean and no pre-existing regressors. Selection here
// precedes the module's LP, Easter and mean significance updates.
struct TramoCalendarFCandidate {
    double residual_ssq;
    int degrees_of_freedom;  // Before subtracting the two free MA parameters.
    double mean_square;
    double lower_probability;
};

enum class TramoTradingDays { None, WorkingDays, TradingDays };

struct TramoCalendarFTestResult {
    double base_ssq;
    TramoCalendarFCandidate working_days;
    TramoCalendarFCandidate trading_days;
    TramoTradingDays selected;
};

TramoCalendarFTestResult tramo_calendar_f_test(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision);

} // namespace kalos
#endif
