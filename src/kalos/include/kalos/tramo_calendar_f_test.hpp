#ifndef KALOS_TRAMO_CALENDAR_F_TEST_HPP
#define KALOS_TRAMO_CALENDAR_F_TEST_HPP

#include "kalos/tramo_differencing_battery.hpp"

#include <array>
#include <optional>
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

struct TramoStockCalendar {
    TramoTradingDays trading_days = TramoTradingDays::None;
    bool leap_year = false;
    bool easter = false;
};

struct TramoFittedCalendar {
    TramoStockCalendar selection;
    std::vector<double> coefficients; // TD contrasts, leap year, then Easter
};

struct calendar_domain;
std::vector<std::vector<double>> tramo_stock_calendar_columns(const TramoStockCalendar& selection,
    const calendar_domain& domain);

std::vector<double> tramo_stock_calendar_effect(const TramoFittedCalendar& fit, const calendar_domain& domain);

// Original-series scale of the public cal and cal_f getters.
std::vector<double> tramo_stock_calendar_effect_original_scale(const TramoFittedCalendar& fit,
    const calendar_domain& domain, bool log);

struct TramoCalendarFTestResult {
    double base_ssq;
    TramoCalendarFCandidate working_days;
    TramoCalendarFCandidate trading_days;
    TramoTradingDays selected;
};

TramoCalendarFTestResult tramo_calendar_f_test(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision);

// Initial TRfull regression module on a monthly airline description with a
// mean, no existing regressors, and no length-of-period preadjustment.
struct TramoFullCalendarResult {
    TramoCalendarFTestResult comparison;
    double mean_t;
    std::optional<double> leap_year_t;
    double easter_t;
    bool keep_mean;
    bool keep_leap_year;
    bool keep_easter;
    // The module marks changes to mean, LP and Easter in its processing status.
    bool module_changed() const { return !keep_mean || keep_leap_year || keep_easter; }
};

TramoFullCalendarResult tramo_full_initial_calendar(const std::vector<double>& series,
    int start_year, int start_month, bool log);

// When the seasonality pretest clears seasonal differencing and seasonal MA,
// all calendar candidates inherit the resulting nonseasonal IMA(1,1) model.
TramoFullCalendarResult tramo_full_initial_calendar(const std::vector<double>& series,
    int start_year, int start_month, bool log, bool seasonal);

// The late TradingDaysController checks the concentrated-likelihood residuals
// against a working-days contrast on their actual right-aligned monthly domain.
struct TramoResidualTradingDays {
    double statistic;
    double upper_probability;
    bool needs_processing;
};

// A missing result means OLS could not be evaluated; the controller skips the refit.
std::optional<TramoResidualTradingDays> tramo_residual_trading_days(
    const std::vector<double>& residuals, int estimation_start_year, int estimation_start_month,
    int estimation_length, double significance);

// Late TradingDaysController for a monthly airline model on its estimation
// scale. The initial model has an optional mean. Every candidate refit is
// recorded, including candidates rejected by the comparison.
struct TramoCalendarAirlineFit {
    bool mean;
    bool working_days;
    std::vector<double> ma_parameters;
    std::vector<double> coefficients;
    std::vector<double> residuals;
    double residual_ssq;
    std::optional<double> mean_t;
    BatteryPanel statistics;
    // Automatically detected (position, AO/LS/TC factory index), in regression order.
    std::vector<std::array<int, 2>> outliers;
};

struct TramoTradingDaysControllerResult {
    TramoCalendarAirlineFit original;
    std::optional<TramoResidualTradingDays> residual_test;
    std::vector<TramoCalendarAirlineFit> candidates;
    bool changed;
    const TramoCalendarAirlineFit& selected() const { return changed ? candidates.back() : original; }
};

TramoTradingDaysControllerResult tramo_airline_trading_days_controller(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision);

// Detect AO, LS and TC variables on each fitted description when a critical value is supplied.
TramoTradingDaysControllerResult tramo_airline_trading_days_controller(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision,
    std::optional<double> outlier_critical_value);

} // namespace kalos
#endif
