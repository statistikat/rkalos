#ifndef KALOS_CALENDAR_GENERATION_HPP
#define KALOS_CALENDAR_GENERATION_HPP

#include <string>
#include <vector>

namespace kalos {

// Item 21-I.1: the calendar-generation mirror. This unit reproduces, column
// for column and bit for bit, what the LIVE pinned jdplus 3.9.0 factories
// emit as regression design data (jdplus.toolkit.base.core.modelling.
// regression.{GenericTradingDaysFactory, LPFactory, EasterFactory,
// JulianEasterFactory, RampFactory, TsVariablesFactory}, driven through
// Regression.matrix over a TsDomain): the day-count table tdCount, the
// meanDays normalizers, the day-clustering contrast columns, the
// LengthOfPeriod/LeapYear columns, the Gregorian and Julian Easter computus
// and windows, ramps, user-variable placement, and the pinned name and
// description conventions. It is the design-algebra layer the 23-A builder
// registered variables for (deferred there per register D10); NO detection
// decisions live here (the {none, 1, 8, 15} candidate batteries are
// 21-I.2/I.3).
//
// Oracle: tmp_opencode/i21probe double captures (validation/evidence/native/
// calendar-generation/, byte-identical across two OpenJDK 21 JVMs; primary
// 4128-line gen.out + 5904-line envelope-shape gen2.out). Every machine was
// replayed numerically bit-exact (3960 + 5672 matrix rows, computus,
// probability tables) before this mirror was written.
//
// Envelope (fail-closed std::invalid_argument, never a silent
// reinterpretation):
//   - calendar_domain: frequency 12, 4 or 1 with an aligned first month,
//     positive length, year in [1, 9999] (the corpus never leaves it);
//   - trading-day and Easter-family functions accept monthly and quarterly
//     domains - the faithful jar gates (GenericTradingDaysFactory.fill and
//     EasterFactory.fill return false for annualFrequency <= 3, and only 12
//     and 4 pass; JulianEasterFactory rejects annualFrequency < 3). Where the
//     jar would silently leave the (zero-initialized) column empty, this API
//     REFUSES instead of emitting data the pinned path never demonstrates;
//   - LeapYear / LengthOfPeriod / ramp / user columns accept {12, 4, 1};
//   - Easter durations [1, 28] (the Julian factory indexes 28-entry cycle
//     tables and throws ArrayIndexOutOfBounds beyond 28; the Gregorian
//     envelope adopts the same bound, and PreComputed indexes 25-entry mean
//     tables so it requires [1, 25]); endPosition is faithful: only -1 and
//     +1 are special, any other value behaves as 0;
//   - user data must carry the domain's own frequency (mixed-unit alignment
//     arithmetic is unprobed, hence refused).
//
// Faithful quirk mirrors (all oracle-pinned, do not "fix"): the reversed
// ramp (start after end) is the ALL-ZERO column, not a mirrored staircase;
// the LeapYear +0.75 strikes every 4 periods from the FIRST leap year
// (never skipping century non-leap years), and the search for that first
// leap year is bounds-free (TsDomain.get is pure arithmetic), so a range
// without one still completes with only the -0.25 cells set; GenericTradingDays.raw
// and .contrasts produce BIT-IDENTICAL columns (GenericTradingDaysVariable
// discards the contrast flag and the factory always runs the contrast
// machine); the Gregorian Easter factory writes month-pair shares with SET,
// the Julian factory accumulates with ADD (and on a quarterly domain its
// March and April contributions land on the same cell and add).

// TsDomain mirror. `month` is the first month of the first period, aligned
// with the frequency: 1..12 (monthly), 1/4/7/10 (quarterly), 1 (annual).
// The epoch convention is the jdplus DEFAULT_EPOCH (position arithmetic in
// months); every captured domain is epoch-normal.
struct calendar_domain {
    int year = 0;
    int month = 1;
    int frequency = 12;
    int length = 0;
};

// DayClustering mirror (the 9 clusterings the pinned generator reaches;
// TD2d's TRADING-DAYS TYPE is TDuser - the jar tag/nameOf print "TDuser").
enum class day_clustering { td2, td2c, td2d, td3, td3c, td4, td4c, td6, td7 };

// LengthOfPeriodType mirror. None keeps dim 1 but fills nothing: through
// Regression.matrix the column is the zero column.
enum class calendar_lop { none, leap_year, length_of_period };

// EasterVariable.Correction mirror (builder default is simple).
enum class easter_correction { none, simple, pre_computed, theoretical };

// GenericTradingDaysFactory.tdCount: per-cell weekday counts, rows
// Monday..Sunday, columns = domain periods. The quarterly machine counts
// every day whose date falls inside each period span.
std::vector<std::vector<int>> td_count(const calendar_domain& domain);

// GenericTradingDaysFactory.meanDays: mean_days_yearly(12) is the pinned
// MDAYS/7 twelve-vector; mean_days cycles it from the domain start for
// monthly domains and sums the per-cell member months for lower
// frequencies (captured bit-exact at 12, 4 and 1).
std::vector<double> mean_days(const calendar_domain& domain);
std::vector<double> mean_days_yearly(int frequency);

// GenericTradingDays contrast columns (columns = groups minus the reference
// group, ascending group ids, reference = lowest id, x = refsum/|ref|
// divided FIRST, col = sum_g - |g|*x). The raw(GenericTradingDays) path is
// bit-identical (see header note), so this single machine serves both.
std::vector<std::vector<double>> trading_days_contrast(day_clustering clustering,
                                                       const calendar_domain& domain);

// LPFactory: one design column. none yields the zero column; leap_year and
// length_of_period follow the pinned machines (see header notes).
std::vector<double> length_of_period_column(calendar_lop type,
                                            const calendar_domain& domain);

// EasterFactory (Gregorian): the dim-1 design column, duration [1, 28],
// endPosition faithful (-1/+1 special, others behave as 0).
std::vector<double> easter_column(int duration, easter_correction correction,
                                  int end_position, const calendar_domain& domain);

// JulianEasterFactory: dim-1 design column, duration [1, 28] (28-entry
// cycle tables), gregorianDates selects the JDN-converted date.
std::vector<double> julian_easter_column(int duration, bool gregorian_dates,
                                         const calendar_domain& domain);

// Calendar.Easter computus mirror. easter_date is the single Gregorian
// computus the factory uses (easter == easter2 pinned 1985..2015);
// julian_easter_date converts through the jar JDN pair. prob_julian_easter
// guards [0, 43) (JD[43]/1064), prob_easter guards [0, 35) with the three
// pinned branches (d<6 ramp, 6..27 constant, 28..34 tail).
struct calendar_date {
    int year = 0, month = 0, day = 0;
    bool operator==(const calendar_date& o) const {
        return year == o.year && month == o.month && day == o.day;
    }
};
calendar_date easter_date(int year);
calendar_date julian_easter_date(int year, bool gregorian_dates);
double prob_julian_easter(int day);
double prob_easter(int day);

// RampFactory: ramp from (start_year, start_month) to (end_year,
// end_month), endpoints clamped through TsPeriod.of truncation. The column
// is -1 at/before the start position, linear ((i - i1)/n - 1) strictly
// between, 0 at/after the end; end <= start gives the all-zero column and a
// start past the last cell gives all -1.
std::vector<double> ramp_column(int start_year, int start_month, int end_year,
                                int end_month, const calendar_domain& domain);

// TsVariablesFactory (user data, same frequency as the domain): the data is
// placed at its CALENDAR positions inside the domain, out-of-domain cells
// keep 0.0, and a data range not overlapping the domain gives the zero
// column (the pinned "too short" warning path).
std::vector<double> user_column(const std::vector<double>& data, int data_year,
                                int data_month, const calendar_domain& domain);

// Pinned name/description conventions (oracle-echoed strings).
// trading_days_tag: "TD2".."TD7" plus "TDuser" for td2d;
// trading_days_name: nameOf/description ("td2".."td7", "tduser");
// trading_days_default_name(k): the pinned ITradingDaysVariable.defaultName
// (always "td"); trading_days_variable_name: ITradingDaysVariable.NAME
// ("td"); length_of_period_variable_name: "lp"; easter_description:
// "easter [d]"; julian_easter_description: "julian easter [d]";
// length_of_period_description: None/LengthOfPeriod "Length of period",
// LeapYear "Leap year"; ramp_description: "ramp 1992-03 / 1994-06" (monthly
// precision, the captured format - other domain frequencies are refused).
const char* trading_days_tag(day_clustering clustering);
const char* trading_days_name(day_clustering clustering);
const char* trading_days_default_name(int coefficient_count);
const char* trading_days_variable_name();
const char* length_of_period_variable_name();
std::string easter_description(int duration);
std::string julian_easter_description(int duration);
const char* length_of_period_description(calendar_lop type);
std::string ramp_description(int start_year, int start_month, int end_year,
                             int end_month, int frequency);

// TradingDaysType pinned table (nvar, names/contrastNames as the jar
// Arrays.toString renders them - "null" where the jar prints null - and the
// slash description; TDuser renders nvar -1). Index: NONE, TD7, TD4, TD3c,
// TD3, TD2, TD2c, TD2d, TD4c, TD6, TDuser (the jar values() order).
struct trading_days_type_info {
    const char* tag;
    int nvar;
    const char* names;
    const char* contrast_names;
    const char* description;
};
const std::vector<trading_days_type_info>& trading_days_types();

}  // namespace kalos

#endif
