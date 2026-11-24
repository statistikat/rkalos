#ifndef KALOS_RESIDUAL_TRADING_DAYS_TESTS_HPP
#define KALOS_RESIDUAL_TRADING_DAYS_TESTS_HPP

#include <vector>

namespace kalos {

// M7 (item 27 cutover gap): the generic-SA trading-day residual diagnostics,
// the four public scalars diagnostics.td-sa-all, diagnostics.td-sa-last,
// diagnostics.td-i-all and diagnostics.td-i-last. The JVM machine is
// jdplus.sa.base.core.diagnostics.ResidualTradingDaysTests: for each of the
// seasonally-adjusted and irregular components the test td2(series, flast,
// false) feeds TradingDaysTest.maTest, which fits the monthly ARIMA
// (0,1,1)(0,0,0)[12] (one regular difference, one free MA(1) coefficient)
// with the six GenericTradingDays TD7 contrast columns as the only design
// (no mean) through RegSarimaComputer.PROCESSOR (the pinned
// singleton at precision 1e-9, HannanRissanen initializer, Multiple warm
// starts) and then evaluates the pinned JointTest over the concentrated
// likelihood with hyperParametersCount(1): F = (ssq(b)/6)/(ssq/(dim-6-1))
// with the F(6, dim-7) upper probability. The pinned option records flast =
// 0 (validation/evidence/native/residual-trading-days/public-x86_64.txt
// "case <id> <flast>" lines), where td2 tests the whole component, so the
// "-all" and "-last" public names of a component are computed by the one
// same fit and agree bit for bit; this entry reproduces that machine for
// that recorded option and surfaces both names per component.
//
// The replay tests/test_residual_trading_days.cpp bit-compares all four
// names on the 24 double-JVM capture cases (192 doubles).
struct TradingDaysFTest {
    double statistic;
    double p_value;
};

struct ResidualTradingDaysTests {
    TradingDaysFTest sa_all;        // diagnostics.td-sa-all
    TradingDaysFTest sa_last;       // diagnostics.td-sa-last (flast = 0: same fit)
    TradingDaysFTest irregular_all; // diagnostics.td-i-all
    TradingDaysFTest irregular_last;// diagnostics.td-i-last (flast = 0: same fit)
};

// The components run on their own (extraction) domains, given here by their
// start year and month; the frequency is monthly. Both series must be finite
// and long enough for positive test degrees of freedom (more than 20
// observations each).
ResidualTradingDaysTests residual_trading_days_tests(const std::vector<double>& sa,
    int sa_start_year, int sa_start_month, const std::vector<double>& irregular,
    int irregular_start_year, int irregular_start_month);

} // namespace kalos

#endif
