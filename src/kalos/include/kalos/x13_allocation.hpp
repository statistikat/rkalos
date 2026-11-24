#ifndef KALOS_X13_ALLOCATION_HPP
#define KALOS_X13_ALLOCATION_HPP
// Item 22c: the deterministic-effect allocation staging of the pinned X-13
// pipeline — the fixed-coefficient preadjust (X13Kernel.initialStep) and the
// finals restore (X13Kernel.finals) — decoded from the pinned rjd3 3.9.0
// bytecode and qualified bit-exactly against the pinned kernel on the
// validation/evidence/native/x13-allocation corpus. The kernel gives no
// fixed-coefficient entry point (every X13Kernel run estimates the ARIMA
// part), so the contract here is the exact seam the capture proves: given
// the COMPLETE input series, the pinned model-description variable order
// with each variable's pinned regression columns and its coefficient values
// as INPUTS, the staging tables are pure algebra and this header reproduces
// them; everything upstream of the seam (coefficient and ARIMA estimation,
// outlier detection, trading-day tests) and everything the fixed-model path
// with zero horizons and LengthOfPeriodType.None never reaches (missing-row
// interpolation through MissingValueEstimation, the length-of-period
// transform, benchmarking and diagnostics) is outside this contract. The
// nonzero-horizon ARIMA extension route through ExactArimaForecasts
// (the a1a/a1b extensions and the split finals) landed with item 21-H.4 and
// is qualified by docs/contracts/arima-forecasts-contract.md and the
// automatic-pipeline horizon captures. See docs/contracts/x13-allocation-contract.md and
// docs/contracts/x13-allocation-scope-22c.md for the scoping memo.
//
// Pinned semantics reflected here (each pinned bit-exact by a capture
// check): every a-table is backTransform'ed for storage (exp under a log
// model), out = input minus (additive) or divided by (log) the exp'd sum of
// the buckets, the bucket membership is the class predicate AND an EXACT
// regeffect-attribute match (a missing attribute is the Undefined bucket, so
// an outlier outside Trend/Seasonal/Irregular contributes to no table), the
// days-related and moving-holiday buckets ignore the regeffect entirely, the
// per-column accumulation is the pinned DataBlock.addAY rule (zero
// coefficient skips the column, one adds it undamped, otherwise x + a*y with
// plain multiply-then-add), the bucket sums fold left-associatively per row
// in the pinned order, and the restore uses the pinned truth table op
// (divide only for multiplicative and pseudo_additive), invOp (multiply for
// multiplicative, LOG-ADDITIVE and pseudo_additive — the pinned
// DecompositionMode.isMultiplicative() includes LogAdditive — add only for
// additive) and mean (0 only for additive), with the c17 correction rule
// y[i] exactly where c17[i] == 0 (X11Utility.correctSeries).
#include "arima_forecasts.hpp"
#include "x11_decomposition.hpp"
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
namespace kalos {
// jdplus.sa.base.api.ComponentType in pinned declaration order. The
// regeffect attribute of a model variable names exactly one of these.
enum class X13Component {
    undefined,
    series,
    trend,
    seasonal,
    seasonally_adjusted,
    irregular,
    calendar_effect
};
// The pinned ModellingUtility class predicates that partition the model
// variables (isOutlier, isDaysRelated, isMovingHoliday, isUser). The
// pinned capture proves the classes mutually exclusive for spec-built
// models; this contract takes them as given per variable.
enum class X13VariableClass { outlier, days_related, moving_holiday, user };
// One model variable: its pinned class, its EXACT regeffect bucket, and one
// coefficient per regression column (the values the pinned run used, fixed
// or estimated — this kernel consumes them as inputs). The columns are the
// pinned Regression.matrix columns over the input calendar.
struct X13AllocationVariable {
    X13VariableClass cls=X13VariableClass::user;
    X13Component effect=X13Component::undefined;
    std::vector<double> coefficient;
    std::vector<std::vector<double>> columns;
};
struct X13PreadjustInput {
    std::vector<double> input; // complete; strictly positive when log
    std::vector<X13AllocationVariable> variables; // pinned model order
    bool log=false;
    // The nonzero-horizon extension of item 21-H.4 (the pinned initialStep
    // branch): with backcast/forecast nonzero the variable columns span the
    // EXTENDED calendar (length n + backcast + forecast, the input rows at
    // [backcast, backcast + n)) — the pinned factory fills evaluated over it
    // (the zero-ended LS/SO patterns extend BACKWARD across the position, the
    // TC ramp runs into the forecast window) — and `arima`, the landing
    // stochastic component on the JVM coefficient surface with the mean
    // correction flag, drives the ExactArimaForecasts route on
    // r = transformed(input) - v27 over the input window.
    int backcast=0,forecast=0;
    std::optional<ArimaForecastModel> arima; // required when a horizon is nonzero
    // Optional factors on the extended monthly domain of a log/LeapYear model.
    // The model response uses the forward factors and the days-related effect,
    // combined deterministic effect and original-series forecasts use the reverse.
    std::vector<double> leap_forward, leap_reverse;
};
// The stored (backTransform'ed) preadjustment tables of X13Kernel.
// initialStep: a1 is the interpolated series (= the complete input) on the
// INPUT calendar, a6 the days-related effect, a7 the moving-holiday
// effect, a8t/a8s/a8i the trend/seasonal/irregular outlier+user buckets,
// a8 their pinned ((t+s)+i) fold, a9cal/a9sa/a9ser/a9u the user calendar,
// seasonally-adjusted, series and undefined buckets, a9 the pinned
// ((cal+sa)+ser) fold, and out the series fed to X-11. Under zero horizons
// all tables share the input length and calendar; under the item 21-H.4
// extension a6..a9u and out span the extended calendar (n + backcast +
// forecast rows), a1 stays the input, and a1b/a1a carry the ARIMA-extended
// raw series on the backcast/forecast windows (backTransform'ed, the pinned
// backcasts + v27 in the back window, forecasts + v27 in the forecast one).
struct X13Preadjustment {
    std::vector<double> a1,a6,a7,a8,a8t,a8s,a8i,a9,a9u,a9cal,a9sa,a9ser,out;
    std::vector<double> a1a,a1b; // empty under zero horizons
    std::string status; // "available", "unsupported:<reason>" or "failed:arima_forecasts"
};
// Fixed-coefficient preadjustment. Throws std::invalid_argument on shape
// violations (missing or non-positive-under-log observations, ragged or
// mismatched columns, non-finite coefficients, a nonzero horizon without a
// stochastic component); a nonzero horizon with a stochastic component runs
// the pinned ExactArimaForecasts route (whose own rejection answers
// "failed:arima_forecasts").
X13Preadjustment x13_preadjust(const X13PreadjustInput& input);

// The restore inputs: the available preadjustment above, the X-11 tables
// D10..D13 and C17 exactly as the pinned X11Results stores them for the same
// calendar (under the pinned LogAdditive results-builder stage — exp plus the
// default ratio bias correction, see the contract — D10..D13 are the exp'd
// and rescaled tables), the effective decomposition mode the pinned
// X13Kernel.updateSpec produced for the X-11 run (additive for non-log
// models, multiplicative for log models unless log_additive or
// pseudo_additive was requested — the forcing rule is part of the contract
// and the fixtures pin it), and — item 21-H.4 — the backcast/forecast row
// counts of the extended calendar the preadjustment and X-11 tables span
// (zero under the zero-horizon contract; the pinned finals then slice the
// a/b extension tables and fit the finals to the a1 calendar).
struct X13RestoreInput {
    X13Preadjustment pre;
    std::vector<double> d10,d11,d12,d13,c17;
    X11DecompositionMode mode=X11DecompositionMode::additive;
    int backcast=0,forecast=0;
};
// The pinned X13Finals: d11final (seasonally adjusted with the trend and
// irregular effects restored), d12final and d13final, d16 (original with the
// seasonal-side effects restored), d18 (calendar-effect-only composite), and
// the finals adjustment tables e1/e2/e3/e11 — all fitted to the a1 calendar.
// Under the item 21-H.4 extension the d11a/d12a/d16a/d18a tables carry the
// forecast-window slices (the pinned fitToDomain over the a1a domain) and
// d11b/... the backcast-window ones; they stay empty under zero horizons.
struct X13Finals {
    std::vector<double> d11final,d12final,d13final,d16,d18,e1,e2,e3,e11;
    std::vector<double> d11a,d12a,d16a,d18a,d11b,d12b,d16b,d18b;
    std::string status; // "available", or "unsupported:<reason>"
};
// Finals restore. Throws std::invalid_argument on shape violations (tables
// of differing lengths, an unavailable preadjustment, a calendar mismatch
// against the horizon counts).
X13Finals x13_restore(const X13RestoreInput& input);
} // namespace kalos
#endif
