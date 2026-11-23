#include "kalos/tramo_calendar_f_test.hpp"
#include "kalos/calendar_generation.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/hannan_rissanen.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace kalos {
namespace {
using Columns = std::vector<std::vector<double>>;

RegressionDesign make_design(const Columns& columns) {
    RegressionDesign design{};
    design.columns = columns.size();
    design.fixed.resize(columns.size());
    design.coefficients.resize(columns.size());
    for (const auto& column : columns)
        design.values.insert(design.values.end(), column.begin(), column.end());
    return design;
}

// RegArimaUtility.processor(true, precision): Hannan–Rissanen on the
// differenced linear model, then one warm LM pass. The final likelihood uses
// DEFAULT_COMPUTER, whose QR layout differs from the optimization objective.
double candidate_ssq(const std::vector<double>& series, const std::vector<double>& differences,
                     const std::vector<double>& filter, Columns columns, bool mean, double precision) {
    Columns filtered;
    if (mean) filtered.emplace_back(differences.size(), 1.0);
    for (const auto& column : columns)
        filtered.push_back(detail::pinned_backfilter_apply(filter, column));
    const auto hr = hannan_rissanen_start(differences, make_design(filtered), {0, 0, 1, 1, 12});

    SarimaEstimationRequest request{};
    request.series.values = series;
    request.series.missing.resize(series.size());
    request.period = 12;
    request.d = request.seasonal_d = 1;
    if (mean) columns.insert(columns.begin(), std::vector<double>(series.size(), 1.0));
    request.design = make_design(columns);
    for (const int block : {1, 3}) {
        request.blocks[block].initial.resize(1);
        request.blocks[block].fixed.resize(1);
    }
    request.control = {precision, 1e-4, true, false, true, true,
        EstimationStartingPoint::HannanRissanen, EstimationSolver::PinnedLevenbergMarquardt};
    std::array<std::vector<double>, 4> start{{{}, hr.coefficients[2], {}, hr.coefficients[3]}};
    const auto fit = optimize_sarima_from(request, start);
    if (!std::isfinite(fit.log_likelihood))
        throw std::runtime_error("TRAMO calendar F-test candidate estimation failed");
    const std::vector<double> ma{fit.coefficients[1].front(), fit.coefficients[3].front()};
    const auto likelihood = filtered.empty()
        ? detail::pinned_deviances_evaluate(differences, 0, 0, 1, 1, 12, {}, ma)
        : detail::pinned_deviances_evaluate_multi(differences, 0, 0, 1, 1, 12,
              filtered, ma, detail::PinnedComputerShape::Fixed);
    if (!likelihood || !(likelihood->ssq_err > 0.0))
        throw std::runtime_error("TRAMO calendar F-test candidate has no usable residual variance");
    return likelihood->ssq_err;
}

TramoCalendarFCandidate compare(double base, double residual_ssq, int df, int contrasts) {
    const int denominator_df = df - 2;
    const double mean_square = residual_ssq / denominator_df;
    const double f = (base - residual_ssq) / (mean_square * contrasts);
    const double probability = f <= 0.0 ? 0.0 : 1.0 - detail::f_upper_p(f, contrasts, denominator_df);
    if (!std::isfinite(probability))
        throw std::runtime_error("TRAMO calendar F-test probability is undefined");
    return {residual_ssq, df, mean_square, probability};
}
} // namespace

TramoCalendarFTestResult tramo_calendar_f_test(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision) {
    if (!(significance >= 0 && significance <= 1) || !(precision > 0) || !std::isfinite(precision))
        throw std::invalid_argument("TRAMO calendar F-test requires a probability and positive finite precision");
    // Thirteen differences, six TD contrasts, an optional mean, two fitted MA coefficients.
    if (series.size() <= static_cast<std::size_t>(21 + mean) ||
        series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("TRAMO calendar F-test requires positive residual degrees of freedom");
    const calendar_domain domain{start_year, start_month, 12, static_cast<int>(series.size())};
    const auto working = trading_days_contrast(day_clustering::td2, domain);
    const auto trading = trading_days_contrast(day_clustering::td7, domain);
    const auto filter = detail::pinned_differencing_polynomial(12, 1, 1);
    const auto differences = detail::pinned_backfilter_apply(filter, series);
    const int df = static_cast<int>(differences.size()) - mean;
    const double base = candidate_ssq(series, differences, filter, {}, mean, precision);
    const auto one = compare(base, candidate_ssq(series, differences, filter, working, mean, precision), df-1, 1);
    const auto six = compare(base, candidate_ssq(series, differences, filter, trading, mean, precision), df-6, 6);
    const double threshold = 1.0 - significance;
    TramoTradingDays selected = TramoTradingDays::None;
    if (six.lower_probability > one.lower_probability && six.lower_probability > threshold)
        selected = TramoTradingDays::TradingDays;
    else if (one.lower_probability >= threshold)
        selected = TramoTradingDays::WorkingDays;
    return {base, one, six, selected};
}
} // namespace kalos
