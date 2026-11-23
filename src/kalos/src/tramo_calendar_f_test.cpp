#include "kalos/tramo_calendar_f_test.hpp"
#include "kalos/calendar_generation.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/tramo_regression_tests.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/hannan_rissanen.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

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
                     const std::vector<double>& filter, Columns columns, bool mean, double precision,
                     std::vector<double>* tstats = nullptr) {
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
    detail::PinnedRegressionDiagnostics diagnostics;
    const auto likelihood = filtered.empty()
        ? detail::pinned_deviances_evaluate(differences, 0, 0, 1, 1, 12, {}, ma)
        : detail::pinned_deviances_evaluate_multi(differences, 0, 0, 1, 1, 12,
              filtered, ma, detail::PinnedComputerShape::Fixed, tstats ? &diagnostics : nullptr);
    if (!likelihood || !(likelihood->ssq_err > 0.0))
        throw std::runtime_error("TRAMO calendar F-test candidate has no usable residual variance");
    if (tstats) {
        std::vector<double> diagonal(filtered.size());
        for (std::size_t j = 0; j < diagonal.size(); ++j)
            diagonal[j] = diagnostics.unscaled_covariance.at(j * diagonal.size() + j);
        *tstats = tramo_test_tstats({likelihood->b, diagonal, likelihood->ssq_err,
                                     static_cast<int>(differences.size()), 2}, true);
    }
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
TramoCalendarFTestResult compare_calendar_models(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision,
    std::vector<double>* selected_tstats) {
    const bool full = selected_tstats != nullptr;
    if (!(significance >= 0 && significance <= 1) || !(precision > 0) || !std::isfinite(precision))
        throw std::invalid_argument("TRAMO calendar F-test requires a probability and positive finite precision");
    // The largest candidate includes 13 differencing losses, six TD contrasts,
    // two MA coefficients, any mean, and in TRfull a leap-year and Easter term.
    if (series.size() <= static_cast<std::size_t>(21 + mean + (full ? 2 : 0)) ||
        series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("TRAMO calendar F-test requires positive residual degrees of freedom");
    const calendar_domain domain{start_year, start_month, 12, static_cast<int>(series.size())};
    const auto working = trading_days_contrast(day_clustering::td2, domain);
    const auto trading = trading_days_contrast(day_clustering::td7, domain);
    const auto filter = detail::pinned_differencing_polynomial(12, 1, 1);
    const auto differences = detail::pinned_backfilter_apply(filter, series);
    const auto easter = full ? easter_column(6, easter_correction::simple, 0, domain) : std::vector<double>{};
    auto estimate = [&](Columns columns, std::vector<double>* stats = nullptr) {
        if (full) columns.push_back(easter);
        return candidate_ssq(series, differences, filter, std::move(columns), mean, precision, stats);
    };
    const int df = static_cast<int>(differences.size()) - mean - full;
    const double base = estimate({}, selected_tstats);
    const auto one = compare(base, estimate(working), df-1, 1);
    const auto six = compare(base, estimate(trading), df-6, 6);
    const double threshold = 1.0 - significance;
    TramoTradingDays selected = TramoTradingDays::None;
    if (six.lower_probability > one.lower_probability && six.lower_probability > threshold)
        selected = TramoTradingDays::TradingDays;
    else if (one.lower_probability >= threshold)
        selected = TramoTradingDays::WorkingDays;
    if (full && selected != TramoTradingDays::None) {
        auto columns = selected == TramoTradingDays::WorkingDays ? working : trading;
        columns.push_back(length_of_period_column(calendar_lop::leap_year, domain));
        estimate(std::move(columns), selected_tstats);
    }
    return {base, one, six, selected};
}
} // namespace

TramoCalendarFTestResult tramo_calendar_f_test(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision) {
    return compare_calendar_models(series, start_year, start_month, mean, significance, precision, nullptr);
}

TramoFullCalendarResult tramo_full_initial_calendar(const std::vector<double>& series,
    int start_year, int start_month, bool log) {
    auto estimation_series = series;
    if (log) {
        for (auto& value : estimation_series) {
            if (!(value > 0.0)) throw std::invalid_argument("TRfull log transform requires positive observations");
            value = detail::fdlibm_log(value);
        }
    }
    std::vector<double> tstats;
    const auto comparison = compare_calendar_models(estimation_series, start_year, start_month,
                                                     true, 0.01, 1e-5, &tstats);
    // Test-model columns are mean, any TD contrasts, LP (when TD selected), Easter.
    const auto lp = comparison.selected == TramoTradingDays::None ? std::nullopt :
        std::optional<double>(tstats.at(tstats.size() - 2));
    return {comparison, tstats.front(), lp, tstats.back(), std::abs(tstats.front()) > 1.96,
            lp && std::abs(*lp) > 2.0, std::abs(tstats.back()) > 2.2};
}
} // namespace kalos
