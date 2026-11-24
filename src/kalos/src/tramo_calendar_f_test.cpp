#include "kalos/tramo_calendar_f_test.hpp"
#include "kalos/calendar_generation.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/tramo_regression_tests.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/hannan_rissanen.hpp"
#include "kalos/outlier_scoring.hpp"
#include "kalos/tramo_outlier_rounds.hpp"

#include <algorithm>
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

struct CalendarLikelihood {
    std::vector<double> ma;
    detail::PinnedMachinePoint point;
    detail::PinnedRegressionDiagnostics regression;
    std::vector<double> covariance_diagonal;
    int dimension;

    std::vector<double> tstats() const {
        return tramo_test_tstats({point.b, covariance_diagonal, point.ssq_err, dimension, static_cast<int>(ma.size())}, true);
    }
};

// RegArimaUtility.processor(true, precision): Hannan–Rissanen on the
// differenced linear model, then one warm LM pass. The final likelihood uses
// DEFAULT_COMPUTER, whose QR layout differs from the optimization objective.
CalendarLikelihood fit_calendar_likelihood(const std::vector<double>& series,
    const std::vector<double>& differences, const std::vector<double>& filter, Columns columns,
    bool mean, double precision, bool collect_diagnostics, bool seasonal = true) {
    Columns filtered;
    if (mean) filtered.emplace_back(differences.size(), 1.0);
    for (const auto& column : columns)
        filtered.push_back(detail::pinned_backfilter_apply(filter, column));
    const auto hr = hannan_rissanen_start(differences, make_design(filtered), {0, 0, 1, seasonal ? 1 : 0, 12});

    SarimaEstimationRequest request{};
    request.series.values = series;
    request.series.missing.resize(series.size());
    request.period = 12;
    request.d = 1;
    request.seasonal_d = seasonal;
    if (mean) columns.insert(columns.begin(), std::vector<double>(series.size(), 1.0));
    request.design = make_design(columns);
    for (int block = 1; block <= (seasonal ? 3 : 1); block += 2) {
        request.blocks[block].initial.resize(1);
        request.blocks[block].fixed.resize(1);
    }
    request.control = {precision, 1e-4, true, false, true, true,
        EstimationStartingPoint::HannanRissanen, EstimationSolver::PinnedLevenbergMarquardt};
    std::array<std::vector<double>, 4> start{{{}, hr.coefficients[2], {}, hr.coefficients[3]}};
    const auto fit = optimize_sarima_from(request, start);
    if (!std::isfinite(fit.log_likelihood))
        throw std::runtime_error("TRAMO calendar F-test candidate estimation failed");
    auto ma = fit.coefficients[1];
    ma.insert(ma.end(), fit.coefficients[3].begin(), fit.coefficients[3].end());
    detail::PinnedRegressionDiagnostics diagnostics;
    auto likelihood = filtered.empty() && !collect_diagnostics
        ? detail::pinned_deviances_evaluate(differences, 0, 0, 1, seasonal, 12, {}, ma)
        : detail::pinned_deviances_evaluate_multi(differences, 0, 0, 1, seasonal, 12,
            filtered, ma, detail::PinnedComputerShape::Fixed, collect_diagnostics ? &diagnostics : nullptr);
    if (!likelihood || !(likelihood->ssq_err > 0.0))
        throw std::runtime_error("TRAMO calendar candidate has no usable residual variance");
    std::vector<double> diagonal;
    if (collect_diagnostics) {
        const std::size_t columns_count = filtered.size();
        for (std::size_t j = 0; j < columns_count; ++j)
            diagonal.push_back(diagnostics.unscaled_covariance.at(j * columns_count + j));
    }
    return {std::move(ma), std::move(*likelihood), std::move(diagnostics), std::move(diagonal),
            static_cast<int>(differences.size())};
}

TramoCalendarAirlineFit controller_fit(CalendarLikelihood fit, int observations, bool mean, bool working_days,
    std::vector<std::array<int, 2>> outliers) {
    const auto t = fit.tstats();
    TramoModelStatisticsRequest stats;
    stats.residuals = std::move(fit.regression.qr_residuals);
    stats.coefficients = std::move(fit.point.b);
    stats.cov_diagonal = std::move(fit.covariance_diagonal);
    stats.ssq = fit.point.ssq_err;
    stats.dim = fit.dimension;
    stats.degrees_of_freedom = stats.dim - static_cast<int>(stats.coefficients.size());
    stats.sigma2 = stats.ssq / stats.dim;
    stats.free_parameters = 2;
    stats.estimation_domain = observations;
    stats.ami_outliers = static_cast<int>(outliers.size());
    const auto panel = tramo_model_comparison_panel(stats);
    if (!panel) throw std::runtime_error("TRAMO calendar candidate has no comparison panel");
    return {mean, working_days, std::move(fit.ma), std::move(stats.coefficients), std::move(stats.residuals),
            stats.ssq, mean ? std::optional<double>(t.front()) : std::nullopt, *panel, std::move(outliers)};
}

TramoCalendarFCandidate compare(double base, double residual_ssq, int df, int contrasts, int fitted_parameters) {
    const int denominator_df = df - fitted_parameters;
    const double mean_square = residual_ssq / denominator_df;
    const double f = (base - residual_ssq) / (mean_square * contrasts);
    const double probability = f <= 0.0 ? 0.0 : 1.0 - detail::f_upper_p(f, contrasts, denominator_df);
    if (!std::isfinite(probability))
        throw std::runtime_error("TRAMO calendar F-test probability is undefined");
    return {residual_ssq, df, mean_square, probability};
}
TramoCalendarFTestResult compare_calendar_models(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision,
    std::vector<double>* selected_tstats, bool seasonal = true) {
    const bool full = selected_tstats != nullptr;
    if (!(significance >= 0 && significance <= 1) || !(precision > 0) || !std::isfinite(precision))
        throw std::invalid_argument("TRAMO calendar F-test requires a probability and positive finite precision");
    // The largest candidate includes 13 differencing losses, six TD contrasts,
    // two MA coefficients, any mean, and in TRfull a leap-year and Easter term.
    if (series.size() <= static_cast<std::size_t>(8 + 13 * seasonal + mean + (full ? 2 : 0)) ||
        series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("TRAMO calendar F-test requires positive residual degrees of freedom");
    const calendar_domain domain{start_year, start_month, 12, static_cast<int>(series.size())};
    const auto working = trading_days_contrast(day_clustering::td2, domain);
    const auto trading = trading_days_contrast(day_clustering::td7, domain);
    const auto filter = detail::pinned_differencing_polynomial(12, 1, seasonal);
    const auto differences = detail::pinned_backfilter_apply(filter, series);
    const auto easter = full ? easter_column(6, easter_correction::simple, 0, domain) : std::vector<double>{};
    auto estimate = [&](Columns columns, std::vector<double>* stats = nullptr) {
        if (full) columns.push_back(easter);
        const auto fit = fit_calendar_likelihood(series, differences, filter, std::move(columns),
                                                mean, precision, stats != nullptr, seasonal);
        if (stats) *stats = fit.tstats();
        return fit.point.ssq_err;
    };
    const int df = static_cast<int>(differences.size()) - mean - full;
    const double base = estimate({}, selected_tstats);
    const auto one = compare(base, estimate(working), df-1, 1, 1 + seasonal);
    const auto six = compare(base, estimate(trading), df-6, 6, 1 + seasonal);
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

std::vector<std::vector<double>> tramo_stock_calendar_columns(const TramoStockCalendar& selection,
    const calendar_domain& domain) {
    Columns columns;
    switch (selection.trading_days) {
    case TramoTradingDays::None: break;
    case TramoTradingDays::WorkingDays: columns = trading_days_contrast(day_clustering::td2, domain); break;
    case TramoTradingDays::TradingDays: columns = trading_days_contrast(day_clustering::td7, domain); break;
    default: throw std::invalid_argument("unknown TRAMO trading-days selection");
    }
    if (selection.leap_year) columns.push_back(length_of_period_column(calendar_lop::leap_year, domain));
    if (selection.easter) columns.push_back(easter_column(6, easter_correction::simple, 0, domain));
    return columns;
}

std::vector<double> tramo_stock_calendar_effect(const TramoFittedCalendar& fit, const calendar_domain& domain) {
    const auto columns = tramo_stock_calendar_columns(fit.selection, domain);
    if (columns.size() != fit.coefficients.size() || domain.length < 0)
        throw std::invalid_argument("TRAMO calendar coefficient count or domain is invalid");
    std::vector<double> effect(domain.length, 0.0);
    for (std::size_t j = 0; j < columns.size(); ++j) {
        const double coefficient = fit.coefficients[j];
        if (!std::isfinite(coefficient)) throw std::invalid_argument("non-finite calendar coefficient");
        for (std::size_t i = 0; i < effect.size(); ++i)
            effect[i] += coefficient * columns[j][i];
    }
    return effect;
}

TramoCalendarFTestResult tramo_calendar_f_test(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision) {
    return compare_calendar_models(series, start_year, start_month, mean, significance, precision, nullptr);
}

TramoFullCalendarResult tramo_full_initial_calendar(const std::vector<double>& series,
    int start_year, int start_month, bool log) {
    return tramo_full_initial_calendar(series, start_year, start_month, log, true);
}

TramoFullCalendarResult tramo_full_initial_calendar(const std::vector<double>& series,
    int start_year, int start_month, bool log, bool seasonal) {
    auto estimation_series = series;
    if (log) {
        for (auto& value : estimation_series) {
            if (!(value > 0.0)) throw std::invalid_argument("TRfull log transform requires positive observations");
            value = detail::fdlibm_log(value);
        }
    }
    std::vector<double> tstats;
    const auto comparison = compare_calendar_models(estimation_series, start_year, start_month,
                                                     true, 0.01, 1e-5, &tstats, seasonal);
    // Test-model columns are mean, any TD contrasts, LP (when TD selected), Easter.
    const auto lp = comparison.selected == TramoTradingDays::None ? std::nullopt :
        std::optional<double>(tstats.at(tstats.size() - 2));
    return {comparison, tstats.front(), lp, tstats.back(), std::abs(tstats.front()) > 1.96,
            lp && std::abs(*lp) > 2.0, std::abs(tstats.back()) > 2.2};
}
std::optional<TramoResidualTradingDays> tramo_residual_trading_days(
    const std::vector<double>& residuals, int estimation_start_year, int estimation_start_month,
    int estimation_length, double significance) {
    if (estimation_start_month < 1 || estimation_start_month > 12 ||
        estimation_start_year < 1 || estimation_start_year > 9999 ||
        estimation_length < 0 || residuals.size() > static_cast<std::size_t>(estimation_length))
        throw std::invalid_argument("TRAMO residual trading-day test requires a valid monthly estimation domain");
    if (!(significance >= 0 && significance <= 1))
        throw std::invalid_argument("TRAMO residual trading-day significance must be in [0,1]");
    if (residuals.size() < 2) return std::nullopt;
    const long long residual_start = static_cast<long long>(estimation_start_year) * 12 +
        estimation_start_month - 1 + estimation_length - static_cast<int>(residuals.size());
    const calendar_domain domain{static_cast<int>(residual_start / 12),
        static_cast<int>(residual_start % 12) + 1, 12, static_cast<int>(residuals.size())};
    const auto columns = trading_days_contrast(day_clustering::td2, domain);
    double total_ssq = 0.0;
    for (const double value : residuals) {
        if (!std::isfinite(value))
            throw std::invalid_argument("TRAMO residual trading-day test requires finite residuals");
        total_ssq += value * value;
    }
    const auto fit = detail::pinned_ols_fit(residuals, columns);
    if (!fit || fit->rank != 1) return std::nullopt;
    const double error_ssq = std::min(fit->residual_ssq, total_ssq);
    const double error_df = static_cast<double>(residuals.size() - 1);
    const double statistic = (total_ssq - error_ssq) / (error_ssq / error_df);
    const double p = detail::f_upper_p(statistic, 1.0, error_df);
    return TramoResidualTradingDays{statistic, p, p < significance};
}

TramoTradingDaysControllerResult tramo_airline_trading_days_controller(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision,
    std::optional<double> outlier_critical_value) {
    if (!(significance >= 0 && significance <= 1) || !(precision > 0) || !std::isfinite(precision))
        throw std::invalid_argument("TRAMO trading-day controller requires a probability and positive precision");
    if (series.size() < 20 || series.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("TRAMO airline controller requires at least 20 monthly observations");
    if (outlier_critical_value && (!std::isfinite(*outlier_critical_value) || *outlier_critical_value <= 0.0))
        throw std::invalid_argument("TRAMO outlier detection requires a positive finite critical value");
    const calendar_domain domain{start_year, start_month, 12, static_cast<int>(series.size())};
    const auto td = trading_days_contrast(day_clustering::td2, domain);
    const auto filter = detail::pinned_differencing_polynomial(12, 1, 1);
    const auto differences = detail::pinned_backfilter_apply(filter, series);
    auto estimate = [&](bool fit_mean, Columns columns) {
        const bool working_days = !columns.empty();
        std::vector<std::array<int, 2>> outliers;
        if (outlier_critical_value) {
            TramoOutlierScanRequest scan{};
            scan.y = series; // already on the estimation scale
            scan.frequency = 12;
            scan.mean = fit_mean;
            scan.d = scan.bd = scan.q = scan.bq = 1;
            scan.ao = scan.ls = scan.tc = true;
            scan.tc_rate = 0.7;
            scan.precision = precision;
            scan.cv = *outlier_critical_value;
            scan.max_round = 30;
            scan.max_outliers = 50;
            scan.upper_bound = domain.length;
            scan.calendar_regressors = columns;
            const auto detected = tramo_outlier_rounds(scan);
            // ModelEstimator runs calc() even when the detector reports Failed.
            if (detected.status == TramoOutlierRoundsResult::Status::Ok) outliers = detected.outliers;
            std::stable_sort(outliers.begin(), outliers.end(),
                [](const auto& a, const auto& b) { return a[0] < b[0]; });
            for (const auto& found : outliers) {
                columns.emplace_back(series.size(), 0.0);
                fill_outlier_column(static_cast<OutlierType>(found[1]), found[0], 12, scan.tc_rate, columns.back());
            }
        }
        auto fit = fit_calendar_likelihood(series, differences, filter, columns, fit_mean, precision, true);
        return controller_fit(std::move(fit), domain.length, fit_mean, working_days, std::move(outliers));
    };
    TramoTradingDaysControllerResult result{};
    result.original = estimate(mean, {});
    result.residual_test = tramo_residual_trading_days(result.original.residuals, start_year,
        start_month, domain.length, significance);
    if (!result.residual_test || !result.residual_test->needs_processing) return result;
    result.candidates.push_back(estimate(mean, td));
    const auto t = result.candidates.front().mean_t;
    if (t && !(std::abs(*t) >= 1.96))
        result.candidates.push_back(estimate(false, td));
    result.changed = tramo_candidate_preferred(result.original.statistics, result.candidates.back().statistics);
    return result;
}
TramoTradingDaysControllerResult tramo_airline_trading_days_controller(const std::vector<double>& series,
    int start_year, int start_month, bool mean, double significance, double precision) {
    return tramo_airline_trading_days_controller(series, start_year, start_month, mean,
        significance, precision, std::nullopt);
}
} // namespace kalos
