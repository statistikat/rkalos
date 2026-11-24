#include "kalos/combined_seasonality.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/seats_extract_wk.hpp"
#include "kalos/x13_fit.hpp"
#include "kalos/tramo_transform_decision.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "mlj_engine.hpp"
#include "pinned_cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace kalos {
namespace {

struct FTest { double value, p; };

FTest anova(double model_ss, int model_df, double residual_ss, int residual_df) {
    const double f = (model_ss / model_df) * (residual_df / std::max(residual_ss, 0.0));
    return {f, detail::f_upper_p(f, model_df, residual_df)};
}

FTest stable_seasons(const std::vector<double>& si, int period) {
    const int n = static_cast<int>(si.size());
    const double mean = std::accumulate(si.begin(), si.end(), 0.0) / n;
    double total = 0, between = 0;
    for (double value : si) total += (value - mean) * (value - mean);
    for (int position = 0; position < period; ++position) {
        double sum = 0;
        int count = 0;
        for (std::size_t j = position; j < si.size(); j += period) {
            sum += si[j];
            ++count;
        }
        const double delta = sum / count - mean;
        between += delta * delta * count;
    }
    return anova(between, period - 1, total - between, n - period);
}

FTest moving_seasons(const std::vector<double>& si, int period, int first, double center) {
    const int skip = (period - (first - 1)) % period;
    const int years = (static_cast<int>(si.size()) - skip) / period;
    std::vector<double> seasons(period), annual(years);
    double mean = 0;
    for (int y = 0; y < years; ++y) {
        for (int s = 0; s < period; ++s) {
            const double value = std::abs(si[skip + y * period + s] - center);
            seasons[s] += value;
            annual[y] += value;
        }
        mean += annual[y];
    }
    mean /= years * period;
    double total = 0, between_years = 0, between_seasons = 0;
    for (int t = 0; t < years * period; ++t) {
        const double delta = std::abs(si[skip + t] - center) - mean;
        total += delta * delta;
    }
    for (double sum : annual) {
        const double delta = sum / period - mean;
        between_years += delta * delta;
    }
    for (double sum : seasons) {
        const double delta = sum / years - mean;
        between_seasons += delta * delta;
    }
    between_years *= period;
    between_seasons *= years;
    return anova(between_years, years - 1, total - between_seasons - between_years, (years - 1) * (period - 1));
}

// Rank denominators and offsets in the public statistic are JVM signed int products.
double signed_int_product(std::uint32_t a, std::uint32_t b) {
    const std::uint32_t bits = a * b;
    return bits <= 0x7fffffffu ? static_cast<double>(bits) : static_cast<double>(bits) - 4294967296.0;
}

// Kruskal-Wallis with mean ranks for ties and the public uncorrected statistic.
double seasonal_ranks(const std::vector<double>& si, int period) {
    std::vector<std::size_t> order(si.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](auto i, auto j) { return si[i] < si[j]; });
    std::vector<double> ranks(period);
    std::vector<int> counts(period);
    for (std::size_t i = 0; i < order.size();) {
        std::size_t end = i + 1;
        while (end < order.size() && si[order[end]] == si[order[i]]) ++end;
        const double rank = double(i) + 0.5 * double(end - i + 1);
        for (; i < end; ++i) {
            const auto season = order[i] % period;
            ranks[season] += rank;
            ++counts[season];
        }
    }
    double sum = 0;
    for (int j = 0; j < period; ++j) sum += ranks[j] * ranks[j] / counts[j];
    const auto n = static_cast<std::uint32_t>(si.size());
    return 12 * sum / signed_int_product(n, n + 1u) - signed_int_product(3u, n + 1u);
}

}  // namespace

std::optional<CombinedSiSeasonality> combined_si_seasonality(const std::vector<double>& si, int frequency,
                                                             int first_period, bool multiplicative) {
    if ((frequency != 4 && frequency != 12) || first_period < 1 || first_period > frequency ||
        si.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("combined SI seasonality requires a valid monthly or quarterly domain");
    if (!std::all_of(si.begin(), si.end(), [](double x) { return std::isfinite(x); }))
        throw std::invalid_argument("combined SI seasonality requires complete finite observations");
    if (si.size() < static_cast<std::size_t>(3 * frequency)) return std::nullopt;
    const auto stable = stable_seasons(si, frequency);
    const auto moving = moving_seasons(si, frequency, first_period, multiplicative ? 1.0 : 0.0);
    const double rank = seasonal_ranks(si, frequency);
    const double rank_p = detail::chi2_upper_p(frequency - 1, rank);
    const double t1 = std::min(7 / stable.value, 9.0);
    const double t2 = stable.value == 0 ? 9.0 : std::min(3 * moving.value / stable.value, 9.0);
    std::string verdict = "Present";
    if (stable.p >= 0.001 || (moving.p < 0.05 && std::sqrt((t1 + t2) / 2) >= 1)) verdict = "None";
    else if (t1 >= 1 || t2 >= 1 || rank_p >= 0.001) verdict = "ProbablyNone";
    return CombinedSiSeasonality{verdict, stable.value, stable.p, moving.value, moving.p, rank, rank_p};
}

std::optional<CombinedSiSeasonality> seats_si_seasonality(const std::vector<seats_component_series>& wk_final,
                                                          int frequency, bool multiplicative) {
    if (frequency != 4 && frequency != 12) throw std::invalid_argument("invalid WK seasonality frequency");
    const seats_component_series* seasonal = nullptr;
    const seats_component_series* irregular = nullptr;
    for (const auto& row : wk_final) {
        if (row.information != seats_information::value) continue;
        if (row.component == seats_component::seasonal) seasonal = &row;
        if (row.component == seats_component::irregular) irregular = &row;
    }
    const auto* domain = seasonal ? seasonal : irregular;
    if (!domain) return std::nullopt;
    const int width = 12 / frequency;
    if (domain->start_month < 1 || domain->start_month > 12 || (domain->start_month - 1) % width)
        throw std::invalid_argument("invalid WK SI start month");
    auto values = domain->values;
    if (seasonal && irregular) {
        if (seasonal->start_year != irregular->start_year || seasonal->start_month != irregular->start_month ||
            values.size() != irregular->values.size()) throw std::invalid_argument("WK SI domains disagree");
        for (std::size_t i = 0; i < values.size(); ++i)
            values[i] = multiplicative ? values[i] * irregular->values[i] : values[i] + irregular->values[i];
    }
    return combined_si_seasonality(values, frequency, 1 + (domain->start_month - 1) / width, multiplicative);
}

std::optional<CombinedSiSeasonality> x13_si_seasonality(const X13FitResult& fit) {
    if (fit.status != "available") return std::nullopt;
    const auto& d8 = fit.selection.d.d8;
    if (fit.backcast < 0 || fit.forecast < 0 || d8.missing.size() != d8.values.size() ||
        static_cast<std::size_t>(fit.backcast) + static_cast<std::size_t>(fit.forecast) > d8.values.size())
        throw std::invalid_argument("invalid extended X13 D8 domain");
    const auto n = d8.values.size() - fit.forecast - fit.backcast;
    std::vector<double> si(n);
    const bool log_additive = fit.mode_effective == X11DecompositionMode::log_additive;
    for (std::size_t j = 0; j < n; ++j) {
        const auto t = j + fit.backcast;
        if (d8.missing[t]) return std::nullopt;
        si[j] = log_additive ? detail::fdlibm_exp(d8.values[t]) : d8.values[t];
    }
    // The public decomposition mode classifies LogAdditive as multiplicative
    // and PseudoAdditive as additive for this test's center.
    return combined_si_seasonality(si, fit.frequency, fit.first_period,
                                   log_additive || fit.mode_effective == X11DecompositionMode::multiplicative);
}

AdjustedSeasonality monthly_adjusted_seasonality(const std::vector<double>& series,
    int differencing_order, bool mean) {
    if (differencing_order < 0 || differencing_order > 3 ||
        !std::all_of(series.begin(), series.end(), [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument("adjusted seasonality requires finite observations and a fitted differencing order");
    if (series.size() <= static_cast<std::size_t>(differencing_order)) return {};

    auto data = series;
    for (int order = 0; order < differencing_order; ++order) {
        for (std::size_t i = 0; i + 1 < data.size(); ++i) data[i] = data[i+1] - data[i];
        data.pop_back();
    }
    if (mean) {
        const double average = std::accumulate(data.begin(), data.end(), 0.) / data.size();
        for (double& value : data) value -= average;
    }

    AdjustedSeasonality result;
    const auto qs = monthly_qs_test(data);
    if (qs.valid) result.qs = AdjustedSeasonalityTest{qs.value, qs.p_value};
    // FTest.processAr regresses on a constant, one AR lag and eleven
    // periodic contrasts. Phase zero starts at the first response row.
    constexpr int months = 12, contrasts = months - 1, columns = contrasts + 2;
    const auto rows = data.size() - 1;
    if (rows < columns) return result;
    const std::vector<double> response(data.begin()+1, data.end());
    std::vector<std::vector<double>> design(columns, std::vector<double>(rows, 0.));
    std::fill(design[0].begin(), design[0].end(), 1.);
    std::copy_n(data.begin(), rows, design[1].begin());
    for (std::size_t i = 0; i < rows; ++i) {
        const int phase = static_cast<int>(i % months);
        for (int j = 0; j < contrasts; ++j)
            design[2+j][i] = phase == j ? 1. : phase == contrasts ? -1. : 0.;
    }
    std::vector<double> covariance;
    auto fit = detail::pinned_ols_fit(response, design, &covariance);
    if (!fit || fit->rank != columns) return result;
    if (rows == columns) {
        const double undefined = std::numeric_limits<double>::quiet_NaN();
        result.seasonal_f = AdjustedSeasonalityTest{undefined, undefined};
        return result;
    }
    std::vector<double> block(contrasts * contrasts);
    for (int col = 0; col < contrasts; ++col)
        for (int row = 0; row < contrasts; ++row)
            block[col * contrasts + row] = covariance[(row + 2) * columns + col + 2];
    if (!detail::pinned_lcholesky(block, contrasts)) return result;
    std::vector<double> beta(fit->coefficients.begin()+2, fit->coefficients.end());
    mlj::solve_lx(block, contrasts, beta);
    const double df = static_cast<double>(rows - columns);
    const double f = (mlj::ssq(beta.data(), beta.size()) / contrasts) / (fit->residual_ssq / df);
    result.seasonal_f = AdjustedSeasonalityTest{f, detail::f_upper_p(f, contrasts, df)};
    return result;
}

AdjustedSeasonality seats_adjusted_seasonality(const std::vector<seats_component_series>& wk_initial,
    int differencing_order, bool mean) {
    for (const auto& item : wk_initial) {
        if (item.component == seats_component::seasonally_adjusted && item.information == seats_information::value)
            return monthly_adjusted_seasonality(item.values, differencing_order, mean);
    }
    return {};
}

AdjustedSeasonality x13_adjusted_seasonality(const X13FitResult& fit, int differencing_order, bool mean) {
    if (fit.status != "available" || fit.frequency != 12) return {};
    if (fit.backcast < 0 || fit.forecast < 0 ||
        static_cast<std::size_t>(fit.backcast) + static_cast<std::size_t>(fit.forecast) > fit.d11.size())
        throw std::invalid_argument("invalid extended X-13 D11 domain");
    std::vector<double> observed(fit.d11.begin() + fit.backcast, fit.d11.end() - fit.forecast);
    if (fit.mode_effective == X11DecompositionMode::multiplicative ||
        fit.mode_effective == X11DecompositionMode::log_additive) {
        for (auto& value : observed) {
            if (!(value > 0.)) return {};
            value = detail::fdlibm_log(value);
        }
    }
    return monthly_adjusted_seasonality(observed, differencing_order, mean);
}

}  // namespace kalos
