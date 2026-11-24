#include "kalos/combined_seasonality.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/seats_extract_wk.hpp"
#include "kalos/x13_fit.hpp"

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

}  // namespace kalos
