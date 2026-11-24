// Monthly spectral-peak classifications used by TRAMO SeasonalityTests.
#include "kalos/tramo_differencing_battery.hpp"
#include "kalos/tramo_transform_decision.hpp"
#include "kalos/detail/hostlineage.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/detail/special_functions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace kalos {
namespace {

constexpr double pi = 0x1.921fb54442d18p+1;

std::vector<double> tapered_autocovariances(const std::vector<double>& values, int window) {
    const double n = static_cast<double>(values.size());
    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
    std::vector<double> centered(values);
    for (double& x : centered) x -= mean;

    std::vector<double> covariance(window);
    for (int lag = 0; lag < window; ++lag) {
        double product = 0.0;
        for (std::size_t i = lag; i < centered.size(); ++i) product += centered[i - lag] * centered[i];
        covariance[lag] = product / n;
    }
    const double variance = covariance[0];
    covariance[0] = 1.0;
    for (int lag = 1; lag < window; ++lag) {
        const double weight = 0.5 * (1.0 + detail::hotspot_spectrum_cos(pi * (static_cast<double>(lag) / window)));
        covariance[lag] *= weight / variance;
    }
    return covariance;
}

std::array<double, 4> tukey_degrees_of_freedom(std::size_t n, int window) {
    if (window == 79) return {6.35251, 19.6308, 2.29316, 6.55412};
    // Monthly long-window calibration of the two interior and the two
    // Nyquist F-test degrees of freedom, as functions of sample length.
    constexpr double coefficients[4][3] = {
        {0.54630, 2.93030, 2.2042}, {1.1329, 7.6924, 10.8795},
        {-0.3492, 1.533, 2.7696}, {0.9829, 3.8217, 6.9345}
    };
    std::array<double, 4> df{};
    for (std::size_t k = 0; k < df.size(); ++k) {
        const auto& c = coefficients[k];
        df[k] = c[0] + c[1] * (n / 100.0) + c[2] * (100.0 / n);
    }
    return df;
}

} // namespace

std::vector<PeakTukey> tramo_monthly_tukey_peaks(const std::vector<double>& values) {
    if (!std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); }))
        throw std::invalid_argument("monthly Tukey spectrum requires finite data");
    if (values.size() < 80) return {};
    const int window = values.size() >= 120 ? 112 : 79;
    const auto covariance = tapered_autocovariances(values, window);
    const auto density = [&](int index) {
        const double omega = index * (2.0 * pi) / window;
        double sum = 1.0;
        for (int lag = 1; lag < window; ++lag)
            sum += 2.0 * covariance[lag] * detail::hotspot_spectrum_cos(lag * omega);
        return sum < 0.0 ? 0.0 : sum;
    };
    const auto df = tukey_degrees_of_freedom(values.size(), window);
    const auto classify = [](double ratio, double numerator_df, double denominator_df) {
        const double probability = ratio > 0.0 ? 1.0 - detail::f_upper_p(ratio, numerator_df, denominator_df) : 0.0;
        if (probability > 0.99) return PeakTukey::T;
        if (probability > 0.9) return PeakTukey::t;
        return PeakTukey::none;
    };
    const std::array<int, 5> frequencies = window == 112 ? std::array<int, 5>{9, 19, 28, 37, 47}
                                                        : std::array<int, 5>{7, 13, 20, 26, 33};
    std::vector<PeakTukey> peaks;
    for (int f : frequencies) {
        const double ratio = 2.0 * density(f) / (density(f + 1) + density(f - 1));
        peaks.push_back(classify(ratio, df[0], df[1]));
    }
    const int nyquist = window / 2;
    const double below = density(nyquist - 1), at = density(nyquist);
    peaks.push_back(classify(below == 0.0 ? at : at / below, df[2], df[3]));
    return peaks;
}

std::vector<PeakAr> tramo_monthly_ar_peaks(const std::vector<double>& values) {
    if (!std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); }))
        throw std::invalid_argument("monthly AR spectrum requires finite data");
    constexpr int lags = 30; // monthly AutoRegressiveSpectrumTest default
    if (values.size() < 2 * lags) return {};
    const std::size_t n = values.size();
    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
    std::vector<double> centered(values);
    for (double& x : centered) x -= mean;
    const std::vector<double> response(centered.begin() + lags, centered.end());
    std::vector<std::vector<double>> columns;
    for (int lag = 1; lag <= lags; ++lag)
        columns.emplace_back(centered.begin() + lags - lag, centered.end() - lag);
    const auto fit = detail::pinned_fast_ols_fit(response, columns);
    if (!fit) return {};
    const double variance = fit->residual_ssq / response.size();

    constexpr int intervals = 60;
    constexpr double step = pi / intervals;
    std::array<double, intervals + 1> frequency;
    for (int i = 0; i <= intervals; ++i) frequency[i] = step * i;
    const double month = 365.25 / 12;
    const double trading_day = (2.0 * pi / 7) * (month - 7.0 * std::floor(month / 7));
    const int offset = static_cast<int>(trading_day / step + 0.5);
    for (int i = -1; i <= 1; ++i) frequency[offset + i] = trading_day + i * step;
    std::array<double, intervals + 1> spectrum;
    for (int i = 0; i <= intervals; ++i) {
        double real = 1.0, imaginary = 0.0;
        for (int j = 0; j < lags; ++j) {
            const double angle = (j + 1) * frequency[i];
            real -= fit->coefficients[j] * detail::hotspot_spectrum_cos(angle);
            imaginary -= fit->coefficients[j] * std::sin(angle);
        }
        const double power = variance / (real * real + imaginary * imaginary);
        spectrum[i] = power <= 0.0 ? 0.0 : 10.0 * std::log10(power);
    }
    auto ordered = spectrum;
    std::sort(ordered.begin(), ordered.end());
    const double median = ordered[intervals / 2], range = ordered.back() - ordered.front();
    std::vector<PeakAr> peaks;
    for (int i = 10; i <= intervals; i += 10) {
        const double height = spectrum[i];
        PeakAr peak = PeakAr::none;
        if (height > median) {
            // The JVM uses the one-sided calibration for a nonpositive
            // interior ordinate too; Nyquist always uses its left neighbor.
            const bool symmetric = i < intervals && height > 0.0;
            double rise = height - spectrum[i - 1];
            if (i < intervals) {
                const double right = height - spectrum[i + 1];
                rise = symmetric ? std::min(rise, right) : right;
            }
            const double ratio = rise / range;
            // The calibration table maps .8 + .0025 * count to significance.
            // Its 41st and 77th values are the first to exceed .9 and .99.
            if (ratio >= (symmetric ? 0.0968 : 0.2243)) peak = PeakAr::A;
            else if (ratio >= (symmetric ? 0.0353 : 0.1224)) peak = PeakAr::a;
        }
        peaks.push_back(peak);
    }
    return peaks;
}

BatteryTests tramo_monthly_seasonal_tests(const std::vector<double>& values, TramoSeasonalTestInput input) {
    TramoSeasonalityStatistics stats;
    std::vector<double> differenced;
    if (input == TramoSeasonalTestInput::linearized) {
        auto report = tramo_seasonality(values);
        stats = {report.qs, report.friedman};
        differenced = std::move(report.differenced);
    } else {
        stats = tramo_seasonality_statistics(values);
    }
    const auto& sequence = input == TramoSeasonalTestInput::linearized ? differenced : values;
    BatteryTests panel;
    panel.qs_v = stats.qs.value;
    panel.qs_p = stats.qs.p_value;
    panel.np_v = stats.friedman.value;
    panel.np_p = stats.friedman.p_value;
    panel.differenced_length = static_cast<int>(sequence.size());
    const auto ar = tramo_monthly_ar_peaks(sequence);
    const auto tukey = tramo_monthly_tukey_peaks(sequence);
    panel.ar_null = ar.empty();
    panel.tukey_null = tukey.empty();
    if (!panel.ar_null && !panel.tukey_null) {
        for (std::size_t i = 0; i < ar.size(); ++i) panel.peaks.push_back({ar[i], tukey[i]});
    }
    return panel;
}

} // namespace kalos
