// Monthly Tukey spectral-peak classification used by TRAMO SeasonalityTests.
#include "kalos/tramo_differencing_battery.hpp"
#include "kalos/detail/hostlineage.hpp"
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
        for (std::size_t i = lag; i < centered.size(); ++i) product += centered[i-lag] * centered[i];
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

std::array<double,4> tukey_degrees_of_freedom(std::size_t n, int window) {
    if (window == 79) return {6.35251, 19.6308, 2.29316, 6.55412};
    // Monthly long-window calibration of the two interior and the two
    // Nyquist F-test degrees of freedom, as functions of sample length.
    constexpr double coefficients[4][3] = {
        {0.54630, 2.93030, 2.2042}, {1.1329, 7.6924, 10.8795},
        {-0.3492, 1.533, 2.7696}, {0.9829, 3.8217, 6.9345}
    };
    std::array<double,4> df{};
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
    const std::array<int,5> frequencies = window == 112 ? std::array<int,5>{9, 19, 28, 37, 47}
                                                        : std::array<int,5>{7, 13, 20, 26, 33};
    std::vector<PeakTukey> peaks;
    for (int f : frequencies) {
        const double ratio = 2.0 * density(f) / (density(f+1) + density(f-1));
        peaks.push_back(classify(ratio, df[0], df[1]));
    }
    const int nyquist = window / 2;
    const double below = density(nyquist-1), at = density(nyquist);
    peaks.push_back(classify(below == 0.0 ? at : at / below, df[2], df[3]));
    return peaks;
}

} // namespace kalos
