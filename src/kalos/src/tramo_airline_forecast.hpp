#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace kalos::detail {

// Conditional distribution of the last 13 innovations of the differenced
// monthly airline model. Each observed response or design column uses the
// same gains, with its own conditional mean. Future shocks are independent.
class AirlineForecastBoundary {
public:
    using State = std::array<double, 13>;

    AirlineForecastBoundary(double theta, double seasonal_theta)
        : theta_(theta), seasonal_(seasonal_theta), interaction_(theta * seasonal_theta) {
        for (int i = 0; i < 13; ++i) covariance_[i][i] = 1.0;
    }

    double advance() {
        State cross{};
        for (int i = 0; i < 13; ++i)
            cross[i] = theta_ * covariance_[i][0] + seasonal_ * covariance_[i][11]
                     + interaction_ * covariance_[i][12];
        const double variance = 1.0 + theta_ * cross[0] + seasonal_ * cross[11] + interaction_ * cross[12];
        if (!(variance > 0.0) || !std::isfinite(variance)) return variance;
        gain_[0] = 1.0;
        std::copy_n(cross.begin(), 12, gain_.begin() + 1);
        for (int i = 12; i >= 0; --i) {
            for (int j = 12; j >= 0; --j) {
                const double prior = i && j ? covariance_[i - 1][j - 1] : double(i == 0 && j == 0);
                covariance_[i][j] = prior - gain_[i] * gain_[j] / variance;
            }
        }
        return variance;
    }

    double observe(State& state, double value, double variance) const {
        const double innovation = value - theta_ * state[0] - seasonal_ * state[11] - interaction_ * state[12];
        for (int i = 12; i >= 0; --i)
            state[i] = (i ? state[i - 1] : 0.0) + gain_[i] * innovation / variance;
        return innovation;
    }

    // The impulse response of 1 / ((1-B)(1-B^12)).
    static double integrated_impulse(int lag) { return lag < 0 ? 0.0 : double(1 + lag / 12); }

    double forecast_variance(int horizon, State& weights, double& cumulative_variance) const {
        const auto g = integrated_impulse;
        const double psi = g(horizon - 1) + theta_ * g(horizon - 2) + seasonal_ * g(horizon - 13)
                         + interaction_ * g(horizon - 14);
        cumulative_variance += psi * psi;
        for (int i = 0; i < 13; ++i) {
            if (i == 0) weights[i] += theta_ * g(horizon - 1);
            if (i < 12) weights[i] += seasonal_ * g(horizon - 12 + i);
            weights[i] += interaction_ * g(horizon - 13 + i);
        }
        double process = cumulative_variance;
        for (int i = 0; i < 13; ++i)
            for (int j = 0; j < 13; ++j) process += weights[i] * covariance_[i][j] * weights[j];
        return process;
    }

private:
    double theta_, seasonal_, interaction_;
    std::array<State, 13> covariance_{};
    State gain_{};
};

}  // namespace kalos::detail
