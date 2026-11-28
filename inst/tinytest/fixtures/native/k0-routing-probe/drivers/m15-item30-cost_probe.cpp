// Item-29 probe slice: cost attribution for the exact diffuse path and a
// frozen-gain tail prototype. Scratch probe; read-only against the library.
//
// usage: cost_probe cases.csv phase
//   phase = eval | init | nsweep | tail | reg
// eval : median wall time of sarima_innovations per call (the per-evaluation
//        covariance machine), R repetitions.
// init : stationary_arma + stationary_initialization + companion assembly +
//        FiniteStateFactor construction (the pre-loop phase replicated from
//        sarima_innovations, public API only).
// nsweep: same over truncated series lengths -> per-step slope vs intercept.
// tail : frozen-gain scalar recursion over the post-transient tail, using the
//        converged (F-bar, K-bar) read from the save_gains run -> prototype
//        cost of the steady-state path per tail step.
// reg  : median wall time of sarima_regression_likelihood per call (the
//        objective actually evaluated by the fixed-order ML path with x).
#include "kalos/diffuse_state_space.hpp"
#include "kalos/exact_history.hpp"
#include "kalos/finite_state_factor.hpp"
#include "kalos/sparse_history.hpp"
#include "kalos/regression_likelihood.hpp"
#include "kalos/state_space.hpp"
#include "kalos/numerics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using clock_type = std::chrono::steady_clock;

std::vector<double> split_list(const std::string& text, char sep = ';') {
    std::vector<double> out;
    std::stringstream in(text);
    std::string item;
    while (std::getline(in, item, sep)) {
        if (!item.empty()) out.push_back(std::stod(item));
    }
    return out;
}

std::vector<std::string> split_fields(const std::string& line, char sep) {
    std::vector<std::string> out;
    std::stringstream in(line);
    std::string item;
    while (std::getline(in, item, sep)) out.push_back(item);
    return out;
}

void load_series(const std::string& path, std::vector<double>& values,
                 std::vector<std::uint8_t>& missing) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        auto fields = split_fields(line, ',');
        if (fields.empty() || fields[0].empty() || fields[0] == "NA") {
            values.push_back(0);
            missing.push_back(1);
            continue;
        }
        values.push_back(std::stod(fields[0]));
        missing.push_back(fields.size() > 1 && fields[1] == "1" ? 1 : 0);
    }
}

void load_design(const std::string& path, std::vector<double>& values,
                 std::size_t& columns) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    std::vector<std::vector<double>> rows;
    columns = 0;
    while (std::getline(in, line)) {
        auto fields = split_fields(line, ',');
        if (fields.size() > columns) columns = fields.size();
        rows.push_back({});
        for (const auto& field : fields) rows.back().push_back(std::stod(field));
    }
    values.assign(columns * rows.size(), 0);
    for (std::size_t t = 0; t < rows.size(); ++t)
        for (std::size_t j = 0; j < rows[t].size(); ++j) values[j * rows.size() + t] = rows[t][j];
}

double median_of(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: cost_probe cases.csv phase [reps]\n";
        return 2;
    }
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    setenv("OMP_NUM_THREADS", "1", 1);
    setenv("KALOS_NUM_THREADS", "1", 1);
    const std::string phase = argv[2];
    const int reps = argc > 4 ? std::stoi(argv[3]) : 9;
    std::ifstream cases(argv[1]);
    std::string line;
    std::getline(cases, line);
    while (std::getline(cases, line)) {
        if (line.empty()) continue;
        auto f = split_fields(line, ',');
        if (f.size() < 11) continue;
        const std::string id = f[0];
        const int freq = std::stoi(f[1]);
        const int d = std::stoi(f[2]);
        const int seasonal_d = std::stoi(f[3]);
        const auto ar = split_list(f[4], ':');
        const auto ma = split_list(f[5], ':');
        const auto sar = split_list(f[6], ':');
        const auto sma = split_list(f[7], ':');
        std::vector<double> series_values;
        std::vector<std::uint8_t> series_missing;
        load_series(f[8], series_values, series_missing);
        std::vector<double> design;
        std::size_t design_cols = 0;
        if (!f[9].empty()) load_design(f[9], design, design_cols);
        if (f.size() > 11 && !f[11].empty()) design_cols = static_cast<std::size_t>(std::stoul(f[11]));
        for (const auto index : split_list(f[10])) {
            const auto i = static_cast<std::size_t>(index) - 1;
            if (i < series_missing.size()) series_missing[i] = 1;
        }

        if (phase == "nsweep") {
            for (const std::size_t n : {24, 48, 96, 120, 144, 240}) {
                if (n > series_values.size()) continue;
                const kalos::MaskedSeries series{
                    std::vector<double>(series_values.begin(), series_values.begin() +
                        static_cast<std::ptrdiff_t>(n)),
                    std::vector<std::uint8_t>(series_missing.begin(), series_missing.begin() +
                        static_cast<std::ptrdiff_t>(n))};
                std::vector<double> times;
                for (int r = 0; r < std::max(3, reps / 2); ++r) {
                    const auto t0 = clock_type::now();
                    const auto lik = kalos::sarima_innovations(series, ar, ma, sar, sma,
                                                               freq, d, seasonal_d);
                    const auto dt = std::chrono::duration<double, std::milli>(
                        clock_type::now() - t0).count();
                    times.push_back(dt);
                    if (lik.likelihood.observed_count == 0) std::cerr << "x";
                }
                std::printf("NSWEEP %s n=%zu median_ms=%.6f\n", id.c_str(), n,
                            median_of(times));
            }
            continue;
        }

        const kalos::MaskedSeries series{series_values, series_missing};

        if (phase == "eval") {
            std::vector<double> times;
            for (int r = 0; r < reps; ++r) {
                const auto t0 = clock_type::now();
                const auto lik = kalos::sarima_innovations(series, ar, ma, sar, sma,
                                                           freq, d, seasonal_d);
                times.push_back(std::chrono::duration<double, std::milli>(
                    clock_type::now() - t0).count());
                if (lik.likelihood.observed_count == 0) std::cerr << "x";
            }
            std::printf("EVAL %s n=%zu median_ms=%.6f\n", id.c_str(),
                        series_values.size(), median_of(times));
        } else if (phase == "init") {
            // Replicates (public API only) everything sarima_innovations does
            // before the time loop: stationary representation, stationary
            // initialization, companion assembly, finite-factor construction.
            std::vector<double> times;
            for (int r = 0; r < reps; ++r) {
                const auto t0 = clock_type::now();
                auto arma = kalos::stationary_arma(ar, ma, sar, sma, freq);
                const auto init = kalos::stationary_initialization(arma.state_space);
                const auto diff = kalos::difference_polynomial(d, seasonal_d, freq);
                const auto r_ = arma.state_space.transition.rows();
                const auto k = diff.size() - 1;
                const auto m = r_ + k;
                kalos::DenseMatrix transition(m, m), loading(m, 1), initial_finite(m, m);
                std::vector<double> z(m, 0);
                for (std::size_t i = 0; i < r_; ++i) {
                    z[i] = arma.state_space.observation[i];
                    loading(i, 0) = arma.state_space.noise_loading(i, 0);
                    for (std::size_t j = 0; j < r_; ++j) {
                        transition(i, j) = arma.state_space.transition(i, j);
                        initial_finite(i, j) = init.covariance(i, j);
                    }
                }
                std::vector<double> basis(m * k, 0);
                if (k > 0) {
                    for (std::size_t i = 0; i < k; ++i) {
                        z[r_ + i] = -diff[i + 1];
                        basis[r_ + i + m * i] = 1;
                    }
                    for (std::size_t j = 0; j < m; ++j) transition(r_, j) = z[j];
                    for (std::size_t i = 1; i < k; ++i)
                        transition(r_ + i, r_ + i - 1) = 1;
                }
                kalos::detail::CountsHistoryCertificate history(d, seasonal_d, freq);
                kalos::detail::SparseHistoryBasis stable_history(d, seasonal_d, freq);
                std::unique_ptr<kalos::detail::FiniteStateFactor> factor;
                if (k > 0)
                    factor = std::make_unique<kalos::detail::FiniteStateFactor>(initial_finite);
                times.push_back(std::chrono::duration<double, std::milli>(
                    clock_type::now() - t0).count());
                if (factor && factor->dimension() == 0) std::cerr << "x";
            }
            std::printf("INIT %s median_ms=%.6f\n", id.c_str(), median_of(times));
        } else if (phase == "reg") {
            if (design.empty()) continue;
            kalos::RegressionDesign design_spec;
            design_spec.columns = design_cols;
            design_spec.values = design;
            design_spec.fixed = std::vector<std::uint8_t>(design_cols, 0);
            design_spec.coefficients = std::vector<double>(design_cols, 0);
            std::vector<double> times;
            for (int r = 0; r < reps; ++r) {
                const auto t0 = clock_type::now();
                const auto reg = kalos::sarima_regression_likelihood(
                    series, design_spec, ar, ma, sar, sma, freq, d, seasonal_d);
                times.push_back(std::chrono::duration<double, std::milli>(
                    clock_type::now() - t0).count());
                if (reg.status.empty()) std::cerr << "x";
            }
            std::printf("REG %s n=%zu median_ms=%.6f\n", id.c_str(),
                        series_values.size(), median_of(times));
        } else if (phase == "tail") {
            // One exact run with gains; F-bar and the final gain column define
            // the frozen recursion; time the frozen tail only.
            const auto lik = kalos::sarima_innovations(series, ar, ma, sar, sma,
                                                       freq, d, seasonal_d, 0.0, 1.0,
                                                       {}, true);
            const auto& l = lik.likelihood;
            const auto m = lik.model.transition.rows();
            const auto n = series_values.size();
            double fbar = 0;
            std::size_t tail_start = 0;
            for (std::size_t t = 0; t < n; ++t) {
                if (l.update_kind[t] == 2) {
                    fbar = l.finite_variances[t];
                    tail_start = t;
                }
            }
            // steady switch: last index where |F-fbar|/fbar > 1e-14 among finite steps
            std::size_t tstar = 0;
            for (std::size_t t = 0; t < n; ++t) {
                if (l.update_kind[t] != 2) continue;
                if (std::abs(l.finite_variances[t] - fbar) / fbar > 1e-14) tstar = t + 1;
            }
            std::vector<double> kbar(m, 0);
            for (std::size_t i = 0; i < m; ++i) kbar[i] = l.observation_gains[i + m * (n - 1)];
            const auto& z = lik.model.observation;
            const auto& T = lik.model.transition;
            // Warm the state to the exact filter's final state (cost of warm
            // state is part of the exact prefix in a real implementation).
            auto state = l.final_state;
            std::vector<double> times;
            for (int r = 0; r < reps * 20; ++r) {
                const auto t0 = clock_type::now();
                double ssq = 0, slog = 0;
                for (std::size_t t = tstar; t < n; ++t) {
                    if (series_missing[t]) continue; // frozen path must handle gaps
                    double pred = 0;
                    for (std::size_t i = 0; i < m; ++i) pred += z[i] * state[i];
                    const double nu = series_values[t] - pred;
                    for (std::size_t i = 0; i < m; ++i)
                        state[i] += kbar[i] * (nu / fbar);
                    std::vector<double> next(m, 0);
                    for (std::size_t i = 0; i < m; ++i) {
                        double value = 0;
                        for (std::size_t j = 0; j < m; ++j)
                            if (T(i, j) != 0) value += T(i, j) * state[j];
                        next[i] = value;
                    }
                    state = std::move(next);
                    ssq += nu * nu / fbar;
                    slog += std::log(fbar);
                }
                times.push_back(std::chrono::duration<double, std::micro>(
                    clock_type::now() - t0).count());
                if (ssq < 0) std::cerr << "x";
            }
            std::printf("TAIL %s n=%zu tstar=%zu tail=%zu fbar=%.17g tail_us=%.3f perstep_us=%.4f\n",
                        id.c_str(), n, tstar, n - tstar, fbar, median_of(times),
                        median_of(times) / static_cast<double>(n - tstar));
        }
    }
    return 0;
}
