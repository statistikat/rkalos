// Item-29 probe slice: exact-diffuse-path golden generator + transient table.
// Read-only against the kalos library (ec17543 class). Scratch probe, not
// production code; linked against the probe-build libkalos.a.
//
// usage: filter_probe cases.csv outdir
// cases.csv columns: id,freq,d,seasonal_d,ar,ma,sar,sma,series,x,missing
//   ar/ma/sar/sma: semicolon-separated lists (may be empty)
//   series,x: csv paths (header line, value in column 1, missing flag in
//             column 2 when present); x empty = no regressor
//   missing: semicolon-separated 1-based indices into the series ("" allowed)
#include "kalos/diffuse_state_space.hpp"
#include "kalos/regression_likelihood.hpp"
#include "kalos/state_space.hpp"
#include "kalos/numerics.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {

std::vector<double> split_list(const std::string& text, char sep) {
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
    std::getline(in, line); // header
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

// Multi-column design reader: every field after the header is a finite
// design value; values come out column-major with ncols columns.
void load_design(const std::string& path, std::vector<double>& values,
                 std::size_t& columns) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line); // header
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

std::string g(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.17g", value);
    return buffer;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: filter_probe cases.csv outdir\n";
        return 2;
    }
    // Match the library contract: single-threaded, no FMA contraction at the
    // BLAS boundary is a build property; here pin one worker.
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    setenv("OMP_NUM_THREADS", "1", 1);
    setenv("KALOS_NUM_THREADS", "1", 1);
    std::ifstream cases(argv[1]);
    std::string line;
    std::getline(cases, line); // header
    while (std::getline(cases, line)) {
        if (line.empty()) continue;
        auto f = split_fields(line, ',');
        if (f.size() < 11) { std::cerr << "bad case line: " << line << "\n"; return 2; }
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
        const auto extra_missing = split_list(f[10], ';');
        for (const auto index : extra_missing) {
            const auto i = static_cast<std::size_t>(index) - 1;
            if (i < series_missing.size()) series_missing[i] = 1;
        }
        const kalos::MaskedSeries series{series_values, series_missing};

        std::ofstream out(std::string(argv[2]) + "/" + id + ".txt");
        out << "CASE " << id << " n=" << series_values.size() << " freq=" << freq
            << " d=" << d << " D=" << seasonal_d << "\n";

        // 1. Bare diffuse-filter likelihood at variance 1 (the per-evaluation
        //    covariance machine of the fixed-order ML path).
        auto lik = kalos::sarima_innovations(series, ar, ma, sar, sma, freq, d,
                                             seasonal_d, 0.0, 1.0, {}, false);
        const auto& l = lik.likelihood;
        out << "FILTER m_state=" << lik.model.transition.rows()
            << " r_stat=" << lik.model.transition.rows() - l.initial_diffuse_rank
            << " k_hist=" << l.initial_diffuse_rank
            << " observed=" << l.observed_count << " diffuse=" << l.diffuse_count
            << " finite=" << l.finite_count
            << " last_diffuse=" << l.last_diffuse_index
            << " remaining=" << l.remaining_diffuse_rank << "\n";
        out << "SUMS ssq=" << g(l.sum_squared_standardized)
            << " slogF=" << g(l.sum_log_finite_variances)
            << " slogD=" << g(l.sum_log_diffuse_variances)
            << " fixedll=" << g(l.fixed_log_likelihood)
            << " profiledvar=" << g(l.profiled_variance)
            << " profiledll=" << g(l.profiled_log_likelihood)
            << " status=" << l.likelihood_status << "," << l.profile_status << "\n";
        out << "NUMERICS rankchecks=" << l.numerics.propagation_rank_checks
            << " maxproj=" << g(l.numerics.maximum_projection_residual)
            << " exact_hist=" << (l.numerics.exact_history_certification ? 1 : 0)
            << " exactzero=" << l.numerics.exact_history_zero_count
            << " exactpivot=" << l.numerics.exact_history_pivot_count
            << " stable_hist=" << (l.numerics.stable_history_propagation ? 1 : 0)
            << " histrebuild=" << l.numerics.history_reconstructions
            << " histframes=" << l.numerics.history_frame_rebuilds
            << " histnormres=" << g(l.numerics.maximum_history_normalization_residual)
            << " sqrtfinite=" << (l.numerics.square_root_finite_covariance ? 1 : 0)
            << " compressions=" << l.numerics.finite_factor_compressions
            << " clipped=" << l.numerics.finite_factor_clipped_eigenvalues
            << " initres=" << g(l.numerics.finite_factor_initialization_residual) << "\n";
        out << "FTRACE";
        for (std::size_t t = 0; t < l.finite_variances.size(); ++t) {
            out << ' ' << g(l.finite_variances[t]);
        }
        out << "\n";
        out << "UKIND";
        for (const auto kind : l.update_kind) out << ' ' << kind;
        out << "\n";

        // Transient table: F-bar = last observed finite variance; t*(tau) =
        // first index from which all observed finite variances stay within
        // tau (relative) of F-bar, ignoring missing times.
        double fbar = 0;
        for (std::size_t t = 0; t < l.finite_variances.size(); ++t) {
            if (l.update_kind[t] == 2) fbar = l.finite_variances[t];
        }
        out << "FBAR " << g(fbar) << "\n";
        for (const double tau : {1e-10, 1e-12, 1e-14}) {
            std::size_t tstar = 0;
            for (std::size_t t = 0; t < l.finite_variances.size(); ++t) {
                if (l.update_kind[t] != 2) continue;
                const double rel = std::abs(l.finite_variances[t] - fbar) / fbar;
                if (rel > tau) tstar = t + 1;
            }
            out << "TSTAR tau=" << tau << " tstar=" << tstar << "\n";
        }

        // 2. Regression-concentrated likelihood surface (the actual objective
        //    the fixed-order ML path evaluates per step), design = x columns.
        if (!design.empty()) {
            kalos::RegressionDesign design_spec;
            design_spec.columns = design_cols;
            design_spec.values = design;
            design_spec.fixed = std::vector<std::uint8_t>(design_cols, 0);
            design_spec.coefficients = std::vector<double>(design_cols, 0);
            const auto reg = kalos::sarima_regression_likelihood(
                series, design_spec, ar, ma, sar, sma, freq, d, seasonal_d);
            out << "REG status=" << reg.status
                << " rank=" << reg.rank << " residdf=" << reg.residual_df
                << " finite=" << reg.finite_count
                << " ssq=" << g(reg.sum_squared_standardized)
                << " slogF=" << g(reg.sum_log_finite_variances)
                << " slogD=" << g(reg.sum_log_diffuse_variances)
                << " fixedll=" << g(reg.fixed_log_likelihood)
                << " profiledll=" << g(reg.profiled_log_likelihood)
                << " ztol=" << g(reg.whitening_zero_tolerance) << "\n";
            out << "REGCOEF";
            for (const auto value : reg.coefficients) out << ' ' << g(value);
            out << "\nREGNORM";
            for (const auto value : reg.whitening_norm_ratio) out << ' ' << g(value);
            out << "\nREGZERO";
            for (const auto flag : reg.numerically_zero_columns) out << ' ' << static_cast<int>(flag);
            out << "\n";
        }
        out.close();
        std::cout << "OK " << id << "\n";
    }
    return 0;
}
