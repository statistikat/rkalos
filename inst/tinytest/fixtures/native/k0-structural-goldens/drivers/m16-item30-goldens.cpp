// Item-30 adoption slice (option (ii), maintainer ruling 2026-11-28): byte-identity
// oracle capture. Dumps the full named DiffuseLikelihood (%.17g, the frozen
// dossier §7 dump protocol of m15-item30-route_probe.cpp) for every frozen
// k = 0 shape, the full item-29 corpus correctness set (the non-k0 controls,
// including the D4-armed flagship and the small shapes), and the k = 0 gap /
// mean != 0 / near-boundary-root / zero-horizon / refusal arms the dossier §7
// coverage row assigns to the implementation slice. Captured at HEAD BEFORE
// the structural k = 0 steps ship; the test test_diffuse_k0_goldens replays
// every line bit-for-bit after the change.
//
// usage: m16-item30-goldens <corpus-dir>
//   corpus-dir = validation/evidence/native/diffuse-steady-state-probe
#include "kalos/diffuse_state_space.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace kalos;

struct Case {
    std::string id;
    int freq, d, seasonal_d;
    std::vector<double> ar, ma, sar, sma;
    std::vector<double> values;
    std::vector<std::uint8_t> missing;
};

std::vector<double> split_list(const std::string& text, char sep = ';') {
    std::vector<double> out;
    std::stringstream in(text);
    std::string item;
    while (std::getline(in, item, sep)) if (!item.empty()) out.push_back(std::stod(item));
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

Case load_case(const std::string& corpus_dir, const std::string& want) {
    std::ifstream cases(corpus_dir + "/cases.csv");
    std::string line;
    std::getline(cases, line);
    while (std::getline(cases, line)) {
        if (line.empty()) continue;
        auto f = split_fields(line, ',');
        if (f.size() < 11 || f[0] != want) continue;
        Case c;
        c.id = f[0];
        c.freq = std::stoi(f[1]);
        c.d = std::stoi(f[2]);
        c.seasonal_d = std::stoi(f[3]);
        c.ar = split_list(f[4], ':');
        c.ma = split_list(f[5], ':');
        c.sar = split_list(f[6], ':');
        c.sma = split_list(f[7], ':');
        load_series(corpus_dir + "/" + f[8], c.values, c.missing);
        for (const auto index : split_list(f[10])) {
            const auto i = static_cast<std::size_t>(index) - 1;
            if (i < c.missing.size()) c.missing[i] = 1;
        }
        return c;
    }
    throw std::runtime_error("case not found: " + want);
}

std::vector<std::string> case_ids() {
    return {"real-airline-noreg", "real-airline-x", "real-airline-gaps-x",
            "real-airline-redundant-x", "real-airline-zerocol-x", "syn-airline12",
            "syn-stat-arma21", "syn-stat-sarma11", "syn-airline4", "syn-mixed4",
            "nearur-airline", "nearni-airline", "d2-arma12", "slow-seasonal",
            "leadgap-airline-x"};
}

void print_g(const char* tag, double v) {
    std::printf("%s %.17g\n", tag, v);
}

// The frozen dossier §7 dump protocol (m15-item30-route_probe.cpp dump).
void dump_likelihood(const DiffuseLikelihood& l) {
    std::printf("STATUS likelihood=%s profile=%s\n", l.likelihood_status.c_str(), l.profile_status.c_str());
    std::printf("COUNTS observed=%zu diffuse=%zu finite=%zu initrank=%zu remain=%zu lastdiff=%zu\n",
                l.observed_count, l.diffuse_count, l.finite_count, l.initial_diffuse_rank,
                l.remaining_diffuse_rank, l.last_diffuse_index);
    print_g("SUM ssq", l.sum_squared_standardized);
    print_g("SUM logfinite", l.sum_log_finite_variances);
    print_g("SUM logdiffuse", l.sum_log_diffuse_variances);
    print_g("SUM fixed_ll", l.fixed_log_likelihood);
    print_g("SUM profiled_variance", l.profiled_variance);
    print_g("SUM profiled_ll", l.profiled_log_likelihood);
    const auto& n = l.numerics;
    std::printf("NUMERICS rankchecks=%zu exactcert=%d zerocount=%zu pivotcount=%zu stable=%d rebuilds=%zu sqr=%d compressions=%zu clipped=%zu routes=%zu route=%zu poststeps=%zu postrebuilds=%zu\n",
                n.propagation_rank_checks, n.exact_history_certification ? 1 : 0,
                n.exact_history_zero_count, n.exact_history_pivot_count,
                n.stable_history_propagation ? 1 : 0, n.history_reconstructions,
                n.square_root_finite_covariance ? 1 : 0, n.finite_factor_compressions,
                n.finite_factor_clipped_eigenvalues, n.history_frame_rebuilds,
                n.post_pivot_route, n.post_pivot_certified_steps, n.post_pivot_rebuilds);
    std::printf("RES projection=%.17g historynorm=%.17g factorinit=%.17g postres=%.17g postrebuildres=%.17g\n",
                n.maximum_projection_residual, n.maximum_history_normalization_residual,
                n.finite_factor_initialization_residual, n.post_pivot_certificate_residual,
                n.post_pivot_rebuild_residual);
    std::printf("UKIND");
    for (const auto kind : l.update_kind) std::printf(" %d", kind);
    std::printf("\n");
    for (std::size_t t = 0; t < l.predictions.size(); ++t) print_g("PRED", l.predictions[t]);
    for (std::size_t t = 0; t < l.finite_variances.size(); ++t) print_g("FVAR", l.finite_variances[t]);
    for (std::size_t t = 0; t < l.innovations.values.size(); ++t) print_g("INNO", l.innovations.values[t]);
    for (std::size_t t = 0; t < l.standardized_innovations.values.size(); ++t) print_g("STDI", l.standardized_innovations.values[t]);
    for (const auto v : l.final_state) print_g("FSTATE", v);
    const auto print_matrix = [](const char* tag, const DenseMatrix& p) {
        for (const auto v : p.values()) print_g(tag, v);
    };
    print_matrix("INITCOV", l.initial_finite_covariance);
    print_matrix("INITDIFF", l.initial_diffuse_covariance);
    print_matrix("FINALCOV", l.final_finite_covariance);
    print_matrix("FINALDIFF", l.final_diffuse_covariance);
}

void emit_case(const Case& c, double mean) {
    std::printf("ARM %s\n", c.id.c_str());
    std::printf("CFG freq=%d d=%d D=%d mean=%.17g n=%zu\n", c.freq, c.d, c.seasonal_d,
                mean, c.values.size());
    const MaskedSeries series{c.values, c.missing};
    try {
        const auto lik = sarima_innovations(series, c.ar, c.ma, c.sar, c.sma,
                                            c.freq, c.d, c.seasonal_d, mean);
        dump_likelihood(lik.likelihood);
    } catch (const std::exception& error) {
        std::printf("THROW %s\n", error.what());
    }
}

void emit_refusal(const char* id, const std::vector<double>& ar, const std::vector<double>& ma,
                  const std::vector<double>& sar, const std::vector<double>& sma,
                  int freq, int d, int seasonal_d, double mean,
                  const std::vector<double>& values, const std::vector<std::uint8_t>& missing) {
    std::printf("ARM %s\n", id);
    std::printf("CFG freq=%d d=%d D=%d mean=%.17g n=%zu\n", freq, d, seasonal_d, mean,
                values.size());
    const MaskedSeries series{values, missing};
    try {
        const auto lik = sarima_innovations(series, ar, ma, sar, sma, freq, d, seasonal_d, mean);
        std::printf("NO-THROW observed=%zu\n", lik.likelihood.observed_count);
    } catch (const std::exception& error) {
        std::printf("THROW %s\n", error.what());
    }
}

} // namespace

int main(int argc, char** argv) try {
    if (argc < 2) {
        std::cerr << "usage: m16-item30-goldens <corpus-dir>\n";
        return 2;
    }
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    setenv("OMP_NUM_THREADS", "1", 1);
    setenv("KALOS_NUM_THREADS", "1", 1);
    const std::string corpus = argv[1];
    std::printf("# kalos item-30 byte-identity oracle\n");
    std::printf("# dump protocol: docs/item30-k0-routing-probe-2026-11-28.md S7\n");
    for (const auto& id : case_ids()) {
        const auto c = load_case(corpus, id);
        emit_case(c, 0.0);
    }
    // ---- implementation-slice coverage arms (dossier S7 row 2) ----
    {
        auto c = load_case(corpus, "syn-stat-sarma11");
        c.id = "k0-sarma11-gap";
        for (const auto index : {18U, 42U, 43U, 101U, 200U}) c.missing[index] = 1;
        emit_case(c, 0.0);
    }
    {
        auto c = load_case(corpus, "syn-stat-arma21");
        c.id = "k0-arma21-mean";
        emit_case(c, 1.25);
    }
    {
        auto c = load_case(corpus, "syn-stat-arma21");
        c.id = "k0-nearboundary";
        c.ar = {0.999};
        c.ma = {-0.3};
        emit_case(c, 0.0);
    }
    {
        auto c = load_case(corpus, "syn-stat-sarma11");
        c.id = "k0-nearboundary-seasonal";
        c.ar = {0.5};
        c.ma = {-0.3};
        c.sar = {0.999};
        c.sma = {-0.2};
        emit_case(c, 0.0);
    }
    {
        const auto c = load_case(corpus, "syn-stat-sarma11");
        std::printf("ARM k0-zerohorizon\n");
        std::printf("CFG freq=%d d=%d D=%d mean=0 n=0\n", c.freq, c.d, c.seasonal_d);
        const MaskedSeries series{{}, {}};
        const auto lik = sarima_innovations(series, c.ar, c.ma, c.sar, c.sma,
                                            c.freq, c.d, c.seasonal_d);
        dump_likelihood(lik.likelihood);
    }
    {
        const auto c = load_case(corpus, "syn-stat-arma21");
        emit_refusal("k0-refuse-period", c.ar, c.ma, c.sar, c.sma, 6, 0, 0, 0.0, c.values, c.missing);
        emit_refusal("k0-refuse-dorder", c.ar, c.ma, c.sar, c.sma, 12, 3, 0, 0.0, c.values, c.missing);
        emit_refusal("k0-refuse-dmean", c.ar, c.ma, c.sar, c.sma, 12, 1, 0, 0.5, c.values, c.missing);
        emit_refusal("k0-refuse-nanmean", c.ar, c.ma, c.sar, c.sma, 12, 0, 0,
                     std::numeric_limits<double>::quiet_NaN(), c.values, c.missing);
    }
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
}
