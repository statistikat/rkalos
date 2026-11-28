// item-30 probe (roadmap item 30, probe slice): k = 0 stationary-factorization
// routing. Measurement-only; links the HEAD libkalos.a plus the namespaced
// filter copy (m15-item30-filter_copy.cpp). No core code is touched.
//
// usage: m15-item30-route_probe <command> ...
//   ab     cases.csv case-id reps        interleaved A/B wall timings per route
//   nsweep cases.csv case-id reps        per-length medians per route (slopes)
//   instr  cases.csv case-id route reps  TSC-section decomposition (instrumented)
//   dump   cases.csv case-id route       canonical %.17g digest of the likelihood
//   sample cases.csv case-id reps        SIGPROF self-sampling over library eval
// routes: lib copy0 fac1 fac2 struct3
#include "kalos/diffuse_state_space.hpp"
#include "kalos/regression_likelihood.hpp"
#include "m15-item30-probe_timing.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <ucontext.h>
#include <unordered_map>
#include <vector>

namespace probe30 {

kalos::SarimaLikelihood sarima_innovations_route(const kalos::MaskedSeries&,
    const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    const std::vector<double>&, int, int, int, double, double, const kalos::Interrupt&,
    bool, int);
kalos::SarimaLikelihood sarima_innovations_copy(const kalos::MaskedSeries&,
    const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    const std::vector<double>&, int, int, int);

void probe_calibrate_tsc() {
    using clock = std::chrono::steady_clock;
    const auto c0 = clock::now();
    const auto t0 = probe_now();
    while (std::chrono::duration<double, std::milli>(clock::now() - c0).count() < 50.0) {
        __builtin_ia32_rdtsc();
    }
    const auto t1 = probe_now();
    const auto c1 = clock::now();
    const double ns = std::chrono::duration<double, std::nano>(c1 - c0).count();
    g_tsc_ns = ns / static_cast<double>(t1 - t0);
}

} // namespace probe30

namespace {

using clock_type = std::chrono::steady_clock;

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

Case load_case(const std::string& csv_path, const std::string& want) {
    std::ifstream cases(csv_path);
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
        load_series(f[8], c.values, c.missing);
        for (const auto index : split_list(f[10])) {
            const auto i = static_cast<std::size_t>(index) - 1;
            if (i < c.missing.size()) c.missing[i] = 1;
        }
        return c;
    }
    throw std::runtime_error("case not found: " + want);
}

double median_of(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

enum Route { kLib, kCopy0, kFac1, kFac2, kStruct3, kStat, kRouteCount };
const char* route_name[] = {"lib", "copy0", "fac1", "fac2", "struct3", "stat"};
const char* route_spec[] = {"lib", "copy0", "fac1", "fac2", "struct3", "stat"};

int parse_route(const std::string& text) {
    for (int i = 0; i < kRouteCount; ++i) if (text == route_spec[i]) return i;
    throw std::runtime_error("unknown route: " + text);
}

kalos::MaskedSeries make_series(const Case& c, std::size_t n) {
    return kalos::MaskedSeries{
        std::vector<double>(c.values.begin(), c.values.begin() + static_cast<std::ptrdiff_t>(n)),
        std::vector<std::uint8_t>(c.missing.begin(), c.missing.begin() + static_cast<std::ptrdiff_t>(n))};
}

double run_stat(const Case& c, std::size_t n) {
    // option-(iv) falsification: dispatch the k=0 shape to the existing
    // stationary_innovations primitive (public API, libkalos build).
    auto arma = kalos::stationary_arma(c.ar, c.ma, c.sar, c.sma, c.freq);
    const auto series = make_series(c, n);
    const auto lik = kalos::stationary_innovations(series, arma.state_space);
    return lik.sum_squared_standardized;
}

kalos::SarimaLikelihood run_route(const Case& c, int route, std::size_t n) {
    const auto series = make_series(c, n);
    switch (route) {
        case kLib: return kalos::sarima_innovations(series, c.ar, c.ma, c.sar, c.sma,
                                                    c.freq, c.d, c.seasonal_d);
        case kCopy0: return probe30::sarima_innovations_copy(series, c.ar, c.ma, c.sar, c.sma,
                                                             c.freq, c.d, c.seasonal_d);
        case kFac1: return probe30::sarima_innovations_route(series, c.ar, c.ma, c.sar, c.sma,
                                                             c.freq, c.d, c.seasonal_d, 0.0, 1.0, {}, false, 1);
        case kFac2: return probe30::sarima_innovations_route(series, c.ar, c.ma, c.sar, c.sma,
                                                             c.freq, c.d, c.seasonal_d, 0.0, 1.0, {}, false, 2);
        default: return probe30::sarima_innovations_route(series, c.ar, c.ma, c.sar, c.sma,
                                                          c.freq, c.d, c.seasonal_d, 0.0, 1.0, {}, false, 3);
    }
}

void print_g(const char* tag, double v) {
    std::printf("%s %.17g\n", tag, v);
}

void dump_likelihood(const kalos::DiffuseLikelihood& l) {
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
    const auto print_matrix = [](const char* tag, const kalos::DenseMatrix& p) {
        for (const auto v : p.values()) print_g(tag, v);
    };
    print_matrix("INITCOV", l.initial_finite_covariance);
    print_matrix("INITDIFF", l.initial_diffuse_covariance);
    print_matrix("FINALCOV", l.final_finite_covariance);
    print_matrix("FINALDIFF", l.final_diffuse_covariance);
}

// ---------------- SIGPROF self-sampling (no perf, no ptrace) ----------------
std::vector<uintptr_t> g_pcs;
uintptr_t g_overflow = 0;

void prof_handler(int, siginfo_t*, void* ucontext) {
    auto* ctx = static_cast<ucontext_t*>(ucontext);
    const uintptr_t pc = static_cast<uintptr_t>(ctx->uc_mcontext.gregs[REG_RIP]);
    if (g_pcs.size() < 4000000) g_pcs.push_back(pc); else ++g_overflow;
}

} // namespace

int main(int argc, char** argv) try {
    if (argc < 4) {
        std::cerr << "usage: m15-item30-route_probe <ab|nsweep|instr|dump|sample> cases.csv case-id [args]\n";
        return 2;
    }
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    setenv("OMP_NUM_THREADS", "1", 1);
    setenv("KALOS_NUM_THREADS", "1", 1);
    const std::string command = argv[1];
    const auto c = load_case(argv[2], argv[3]);
    std::printf("CASE %s n=%zu d=%d D=%d period=%d\n", c.id.c_str(), c.values.size(),
                c.d, c.seasonal_d, c.freq);

    if (command == "dump") {
        const int route = parse_route(argv[4]);
        const auto lik = run_route(c, route, c.values.size());
        dump_likelihood(lik.likelihood);
        return 0;
    }

    if (command == "ab") {
        const int reps = argc > 4 ? std::stoi(argv[4]) : 9;
        if (c.d != 0 || c.seasonal_d != 0) { std::cerr << "ab is for k=0 cases\n"; return 2; }
        // Warm every route once (first-touch, lazy LAPACK workspaces).
        for (int r = 0; r < kRouteCount; ++r) run_route(c, r, c.values.size());
        std::vector<std::vector<double>> times(kRouteCount);
        std::vector<int> order = {kLib, kCopy0, kFac1, kFac2, kStruct3, kStat};
        for (int rep = 0; rep < reps; ++rep) {
            std::rotate(order.begin(), order.begin() + rep % kRouteCount, order.end());
            std::vector<double> row(kRouteCount);
            for (const int route : order) {
                const auto t0 = clock_type::now();
                if (route == kStat) {
                    if (run_stat(c, c.values.size()) == 0) std::cerr << "x";
                } else {
                    const auto lik = run_route(c, route, c.values.size());
                    if (lik.likelihood.observed_count == 0) std::cerr << "x";
                }
                row[route] = std::chrono::duration<double, std::milli>(clock_type::now() - t0).count();
            }
            for (int route = 0; route < kRouteCount; ++route) times[route].push_back(row[route]);
            std::printf("REP %d", rep);
            for (int route = 0; route < kRouteCount; ++route) std::printf(" %s=%.6f", route_name[route], row[route]);
            std::printf("\n");
        }
        std::printf("MEDIAN %s", c.id.c_str());
        for (int route = 0; route < kRouteCount; ++route) std::printf(" %s=%.6f", route_name[route], median_of(times[route]));
        std::printf("\n");
        return 0;
    }

    if (command == "nsweep") {
        const int reps = argc > 4 ? std::stoi(argv[4]) : 5;
        for (int route = 0; route < kRouteCount; ++route) run_route(c, route, c.values.size());
        for (const std::size_t n : {24, 48, 96, 120, 144, 240}) {
            if (n > c.values.size()) continue;
            std::printf("NSWEEP %s n=%zu", c.id.c_str(), n);
            for (int route = 0; route < kStat; ++route) {
                std::vector<double> times;
                for (int r = 0; r < reps; ++r) {
                    const auto t0 = clock_type::now();
                    const auto lik = run_route(c, route, n);
                    times.push_back(std::chrono::duration<double, std::milli>(clock_type::now() - t0).count());
                    if (lik.likelihood.observed_count == 0) std::cerr << "x";
                }
                std::printf(" %s=%.6f", route_name[route], median_of(times));
            }
            std::printf("\n");
        }
        return 0;
    }

    if (command == "instr") {
        const int route = parse_route(argv[4]);
        const int reps = argc > 5 ? std::stoi(argv[5]) : 9;
        probe30::probe_calibrate_tsc();
        std::fprintf(stderr, "TSC %.6f ns/tick\n", probe30::g_tsc_ns);
        run_route(c, route, c.values.size());
        probe30::g_probe.on = true;
        std::vector<double> wall;
        for (int r = 0; r < reps; ++r) {
            probe30::probe_reset();
            const auto t0 = clock_type::now();
            const auto lik = run_route(c, route, c.values.size());
            wall.push_back(std::chrono::duration<double, std::milli>(clock_type::now() - t0).count());
            if (r == 0) {
                const auto& p = probe30::g_probe;
                std::printf("SECTIONS arma_ms=%.4f statinit_ms=%.4f assembly_ms=%.4f obs_ms=%.4f upd_ms=%.4f pred_ms=%.4f obs_calls=%llu upd_calls=%llu pred_calls=%llu\n",
                            p.arma_ns / 1e6, p.stat_init_ns / 1e6, p.assembly_ns / 1e6,
                            p.obs_ns / 1e6, p.upd_ns / 1e6, p.pred_ns / 1e6,
                            p.obs_calls, p.upd_calls, p.pred_calls);
            }
            if (lik.likelihood.observed_count == 0) std::cerr << "x";
        }
        probe30::g_probe.on = false;
        std::printf("INSTRWALL %s route=%s median_ms=%.6f\n", c.id.c_str(), route_spec[route], median_of(wall));
        return 0;
    }

    if (command == "sample") {
        const int reps = argc > 4 ? std::stoi(argv[4]) : 120;
        run_route(c, kLib, c.values.size());
        struct sigaction action{};
        action.sa_sigaction = prof_handler;
        action.sa_flags = SA_SIGINFO | SA_RESTART;
        sigaction(SIGPROF, &action, nullptr);
        struct itimerval timer{};
        timer.it_interval.tv_usec = 1000;
        timer.it_value.tv_usec = 1000;
        setitimer(ITIMER_PROF, &timer, nullptr);
        double sink = 0;
        for (int r = 0; r < reps; ++r) {
            const auto lik = run_route(c, kLib, c.values.size());
            sink += lik.likelihood.sum_squared_standardized;
        }
        timer.it_interval.tv_usec = 0;
        timer.it_value.tv_usec = 0;
        setitimer(ITIMER_PROF, &timer, nullptr);
        std::fprintf(stderr, "sink %.17g overflow=%llu\n", sink, static_cast<unsigned long long>(g_overflow));
        std::sort(g_pcs.begin(), g_pcs.end());
        for (std::size_t i = 0; i < g_pcs.size();) {
            std::size_t j = i;
            while (j < g_pcs.size() && g_pcs[j] == g_pcs[i]) ++j;
            std::printf("%p %zu\n", reinterpret_cast<void*>(g_pcs[i]), j - i);
            i = j;
        }
        std::fprintf(stderr, "SAMPLES %zu distinct\n", g_pcs.size());
        return 0;
    }

    std::cerr << "unknown command " << command << "\n";
    return 2;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
}
