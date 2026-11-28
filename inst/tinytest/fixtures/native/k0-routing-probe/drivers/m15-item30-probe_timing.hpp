#ifndef M15_ITEM30_PROBE_TIMING_HPP
#define M15_ITEM30_PROBE_TIMING_HPP
// item-30 probe: shared instrumentation declarations (measurement-only
// scratch; never shipped). Shared by m15-item30-filter_copy.cpp and
// m15-item30-route_probe.cpp.
#include <chrono>
#include <x86intrin.h>

namespace probe30 {

struct ProbeTiming {
    double arma_ns = 0, stat_init_ns = 0, assembly_ns = 0, factor_ns = 0;
    double loop_ns = 0, obs_ns = 0, upd_ns = 0, pred_ns = 0, diffuse_ns = 0;
    unsigned long long obs_calls = 0, upd_calls = 0, pred_calls = 0, diff_calls = 0;
    bool on = false;
};
inline ProbeTiming g_probe{};
inline double g_tsc_ns = 0.2712;  // calibrated by the driver at startup
inline void probe_reset() { const bool was_on = g_probe.on; g_probe = ProbeTiming{}; g_probe.on = was_on; }

inline unsigned long long probe_now() { return __rdtsc(); }

struct ProbeScope {
    unsigned long long start;
    double& sink;
    unsigned long long& calls;
    explicit ProbeScope(double& target, unsigned long long& counter)
        : start(g_probe.on ? probe_now() : 0), sink(target), calls(counter) {}
    ~ProbeScope() { if (g_probe.on) { sink += static_cast<double>(probe_now() - start) * g_tsc_ns; ++calls; } }
};

void probe_calibrate_tsc();  // defined in the driver (chrono vs TSC)

} // namespace probe30
#endif
