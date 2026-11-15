#include "kalos/x11_sigma.hpp"
#include "kalos/parallel.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace kalos {
namespace {
// Independent reimplementation of the pinned default extreme-value weighting
// stage (item 22b.3a); no upstream source is translated. The pinned reference
// rules, thresholds and pass structure are documented in
// docs/contracts/x11-sigma-contract.md and qualified by
// validation/fixtures/native/x11-sigma.
constexpr std::size_t cycles = 5;                    // Five-year moving dispersion.
constexpr double dispersion_cutoff = 1.0e-5;         // Legacy full-weight dispersion rule.
constexpr double threshold_slack = 1.0e-15;          // Legacy threshold comparison slack.
void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }

// One-window moving RMS over non-missing observations; the missing count is
// excluded from the denominator exactly as the pinned single-window rule.
double window_rms(const std::vector<double>& values,const std::vector<std::uint8_t>& missing,
    std::size_t start,std::size_t length,double center) {
    double sum = 0; std::size_t usable = 0;
    for(std::size_t i = 0; i < length; ++i) {
        const auto row = start + i;
        if(missing[row]) continue;
        ++usable;
        const auto deviation = values[row] - center;
        sum += deviation * deviation;
    }
    if(usable == 0) return std::nan(""); // Pinned rule: no usable row leaves the dispersion unavailable.
    // No range throw: the pinned reference accumulates in IEEE doubles, so an
    // out-of-range window yields +Inf dispersion and detection then keeps
    // every weight at one because every deviation stays below the infinite
    // thresholds. The same binary64 arithmetic reproduces that on any platform.
    return std::sqrt(sum / static_cast<double>(usable));
}

std::size_t window_usable(const std::vector<std::uint8_t>& missing,
    std::size_t start,std::size_t length) {
    std::size_t usable = 0;
    for(std::size_t i = 0; i < length; ++i) if(!missing[start + i]) ++usable;
    return usable;
}

// Dispersion windows for the effective span [lo,hi) of the full series.
// Indices in the result are analysis-span indices; `year_of_first` is the
// calendar-year index of span row zero within the reported window list.
struct Span {
    std::size_t lo=0,hi=0;
    std::int64_t first_year=0; // Calendar-year index of row lo.
    std::int64_t year_count=0; // Calendar years touched by the full span.
};

std::vector<X11SigmaWindow> dispersion_windows(const std::vector<double>& values,
    const std::vector<std::uint8_t>& missing,std::size_t frequency,std::int64_t start,
    const Span& span,double center,const Interrupt& interrupt) {
    const auto len = static_cast<std::int64_t>(span.hi - span.lo);
    const auto phase = static_cast<std::int64_t>((start + static_cast<std::int64_t>(span.lo)) % frequency);
    const auto head = phase == 0 ? std::int64_t{0} : static_cast<std::int64_t>(frequency) - phase;
    auto complete = len == 0 ? std::int64_t{0} : 1 + (phase + len - 1) / static_cast<std::int64_t>(frequency);
    if(head > 0) --complete;
    const auto tail = (phase + len) % static_cast<std::int64_t>(frequency);
    if(tail > 0) --complete;
    std::vector<X11SigmaWindow> windows(static_cast<std::size_t>(span.year_count));
    if(len == 0) return windows;
    if(complete < static_cast<std::int64_t>(cycles)) {
        const double dispersion = window_rms(values,missing,span.lo,static_cast<std::size_t>(len),center);
        X11SigmaWindow global;
        global.dispersion = dispersion; global.start = span.lo;
        global.length = static_cast<std::size_t>(len);
        global.usable = window_usable(missing,span.lo,static_cast<std::size_t>(len));
        global.kind = X11SigmaWindowKind::global_single;
        // The pinned rule keeps a single-entry array and clamps every year
        // onto it; the caller reproduces that clamp when reading windows.
        windows[0] = global;
        windows.resize(1);
        return windows;
    }
    const auto y0 = span.first_year + 2 + (head > 0 ? 1 : 0); // First centered year.
    const auto y1 = y0 + complete - static_cast<std::int64_t>(cycles); // Last centered year.
    std::vector<char> done(static_cast<std::size_t>(span.year_count), 0);
    const auto set_window = [&](std::int64_t year,std::size_t begin,std::size_t size,
        X11SigmaWindowKind kind) {
        auto& entry = windows[static_cast<std::size_t>(year)];
        entry.dispersion = window_rms(values,missing,begin,size,center);
        entry.start = begin; entry.length = size;
        entry.usable = window_usable(missing,begin,size);
        entry.kind = kind;
    };
    // Centered five-calendar-year windows; each is independent and the inner
    // accumulation order stays serial, so results match the serial build.
    const auto central = static_cast<std::size_t>(y1 - y0 + 1);
    parallel_for(central,1,interrupt,[&](std::size_t index) {
        const auto at = static_cast<std::size_t>(head) + span.lo + index * frequency;
        set_window(y0 + static_cast<std::int64_t>(index),at,frequency*cycles,
            X11SigmaWindowKind::central);
        done[static_cast<std::size_t>(y0 + static_cast<std::int64_t>(index))] = 1;
    });
    X11SigmaWindow leading;
    if(head > 0) {
        leading.dispersion = window_rms(values,missing,span.lo,
            static_cast<std::size_t>(head) + frequency*cycles,center);
        leading.start = span.lo;
        leading.length = static_cast<std::size_t>(head) + frequency*cycles;
        leading.usable = window_usable(missing,leading.start,leading.length);
        leading.kind = X11SigmaWindowKind::leading_extension;
    } else {
        leading = windows[static_cast<std::size_t>(y0)];
        leading.kind = X11SigmaWindowKind::boundary_reuse;
    }
    for(std::int64_t year = 0; year < y0; ++year) {
        windows[static_cast<std::size_t>(year)] = leading;
        done[static_cast<std::size_t>(year)] = 1;
    }
    X11SigmaWindow trailing;
    if(tail > 0) {
        const auto begin = static_cast<std::size_t>(head) + span.lo + (central - 1) * frequency;
        trailing.dispersion = window_rms(values,missing,begin,span.hi - begin,center);
        trailing.start = begin; trailing.length = span.hi - begin;
        trailing.usable = window_usable(missing,begin,span.hi - begin);
        trailing.kind = X11SigmaWindowKind::trailing_extension;
    } else {
        trailing = windows[static_cast<std::size_t>(y1)];
        trailing.kind = X11SigmaWindowKind::boundary_reuse;
    }
    for(std::int64_t year = y1 + 1; year < span.year_count; ++year) {
        windows[static_cast<std::size_t>(year)] = trailing;
        done[static_cast<std::size_t>(year)] = 1;
    }
    for(const char entry : done) if(!entry)
        throw std::logic_error("X-11 dispersion windows left an uncovered year");
    return windows;
}

// One detection pass over the full span against per-year dispersions. The
// pinned year boundaries advance from `period - start` in steps of period.
void detect(const std::vector<double>& values,const std::vector<std::uint8_t>& missing,
    std::size_t frequency,std::int64_t start,const std::vector<X11SigmaWindow>& windows,
    std::int64_t year_count,X11SigmaPolicy policy,double lower_sigma,double upper_sigma,
    double center,std::vector<double>& weights,std::vector<X11SigmaWeightState>& states,
    const Interrupt& interrupt) {
    std::int64_t yStart = 0, yEnd = static_cast<std::int64_t>(frequency) - start;
    for(std::int64_t year = 0; year < year_count; ++year) {
        yEnd = std::min(yEnd, static_cast<std::int64_t>(values.size()));
        const auto& window = windows[std::min<std::size_t>(
            static_cast<std::size_t>(year), windows.size() - 1)];
        const auto dispersion = window.dispersion;
        const bool exempt = policy == X11SigmaPolicy::legacy &&
            !std::isnan(dispersion) && std::abs(dispersion) < dispersion_cutoff;
        const bool without = std::isnan(dispersion);
        const auto lower = dispersion * lower_sigma, upper = dispersion * upper_sigma;
        for(auto row = static_cast<std::size_t>(yStart); row < static_cast<std::size_t>(yEnd); ++row) {
            if((row & 1023U) == 0) checkpoint(interrupt);
            if(missing[row]) { weights[row] = 1; states[row] = X11SigmaWeightState::missing_input; continue; }
            if(without) { weights[row] = 1; states[row] = X11SigmaWeightState::full_without_dispersion; continue; }
            if(exempt) { weights[row] = 1; states[row] = X11SigmaWeightState::full_below_dispersion_cutoff; continue; }
            const auto deviation = std::abs(values[row] - center);
            if(policy == X11SigmaPolicy::legacy) {
                if(deviation - upper > threshold_slack) { weights[row] = 0; states[row] = X11SigmaWeightState::zero; }
                else if(deviation - lower > threshold_slack) {
                    const auto value = (upper - deviation) / (upper - lower);
                    weights[row] = value;
                    states[row] = value <= 0 ? X11SigmaWeightState::zero : X11SigmaWeightState::ramp;
                } else { weights[row] = 1; states[row] = X11SigmaWeightState::full; }
            } else if(dispersion == 0) {
                // Scale limit of the ramp: full weight only for an exactly
                // zero deviation, zero weight for any nonzero deviation.
                if(deviation == 0) { weights[row] = 1; states[row] = X11SigmaWeightState::full; }
                else { weights[row] = 0; states[row] = X11SigmaWeightState::zero; }
            } else if(deviation > upper) { weights[row] = 0; states[row] = X11SigmaWeightState::zero; }
            else if(deviation > lower) {
                const auto value = (upper - deviation) / (upper - lower);
                weights[row] = value;
                states[row] = value <= 0 ? X11SigmaWeightState::zero : X11SigmaWeightState::ramp;
            } else { weights[row] = 1; states[row] = X11SigmaWeightState::full; }
        }
        yStart = yEnd;
        yEnd = std::min(yEnd + static_cast<std::int64_t>(frequency), static_cast<std::int64_t>(values.size()));
    }
}
} // namespace

X11SigmaWeights x11_sigma_weights(const MaskedSeries& series,int frequency,int start_period,
    X11SigmaPolicy policy,int backcasts,int forecasts,bool exclude_projections,
    double lower_sigma,double upper_sigma,bool multiplicative,const Interrupt& interrupt) {
    series.validate();
    if(frequency != 4 && frequency != 12)
        throw std::invalid_argument("X-11 sigma weights support monthly and quarterly series");
    if(start_period < 1 || start_period > frequency)
        throw std::invalid_argument("X-11 sigma start period must be a calendar period of the frequency");
    if(!(upper_sigma > lower_sigma) || !(lower_sigma > 0.5))
        throw std::invalid_argument("X-11 sigma limits must satisfy upper > lower > 0.5");
    if(backcasts < 0 || forecasts < 0)
        throw std::invalid_argument("X-11 sigma projection horizons must be nonnegative");
    const auto n = series.values.size();
    if(n == 0) throw std::invalid_argument("X-11 sigma weights require at least one observation");
    if(n < static_cast<std::size_t>(frequency))
        // The pinned reference has no defined behavior below one full calendar
        // period (its first-year boundary walk runs off the array); undefined
        // behavior is not a result to reproduce.
        throw std::invalid_argument("X-11 sigma weights require at least one full calendar period");
    if(exclude_projections &&
        static_cast<std::size_t>(backcasts) + static_cast<std::size_t>(forecasts) > n)
        // The pinned span trim rejects a range whose end precedes its start.
        throw std::invalid_argument("X-11 sigma projection horizons exceed the series");
    checkpoint(interrupt);
    const auto center = multiplicative ? 1.0 : 0.0;
    const auto start = static_cast<std::int64_t>(start_period - 1);
    X11SigmaWeights result;
    result.weights.assign(n, 1);
    result.states.assign(n, X11SigmaWeightState::full);
    Span span; span.hi = n;
    if(exclude_projections) {
        span.lo = static_cast<std::size_t>(backcasts);
        span.hi = n - static_cast<std::size_t>(forecasts);
    }
    span.first_year = static_cast<std::int64_t>(start + static_cast<std::int64_t>(span.lo)) / frequency;
    span.year_count = 1 + (start + static_cast<std::int64_t>(n) - 1) / frequency;
    if(span.lo >= span.hi) {
        // The pinned rule gives an empty span one unavailable dispersion and
        // detection then keeps every weight at one; the core reports that
        // state explicitly instead of emitting NaN windows.
        result.status = "no_effective_observations";
        result.passes = 1;
        result.first_pass_weights = result.weights;
        for(std::size_t row = 0; row < n; ++row)
            result.states[row] = series.missing[row] ? X11SigmaWeightState::missing_input
                : X11SigmaWeightState::full_without_dispersion;
        return result;
    }
    auto windows = dispersion_windows(series.values,series.missing,
        static_cast<std::size_t>(frequency),start,span,center,interrupt);
    std::vector<double> first(n, 1);
    std::vector<X11SigmaWeightState> first_states(n, X11SigmaWeightState::full);
    detect(series.values,series.missing,static_cast<std::size_t>(frequency),start,windows,
        span.year_count,policy,lower_sigma,upper_sigma,center,first,first_states,interrupt);
    result.first_pass_weights = first;
    result.first_pass_windows = windows;
    bool any = false;
    for(const auto value : first) if(value < 1) { any = true; break; }
    if(!any) {
        result.weights = std::move(first);
        result.states = std::move(first_states);
        result.windows = std::move(windows);
        result.passes = 1;
        result.status = "available";
        return result;
    }
    // Second pass: the pinned rule removes only exact zero weights and then
    // recomputes dispersions while still detecting on the original values.
    auto excluded = series.missing;
    for(std::size_t row = 0; row < n; ++row) if(first[row] == 0) {
        excluded[row] = 1;
        result.first_pass_zeroed.push_back(row);
    }
    auto second_windows = dispersion_windows(series.values,excluded,
        static_cast<std::size_t>(frequency),start,span,center,interrupt);
    detect(series.values,series.missing,static_cast<std::size_t>(frequency),start,second_windows,
        span.year_count,policy,lower_sigma,upper_sigma,center,result.weights,result.states,interrupt);
    result.windows = std::move(second_windows);
    result.passes = 2;
    result.status = "available";
    checkpoint(interrupt);
    return result;
}
} // namespace kalos
