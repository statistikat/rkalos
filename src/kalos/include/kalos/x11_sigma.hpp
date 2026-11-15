#ifndef KALOS_X11_SIGMA_HPP
#define KALOS_X11_SIGMA_HPP
// Item 22b.3a: moving five-year dispersion windows and two-pass sigma weights
// for the default X-11 extreme-value stage. Two rules are delivered and both
// are qualified; the policy must be stated by the caller because the default
// choice remains an open maintainer decision (see docs/contracts/x11-sigma-
// contract.md). Replacement values and correction factors are separate work
// (22b.3b) and are not produced here.
#include "numerics.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
namespace kalos {
// `legacy` reproduces the pinned JDemetra+ default corrector exactly: years
// whose dispersion magnitude is below 1e-5 keep full weight and threshold
// comparisons use a 1e-15 absolute slack. `accurate` applies the scale-limit
// rule: no dispersion cutoff, no slack, and a zero dispersion grants full
// weight only to exactly zero deviations. Neither is a promoted default.
enum class X11SigmaPolicy { legacy, accurate };
enum class X11SigmaWindowKind {
    global_single,      // Fewer than five complete years: one dispersion over the span.
    leading_extension,  // Wider head window shared by all years before the first centered one.
    central,            // Exactly five complete calendar years centered on this year.
    trailing_extension, // Window from the last centered start through the span end.
    boundary_reuse      // Year inherits the adjacent centered year's dispersion value.
};
enum class X11SigmaWeightState {
    full,                           // Deviation within the lower threshold: weight one.
    ramp,                           // Strictly between thresholds: weight in (0,1).
    zero,                           // Deviation beyond the upper threshold: weight zero.
    full_below_dispersion_cutoff,   // Legacy-only exempt year: |dispersion| < 1e-5.
    full_without_dispersion,        // Dispersion is unavailable (NaN): weight stays one.
    missing_input                   // Masked observation: never flagged, weight one.
};
struct X11SigmaWindow {
    double dispersion=0; // Moving RMS of (value - 0 or 1) over the window.
    std::size_t start=0; // First analysis-span index used.
    std::size_t length=0; // Observation count of the window.
    std::size_t usable=0; // Non-missing observations inside the window.
    X11SigmaWindowKind kind=X11SigmaWindowKind::central;
};
struct X11SigmaWeights {
    std::vector<double> weights;   // One per observation; missing inputs keep one.
    std::vector<X11SigmaWeightState> states;
    // One entry per calendar year of the analysis span, in order. With fewer
    // than five complete years a single global entry is retained and reused.
    std::vector<X11SigmaWindow> windows;
    std::vector<double> first_pass_weights;   // Equal to weights after one pass.
    std::vector<X11SigmaWindow> first_pass_windows;
    std::vector<std::size_t> first_pass_zeroed; // Rows excluded before pass two.
    int passes=0;
    std::string status; // "available" or "no_effective_observations".
};
// `series` is the current irregular-type sequence (for example B13 semantics);
// frequency is 4 or 12 and start_period is the 1-based calendar period of the
// first row. Projections stay part of the span unless `exclude_projections`
// drops exactly `backcasts` leading and `forecasts` trailing rows from the
// dispersion analysis; weights are always reported for every row. Requires
// upper > lower > 0.5 as in the pinned sigma-limit setter. Legacy numeric
// equality holds operation for operation with the pinned reference.
X11SigmaWeights x11_sigma_weights(const MaskedSeries& series,int frequency,int start_period,
    X11SigmaPolicy policy,int backcasts=0,int forecasts=0,bool exclude_projections=false,
    double lower_sigma=1.5,double upper_sigma=2.5,bool multiplicative=false,
    const Interrupt& interrupt={});
} // namespace kalos
#endif
