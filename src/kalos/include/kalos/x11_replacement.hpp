#ifndef KALOS_X11_REPLACEMENT_HPP
#define KALOS_X11_REPLACEMENT_HPP
// Item 22b.3b: same-phase replacement values, neighbor provenance, endpoint
// extension, insufficient-neighbor fallback, correction factors and the
// applyCorrections identity of the default X-11 extreme-value stage. Consumes
// the final observation weights qualified by item 22b.3a
// (x11_sigma.hpp/docs/contracts/x11-sigma-contract.md) and adds the
// replacement pass on top; see docs/contracts/x11-replacement-contract.md.
// B17/B20 table placement, the CV correction variant and period/group
// specific alternatives are not delivered or claimed here.
#include "numerics.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
namespace kalos {
// Per-row provenance of the replacement pass.
enum class X11ReplacementKind {
    none,           // Final weight is exactly one: no replacement, corrections cancel.
    four_neighbors, // Blend of the row and four fully weighted same-phase rows.
    period_average, // Fewer than four fully weighted same-phase rows exist.
    outside_span    // Projection row dropped from the replacement span.
};
struct X11ReplacementPlan {
    std::vector<double> replacements; // NaN wherever no replacement exists.
    std::vector<X11ReplacementKind> kinds;
    // Neighbor positions in pinned discovery order (side walks first, then
    // the opposite-side extension); -1 padded and empty for non-blend rows.
    std::vector<std::array<std::int64_t,4>> neighbors;
    std::vector<std::int64_t> neighbor_counts; // Exactly 4 or 0.
    std::vector<double> period_averages; // One per calendar phase of the series.
    std::vector<double> correction_factors; // Zero (additive) or one (multiplicative) defaults.
    std::vector<double> corrected; // Original effective values with replacements applied.
    std::int64_t fallbacks=0; // period_average row count.
    std::string status; // "available".
};
// `weights` are the final item 22b.3a weights (same length as the series, all
// finite; legacy slack may make one slightly negative, which is kept). Same-
// phase rows are rows `i`, `i +/- frequency`: no calendar start enters the
// walk, matching the pinned reference. Masked rows are effective-NaN values
// that keep weight one, so they can be selected as fully weighted neighbors
// and then poison a blend to NaN exactly like the pinned reference; the same
// holds for the period-average fallback. With `exclude_projections` the span
// drops `backcasts` leading and `forecasts` trailing rows before replacement
// (requiring `backcasts + forecasts <= n`; beyond that the pinned reference
// walks outside the output array); outside rows report kind `outside_span`
// and keep a NaN replacement. Correction factors follow the pinned getter:
// additive `x*(1-w)` and multiplicative `x/(1+w*(x-1))` for weights below one
// and the mode default (`0`/`1`) elsewhere; a zero multiplicative denominator
// reproduces the pinned IEEE result. The corrected series is the effective
// value with non-NaN replacements substituted, exactly the pinned
// apply-corrections identity.
X11ReplacementPlan x11_replacements(const MaskedSeries& series,int frequency,
    const std::vector<double>& weights,int backcasts=0,int forecasts=0,
    bool exclude_projections=false,bool multiplicative=false,
    const Interrupt& interrupt={});
} // namespace kalos
#endif
