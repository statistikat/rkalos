#ifndef KALOS_X11_B_PASS_HPP
#define KALOS_X11_B_PASS_HPP
// Item 22b.4 S1: the complete fixed additive B pass of the default X-11
// decomposition, assembled from the qualified primitives of items 22a/22b.2
// (filters, seasonal processing, normalization) and 22b.3a/3b (sigma weights
// with an explicit required policy, replacements and correction factors),
// with the pinned fixed Henderson trend filter (local-polynomial weights and
// mmsre asymmetric endpoints) and the pinned seasonal processing conventions
// reproduced in the internal detail/x11_pinned.hpp machinery shared with the
// C pass (item 22b.4 S2).
// See docs/contracts/x11-bpass-contract.md.
// The C/D passes, automatic filter choice, cross-validation, multiplicative
// and non-seasonal modes, ARIMA-regression inputs and the pinned average
// endpoint fallback for series shorter than the trend filter are not
// delivered or claimed.
#include "x11_filters.hpp"
#include "x11_replacement.hpp"
#include "x11_seasonal_processing.hpp"
#include "x11_sigma.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <vector>
namespace kalos {
// One extreme-value correction round of the B pass: the pinned detection
// series, the 0-based pinned analysis start position, the policy-selected
// sigma weights and the replacement pass computed from them (correction
// history; the round's public tables are read from X11BPass).
struct X11BCorrectionRound {
    MaskedSeries detection;
    std::int64_t analysis_position=0;
    X11SigmaPolicy policy=X11SigmaPolicy::legacy;
    X11SigmaWeights sigma;
    X11ReplacementPlan replacements;
};
// All pinned B tables. B2..B4g live on the central span (length n-2*b2drop),
// every remaining table on the original span (length n), exactly like the
// pinned getters. B4/B9 are correction (replacement) tables: NaN where no
// replacement exists. B17/B20 are the final-round observation weights and
// correction factors. Missing flags mirror NaN values.
struct X11BPass {
    std::size_t b2drop=0;
    int first_period=0; // 1-based echo, checked by the C pass
    MaskedSeries b1,b2,b3,b4,b4a,b4anorm,b4d,b4g,b5,b6,b7,b8,b9,b9g,b10,b11,
        b13,b17,b20;
    // Rounds for the detection series b4d, b9d and b13 in pinned order.
    std::array<X11BCorrectionRound,3> rounds;
    std::vector<int> initial_filters,final_filters; // calendar-period order
    int trend_length=0;
    std::string status; // "available"
};
// Fixed additive B pass on the supplied series (the pinned B1 input).
// `initial_filters` and `final_filters` are one fixed request per calendar
// period (0 stable, N for 3xN), matching the pinned per-phase arrays;
// `trend_length` is the fixed Henderson length used from B7 onward and
// requires `series.values.size() >= trend_length` (below that the pinned
// reference replaces the endpoint region with a mean fallback that is not
// delivered). The sigma policy is a required argument and nothing is
// promoted: the pinned fixture was qualified with X11SigmaPolicy::legacy.
// `first_period` is the 1-based calendar position of the first row.
X11BPass x11_b_pass(const MaskedSeries& series,int frequency,int first_period,
    const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
    int trend_length,X11SigmaPolicy policy,
    double lower_sigma=1.5,double upper_sigma=2.5,const Interrupt& interrupt={});
} // namespace kalos
#endif
