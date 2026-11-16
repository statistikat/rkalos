#ifndef KALOS_X11_C_PASS_HPP
#define KALOS_X11_C_PASS_HPP
// Item 22b.4 S2: the complete fixed additive C pass of the default X-11
// decomposition, composed from the item 22b.4 S1 B pass exactly like the
// pinned X11Kernel: the final sigma weights and correction factors of the B
// pass correct the original series into C1, the final seasonal and Henderson
// stages re-run over it, and the qualified items 22b.3a/3b units perform the
// final extreme-value analysis on C13 with an explicit required policy.
// See docs/contracts/x11-cpass-contract.md.
// The D pass, automatic filter choice, cross-validation, multiplicative and
// non-seasonal modes and ARIMA-regression inputs are not delivered or
// claimed.
#include "x11_b_pass.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
namespace kalos {
// The final extreme-value analysis of the C pass on the C13 detection series:
// the policy-selected sigma weights and the replacement pass computed from
// them (correction history; the delivered C17/C20 tables are its weights and
// correction factors).
struct X11CCorrectionRound {
    MaskedSeries detection;
    std::int64_t analysis_position=0;
    X11SigmaPolicy policy=X11SigmaPolicy::legacy;
    X11SigmaWeights sigma;
    X11ReplacementPlan replacements;
};
// All pinned C tables. C2, C4 and C5a live on the central span (length
// n-2*c2drop), every remaining table on the original span (length n), exactly
// like the pinned getters. C5 is the extended normalization result (full
// span). Missing flags mirror NaN values.
struct X11CPass {
    std::size_t c2drop=0;
    MaskedSeries c1,c2,c4,c5,c5a,c6,c7,c9,c10,c10a,c11,c13,c17,c20;
    // The final-stage analysis on C13 (pinned IExtremeValuesCorrector state).
    X11CCorrectionRound round;
    std::vector<int> initial_filters,final_filters; // calendar-period order
    int trend_length=0;
    std::string status; // "available"
};
// Fixed additive C pass over `series` (the pinned B1 input) and its finished
// fixed additive B pass `b`. The composition mirrors the pinned X11Kernel:
// C1 = series - b.b20 and both passes must agree on frequency, first period,
// the fixed seasonal requests, the Henderson length, the sigma policy and
// the exact series the B pass consumed. As in the B pass the sigma policy is
// a required argument and nothing is promoted; the pinned fixture was
// qualified with X11SigmaPolicy::legacy. `first_period` is the 1-based
// calendar position of the first row.
X11CPass x11_c_pass(const MaskedSeries& series,int frequency,int first_period,
    const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
    int trend_length,X11SigmaPolicy policy,const X11BPass& b,
    double lower_sigma=1.5,double upper_sigma=2.5,const Interrupt& interrupt={});
} // namespace kalos
#endif
