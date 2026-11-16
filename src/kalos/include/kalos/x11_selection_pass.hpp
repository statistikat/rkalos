#ifndef KALOS_X11_SELECTION_PASS_HPP
#define KALOS_X11_SELECTION_PASS_HPP
// Item 22b.4 S4: the pass-level additive X-11 selection flow of the pinned
// X11Kernel — the staged seasonal requests resolved through the qualified
// item 22b.3d seasonal plan, the automatic Henderson trend choices at B7
// (pinned inline ratio >= 1 -> period + 1 branch), C7, D7 and the D final
// step through the qualified item 22b.3c comparison-ratio rule, and the
// global-MSR D-final seasonal decision through the qualified item 22b.3c
// band rule with the pinned attempt loop — composed onto the item 22b.4
// fixed-pass arithmetic. This slice decides nothing new: every selection
// rule is the qualified contract function applied to the pinned I/C ratios
// and MSR tables, and the pinned degenerate choices (a NaN comparison ratio
// taking the longest filter, a NaN MSR table selecting 3x9) stay explicit
// unsupported statuses instead of being reproduced.
// See docs/contracts/x11-selection-pass-contract.md and, for the
// decomposition modes, docs/contracts/x11-multiplicative-contract.md,
// x11-log-additive-contract.md and x11-pseudo-additive-contract.md.
// Cross-validation, non-seasonal modes, the period-specific corrector and
// ARIMA-regression inputs are not delivered or claimed.
#include "x11_d_pass.hpp"
#include "x11_msr_selection.hpp"
#include "x11_seasonal_plan.hpp"
#include <string>
#include <vector>
namespace kalos {
// The comparison-ratio outcome of one automatic Henderson stage. `fixed`
// carries the requested length and (except at the D final stage) no ratio:
// the pinned kernel computes calcICR at a stage only when automatic
// selection engages there. `unsupported` records the pinned ratio that the
// qualified rule refuses (NaN) with no length; the pass then stops with a
// matching status.
enum class X11TrendStageStatus { fixed, selected, unsupported };
struct X11TrendStage {
    X11TrendStageStatus status=X11TrendStageStatus::fixed;
    int length=0;
    double ic_ratio=std::numeric_limits<double>::quiet_NaN();
};
// The global-MSR D-final decision. `engaged` mirrors the pinned
// X11Context.isMSR over the D final stage (the qualified plan's
// requires_msr). The attempt rows and fields follow the qualified item
// 22b.3c selection semantics: evaluations executed (pinned d9msriter), the
// last evaluated global ratio (0 when none ran, the pinned field default),
// the method default flag and its reason. `unsupported` records an
// undefined (NaN) evaluated table, which selects no filter.
struct X11MsrStage {
    bool engaged=false;
    X11MsrIterationPolicy policy=X11MsrIterationPolicy::jdplus;
    X11MsrSelectionStatus status=X11MsrSelectionStatus::incomplete_table;
    X11MsrSelectionReason reason=X11MsrSelectionReason::unavailable;
    int selected_filter=0; // 3, 5 or 9 once decided; 0 while undecided
    bool using_default=false;
    std::size_t iterations=0;
    double global_msr=0.;
    std::vector<X11MsrAttempt> attempts;
};
// The finished selection run. b/c/d carry the tables of the selected state
// (the pinned getSeasFilter and getD9msriter semantics live in
// d.seas_filter and d.msr_iterations); initial_filters/final_filters are
// the resolved calendar-period arrays B, C and the D preliminary stages
// used; the seasonal plans are the qualified per-stage resolutions.
// `status` is "available" or "unsupported:<stage>" naming the first stage
// whose pinned choice the qualified rules refuse; unsupported runs keep
// the stage records and plans but no tables.
struct X11SelectionPass {
    std::string status;
    X11BPass b;
    X11CPass c;
    X11DPass d;
    X11SeasonalPlan initial_plan,final_plan;
    std::vector<int> initial_filters,final_filters;
    X11TrendStage b7,c7,d7,final;
    X11MsrStage msr;
};
// Selection run over `series` (the pinned B1 input, complete: the pinned
// X11Kernel.check refuses missing values). `requests` holds one pinned
// SeasonalFilterOption request per calendar period; `trend_length` is the
// fixed Henderson length or 0 for the pinned automatic selection. The sigma
// policy is required and nothing is promoted. The MSR iteration policy
// defaults to the pinned jdplus loop; the qualified x11_five_year_cap
// alternative stops after five year removals. `backcast`/`forecast` are the
// pinned context horizons of a hand-extended calendar, exactly as in the
// fixed passes. `mode` drives the whole flow (the pinned kernel shares one
// context): the multiplicative mode runs the ratio arithmetic, positivity
// repair and the pinned relative-form MSR table over the ratio-composed
// tables (item 22b.5 M3a) and requires strictly positive complete input;
// log_additive (item 22b.5 M3b-i) runs the additive flow on Math.log of the
// input under the same positivity requirement; pseudo_additive (item
// 22b.5 M3b-ii) runs the pinned X11{B,C,D}StepPseudoAdd pipeline — ratio
// arithmetic with the PseudoAddUtility composition forms, the makePositivity
// trend repair and the removeIrregular-form global-MSR walk with the
// absolute-form MSRE table — with no positivity requirement.
X11SelectionPass x11_selection_pass(const MaskedSeries& series,int frequency,int first_period,
    const std::vector<X11SeasonalRequest>& requests,int trend_length,
    X11SigmaPolicy policy,X11MsrIterationPolicy msr_policy=X11MsrIterationPolicy::jdplus,
    double lower_sigma=1.5,double upper_sigma=2.5,int backcast=0,int forecast=0,
    X11DecompositionMode mode=X11DecompositionMode::additive,
    const Interrupt& interrupt={});
} // namespace kalos
#endif
