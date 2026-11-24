#ifndef KALOS_X13_FIT_HPP
#define KALOS_X13_FIT_HPP
// Item 22c M3-finale: the restricted fixed-model X-13 fit. This is pure
// composition of the qualified pieces — the fixed-coefficient preadjust
// (x13_preadjust, item 22c M3c), the pass-level X-11 selection run over the
// preadjusted series (item 22b.4 S4, all four decomposition modes of item
// 22b.5), and the finals restore (x13_restore, item 22c M3c) — plus the one
// pinned export stage the composition needs: under log_additive the pinned
// X11Kernel results builder (finalResults -> ratioBiasCorrection, the
// X11Spec default bias) rescales the exp'd step tables before X13Kernel.
// finals consumes them; that stage is delivered here as a composition step
// for the Ratio default only (docs/contracts/x13-allocation-scope-22c.md).
//
// The pinned X13Kernel gives no fixed-coefficient entry point (every run
// estimates the ARIMA part and the missing values), so the fit is qualified
// against captures of the full pinned run: validation/evidence/native/
// x13-allocation (zero-horizon seam cases, item 22c M3c) and
// validation/evidence/native/x13-fit (monthly fixed-specification
// end-to-end cases through the x13_fast entry point, jdplus.x13.base.r.
// X13.process). The capture-side identity guards recompute every stored
// table from the pinned pieces and refuse on disagreement, proving the
// fixed path composed here is the path the pinned engine executed.
//
// Deliberate scope, all of it fail-closed (docs/contracts/
// x13-fit-contract.md): coefficients are GIVEN, never estimated; the
// deterministic specification is the x13_allocation vocabulary; the X-11
// requests are fixed seasonal requests per period with a fixed Henderson
// trend length and a REQUIRED sigma policy; the decomposition mode is one
// of the four qualified modes, forced through the pinned X13Kernel.
// updateSpec rule exactly as in the allocation contract. The pinned
// fixed-path horizon domain is exactly {backcast 0, forecast 0}: any
// nonzero requested horizon (negative values included — the pinned kernel
// turns those into |h| * annual frequency) leaves the fixed path through
// the ExactArimaForecasts extension route (a1a/a1b and the extended out
// calendar), which is not delivered, so the fit answers the explicit
// unsupported status instead. Estimation, missing rows, the
// length-of-period transform, benchmarking and diagnostics stay outside.
#include "x11_selection_pass.hpp"
#include "x13_allocation.hpp"
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
namespace kalos {
// One restricted X-13 fit request. `input` is the complete series on the
// calendar (frequency, 1-based first_period); `variables` and `log` are the
// fixed deterministic specification of the allocation contract (pinned
// model order, one coefficient per pinned regression column). The X-11
// block mirrors the pinned X11Spec of the run: one FIXED seasonal request
// per calendar period (x11default and msr requests are outside the fixed
// contract), a fixed Henderson `trend_length` (an odd length >= 5; the
// pinned automatic 0 stays outside the fit), the REQUIRED sigma policy
// (the pinned spec carries no promoted default through this surface), the
// sigma cuttings (the pinned spec defaults 1.5/2.5), the MSR iteration
// policy (inert while every request is fixed), the decomposition mode (the
// pinned X11Spec default is additive), and the requested backcast/forecast
// extension (only the pinned fixed-path domain {0, 0} is delivered).
struct X13FitRequest {
    std::vector<double> input;
    std::vector<X13AllocationVariable> variables; // pinned model order
    bool log=false;
    int frequency=12;
    int first_period=1;
    std::vector<X11SeasonalRequest> seasonal;
    int trend_length=0;
    std::optional<X11SigmaPolicy> sigma; // required; no default
    double lower_sigma=1.5,upper_sigma=2.5;
    X11MsrIterationPolicy msr_policy=X11MsrIterationPolicy::jdplus;
    X11DecompositionMode mode=X11DecompositionMode::additive;
    int backcast=0,forecast=0;
};
// The result of one restricted fit. `status` is "available" or the verbatim
// unsupported reason: the fit boundary's own refusal (non-fixed seasonal
// request, automatic trend length, nonzero requested horizon, missing sigma
// policy) or the selection run's verbatim "unsupported:<stage>" passed
// through with the stage records preserved. Shape violations (an empty or
// ragged request, missing or non-positive-under-log observations, a
// calendar too short for the pinned passes) throw std::invalid_argument
// exactly as the seams below them do. Every provenance field is echoed
// whether or not the run completed; tables stay empty otherwise.
struct X13FitResult {
    std::string status;
    // Provenance: the request as given and the pinned forcing of it.
    int frequency=0,first_period=0;
    bool log=false;
    std::vector<X11SeasonalRequest> seasonal;
    int trend_length=0;
    bool sigma_set=false;
    X11SigmaPolicy sigma=X11SigmaPolicy::legacy;
    double lower_sigma=0.,upper_sigma=0.;
    X11MsrIterationPolicy msr_policy=X11MsrIterationPolicy::jdplus;
    X11DecompositionMode mode_requested=X11DecompositionMode::additive;
    X11DecompositionMode mode_effective=X11DecompositionMode::additive;
    int backcast=0,forecast=0;
    // The composition, in pinned execution order: the stored preadjustment
    // (a-tables and out), the finished X-11 selection run over out (full
    // pass accessibility: b/c/d tables, plans, trend and MSR stages), the
    // X11Results-equivalent tables exactly as finals consumed them, and
    // the pinned finals.
    X13Preadjustment preadjustment;
    X11SelectionPass selection;
    std::vector<double> d10,d11,d12,d13,c17;
    X13Finals finals;
    // The final components on the input calendar: the pinned X13Finals
    // triplet D11final/D12final/D13final (seasonally adjusted, trend and
    // irregular with the deterministic effects restored). The pinned
    // engine stores no separate final-trend or final-irregular table;
    // these are its final-components fields (docs/contracts/
    // x13-fit-contract.md).
    std::vector<double> finalSA,finalTrend,finalIrregular;
};
// The restricted fixed-model X-13 fit. Unsupported-domain requests (see
// the scope header) return the explicit status; the pinned mode forcing
// (pseudo_additive is preserved; otherwise non-log -> additive and log
// -> multiplicative unless log_additive) is applied to the requested mode and echoed.
X13FitResult x13_fit(const X13FitRequest& request);
} // namespace kalos
#endif
