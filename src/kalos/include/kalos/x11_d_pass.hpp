#ifndef KALOS_X11_D_PASS_HPP
#define KALOS_X11_D_PASS_HPP
// Item 22b.4 S3: the complete fixed additive D pass of the default X-11
// decomposition, composed from the item 22b.4 S2 C pass exactly like the
// pinned X11Kernel: the C-pass final correction factors correct the original
// series into D1, the final seasonal and Henderson stages re-run over the
// D tables and the D9 MSRE diagnostic is computed on the original calendar.
// See docs/contracts/x11-dpass-contract.md.
// The MSRE filter selection, the automatic Henderson (the pinned getICRatio
// scalar is inert in fixed mode), cross-validation, the period-specific
// corrector (the pinned D9g is null in the default corrector mode),
// non-seasonal modes and ARIMA-regression inputs are not delivered or
// claimed.
#include "x11_c_pass.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
namespace kalos {
// The pinned D9 MSRE diagnostic table (X11Utility.defaultMsrTable over the
// dropped original calendar at the original-calendar position): per calendar
// phase the scaled mean absolute evolutions of the 7-MA filtered-mean
// residual (irregular) and of the series itself (seasonal), the pinned
// stride difference counts, the per-phase msr ratio and the count-weighted
// global ratio. NaN entries follow the pinned plain-average semantics.
struct X11MsrDiagnostics {
    std::vector<double> mean_irregular_evolutions;
    std::vector<double> mean_seasonal_evolutions;
    std::vector<std::int64_t> counts;
    std::vector<double> msr;
    double global_msr=0.;
};
// All pinned D tables. D2 and D4 live on the central span (length
// n-2*d2drop), every remaining table on the series span (length n), exactly
// like the pinned getters. D11 is the final seasonally adjusted series over
// the reference series and D11bis over the corrected series; D12 is the
// final Henderson trend and D13 = D11 - D12. Missing flags mirror NaN.
struct X11DPass {
    std::size_t d2drop=0;
    int backcast=0,forecast=0; // extended-calendar echo; the original
        // calendar is the row span [backcast, n - forecast)
    MaskedSeries d1,d2,d4,d5,d6,d7,d8,d9,d9bis,d10,d10bis,d11,d11bis,d12,d13;
    X11MsrDiagnostics msr_table; // computed on the original calendar
    std::vector<int> initial_filters,final_filters; // calendar-period order
    std::vector<int> seas_filter; // the seasonal requests the D stage used
    int trend_length=0;
    int d7_henderson_length=0,final_henderson_length=0;
    int msr_iterations=0; // pinned constant in fixed mode: the MSRE filter
                          // selection never engages, so no filter is chosen
    X11DecompositionMode mode=X11DecompositionMode::additive;
    std::string status; // "available"
};
// Fixed D pass over `series` (the pinned reference series, spanning the
// extended calendar when backcast/forecast rows are present) and its
// finished fixed C pass `c`. The composition mirrors the pinned X11Kernel:
// D1 = the pinned ctx.remove of c.c20 from the series and both passes must
// agree on frequency, first period, the fixed seasonal requests, the
// Henderson length, the sigma policy and the exact series the C pass
// consumed. The decomposition mode is inherited from `c`; under the
// multiplicative mode the MSRE diagnostic takes the pinned relative form.
// The sigma policy is a required argument and nothing is promoted; the D
// pass itself performs no extreme-value analysis - its correction round is
// the C pass round carried by `c`. `first_period` is the 1-based calendar
// position of the first row of `series`; `backcast`/`forecast` are the
// pinned context horizons whose dropped span is the original calendar (at
// least one full cycle, used by the MSRE diagnostic exactly like the pinned
// table).
X11DPass x11_d_pass(const MaskedSeries& series,int frequency,int first_period,
    const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
    int trend_length,X11SigmaPolicy policy,const X11CPass& c,
    int backcast=0,int forecast=0,const Interrupt& interrupt={});
} // namespace kalos
#endif
