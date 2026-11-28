#ifndef KALOS_X11_MSTATS_HPP
#define KALOS_X11_MSTATS_HPP
// P4 (X-11 diagnostic battery): the native X13Mstats block. This reproduces
// the pinned jdplus.x13.base.core.x11.Mstatistics battery (the m[] table,
// getM(1..11), getQ(), getQm2()) plus the MDiagnostics "q"/"q2" quality
// grades, computed over a finished X13FitResult exactly as the pinned
// X13Document publishes them on the automatic and estimated-template lanes.
//
// Every formula, loop bound, skip rule and operation order is transcribed
// from the pinned 3.9.0 bytecode (Mstatistics.of -> checkSeries ->
// calcStationaryVariances -> calcSNorm -> calcEvolutions -> calcM1..M11 ->
// getM/getQ/getQm2, with X11Utility.calcAbsMeanVariation(s)/adr and the
// CombinedSeasonality mvalue of M7). The input seam mirrors what the pinned
// Mstatistics consumes from X11Results, X11PreAdjustment and X13Finals:
// O = a1, P = a8 (the whole extended table), Sc/Cc/CIc/Ic = the A1-domain
// slices of D10/D12/D11/D13, Imod/Ome = the finals e3/e1, TD = the finals
// d18 (nulled when constant, the pinned all-equal fold), SI = the A1-domain
// slice of D8 (log-restored through the pinned exp/log pair under
// log_additive, where the step table is stored on the log scale), icr = the
// selection pass's final-stage I/C ratio (AutomaticHenderson.calcICR over
// d11bis) and the D9 record = the fit's D9 filter with the pinned
// X11Results.isD9default() flag, which the pinned X11DStep never writes and
// is therefore always false — the M6-in-q gate reduces to filter == S3X5.
//
// Published quality: the 13 battery words (m1..m11 published through the
// pinned getM clamp to [0, 3], q, qm2) and both grades are replayed against
// validation/evidence/native/x11-mstats (twelve double-JVM cases: the three
// stock RSA3 template fits and the RSA0/RSA1/RSA2c stock estimated fits on
// the x13-refusal-jvm-profiles calendars). The OLS-fed member (M2, with the
// two quality aggregates downstream of it) runs the pinned robust least
// squares (detail::pinned_ols_fit) on the [1 | 0..n-1] mean-correction
// design; the M7 F pair reuses the qualified CombinedSeasonality mirror.
//
// Scope notes. The pinned Mstatistics is refused on the pseudo_additive
// lane (the battery's P/TD/A8-side setup and the PseudoAddUtility forms are
// outside the qualified surface): the answer is the explicit status, never
// a partial battery. Omod/CImod, which the pinned class computes but the
// published battery never consumes, are not reproduced. The tree-level
// `.user_defined.quality.summary` document aggregate has no reproducible
// formula; the grade surface here is exactly the MDiagnostics q/q2 pair.
#include "x11_decomposition.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace kalos {

struct X13FitResult;

// The battery, in pinned publication order. `raw` is the m[] table (the
// values getQ/getQm2 consume, -1.0 for a member the pinned battery leaves
// uncomputed, 0.5*(icr-1) for M3 with its sign kept); `published` is the
// getM(i) view (raw clamped to [0, 3] with the pinned dcmpg/dcmpl semantics,
// so a negative M3 publishes as 0). q_grade/qm2_grade mirror the MDiagnostics
// getDiagnostic outcome ("Good", "Bad", "Severe" or "Undefined" at the pinned
// 1.0/2.0 default thresholds).
struct X11Mstats {
    std::string status; // "available", "unsupported:pseudo_additive" or "failed:trend_regression"
    X11DecompositionMode mode=X11DecompositionMode::additive; // echo of the fit
    std::array<double, 11> raw{{-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.}};
    std::array<bool, 11> used{{false,false,false,false,false,false,false,false,false,false,false}};
    std::array<double, 11> published{{-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.}};
    double q=std::numeric_limits<double>::quiet_NaN();   // getQ()
    double qm2=std::numeric_limits<double>::quiet_NaN(); // getQm2()
    std::string q_grade, qm2_grade;
};

// The battery over a finished fit. A fit whose status is not "available"
// throws std::invalid_argument, as do shape violations (tables that do not
// span the fit calendar, masked D8 rows under an available fit). The
// pseudo_additive mode answers "unsupported:pseudo_additive"; a degenerate
// trend regression (the pinned null LeastSquaresResults, response norm below
// Constants.epsilon) answers "failed:trend_regression" with the empty
// battery. Both refuse-answers leave every word at its untouched default.
X11Mstats x13_mstats(const X13FitResult& fit);

} // namespace kalos
#endif
