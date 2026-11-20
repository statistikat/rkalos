#ifndef KALOS_DETAIL_SPECIAL_FUNCTIONS_HPP
#define KALOS_DETAIL_SPECIAL_FUNCTIONS_HPP

// Operation-for-operation transcription of jdplus.toolkit.base.core.dstats.
// SpecialFunctions (pinned javap dump /tmp/opencode/jdprobe/item21h/
// SF_a02afe.txt, cross-checked in /tmp/opencode/kalos-21h-pinning.md): the
// incomplete-gamma pair behind the Chi-square upper-tail p-values consumed
// by the Ljung-Box statistic in the automatic model controller
// (docs/contracts/arma-diagnostics-contract.md). Math.log and Math.exp route
// through the fdlibm transcription (detail/fdlibm.hpp) so every operation
// matches the JVM operand order; the capture JVM runs with the StubRoutines
// exp/log intrinsics disabled, as for the item 21-D corpus.
//
// The JVM signals invalid arguments and non-convergent searches with
// DStatException; TestsUtility.pvalue catches every such exception and
// stores a NaN p-value. The C++ surface maps each throw site onto
// std::nullopt (chi2_upper_p maps the whole surface onto NaN), and callers
// own the distinction.

#include <optional>

namespace kalos::detail {

// SpecialFunctions.logGamma(a). The JVM throws for a <= 0 or NaN; nullopt
// there. For a > 0 the Lanczos-type kernel is total.
std::optional<double> special_log_gamma(double a);

// SpecialFunctions.IGSeries(a, x): scaled lower incomplete-gamma series,
// sum((n!)/(a+1)..(a+n+1) x^n) * exp(a*log(x) - x - logGamma(a)) form. The
// JVM guards a > 0 && x >= 0 and throws "Too many iterations" after
// m_maxiter = 1000 terms; both map onto nullopt.
std::optional<double> special_ig_series(double a, double x);

// SpecialFunctions.IGContFract(a, x): scaled upper incomplete-gamma
// continued fraction (modified Lentz). Same error surface as ig_series.
std::optional<double> special_ig_cont_fract(double a, double x);

// SpecialFunctions.incGamma(a, x): the regularized lower incomplete gamma;
// the series branch for x < a + 1, one minus the continued fraction
// otherwise (a NaN x takes the continued-fraction branch, as JVM dcmpg
// semantics prescribe). nullopt iff the taken branch fails its guard or
// iteration bound.
std::optional<double> special_inc_gamma(double a, double x);

// Chi2.getProbability(x, Upper) = 1 - incGamma(df / 2, x / 2) (Chi2.class
// pushes (x, df) onto SpecialFunctions.chiSquare(x, df), which guards
// x >= 0 and divides by 2 before calling incGamma). The TestsUtility.pvalue
// exception surface becomes NaN; df <= 0 is rejected earlier by the Chi2
// constructor and is not re-checked here.
double chi2_upper_p(double df, double x);

} // namespace kalos::detail

#endif
