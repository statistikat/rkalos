// Machine-precision transcription of jdplus.toolkit.base.core.dstats.
// SpecialFunctions (logGamma / IGSeries / IGContFract / incGamma / chiSquare)
// and of the Chi2 upper-tail probability used by TestsUtility.pvalue, from
// the pinned bytecode (/tmp/opencode/kalos-21h-pinning.md, section "LjungBox
// p-value machine"). Java guard semantics are reproduced literally: the
// incomplete-gamma pair uses dcmpg tests, so NaN arguments slip past the
// parameter guards and surface as the "Too many iterations" throw, which the
// nullopt return mirrors. Math.log/Math.exp map onto the fdlibm
// transcription, matching the intrinsic-disabled capture JVM.
#include "kalos/detail/special_functions.hpp"

#include "kalos/detail/fdlibm.hpp"

#include <cmath>

namespace kalos::detail {
namespace {

// SpecialFunctions static initializer: c_gammaconst (6 doubles) and
// m_maxiter = 1000. Nothing in the pinned chain calls setMaxIter.
constexpr double kGammaConst[6] = {76.18009172947146, -86.50532032941678,
                                   24.01409824083091, -1.231739572450155,
                                   0.001208650973866179, -5.395239384953e-6};
constexpr int kMaxIter = 1000;
constexpr double kGammaSeed = 1.000000000190015;   // IGSeries/logGamma seed
constexpr double kSqrt2Pi = 2.5066282746310007;    // log'ed per call
constexpr double kTiny = 1.0e-65;                  // Fortran-style tiny
constexpr double kIterFloor = 1.0e-16;             // convergence floors

// logGamma kernel body for a > 0 (guard handled by the caller).
double log_gamma_body(double a) {
  const double t = a + 5.5;
  // series = ((a + 0.5) * log(t)) - t
  double series = (a + 0.5) * fdlibm_log(t) - t;
  double tmp = kGammaSeed;
  for (int j = 0; j < 6; ++j) {
    tmp += kGammaConst[j] / (a + static_cast<double>(j + 1));
  }
  tmp = fdlibm_log(tmp);
  tmp = tmp + fdlibm_log(kSqrt2Pi);
  series = series + tmp;
  return series - fdlibm_log(a);
}

// exp((a * log(x) - x) - logGamma(a)); operand order per bytecode
// (dneg, mul, add, sub, exp) at IGSeries 87-106 and IGContFract 199-220.
std::optional<double> gamma_kernel(double a, double x) {
  const double lg = log_gamma_body(a);
  double t = a * fdlibm_log(x);
  t = t - x;
  t = t - lg;
  return fdlibm_exp(t);
}

}  // namespace

std::optional<double> special_log_gamma(double a) {
  // Guard: dcmpg ifgt - throws unless a > 0.0 (NaN included).
  if (!(a > 0.0)) {
    return std::nullopt;
  }
  return log_gamma_body(a);
}

std::optional<double> special_ig_series(double a, double x) {
  // Guard (dcmpg): throws iff a <= 0 or x < 0; NaN slips through and
  // defeats every convergence test, ending in the max-iteration throw.
  if (a <= 0.0 || x < 0.0) {
    return std::nullopt;
  }
  double del = 1.0 / a;
  double sum = del;
  double ap = a;
  for (int it = 0; it < kMaxIter; ++it) {
    ap = ap + 1.0;
    del = del * (x / ap);
    sum = sum + del;
    if (std::abs(del) < std::abs(sum) * kIterFloor) {
      return sum * gamma_kernel(a, x).value_or(std::nan(""));
    }
  }
  return std::nullopt;  // DStatException "Too many iterations"
}

std::optional<double> special_ig_cont_fract(double a, double x) {
  // Guard: identical dcmpg pattern to ig_series.
  if (a <= 0.0 || x < 0.0) {
    return std::nullopt;
  }
  double b = (x + 1.0) - a;
  double c = 1.0 / kTiny;
  double d = 1.0 / b;
  double h = d;
  bool converged = false;
  for (int i = 1; i < kMaxIter;) {
    const double ana = static_cast<double>(-i) * (static_cast<double>(i) - a);
    b = b + 2.0;
    d = ana * d + b;
    if (std::abs(d) < kTiny) {
      d = kTiny;
    }
    c = b + ana / c;
    if (std::abs(c) < kTiny) {
      c = kTiny;
    }
    d = 1.0 / d;
    const double del = d * c;
    h = h * del;
    // "if |del - 1| > 1e-16 continue; else converged" with JVM dcmpg: NaN
    // compares greater, so a NaN del keeps the loop turning until the
    // iteration bound (the throw), never converges.
    const double gap = std::abs(del - 1.0);
    if (std::isnan(gap) || gap > kIterFloor) {
      ++i;
      continue;
    }
    converged = true;
    break;
  }
  if (!converged) {
    return std::nullopt;  // DStatException "Too many iterations"
  }
  return gamma_kernel(a, x).value_or(std::nan("")) * h;
}

std::optional<double> special_inc_gamma(double a, double x) {
  // Branch test dcmpg ifge: x < a + 1 selects the series; NaN takes the
  // continued fraction.
  if (x < a + 1.0) {
    return special_ig_series(a, x);
  }
  const auto cf = special_ig_cont_fract(a, x);
  if (!cf) {
    return std::nullopt;
  }
  return 1.0 - *cf;
}

double chi2_upper_p(double df, double x) {
  // SpecialFunctions.chiSquare(x, df): dcmpg guard throws iff x < 0; the
  // division order df / 2, x / 2 follows the invocations at chiSquare 18-28.
  if (x < 0.0) {
    return std::nan("");
  }
  const auto lower = special_inc_gamma(df / 2.0, x / 2.0);
  if (!lower) {
    return std::nan("");  // TestsUtility.pvalue catch(Exception) -> NaN
  }
  // getProbability(x, Upper): 1 - lower.
  return 1.0 - *lower;
}

}  // namespace kalos::detail
