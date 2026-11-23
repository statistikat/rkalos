#pragma once

// Shared 24-A/24-B engine: the bit-exact mirror of the pinned jdplus 3.9.0
// validation + decomposition machinery. Extracted verbatim from
// src/tramo_seats_model_validator.cpp at item 24-B so that
// kalos::decompose_ucarima composes on the exact same code the 24-A
// validator runs; the authoritative bytecode notes live in that file and
// in tmp_opencode/jdprobe. Everything here is JVM-bytecode order.

#include "kalos/detail/hostlineage.hpp"
#include "kalos/tramo_seats_model_validator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace kalos {
namespace detail {


using vd = std::vector<double>;

constexpr double kPi = 3.141592653589793;
constexpr double kMax = std::numeric_limits<double>::max();
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();
constexpr double kMinValue = 4.9E-324;  // Double.MIN_VALUE (LU pivot seed)
constexpr double kLuEps = 1e-13;        // LS_Factory.ROBUST_FACTORY (LULinearSystemSolver2 eps)

[[noreturn]] inline void refuse(const std::string& why) {
  throw std::invalid_argument("tramo_seats_model: " + why);
}

// Java comparison flavors: cmpg pushes +1 on NaN, cmpl pushes -1.
inline bool gt_g(double a, double b) { return !(a <= b); }  // dcmpg + ifgt
inline bool ge_g(double a, double b) { return !(a < b); }   // dcmpg + ifge
inline bool lt_g(double a, double b) { return a < b; }      // dcmpg + iflt
inline bool ge_l(double a, double b) { return a >= b; }     // dcmpl + ifge
inline bool lt_l(double a, double b) { return !(a >= b); }  // dcmpl + iflt
inline bool gt_l(double a, double b) { return a > b; }      // dcmpl + ifgt

inline bool is_nan(double v) { return v != v; }

// ---------------------------------------------------------------------------
// Complex (api.math.Complex / ComplexType).
// ---------------------------------------------------------------------------
struct cplx {
  double re = 0.0, im = 0.0;
};

inline cplx c_conj(cplx a) { return {a.re, -a.im}; }
inline cplx c_neg(cplx a) { return {-a.re, -a.im}; }
inline cplx c_minus(cplx a, cplx b) { return {a.re - b.re, a.im - b.im}; }
inline cplx c_times(cplx a, cplx b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
inline cplx c_scale(cplx a, double s) {
  if (s == 1.0) return a;
  if (s == 0.0) return {0.0, 0.0};
  return {a.re * s, a.im * s};
}
inline cplx c_div_d(cplx a, double s) {
  if (s == 1.0) return a;
  return {a.re / s, a.im / s};
}
inline double c_absq(cplx a) { return a.re * a.re + a.im * a.im; }
inline double c_abs(double re, double im) {  // ComplexType.abs
  if (re == 0.0 && im == 0.0) return 0.0;
  double ar = std::fabs(re), ai = std::fabs(im);
  double mx = std::max(ar, ai), mn = std::min(ar, ai);
  if (mn == 0.0) return mx;
  double r = mn / mx;
  return mx * std::sqrt(1.0 + r * r);
}
inline double c_abs(cplx a) { return c_abs(a.re, a.im); }
// Raw libm, intentionally (CI-FDLIBM ledger, 2026-11-21): the JVM Spectrum
// filter evaluations call Math.cos/sin, whose pinned-JVM dispatch is the live
// C2 dcos/dsin stub — IBM-APML lineage, which glibc and the macOS libm both
// derive from (stub-vs-bytecode probe: 337/14753 cos and 399/9507 sin of the
// corpus-reached arguments differ from the FdLibm bytecode, while stub-vs-
// glibc differs at only 19/16; every qualified platform replays the corpora
// bit-exactly with host-libm cos/sin). Routing the vendored fdlibm_cos/
// fdlibm_sin to these seams would BREAK the captured corpora's bit lineage;
// the kernels stay vendored for their golden-proved role (see
// include/kalos/detail/fdlibm.hpp).
inline double c_arg(cplx a) { return kalos::detail::glibc_atan2(a.im, a.re); }

inline cplx c_inv(cplx a) {  // Complex.inv(): branches on |re| >= |im|
  if (!lt_g(std::fabs(a.re), std::fabs(a.im))) {
    double r = a.im / a.re;
    double d = 1.0 / (a.re + r * a.im);
    return {d, d * (-(a.im / a.re))};
  }
  double r = a.re / a.im;
  double d = 1.0 / (a.im + r * a.re);
  return {d * (a.re / a.im), -d};
}

// ---------------------------------------------------------------------------
// Coefficients.ofInternal / getUsedCoefficients (|c| <= 1e-9 trim from the top).
// ---------------------------------------------------------------------------
inline int used_coefficients(const vd& a, double eps) {
  int n = static_cast<int>(a.size());
  while (n > 0) {
    if (gt_g(std::fabs(a[static_cast<std::size_t>(n) - 1]), eps)) break;
    --n;
  }
  return n;
}

inline vd trim_copy(const vd& a) {
  int n = used_coefficients(a, 1e-9);
  if (n == 0) return {0.0};
  if (n == static_cast<int>(a.size())) return a;
  return vd(a.begin(), a.begin() + n);
}

// ---------------------------------------------------------------------------
// Polynomial. The coefficient array is kept EXACTLY as the JVM keeps it: the
// helpers that trim (of/ofInternal/plus/minus-with-flag) do, those that do not
// (times/derivate/divide(double)/valueOf via <init>) do not, because
// degree() is arrayLength - 1 in Java and callers depend on that.
// ---------------------------------------------------------------------------
struct Poly {
  vd c;
  Poly() : c{0.0} {}
  explicit Poly(vd a) : c(std::move(a)) {}
  int degree() const { return static_cast<int>(c.size()) - 1; }
  double get(std::size_t i) const { return c[i]; }
  bool is_zero() const { return c.size() == 1 && c[0] == 0.0; }
  bool is_identity() const { return c.size() == 1 && c[0] == 1.0; }
};

inline Poly p_of(vd a) { return Poly(trim_copy(a)); }   // Polynomial.of / ofInternal
inline Poly p_raw(vd a) { return Poly(std::move(a)); }  // new Polynomial(array)

inline Poly p_one() { return Poly(vd{1.0}); }

inline Poly p_value_of(double scale, const vd& params) {  // Polynomial.valueOf
  if (params.empty()) {
    if (scale == 0.0) return Poly(vd{0.0});
    if (scale == 1.0) return Poly(vd{1.0});
    return Poly(vd{scale});
  }
  vd a(params.size() + 1);
  a[0] = scale;
  for (std::size_t i = 0; i < params.size(); ++i) a[i + 1] = params[i];
  return p_raw(a);
}

// Polynomial.factor(double, int): 1 - x L^n (x == 0 -> ONE).
inline Poly p_factor(double x, int n) {
  if (x == 0.0) return p_one();
  vd a(static_cast<std::size_t>(n) + 1, 0.0);
  a[0] = 1.0;
  a[static_cast<std::size_t>(n)] = -x;
  return p_of(std::move(a));
}

inline vd pl_times(const vd& a, const vd& b) {  // Polynomials.times
  vd r(a.size() + b.size() - 1, 0.0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] == 0.0) continue;
    for (std::size_t j = 0; j < b.size(); ++j) {
      if (b[j] == 0.0) continue;
      r[i + j] += a[i] * b[j];
    }
  }
  return r;
}

inline Poly p_times(const Poly& a, const Poly& b) {  // times(P) / times(P, boolean)
  if (b.is_zero() || a.is_zero()) return Poly(vd{0.0});
  if (b.is_identity()) return a;
  if (a.is_identity()) return b;
  return p_raw(pl_times(a.c, b.c));
}

inline Poly p_times_d(const Poly& a, double d) {  // Polynomial.times(double)
  if (d == 1.0) return a;
  if (d == 0.0 || a.is_zero()) return Poly(vd{0.0});
  vd r(a.c.size());
  for (std::size_t i = 0; i < r.size(); ++i) r[i] = d * a.c[i];
  return p_raw(std::move(r));
}

inline vd pl_minus(const vd& a, const vd& b) {  // Polynomials.minus
  if (a.size() >= b.size()) {
    vd r = a;
    for (std::size_t i = 0; i < b.size(); ++i) r[i] -= b[i];
    return r;
  }
  vd r = b;
  for (std::size_t i = 0; i < a.size(); ++i) r[i] = a[i] - b[i];
  for (std::size_t i = a.size(); i < b.size(); ++i) r[i] = -b[i];
  return r;
}

inline Poly p_minus(const Poly& a, const Poly& b, bool trim) {
  vd r = pl_minus(a.c, b.c);
  return p_raw(trim ? trim_copy(r) : r);
}

inline Poly p_minus_p(const Poly& a, const Poly& b) {  // minus(P): trims iff equal degree
  return p_minus(a, b, a.degree() == b.degree());
}

inline Poly p_divide_d(const Poly& a, double d) { return p_times_d(a, 1.0 / d); }

inline double p_value_at(const Poly& p, double x) {  // Polynomial.evaluateAt (Horner)
  int n = p.degree();
  double v = p.c[static_cast<std::size_t>(n)];
  for (int i = n - 1; i >= 0; --i) v = v * x + p.c[static_cast<std::size_t>(i)];
  return v;
}

inline Poly p_derivate(const Poly& p) {  // Polynomial.derivate
  int n = p.degree();
  if (n == 0) return Poly(vd{0.0});
  vd r(static_cast<std::size_t>(n), 0.0);
  for (int i = 1; i <= n; ++i)
    r[static_cast<std::size_t>(i) - 1] =
        p.c[static_cast<std::size_t>(i)] * static_cast<double>(i);
  return p_raw(std::move(r));
}

// Polynomials.divide(double[],double[]) -> quotient (operands degree-trimmed at
// 1e-9; the quotient is used through Polynomial.ofInternal, i.e. trimmed).
inline vd pl_divide(const vd& a, const vd& b) {
  int n = static_cast<int>(a.size()) - 1;
  while (n >= 0 && !(ge_g(std::fabs(a[static_cast<std::size_t>(n)]), 1e-9))) --n;
  if (n < 0) return {0.0};  // Coefficients.zero()
  int m = static_cast<int>(b.size()) - 1;
  while (m >= 0 && !(ge_g(std::fabs(b[static_cast<std::size_t>(m)]), 1e-9))) --m;
  if (m < 0) return {kNan};  // Coefficients.nan()/negativeInfinity(): unused here
  if (m > n) return {0.0};
  vd ac = a;
  vd q(static_cast<std::size_t>(n - m + 1), 0.0);
  for (int k = n - m; k >= 0; --k) {
    q[static_cast<std::size_t>(k)] =
        ac[static_cast<std::size_t>(m + k)] / b[static_cast<std::size_t>(m)];
    for (int j = m + k - 1; j >= k; --j)
      ac[static_cast<std::size_t>(j)] -=
          q[static_cast<std::size_t>(k)] * b[static_cast<std::size_t>(j - k)];
  }
  return q;
}

inline Poly p_quotient(const Poly& a, const Poly& b) { return p_of(pl_divide(a.c, b.c)); }

// Polynomial.divide(Polynomial, Polynomial) -> Division(remainder, quotient).
struct Division {
  Poly remainder, quotient;
  bool exact() const { return remainder.is_zero(); }
};

inline Division p_divide(const Poly& p, const Poly& q) {
  int n = p.degree(), m = q.degree();
  if (m > n) return Division{p, Poly(vd{0.0})};  // Division(dividend, ZERO)
  std::size_t nn = static_cast<std::size_t>(n) + 1;
  vd ac = p.c, quot(nn, 0.0);
  for (int k = n - m; k >= 0; --k) {
    double v = ac[static_cast<std::size_t>(m + k)] / q.c[static_cast<std::size_t>(m)];
    quot[static_cast<std::size_t>(k)] = v;
    for (int j = m + k - 1; j >= k; --j)
      ac[static_cast<std::size_t>(j)] -= v * q.c[static_cast<std::size_t>(j - k)];
  }
  Poly qq = p_raw(vd(quot.begin(), quot.begin() + static_cast<std::size_t>(n - m + 1)));
  Poly rem = m > 0 ? p_raw(trim_copy(vd(ac.begin(), ac.begin() + static_cast<std::size_t>(m))))
                   : Poly(vd{0.0});
  return Division{rem, std::move(qq)};
}

// ---------------------------------------------------------------------------
// Polynomial.roots() -> FastEigenValuesSolver.factorize: the degree >= 3 qr
// cascade, ported bit-exact at item 24-A.1 (retires the D11 QR refusals).
// Authority: javap -p -c of the pinned 3.9.0 jar
// (tmp_opencode/jdprobe/m524/FastEigenValuesSolver.txt + Rotator.txt +
// Roots_Factory.txt), cross-checked against the v3.9.0 release source
// (ref/FastEigenValuesSolver-v3.9.0.java); Polynomial.roots() dispatches via
// RootsSolver.fastSolver() = Roots_Factory.FASTSOLVER = FastEigenValuesSolver
// ::new (pinned in Roots_Factory.txt). Every operation is a scalar double in
// bytecode order; the matrix products reproduce the PURE-JAVA GEMM loop
// (GEMM.txt apply/addAB + CPointer.set/mul/addAX: per column of C fill 0
// (beta == 0), then ascending inner index ii: s = alpha * B[ii, j], and for
// s != 0 col += s * A(:, ii) elementwise - mul then add, no FMA, and a zero
// scale skips the column entirely (so -0.0 accumulates differently than a
// +0.0 add)). tol = Constants.getEpsilon() = 1.1102230246251565E-16.
// ---------------------------------------------------------------------------

// A JVM exception raised INSIDE DefaultModelDecomposer.decompose's try block
// (bytecode 0-145, catch Exception -> null result). Two QR-cascade paths
// reach it: an it_max-exhausted cascade leaves the roots array null, so the
// consumer (minimizeByRoots / the selectors / fixMa) NPEs on the first null
// entry; a chase of the 2x2 tail block throws ArrayIndexOutOfBoundsException.
// The decomposer converts either into a NULL return (r.decompose and/or
// setVarianceMax null -> the "non decomposable model" event pair). The unit
// mirrors those as this exception, caught at the decompose() seam - they are
// captured DECISIONS, not envelope refusals.
struct jvm_decompose_abort {
  std::string reason;
};

constexpr double kFeEps = 1.1102230246251565E-16;  // Constants.getEpsilon()

inline cplx c_cart(double re, double im) {  // Complex.cart: every singleton branch
                                     // (I/ZERO/NEG_I/ONE/NEG_ONE/NEG_TWO) is
                                     // value-identical to new Complex(re, im)
                                     // except ZERO, which normalizes +-0.0 to
                                     // +0.0 - the only observable difference
  if (re == 0.0) {
    if (im == 1.0) return {0.0, 1.0};
    if (im == 0.0) return {0.0, 0.0};
    if (im == -1.0) return {0.0, -1.0};
  } else if (im == 0.0) {
    if (re == 1.0) return {1.0, 0.0};
    if (re == -1.0) return {-1.0, 0.0};
    if (re == -2.0) return {-2.0, 0.0};
  }
  return {re, im};
}
inline cplx c_cart(double re) { return c_cart(re, 0.0); }

struct fe_rot {
  double c = 1.0, s = 0.0;
};

inline bool fe_simplify(fe_rot& q, double eps) {  // Rotator.simplify(eps)
  if (gt_g(std::fabs(q.s), eps)) return false;  // |s| <= eps (dcmpg)
  q.s = 0.0;
  q.c = gt_l(q.c, 0.0) ? 1.0 : -1.0;            // c > 0 ? 1 : -1 (dcmpl)
  return true;
}

inline void fe_fuse(fe_rot& q1, fe_rot& q2, int type) {  // Rotator.fuse(q1, q2, type)
  double tmp = q1.c * q2.c - q1.s * q2.s;
  if (type == 0) {
    q1.s = q1.s * q2.c + q1.c * q2.s;
    q1.c = tmp;
  } else {
    q2.s = q1.s * q2.c + q1.c * q2.s;
    q2.c = tmp;
  }
}

inline double givens_rotation(double a, double b, fe_rot& gr);
inline void givens_rotation2(double a, double b, fe_rot& gr);

inline void fe_turnover(fe_rot& q1, fe_rot& q2, fe_rot& q3) {  // Rotator.turnover
  double c1 = q1.c, s1 = q1.s;
  double c2 = q2.c, s2 = q2.s;
  double c3 = q3.c, s3 = q3.s;
  double a = s1 * c3 + c1 * c2 * s3;
  double b = s2 * s3;
  fe_rot q;
  double nrm = givens_rotation(a, b, q);
  double c4 = q.c, s4 = q.s;
  a = c1 * c3 - s1 * c2 * s3;
  b = nrm;
  givens_rotation2(a, b, q);
  double c5 = q.c, s5 = q.s;
  double t1 = s1 * s2, t2 = -c1 * s2, t3 = c2;
  a = -t2 * s4 + t3 * c4;
  t2 = t2 * c4 + t3 * s4;
  b = -(-t1 * s5 + t2 * c5);
  givens_rotation2(a, b, q);
  double c6 = q.c, s6 = q.s;
  q1.c = c5;
  q1.s = s5;
  q2.c = c6;
  q2.s = s6;
  q3.c = c4;
  q3.s = s4;
}

inline double givens_rotation(double a, double b, fe_rot& gr) {
  if (b == 0.0) {           // dcmpl ifne: a NaN b falls through to the math
    if (a < 0.0) {          // dcmpg ifge: a NaN takes the +1 branch
      gr.c = -1.0;
      gr.s = 0.0;
      return -a;
    }
    gr.c = 1.0;
    gr.s = 0.0;
    return a;
  }
  double absa = std::fabs(a), absb = std::fabs(b);
  double s, c, r;
  if (absa >= absb) {  // dcmpl iflt: NaN goes to the b-major branch
    s = b / a;
    r = std::sqrt(1.0 + s * s);
    if (a < 0.0) {  // dcmpg
      c = -1.0 / r;
      s *= c;
      r *= -a;
    } else {
      c = 1.0 / r;
      s *= c;
      r *= a;
    }
  } else {
    c = a / b;
    r = std::sqrt(1.0 + c * c);
    if (b < 0.0) {  // dcmpg
      s = -1.0 / r;
      c *= s;
      r *= -b;
    } else {
      s = 1.0 / r;
      c *= s;
      r *= b;
    }
  }
  gr.c = c;
  gr.s = s;
  return r;
}

inline void givens_rotation2(double a, double b, fe_rot& gr) {
  if (b == 0.0) {  // dcmpl ifne; the a < 0 test is dcmpg (NaN -> +1)
    gr.c = a < 0.0 ? -1.0 : 1.0;
    gr.s = 0.0;
    return;
  }
  double absa = std::fabs(a), absb = std::fabs(b);
  double s, c, r;
  if (absa >= absb) {  // dcmpl iflt
    s = b / a;
    r = std::sqrt(1.0 + s * s);
    if (a < 0.0) {  // dcmpg
      c = -1.0 / r;
      s *= c;
    } else {
      c = 1.0 / r;
      s *= c;
    }
  } else {
    c = a / b;
    r = std::sqrt(1.0 + c * c);
    if (b < 0.0) {  // dcmpg
      s = -1.0 / r;
      c *= s;
    } else {
      s = 1.0 / r;
      c *= s;
    }
  }
  gr.c = c;
  gr.s = s;
}

// GEMM.addAB kernel over column-major dense storage (see the block comment).
inline void fe_gemm(double* c0, int ldc, int nr, int nc, const double* A, int lda, int k,
             const double* B, int ldb, double alpha) {
  for (int j = 0; j < nc; ++j) {
    double* col = c0 + static_cast<std::size_t>(j) * ldc;
    for (int i = 0; i < nr; ++i) col[i] = 0.0;  // beta == 0 (all calls)
    for (int ii = 0; ii < k; ++ii) {
      double s = alpha * B[static_cast<std::size_t>(j) * ldb + static_cast<std::size_t>(ii)];
      if (s == 0.0) continue;  // dcmpl ifeq: a zero scale skips the column
      const double* a = A + static_cast<std::size_t>(ii) * lda;
      for (int i = 0; i < nr; ++i) col[i] += s * a[i];
    }
  }
}

// 2x2 eigenValues (bytecode order); D is column-major [m00, m10, m01, m11].
inline void fe_eigen_values(const double* M, cplx* ev, int pos) {
  double m00 = M[0], m10 = M[1], m01 = M[2], m11 = M[3];
  double trace = m00 + m11;
  double detm = m00 * m11 - m10 * m01;
  double disc = trace * trace - 4.0 * detm;
  if (disc < 0.0) {  // dcmpg ifge: NaN -> the real branch
    cplx v = c_cart(trace / 2.0, std::sqrt(-disc) / 2.0);
    ev[pos] = v;
    ev[pos + 1] = c_conj(v);
  } else {
    double sdisc = std::sqrt(disc);
    double qp = std::fabs(trace + sdisc), qm = std::fabs(trace - sdisc);
    if (gt_l(qp, qm)) {  // dcmpl ifle
      double re = (trace + sdisc) / 2.0;
      ev[pos] = c_cart(re);
      ev[pos + 1] = c_cart(detm / re);
    } else if (gt_l(qm, 0.0)) {  // dcmpl ifle
      double re = (trace - sdisc) / 2.0;
      ev[pos] = c_cart(re);
      ev[pos + 1] = c_cart(detm / re);
    } else {
      ev[pos] = cplx{0.0, 0.0};
      ev[pos + 1] = cplx{0.0, 0.0};
    }
  }
}

// The qr cascade itself; solve() returns the roots array exactly as
// factorize()/roots() would expose it. Two JVM paths have no bit-level
// consumer inside the unit because the JVM itself aborts there - both are
// mirrored as jvm_decompose_abort (caught by DefaultModelDecomposer's
// catch(Exception) -> null decision, see decompose()):
//   * it_max (30n) exhaustion: the JVM leaves null roots (consumer NPE);
//   * chase() with stop == 0, which indexes B[-1] (ArrayIndexOutOfBounds).
// Measured on the corpus: q_airline/q_d2 hit exhaustion on the degree-4
// minimizer numerator (JVM and unit agree on the full solver profile,
// stop 2 / itc 120 / chase 120). The bulge() random-shift branch (it_count
// % 25 == 0) is NOT an abort: the JVM is not bit-reproducible there; see the
// substitution note at bulge - the unit pins the observable decisions.
struct fe_solver {
  int n = 0;
  std::vector<cplx> roots;
  std::vector<bool> set;
  std::vector<fe_rot> Q, C, B;
  int start = 0, stop = 0, zero = 0, it_max = 0, it_count = 0, chase_count = 0, tr = 0;
  std::vector<int> its;
  double D[4] = {0.0, 0.0, 0.0, 0.0};
  cplx shifts[2];
  fe_rot B1, B2, B3, tmp, B1b, B2b;

  void initialize(const Poly& p) {
    Q.assign(static_cast<std::size_t>(n), fe_rot{0.0, 1.0});
    Q[static_cast<std::size_t>(n) - 1] = fe_rot{1.0, 0.0};
    C.assign(static_cast<std::size_t>(n), fe_rot{1.0, 0.0});
    B.assign(static_cast<std::size_t>(n), fe_rot{1.0, 0.0});
    double u = n % 2 == 0 ? -1.0 : 1.0, v = -u;
    double pn = p.get(static_cast<std::size_t>(n));
    fe_rot c;
    double r = givens_rotation((v * p.get(0)) / pn, u, c);
    C[static_cast<std::size_t>(n) - 1] = c;
    B[static_cast<std::size_t>(n) - 1] = fe_rot{v * c.s, v * c.c};
    for (int i = n - 1; i > 0; --i) {
      fe_rot cc;
      r = givens_rotation((-p.get(static_cast<std::size_t>(i))) / pn, r, cc);
      C[static_cast<std::size_t>(i) - 1] = cc;
      B[static_cast<std::size_t>(i) - 1] = fe_rot{cc.c, -cc.s};
    }
  }

  void diagonal_block(int k) {
    D[0] = D[1] = D[2] = D[3] = 0.0;  // D.set(0)
    const fe_rot* q = Q.data();
    const fe_rot* c = C.data();
    const fe_rot* b = B.data();
    if (k == 0) {
      double R[4];
      R[0] = -b[0].s / c[0].s;
      double r11 = -b[1].s / c[1].s;
      R[2] = -(b[0].c * b[1].c - r11 * c[0].c * c[1].c) / c[0].s;
      R[1] = 0.0;
      R[3] = r11 * q[1].c;
      double A[4] = {q[0].c, q[0].s, -q[0].s, q[0].c};
      fe_gemm(D, 2, 2, 2, A, 2, 2, R, 2, 1.0);  // GeneralMatrix.aAB_p_bC(1, A, R, 0, D)
    } else {
      std::size_t k0 = static_cast<std::size_t>(k);
      double R[6];  // 3x2 column-major
      double r10 = -b[k0].s / c[k0].s;
      R[1] = r10;
      R[0] = -(b[k0 - 1].c * b[k0].c - r10 * c[k0 - 1].c * c[k0].c) / c[k0 - 1].s;
      double r21 = -b[k0 + 1].s / c[k0 + 1].s;
      R[2] = 0.0;
      R[3] = (b[k0 - 1].c * b[k0].s * b[k0 + 1].c -
              c[k0 - 1].c * (c[k0].c * b[k0].c * b[k0 + 1].c - c[k0 + 1].c * r21) / c[k0].s) /
             c[k0 - 1].s;
      R[4] = -(b[k0].c * b[k0 + 1].c - r21 * c[k0].c * c[k0 + 1].c) / c[k0].s;
      R[5] = r21 * q[k0 + 1].c;
      double A[4];
      double M[4];
      A[0] = q[k0].c;
      A[1] = q[k0].s;
      A[2] = -q[k0].s;
      A[3] = q[k0].c;  // R12 = R.extract(1, 2, 0, 2): rows 1-2, colInc 3
      fe_gemm(M, 2, 2, 2, A, 2, 2, R + 1, 3, 1.0);  // R12.copy(GeneralMatrix.AB(A, R12))
      R[1] = M[0];
      R[2] = M[1];
      R[4] = M[2];
      R[5] = M[3];
      A[0] = q[k0 - 1].c;
      A[1] = q[k0 - 1].s;
      A[2] = -q[k0 - 1].s;
      A[3] = q[k0 - 1].c;  // R01 = R.extract(0, 2, 0, 2) reads the UPDATED row 1
      fe_gemm(M, 2, 2, 2, A, 2, 2, R, 3, 1.0);
      R[0] = M[0];
      R[1] = M[1];
      R[3] = M[2];
      R[4] = M[3];
      D[0] = R[1];  // D.copy(R12)
      D[1] = R[2];
      D[2] = R[4];
      D[3] = R[5];
    }
  }

  void check_deflation() {
    for (int i = stop; i >= 0; --i) {
      if (fe_simplify(Q[static_cast<std::size_t>(i)], kFeEps)) {
        zero = i;
        start = i + 1;
        its[static_cast<std::size_t>(zero)] = it_count;
        it_count = 0;
        return;
      }
    }
  }

  bool update_roots() {
    if (stop == zero) {
      diagonal_block(stop);
      if (stop == 0) {
        roots[0] = c_cart(D[0]);
        roots[1] = c_cart(D[3]);
        set[0] = set[1] = true;
        stop = -1;
      } else {
        roots[static_cast<std::size_t>(stop) + 1] = c_cart(D[3]);
        set[static_cast<std::size_t>(stop) + 1] = true;
        stop--;
        zero = -1;
        start = 0;
      }
      return true;
    } else if (stop - 1 == zero) {
      diagonal_block(stop);
      fe_eigen_values(D, roots.data(), stop);
      set[static_cast<std::size_t>(stop)] = true;
      set[static_cast<std::size_t>(stop) + 1] = true;
      if (stop == 1) {
        diagonal_block(0);
        roots[0] = c_cart(D[0]);
        set[0] = true;
      }
      stop -= 2;
      zero = -1;
      start = 0;
      return true;
    }
    return false;
  }

  void bulge(const char* what, int deg) {
    // JVM branch (pinned javap bulge 0-93): for it_count % 25 == 0 the
    // shifts come from new Random().nextGaussian() twice - a time-based
    // seed. Measured (run32.sh + repeats): the JVM returns DIFFERENT root
    // bits in every process for such polynomials; the corpus enters this
    // branch exactly once - v_simplify's degree-12 stationary AR, JVM run
    // of 30 consecutive bulges (its[] shows it) - and the captured corpus
    // is nevertheless byte-stable because those roots only feed the
    // selector's partition decisions, whose gates sit far from the roots
    // (all |z| = 1.105 against the unit-circle boundary; SECTION 9).
    // Reproducing the JVM bits there is impossible BY DEFINITION; the unit
    // substitutes the eigenvalue shift (the other branch), deterministically
    // reproducing every captured decision. Measured on all corpus and stress
    // polynomials, the eigenvalue branch alone never changes a decision.
    (void)what;
    (void)deg;
    if (it_count % 25 == 0) {
      // SUBSTITUTED SHIFT (see comment): a fixed generic pair, the same
      // trigger as the JVM's random restart. An eigenvalue-shift
      // substitution was measured to STALL on v_simplify's degree-12 AR
      // (the exact-shift step degenerates; that is why the JVM restarts);
      // any generic pair breaks the degeneracy. The pair must never move a
      // decision gate (measured on the corpus: the branch fires only where
      // gates sit far from the roots).
      shifts[0] = c_cart(0.5, 0.5);
      shifts[1] = c_cart(0.5, -0.5);
    } else {
      diagonal_block(stop);
      fe_eigen_values(D, shifts, 0);
    }
    diagonal_block(start);
    double t00 = D[0], t10 = D[1], t01 = D[2], t11 = D[3];
    diagonal_block(start + 1);
    double t21 = D[1];
    double c0 = t00 * t00 + t01 * t10 + shifts[0].re * shifts[1].re -
                shifts[0].im * shifts[1].im - t00 * (shifts[0].re + shifts[1].re);
    double c1 = t10 * (t00 + t11 - shifts[0].re - shifts[1].re);
    double c2 = t10 * t21;
    double r = givens_rotation(c1, c2, B1);
    givens_rotation2(c0, r, B2);
    chase(what, deg);
  }

  void chase(const char* what, int deg) {
    chase_count++;
    if (start == 0) {
      tmp = fe_rot{B2.c, -B2.s};       // tmp.copyInverse(B2)
      B3 = fe_rot{B1.c, -B1.s};        // B3.copyInverse(B1)
      fe_turnover(tmp, B3, Q[0]);
      fe_fuse(B3, Q[1], 1);
      B3 = Q[0];
      Q[0] = tmp;
    } else {
      std::size_t s0 = static_cast<std::size_t>(start);
      tmp.c = B2.c;
      tmp.s = -Q[s0 - 1].c * B2.s;
      B3 = fe_rot{B1.c, -B1.s};
      fe_turnover(tmp, B3, Q[s0]);
      fe_fuse(B3, Q[s0 + 1], 1);
      B3 = Q[s0];
      Q[s0] = tmp;
    }
    for (int i = start; i <= stop - 2; ++i) {
      std::size_t i0 = static_cast<std::size_t>(i);
      if (i < tr - 1) {
        B1b = B1;
        B2b = B2;
        fe_turnover(B[i0 + 1], B[i0 + 2], B1b);
        fe_turnover(B[i0], B[i0 + 1], B2b);
        C[i0] = fe_rot{B[i0].c, -B[i0].s};
        C[i0 + 1] = fe_rot{B[i0 + 1].c, -B[i0 + 1].s};
        C[i0 + 2] = fe_rot{B[i0 + 2].c, -B[i0 + 2].s};
      } else {
        fe_turnover(B[i0 + 1], B[i0 + 2], B1);
        fe_turnover(B[i0], B[i0 + 1], B2);
        fe_turnover(C[i0 + 2], C[i0 + 1], B1);
        fe_turnover(C[i0 + 1], C[i0], B2);
      }
      fe_turnover(Q[i0 + 1], Q[i0 + 2], B1);
      fe_turnover(Q[i0], Q[i0 + 1], B2);
      fe_turnover(B3, B1, B2);
      tmp = B2;  // Rotator.copy
      B2 = B3;
      B3 = B1;
      B1 = tmp;
    }
    std::size_t s0 = static_cast<std::size_t>(stop);
    if (stop < 1)  // JVM indexes B[-1]: ArrayIndexOutOfBoundsException, caught
      throw jvm_decompose_abort{"QR cascade chased a 2x2 tail block (" + std::string(what) +
                                ", degree " + std::to_string(deg) + ")"};
    fe_turnover(B[s0], B[s0 + 1], B1);
    fe_turnover(B[s0 - 1], B[s0], B2);
    fe_turnover(C[s0 + 1], C[s0], B1);
    fe_turnover(C[s0], C[s0 - 1], B2);
    B1.s *= Q[s0 + 1].c;
    fe_fuse(Q[s0], B1, 0);
    fe_turnover(Q[s0 - 1], Q[s0], B2);
    fe_fuse(B3, B2, 0);
    fe_turnover(B[s0], B[s0 + 1], B3);
    fe_turnover(C[s0 + 1], C[s0], B3);
    B3.s *= Q[s0 + 1].c;
    fe_fuse(Q[s0], B3, 0);
    tr -= 2;
  }

  void qr(const Poly& p, const char* what) {
    initialize(p);
    tr = n - 3;
    start = 0;
    stop = n - 2;
    zero = -1;
    it_max = 30 * n;
    it_count = 0;
    chase_count = 0;
    its.assign(static_cast<std::size_t>(n), 0);
    for (int k = 0; k < it_max; ++k) {
      if (stop < 0) {
        return;
      }
      check_deflation();
      if (!update_roots()) {
        it_count++;
        bulge(what, n);
      }
    }
    exhausted = true;  // JVM: roots stay (partly) null; the consumer NPEs
    return;
  }

  bool exhausted = false;

  std::vector<cplx> solve(const Poly& p, const char* what) {
    n = p.degree();
    exhausted = false;
    roots.assign(static_cast<std::size_t>(n), cplx{0.0, 0.0});
    set.assign(static_cast<std::size_t>(n), false);
    qr(p, what);
    if (exhausted)
      throw jvm_decompose_abort{"QR cascade exhausted it_max (" + std::string(what) +
                                ", degree " + std::to_string(n) +
                                "): JVM null roots -> NPE, caught by the decomposer"};
    for (int i = 0; i < n; ++i)
      if (!set[static_cast<std::size_t>(i)])
        refuse(std::string("QR cascade left an unset root (") + what + ", degree " +
               std::to_string(n) + ")");
    return roots;
  }
};

// Polynomial.roots() (fastSolver dispatch). degree 0 -> empty; 1-2 use the
// monic/quadratic closed forms; >= 3 the qr cascade above.
inline std::vector<cplx> p_roots(const Poly& p, const char* what) {
  int n = p.degree();
  if (n >= 3) {
    fe_solver s;
    std::vector<cplx> out = s.solve(p, what);
    return out;
  }
  std::vector<cplx> out;
  if (n <= 0) return out;
  if (n == 1) {  // monic
    out.push_back({-(p.c[0] / p.c[1]), 0.0});
    return out;
  }
  double a = p.c[2], b = p.c[1], cc = p.c[0];
  double d8 = 2.0 * a;
  double disc = b * b - 4.0 * a * cc;
  if (lt_g(disc, 0.0)) {
    double sq = std::sqrt(-disc);
    cplx r0{-b / d8, sq / d8};
    out.push_back(r0);
    out.push_back(c_conj(r0));
  } else {
    double sq = std::sqrt(disc);
    out.push_back({(-b + sq) / d8, 0.0});
    out.push_back({(-b - sq) / d8, 0.0});
  }
  return out;
}

// Polynomial.fromComplexRoots(Complex[], scale).
inline Poly p_from_roots(const std::vector<cplx>& roots, double scale = 1.0) {
  if (roots.empty()) return p_raw(vd{scale});
  std::size_t n = roots.size();
  std::vector<cplx> z(n + 1);
  z[0] = {1.0, 0.0};
  for (std::size_t i = 0; i < n; ++i) {
    z[i + 1] = z[i];  // leading slot (the copy IS the B^(i+1) term)
    for (std::size_t j = i; j >= 1; --j)
      z[j] = c_minus(z[j - 1], c_times(z[j], roots[i]));
    z[0] = c_times(z[0], c_neg(roots[i]));
  }
  double c0 = z[0].re;
  if (c0 == 0.0) refuse("fromComplexRoots with a zero root");
  double k = scale / c0;
  vd a(n + 1);
  for (std::size_t i = 0; i <= n; ++i) a[i] = z[i].re * k;
  a[0] = scale;
  return p_raw(std::move(a));
}

// ComplexUtility.lejaOrder.
inline void leja_order(std::vector<cplx>& z) {
  int n = static_cast<int>(z.size());
  if (n == 0) return;
  double best = 0.0;
  int idx = 0;
  for (int i = 0; i < n; ++i) {
    double s = c_absq(z[static_cast<std::size_t>(i)]);
    if (lt_g(best, s)) {
      best = s;
      idx = i;
    }
  }
  std::swap(z[0], z[static_cast<std::size_t>(idx)]);
  for (int i = 1; i < n; ++i) {
    double sc = 0.0;
    int pick = i;
    for (int j = i; j < n; ++j) {
      double d = 1.0;
      for (int k = 0; k < i; ++k) {
        double re = z[static_cast<std::size_t>(j)].re - z[static_cast<std::size_t>(k)].re;
        double im = z[static_cast<std::size_t>(j)].im - z[static_cast<std::size_t>(k)].im;
        d = d * (re * re + im * im);
      }
      if (lt_g(sc, d)) {
        sc = d;
        pick = j;
      }
    }
    std::swap(z[static_cast<std::size_t>(i)], z[static_cast<std::size_t>(pick)]);
  }
}

// ---------------------------------------------------------------------------
// UnitRoots: D(d) = 1 - L^d, D(d, n) = D(d)^n, S(s, n): 1 at powers k*n, k < s.
// SarimaUtility.differencingFilter(s, d, bd) = D(1, d) * D(s, bd).
// ---------------------------------------------------------------------------
inline Poly ur_d(int d) {
  vd a(static_cast<std::size_t>(d) + 1, 0.0);
  a[0] = 1.0;
  a[static_cast<std::size_t>(d)] = -1.0;
  return p_of(std::move(a));
}

inline Poly ur_d_n(int d, int n) {
  if (n == 0) return p_one();
  Poly base = ur_d(d), r = base;
  for (int i = 1; i < n; ++i) r = p_times(r, base);
  return r;
}

inline Poly ur_s(int s, int n) {
  vd a(static_cast<std::size_t>(n) * (static_cast<std::size_t>(s) - 1) + 1, 0.0);
  for (int i = 0; i < s; ++i) a[static_cast<std::size_t>(i) * n] = 1.0;
  return p_of(std::move(a));
}

inline Poly differencing_filter(int s, int d, int bd) {
  bool have = false;
  Poly acc;
  if (d > 0) {
    acc = ur_d_n(1, d);
    have = true;
  }
  if (bd > 0) {
    Poly sd = ur_d_n(s, bd);
    acc = have ? p_times(acc, sd) : sd;
    have = true;
  }
  return have ? acc : p_one();
}

// ---------------------------------------------------------------------------
// SymmetricFilter: a trimmed polynomial a[0..n] with response
// a[0] + sum_{k>=1} 2 a[k] cos(k x). BackFilter: one-sided coefficients.
// ---------------------------------------------------------------------------
struct SymF {
  vd a;
  SymF() : a{0.0} {}
  explicit SymF(vd v) : a(std::move(v)) {}
  int length() const { return static_cast<int>(a.size()); }
  int upper_bound() const { return static_cast<int>(a.size()) - 1; }
};

inline SymF sym_of(vd raw) { return SymF(trim_copy(raw)); }  // SymmetricFilter.ofInternal

// SymmetricFilter.convolutionOf(IFiniteFilter, double): autocorrelation, scaled.
inline SymF sym_convolution(const vd& w_in, double scale) {
  // BackFilter.weights() indexes polynomial.get(-k) with lower bound -degree,
  // so convolutionOf()'s input is the coefficient vector REVERSED; per-lag
  // products therefore accumulate in descending order (bit-visible on
  // non-palindromic inputs).
  vd w(w_in.rbegin(), w_in.rend());
  std::size_t n = w.size();
  vd r(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i; j < n; ++j) r[j - i] += w[i] * w[j];
  if (scale != 1.0)
    for (std::size_t i = 0; i < n; ++i) r[i] *= scale;
  return sym_of(std::move(r));
}

// The pinned JVM uses the x86-64 HotSpot cosine kernel for frequency
// responses and their derivatives. The separately qualified sine path
// retains its original host lineage.
inline double rfr(const SymF& f, double x) {  // SymmetricFilter.realFrequencyResponse
  std::size_t stop = f.a.size();
  std::size_t k = 0;
  double u = f.a[k++];
  if (k >= stop) return u;
  double t0 = 1.0;
  double cs = kalos::detail::hotspot_spectrum_cos(x);
  double t1 = cs;
  for (;;) {
    u = u + 2.0 * t1 * f.a[k++];
    if (k >= stop) break;
    double tn = 2.0 * cs * t1 - t0;
    t0 = t1;
    t1 = tn;
  }
  return u;
}

// Spectrum$fr.evaluate: f, f', f'' via the angle recurrences.
inline void fr_eval(const SymF& f, double x, double& fv, double& dfv, double& d2fv) {
  std::size_t n = f.a.size();
  fv = f.a[0];
  dfv = 0.0;
  d2fv = 0.0;
  if (n <= 1) return;
  double cs = kalos::detail::hotspot_spectrum_cos(x), sn = kalos::detail::glibc_sin(x);
  double ck = cs, sk = sn;
  for (std::size_t k = 1; k < n; ++k) {
    double ak = f.a[k];
    double tc = 2.0 * ck * ak;
    double ts = 2.0 * sk * ak;
    double cn = ck * cs - sk * sn;
    double sx = sk * cs + ck * sn;
    ck = cn;
    sk = sx;
    fv += tc;
    dfv -= static_cast<double>(k) * ts;
    d2fv -= static_cast<double>(static_cast<int>(k) * static_cast<int>(k)) * tc;
  }
}

// Spectrum$dfr.evaluate (the k^d-weighted one-sided sum).
inline double dfr_eval(const SymF& f, int d, double x) {
  std::size_t n = f.a.size();
  double s = 0.0;
  for (std::size_t k = 1; k < n; ++k) {
    double t = static_cast<double>(k);
    for (int i = 1; i < d; ++i) t *= static_cast<double>(k);
    double trig = (d % 2 == 0) ? kalos::detail::hotspot_spectrum_cos(x * static_cast<double>(k))
                               : kalos::detail::glibc_sin(x * static_cast<double>(k));
    t = t * trig * f.a[k];
    s += t;
  }
  return s;
}

struct Spectrum {
  SymF num, den;
};

inline double spectrum_value(const Spectrum& s, double x) {  // Spectrum.value
  double den = rfr(s.den, x);
  double num = rfr(s.num, x);
  if (gt_l(std::fabs(den), 1e-7)) return num / den;
  if (ge_g(std::fabs(num), 1e-7)) return kNan;
  for (int k = 1; k <= 10; ++k) {
    double dd = dfr_eval(s.den, k, x);
    double nn = dfr_eval(s.num, k, x);
    if (gt_l(std::fabs(dd), 1e-7)) return nn / dd;
    if (!ge_l(std::fabs(nn), 1e-7)) break;
  }
  return kNan;
}

// ---------------------------------------------------------------------------
// SymmetricFrequencyResponse (minimizeByRoots only). m_p =
// Polynomial.of(D2SFR(coefficients)); the transform rows are the Chebyshev T_k
// coefficients (exact integers; ltransform's long table holds the same values).
// ---------------------------------------------------------------------------
struct SFR {
  Poly p;
};

inline SFR sfr_of(const SymF& f) {
  std::size_t n = f.a.size();
  vd b = f.a;
  for (std::size_t k = 1; k < n; ++k) b[k] *= 2.0;  // D2SFR doubling
  if (n > 2) {
    std::vector<vd> t(n);
    t[0] = vd(n, 0.0);
    t[0][0] = 1.0;
    if (n > 1) {
      t[1] = vd(n, 0.0);
      t[1][1] = 1.0;
      for (std::size_t k = 2; k < n; ++k) {
        t[k] = vd(n, 0.0);
        for (std::size_t j = 1; j <= k; ++j) t[k][j] += 2.0 * t[k - 1][j - 1];
        for (std::size_t j = 0; j + 2 <= k; ++j) t[k][j] -= t[k - 2][j];
      }
    }
    vd r(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
      if (b[k] == 0.0) continue;  // Ux: a zero source element is skipped
      // Ux computes x' = U^T x; reading the upper part of the (column-major)
      // transform matrix, U[j][k] = t[k][j] - the cos(k.w) coefficients live
      // at powers j <= k, so the accumulation runs DOWN from k to 0.
      for (std::size_t j = 0; j <= k; ++j) r[j] += t[k][j] * b[k];
    }
    b = std::move(r);
  }
  return SFR{p_of(std::move(b))};
}

inline double sfr_integral(const SFR& s) {  // getIntegral
  double v = s.p.c[0];
  double f = 1.0;
  for (int i = 2; i <= s.p.degree(); i += 2) {
    f *= static_cast<double>(i - 1);
    f /= static_cast<double>(i);
    v += s.p.c[static_cast<std::size_t>(i)] * f;
  }
  return v;
}

// ---------------------------------------------------------------------------
// Spectrum$Minimizer (evaluate / minimizeByRoots / minimizeByGrid + the
// bracketed Newton SpectrumFunction).
// ---------------------------------------------------------------------------
inline double minimizer_evaluate(const Poly& P, const Poly& Q, double x, int depth) {
  if (depth > 64) refuse("minimizer deflation depth exceeded");
  if (lt_g(std::fabs(x), 1e-12)) x = 0.0;
  double p = p_value_at(P, x);
  double q = p_value_at(Q, x);
  if (gt_l(std::fabs(q), 1e-12)) return p / q;
  if (gt_l(std::fabs(p), 1e-12)) return kNan;
  if (x != 0.0) {
    Poly lin = p_of(vd{1.0, -x});  // Polynomial.of({1, -x})
    return minimizer_evaluate(p_quotient(P, lin), p_quotient(Q, lin), x, depth + 1);
  }
  int n = std::min(P.degree(), Q.degree());
  for (int k = 1; k <= n; ++k) {
    double a = P.c[static_cast<std::size_t>(k)];
    double b = Q.c[static_cast<std::size_t>(k)];
    bool ta = gt_l(std::fabs(a), 1e-7);
    bool tb = gt_l(std::fabs(b), 1e-7);
    if (ta && tb) return a / b;
    if (ta) return kNan;
    if (tb) return 0.0;
  }
  return kNan;
}

struct MinResult {
  double min = 0.0;
  double x = 0.0;
};

inline MinResult minimize_by_roots(const Spectrum& s) {
  SFR sn = sfr_of(s.num), sd = sfr_of(s.den);
  double iden = sfr_integral(sd);
  Poly P = p_divide_d(sn.p, iden);
  Poly Q = p_divide_d(sd.p, iden);
  MinResult r;
  if (P.is_zero()) return r;
  if (P.degree() == 0 && Q.degree() == 0) {
    r.min = P.c[0] / Q.c[0];
    return r;
  }
  r.x = 0.0;
  r.min = kMax;
  double v = minimizer_evaluate(P, Q, 0.0, 0);
  if (!is_nan(v)) {
    r.min = v;
    r.x = kPi / 2.0;
  }
  v = minimizer_evaluate(P, Q, 1.0, 0);
  if (!is_nan(v) && lt_g(v, r.min)) {
    r.min = v;
    r.x = 0.0;
  }
  v = minimizer_evaluate(P, Q, -1.0, 0);
  if (!is_nan(v) && lt_g(v, r.min)) {
    r.min = v;
    r.x = kPi;
  }
  Poly w = (Q.degree() > 0) ? p_minus_p(p_times(p_derivate(P), Q), p_times(P, p_derivate(Q)))
                            : p_derivate(P);
  for (const cplx& z : p_roots(w, "minimizer numerator")) {
    if (ge_g(std::fabs(z.im), 1e-7)) continue;
    double re = z.re;
    if (!(gt_l(re, -1.0)) || !(lt_g(re, 1.0))) continue;
    double e = minimizer_evaluate(P, Q, re, 0);
    if (!is_nan(e) && lt_g(e, r.min)) {
      r.min = e;
      r.x = kalos::detail::glibc_acos(re);
    }
  }
  return r;
}

struct FnInstance {
  double pt = 0.0, f = 0.0, df = 0.0, d2f = 0.0;
};

inline FnInstance fn_instance(const Spectrum& s, double pt) {  // SpectrumFunctionInstance
  FnInstance out;
  out.pt = pt;
  double nf, ndf, nd2f, df, ddf, dd2f;
  fr_eval(s.num, pt, nf, ndf, nd2f);
  fr_eval(s.den, pt, df, ddf, dd2f);
  out.f = nf / df;
  out.df = (ndf * df - nf * ddf) / (df * df);
  out.d2f = ((nd2f * df - nf * dd2f) * df - 2.0 * (ndf * df - nf * ddf) * ddf) / (df * df * df);
  return out;
}

// SpectrumFunction.min(double): false means the JVM returned null.
inline bool fn_min(const Spectrum& s, double a, double b, double x0, FnInstance& out) {
  double den = rfr(s.den, x0);
  double x = x0;
  if (lt_g(std::fabs(den), 1e-7)) {
    for (;;) {
      x += 1e-7;
      den = rfr(s.den, x);
      if (gt_g(x, b)) break;
      if (!lt_g(std::fabs(den), 1e-7)) break;
    }
    if (!(x <= b)) return false;
  }
  FnInstance inst = fn_instance(s, x);
  double fmin = inst.f;
  double xn = x, xprev = x;
  int iter = 0;
  for (;;) {
    xprev = xn;
    xn = xn - inst.df / inst.d2f;
    if (is_nan(xn)) break;
    if (lt_g(xn, a))
      xn = a;
    else if (gt_l(xn, b))
      xn = b;
    FnInstance cand = fn_instance(s, xn);
    if (lt_g(cand.f, fmin)) {
      inst = cand;
      fmin = cand.f;
    }
    ++iter;
    if (iter >= 200) break;
    if (!(gt_l(std::fabs(xn - xprev), 1e-12))) break;
  }
  out = inst;
  return true;
}

inline MinResult minimize_by_grid(const Spectrum& s) {
  MinResult r;
  if (s.num.length() == 1 && s.den.length() == 1) {
    r.min = spectrum_value(s, 0.0);
    return r;
  }
  r.x = 0.0;
  r.min = kMax;
  double v = spectrum_value(s, 0.0);
  if (!is_nan(v)) {
    r.min = v;
    r.x = 0.0;
  }
  v = spectrum_value(s, kPi);
  if (!is_nan(v) && lt_g(v, r.min)) {
    r.min = v;
    r.x = kPi;
  }
  int k = s.num.upper_bound() + s.den.upper_bound() - 1;
  double step = kPi / static_cast<double>(k);
  double h = step / 2.0;
  for (int i = 0; i < k; ++i) {
    double x1 = h + step;
    double xl = h;
    double den = rfr(s.den, xl);
    while (!gt_g(den, 0.0) && lt_g(xl, x1)) {
      xl += step / 7.0;
      den = rfr(s.den, xl);
    }
    bool ok = lt_g(xl, x1);
    double xr = x1;
    if (ok) {
      den = rfr(s.den, x1);
      while (!gt_g(den, 0.0) && gt_g(xr, xl)) {
        xr -= step / 7.0;
        den = rfr(s.den, xr);
      }
      ok = gt_g(xr, xl);
    }
    if (ok) {
      FnInstance inst;
      if (fn_min(s, xl, xr, (xl + xr) / 2.0, inst)) {
        if (lt_g(inst.f, r.min)) {
          r.min = inst.f;
          r.x = inst.pt;
        }
      }
    }
    h += step;
  }
  return r;
}

inline MinResult minimize_spectrum(const Spectrum& s) {
  // JVM dispatch: denom.length() <= 6 where SymmetricFilter.length() =
  // 2*degree+1, i.e. the stored half-coefficient count is <= 3.
  if (s.den.length() <= 3) return minimize_by_roots(s);
  return minimize_by_grid(s);
}

// ---------------------------------------------------------------------------
// IRootSelector implementations (toolkit ucarima).
// ---------------------------------------------------------------------------
struct Selector {
  Poly selected, not_selected;
  bool has_selected = true, has_not_selected = true;
  virtual ~Selector() = default;
  virtual bool accept(const cplx& z) const = 0;
  // AbstractRootSelector.select(Polynomial)
  bool select(const Poly& p, const char* what) {
    if (p.degree() == 0) {
      has_selected = false;
      not_selected = p;
      has_not_selected = true;
      return false;
    }
    std::vector<cplx> roots = p_roots(p, what);
    std::vector<cplx> sel, ns;
    for (const cplx& z : roots) (accept(z) ? sel : ns).push_back(z);
    if (sel.empty()) {
      has_selected = false;
      not_selected = p;
      has_not_selected = true;
      return false;
    }
    if (ns.empty()) {
      has_not_selected = false;
      selected = p;
      has_selected = true;
      return true;
    }
    leja_order(sel);
    leja_order(ns);
    Poly s = p_from_roots(sel);
    selected = p_times_d(s, p.get(0) / s.get(0));
    has_selected = true;
    Poly n = p_from_roots(ns);
    not_selected = p_divide_d(n, n.get(0));
    has_not_selected = true;
    return true;
  }
  virtual bool select_unit_roots(const Poly& p) = 0;
};

// TrendCycleSelector: |1/z| bound + low-frequency window on |arg(1/z)|;
// unit roots: greedily extract the (1 - L) factors.
struct TrendCycleSelector : Selector {
  double bound, lfreq;
  explicit TrendCycleSelector(double b) : bound(b), lfreq(kPi / 12.0 - 1e-6) {}
  void set_default_low_freq_threshold(int period) {
    lfreq = kPi / static_cast<double>(period) - 1e-6;
  }
  bool accept(const cplx& z) const override {
    cplx zi = c_inv(z);
    if (z.im == 0.0) return !lt_l(zi.re, bound);
    if (lt_l(c_abs(zi), bound)) return false;
    return !gt_g(std::fabs(c_arg(zi)), lfreq);
  }
  bool select_unit_roots(const Poly& p) override {
    selected = p_one();
    has_selected = true;
    not_selected = p;
    has_not_selected = true;
    Poly d(vd{1.0, -1.0});  // UnitRoots.D1
    for (;;) {
      Division div = p_divide(not_selected, d);
      if (!div.exact()) break;
      selected = p_times(selected, d);
      not_selected = div.quotient;
      if (p.degree() < d.degree()) break;
    }
    return selected.degree() > 0;
  }
};

// SeasonalSelector: roots near 2 pi i / frequency; unit roots: the S(s,1) factors.
struct SeasonalSelector : Selector {
  double k = 0.5, epsphi = 2.0;
  int frequency;
  SeasonalSelector(int freq, double eps) : epsphi(eps), frequency(freq) {}
  void set_k(double v) { k = v; }
  bool accept(const cplx& z) const override {
    if (lt_g(std::fabs(z.im), 1e-6)) return lt_g(1.0 / z.re, -k);
    double w = 6.283185307179586 / static_cast<double>(frequency);
    double a = std::fabs(c_arg(z));
    double tol = epsphi / 180.0 * kPi;
    for (int i = 1; i <= frequency / 2; ++i)
      if (!gt_g(std::fabs(w * static_cast<double>(i) - a), tol)) return true;
    return false;
  }
  bool select_unit_roots(const Poly& p) override {
    selected = p_one();
    has_selected = true;
    not_selected = p;
    has_not_selected = true;
    if (frequency == 1) return false;
    Poly s = ur_s(frequency, 1);
    for (;;) {
      Division div = p_divide(not_selected, s);
      if (!div.exact()) break;
      selected = p_times(selected, s);
      not_selected = div.quotient;
      if (p.degree() < s.degree()) break;
    }
    return selected.degree() > 0;
  }
};

// ---------------------------------------------------------------------------
// ArimaModel (structural part only: stationary AR, unit roots, symmetric MA).
// Two shapes mirror the two constructors: the symmetric form (ar, delta, sma)
// and the one-sided form (ar, delta, ma, var) whose symmetric MA is
// convolutionOf(ma, var). The one-sided MA / variance getters of the symmetric
// form need the 24-B spectral factorization and are not modelled.
// ---------------------------------------------------------------------------
struct Arima {
  Poly ar, delta;
  bool symmetric = true;
  SymF sma;   // symmetric form
  vd ma1;     // one-sided form
  double var = 0.0;

  Poly stationary_ar() const { return ar; }
  Poly non_stationary_ar() const { return delta; }
  Poly full_ar() const { return p_times(ar, delta); }  // ArimaModel.getAr()

  SymF symmetric_ma() const {
    if (symmetric) return sma;
    // AbstractArimaModel.symmetricMa() = convolutionOf(ma, var) at the
    // ArimaModel's OWN innovation variance (0.0 for the checkSpecialCases
    // signal, 1.0 for every model built on the captured path).
    return sym_convolution(ma1, var);
  }
  Spectrum spectrum() const {  // AbstractArimaModel.getSpectrum()
    Spectrum s;
    s.num = symmetric_ma();
    s.den = sym_convolution(full_ar().c, 1.0);  // symmetricAr() = convolutionOf(getAr())
    return s;
  }
};

inline Arima arima_symmetric(Poly ar, Poly delta, SymF sma) {
  Arima a;
  a.ar = std::move(ar);
  a.delta = std::move(delta);
  a.symmetric = true;
  a.sma = std::move(sma);
  return a;
}

inline Arima arima_onesided(Poly ar, Poly delta, Poly ma, double var) {
  Arima a;
  a.ar = std::move(ar);
  a.delta = std::move(delta);
  a.symmetric = false;
  a.ma1 = std::move(ma.c);
  a.var = var;
  return a;
}

// ---------------------------------------------------------------------------
// LULinearSystemSolver2 (robustSolver). The flat working array is column-major
// with column stride n (tsar storage order): index(i, j) = i + j * n.
// Returns false for the JVM's MatrixException paths.
// ---------------------------------------------------------------------------
inline bool lu_solve(int n, int l, vd& a, double eps, std::vector<int>& piv) {
  auto at = [n](std::size_t row, std::size_t col) {
    return row + col * static_cast<std::size_t>(n);
  };
  piv.assign(static_cast<std::size_t>(l), -1);
  std::vector<int>& ipiv = piv;
  for (int k = 0; k < n; ++k) {
    double best = kMinValue;
    int ip = -1;
    for (int i = 0, idx = k; i < n; ++i, idx += n) {
      if (ipiv[static_cast<std::size_t>(i)] != -1) continue;
      double v = a[static_cast<std::size_t>(idx)];
      double av = std::fabs(v);
      if (gt_l(av, best)) {  // dcmpl + ifle: a NaN candidate is not taken
        ip = i;
        best = av;
      }
    }
    if (ip < 0) return false;  // MatrixException()
    ipiv[static_cast<std::size_t>(ip)] = k;
    double piv = a[at(static_cast<std::size_t>(k), static_cast<std::size_t>(ip))];
    if (!ge_g(std::fabs(piv), eps)) return false;  // MatrixException("m_err_sing")
    double d = 1.0 / piv;
    for (int j = 0; j < n; ++j) {
      if (j == k) continue;
      double v = a[at(static_cast<std::size_t>(j), static_cast<std::size_t>(ip))];
      if (v == 0.0) continue;
      double m = d * v;
      for (int c = 0, ic = k, jc = j; c < l; ++c, ic += n, jc += n) {
        if (ipiv[static_cast<std::size_t>(c)] != -1) continue;
        a[static_cast<std::size_t>(jc)] -= m * a[static_cast<std::size_t>(ic)];
      }
    }
    for (int c = 0, ic = k; c < l; ++c, ic += n) {
      if (ipiv[static_cast<std::size_t>(c)] != -1) continue;
      a[static_cast<std::size_t>(ic)] *= d;
    }
  }
  return true;
}

// ---------------------------------------------------------------------------
// RootDecomposer.calc().
// ---------------------------------------------------------------------------
struct DecomposeFail {
  bool threw = false;
};

// One RootDecomposer pass. `fail` is set when the JVM would throw (singular LU
// / out-of-range assembly), which DefaultModelDecomposer turns into null.
inline bool root_decompose(const Arima& model, Selector& sel, Arima& signal, Arima& noise,
                    bool& fail) {
  fail = false;
  // splitRoots()
  Poly m_sar = p_one(), m_nar = p_one(), m_sur = p_one(), m_nur = p_one();
  sel.select(model.stationary_ar(), "stationary AR");
  if (sel.has_selected) m_sar = sel.selected;
  if (sel.has_not_selected) m_nar = sel.not_selected;
  sel.select_unit_roots(model.non_stationary_ar());
  if (sel.has_selected) m_sur = sel.selected;
  if (sel.has_not_selected) m_nur = sel.not_selected;

  SymF m_sma = model.symmetric_ma();
  SymF m_sfds = sym_convolution(p_times(m_sur, m_sar).c, 1.0);
  SymF m_sfdn = sym_convolution(p_times(m_nur, m_nar).c, 1.0);

  if (m_sfds.length() == 1) {  // checkSpecialCases(): the whole MA goes to the noise
    signal = arima_onesided(p_one(), m_sur, p_one(), 0.0);
    noise = arima_symmetric(m_nar, m_nur, m_sma);
    return true;
  }

  vd a = m_sma.a;
  vd b = m_sfds.a;
  vd c = m_sfdn.a;  // m_sfdn is never null in this transcription

  int na = static_cast<int>(a.size()) - 1;
  int nb = static_cast<int>(b.size()) + static_cast<int>(c.size()) - 2;
  int ub = static_cast<int>(b.size()) - 1;
  int uc = static_cast<int>(c.size()) - 1;
  int n1 = ub > 0 ? ub - 1 : 0;
  int n2 = uc > 0 ? uc - 1 : 0;
  if (na >= nb) n2 = na - ub;
  int d1 = n1 + 1;
  int d2 = n2 + 1;
  int dim = d1 + d2;
  if (dim <= 0) {
    fail = true;
    return false;
  }
  std::size_t N = static_cast<std::size_t>(dim);
  int l = dim + 1;
  vd mat(N * static_cast<std::size_t>(l), 0.0);
  auto set = [&](int row, int col, double v) {
    if (row < 0 || col < 0 || row >= dim || col >= dim) {
      fail = true;
      return;
    }
    mat[static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * N] = v;
  };
  auto get = [&](int row, int col) -> double {
    if (row < 0 || col < 0 || row >= dim || col >= dim) {
      fail = true;
      return 0.0;
    }
    return mat[static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * N];
  };
  a[0] /= 2.0;
  b[0] /= 2.0;
  c[0] /= 2.0;
  for (int i = 0; i <= n1 && !fail; ++i)
    for (int j = 0; j <= uc && !fail; ++j) set(i + j, i, c[static_cast<std::size_t>(j)]);
  for (int i = 0; i <= n2 && !fail; ++i)
    for (int j = 0; j <= ub && !fail; ++j) set(i + j, d1 + i, b[static_cast<std::size_t>(j)]);
  for (int i = 0; i <= n1 && !fail; ++i)
    for (int j = 0; j <= uc && !fail; ++j) {
      int k = i - j;
      if (k < 0) k = -k;
      set(k, i, get(k, i) + c[static_cast<std::size_t>(j)]);
    }
  for (int i = 0; i <= n2 && !fail; ++i)
    for (int j = 0; j <= ub && !fail; ++j) {
      int k = i - j;
      if (k < 0) k = -k;
      set(k, d1 + i, get(k, d1 + i) + b[static_cast<std::size_t>(j)]);
    }
  vd rhs(N, 0.0);
  for (int i = 0; i <= na; ++i) {
    if (i >= dim || static_cast<std::size_t>(i) >= a.size()) {
      fail = true;
      break;
    }
    rhs[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i)];
  }
  if (fail) return false;
  for (std::size_t i = 0; i < N; ++i)
    mat[i + N * N] = rhs[i];  // X copied after the n*n block
  std::vector<int> ipiv;
  if (!lu_solve(dim, l, mat, kLuEps, ipiv)) return false;
  // solve() read-back: X(i, c) = a[ipiv[i] + (n + c) * n].
  vd sol(N, 0.0);
  for (int i = 0; i < dim; ++i)
    sol[static_cast<std::size_t>(i)] =
        mat[static_cast<std::size_t>(ipiv[static_cast<std::size_t>(i)]) + N * N];
  vd x1(static_cast<std::size_t>(d1), 0.0), x2(static_cast<std::size_t>(d2), 0.0);
  for (int i = 0; i < d1; ++i) x1[static_cast<std::size_t>(i)] = sol[static_cast<std::size_t>(i)];
  for (int i = 0; i < d2; ++i) x2[static_cast<std::size_t>(i)] = sol[static_cast<std::size_t>(d1 + i)];
  x1[0] *= 2.0;  // RootDecomposer: get(0) doubled back (a[0] was halved)
  x2[0] *= 2.0;
  signal = arima_symmetric(m_sar, m_sur, sym_of(std::move(x1)));
  noise = arima_symmetric(m_nar, m_nur, sym_of(std::move(x2)));
  return true;
}

// ---------------------------------------------------------------------------
// ModelDecomposer.decompose(IArimaModel): comps[i] = stage i's signal, last
// slot = the remaining noise; an exception anywhere -> null (fail).
// ---------------------------------------------------------------------------
inline bool model_decompose(const Arima& model, std::vector<std::unique_ptr<Selector>>& sels,
                     std::vector<Arima>& comps, bool& fail) {
  fail = false;
  comps.clear();
  Arima cur = model;
  for (std::size_t i = 0; i < sels.size(); ++i) {
    Arima sig, noi;
    if (!root_decompose(cur, *sels[i], sig, noi, fail)) return false;
    comps.push_back(sig);
    cur = noi;
  }
  comps.push_back(cur);
  return true;
}

// ---------------------------------------------------------------------------
// UcarimaModel.setVarianceMax(-1, noisy): the decomposability decision. Only
// the minima sum and its sign are mirrored; the components' minus()/plus()
// need the 24-B spectral factorization, so the allocated values are not built.
// ---------------------------------------------------------------------------
struct VarianceMax {
  bool ok = false;    // a UcarimaModel came back (vsum >= 0, or noisy)
  double sum = 0.0;   // the minima sum (setVarianceMax's "variance floor")
  int count = 0;      // components count of the resulting UcarimaModel
};

inline VarianceMax variance_max(const std::vector<Arima>& raw, bool noisy, bool chain_threw) {
  VarianceMax vm;
  if (chain_threw) return vm;  // decompose -> null -> setVarianceMax(null) -> null
  vm.count = static_cast<int>(raw.size()) + 1;  // idx < 0 appends ArimaModel.NULL
  double sum = 0.0;
  for (const Arima& c : raw) {
    MinResult m = minimize_spectrum(c.spectrum());
    if (m.min == 0.0) continue;  // zero-minimum components are skipped
    sum += m.min;
  }
  vm.sum = sum;
  vm.ok = noisy || !lt_g(sum, 0.0);
  return vm;
}

// ---------------------------------------------------------------------------
// SarimaModel + SarimaModel.Builder (parameter vectors define the orders).
// ---------------------------------------------------------------------------
struct Sarima {
  int period = 12, d = 0, bd = 0;
  vd phi, bphi, theta, btheta;

  int p() const { return static_cast<int>(phi.size()); }
  int q() const { return static_cast<int>(theta.size()); }
  int bp() const { return static_cast<int>(bphi.size()); }
  int bq() const { return static_cast<int>(btheta.size()); }

  // SarimaModel.getRegularAR/getRegularMA: valueOf(1, params) (no trim: the
  // constructor stores the raw arrays; ofInternal trims only trailing 1e-9).
  Poly regular_ar() const { return p_value_of(1.0, phi); }
  Poly regular_ma() const { return p_value_of(1.0, theta); }
  // SarimaModel.seasonalAR: sparse at stride `period`; bp == 1 is the
  // factor(-bphi[0], period) special case. For bp >= 2 the pinned bytecode
  // loops `for (i = period, j = 0; i < bphi.length; i += period, ++j)
  // a[i] = bphi[j]` - the bound compares the STRIDE index with the parameter
  // count, so the body never runs at period 4/12 and ofInternal trims the
  // result to ONE: the stationary seasonal AR silently DROPS (measured,
  // Scratch30: SarimaModel(0,1,0)(2,1,2).getStationaryAr() == 1 + phi.B).
  // Reproduced verbatim (24-A.1, pinning notes SECTION 9).
  Poly seasonal_ar() const {
    if (bphi.empty()) return p_one();
    if (bphi.size() == 1) return p_factor(-bphi[0], period);
    vd a(static_cast<std::size_t>(bphi.size()) * static_cast<std::size_t>(period) + 1, 0.0);
    a[0] = 1.0;
    for (int i = period, j = 0; i < static_cast<int>(bphi.size()); i += period, ++j)
      a[static_cast<std::size_t>(i)] = bphi[static_cast<std::size_t>(j)];
    return p_of(std::move(a));
  }
  // SarimaModel.seasonalMA: same structure and the same stride-vs-count loop
  // (the sign is kept for bq == 1's factor(-btheta[0], period)); bq >= 2
  // silently drops to ONE on the pinned JVM (v_cut_bq2's whole-MA-to-noise
  // decomposition is that path; Scratch30 ground truth).
  Poly seasonal_ma() const {
    if (btheta.empty()) return p_one();
    if (btheta.size() == 1) return p_factor(-btheta[0], period);
    vd a(static_cast<std::size_t>(btheta.size()) * static_cast<std::size_t>(period) + 1, 0.0);
    a[0] = 1.0;
    for (int i = period, j = 0; i < static_cast<int>(btheta.size()); i += period, ++j)
      a[static_cast<std::size_t>(i)] = btheta[static_cast<std::size_t>(j)];
    return p_of(std::move(a));
  }
  Poly stationary_ar() const { return p_times(regular_ar(), seasonal_ar()); }
  Poly non_stationary_ar() const { return differencing_filter(period, d, bd); }
  Poly ma() const { return p_times(regular_ma(), seasonal_ma()); }

  // ArimaModel.of(SarimaModel): the one-sided shape. The SarimaModel's own
  // innovation variance is the SarimaSpec default (always 1.0 on the captured
  // path): SeatsModelSpec.innovationVariance lives on the SeatsModel and
  // never reaches the decomposition (golden "model var 1.0" under a request
  // var of 2.0, variance-free r sfd facts).
  Arima to_arima() const {
    return arima_onesided(stationary_ar(), non_stationary_ar(), ma(), 1.0);
  }
};

// SarimaModel.Builder.adjust(double[]): drop trailing |v| < 1e-6.
inline void builder_adjust(vd& v) {
  std::size_t n = v.size();
  while (n > 0 && !(ge_g(std::fabs(v[n - 1]), 1e-6))) --n;
  v.resize(n);
}

// ---------------------------------------------------------------------------
// jdplus.tramoseats.base.core.seats.DefaultModelValidator.
// ---------------------------------------------------------------------------
struct Validator {
  double xl = 0.95, eps = 1e-4;
  Sarima model;
  bool changed = false;

  // canSimplify: non-empty vector with |last| < 1e-6.
  bool can_simplify(const vd& v) const {
    return !v.empty() && lt_g(std::fabs(v[v.size() - 1]), 1e-6);
  }

  // simplifyModel: ONE pass of toBuilder().adjustOrders(true).build(), which
  // trims trailing near-zero parameters from all four vectors at once.
  bool simplify_model() {
    if (can_simplify(model.phi) || can_simplify(model.bphi) || can_simplify(model.theta) ||
        can_simplify(model.btheta)) {
      builder_adjust(model.phi);
      builder_adjust(model.bphi);
      builder_adjust(model.theta);
      builder_adjust(model.btheta);
      return true;
    }
    return false;
  }

  // stabilizeMA(cutOff, theta): the XL projection (unit-degree clamp branch +
  // root inversion re-scaled to |z| = 1/cutOff).
  bool stabilize_ma(double cut, vd& th) {
    if (th.size() == 1) {
      double v = th[0];
      if (lt_g(v, -cut)) {
        th[0] = -cut;
        return true;
      }
      if (gt_l(v, cut)) {
        th[0] = cut;
        return true;
      }
      return false;
    }
    if (th.empty()) return false;
    Poly p = p_value_of(1.0, th);
    std::vector<cplx> roots = p_roots(p, "stabilizeMA");
    bool any = false;
    for (std::size_t i = 0; i < roots.size(); ++i) {
      double m = 1.0 / c_abs(roots[i]);
      if (gt_g(m, cut)) {
        any = true;
        roots[i] = c_scale(roots[i], m / cut);
      }
    }
    if (!any) return false;
    Poly q = p_from_roots(roots);
    q = p_divide_d(q, q.get(0));
    for (std::size_t i = 0; i < th.size(); ++i) th[i] = q.get(i + 1);
    return true;
  }

  // fixMaUnitRoots(): only reached when xl >= 1. The deg-1 quirk (an outside
  // value is sent to -1, never to +1) is a faithful transcription, not a bug.
  bool fix_ma_unit_roots() {
    bool any = false;
    double lim = 1.0 - eps;
    vd& th = model.theta;
    vd& bth = model.btheta;
    if (!bth.empty()) {
      double pw = kalos::detail::glibc_pow(lim, static_cast<double>(model.period));
      double v = bth[0];
      if (lt_g(v, -pw)) {
        bth[0] = -1.0;
        any = true;
      } else if (gt_l(v, pw)) {
        bth[0] = 1.0;
        any = true;
      }
    }
    if (th.size() == 1) {
      double v = th[0];
      if (lt_g(v, -lim)) {
        th[0] = -1.0;
        any = true;
      } else if (gt_l(v, lim)) {
        th[0] = -1.0;
        any = true;
      }
    } else if (th.size() > 1) {
      Poly p = p_value_of(1.0, th);
      std::vector<cplx> roots = p_roots(p, "fixMaUnitRoots");
      bool sub = false;
      for (std::size_t i = 0; i < roots.size(); ++i) {
        double m = c_abs(roots[i]);
        if (lt_g(m, lim)) {
          sub = true;
          roots[i] = c_div_d(roots[i], m);
        }
      }
      if (sub) {
        Poly q = p_from_roots(roots);
        for (std::size_t i = 0; i < th.size(); ++i) th[i] = q.get(i + 1);
        any = true;
      }
    }
    return any;
  }

  bool change_ma() {
    if (!ge_g(xl, 1.0)) {
      bool any = false;
      if (stabilize_ma(xl, model.theta)) any = true;
      if (stabilize_ma(xl, model.btheta)) any = true;
      return any;
    }
    return fix_ma_unit_roots();
  }

  // validate(): true when the model was left untouched. changeAR(false) is a
  // no-op in 3.9.0 (its body is commented out in the released class).
  bool validate(const Sarima& m) {
    model = m;
    bool a = simplify_model();
    bool b = change_ma();
    bool c = false;
    changed = a || b || c;
    return !changed;
  }
};

// ---------------------------------------------------------------------------
// DefaultModelDecomposer.decompose wiring (the two selectors exactly as
// SeatsToolkit.of wires them from the DecompositionSpec knobs).
// ---------------------------------------------------------------------------
struct EngineFacts {
  bool chain_threw = false;  // ModelDecomposer.decompose threw -> um null
  std::vector<Arima> raw;    // the raw (pre-allocation) components
  VarianceMax vm;
};

inline EngineFacts decompose(const Sarima& m1, double trend_boundary,
                      double seas_tolerance, double seas_boundary, double seas_boundary_pi,
                      bool noisy) {
  EngineFacts ef;
  Arima model = m1.to_arima();
  std::vector<std::unique_ptr<Selector>> sels;
  std::unique_ptr<TrendCycleSelector> tcs(new TrendCycleSelector(trend_boundary));
  tcs->set_default_low_freq_threshold(m1.period);
  std::unique_ptr<SeasonalSelector> ss(new SeasonalSelector(m1.period, seas_tolerance));
  ss->set_k(m1.bd > 0 || m1.bp() > 0 ? seas_boundary : seas_boundary_pi);
  sels.push_back(std::move(tcs));
  sels.push_back(std::move(ss));
  bool fail = false;
  try {
    if (!model_decompose(model, sels, ef.raw, fail)) {
      ef.chain_threw = true;
      ef.raw.clear();
      return ef;
    }
  } catch (const jvm_decompose_abort&) {
    // JVM: null roots NPE / chase AIOOBE inside the decomposer try -> the
    // whole decompose() returns null (r.decompose null, no setVarianceMax).
    ef.chain_threw = true;
    ef.raw.clear();
    return ef;
  }
  try {
    ef.vm = variance_max(ef.raw, noisy, false);
  } catch (const jvm_decompose_abort&) {
    // Same catch in the JVM one step later: setVarianceMax came back null,
    // the raw components themselves stand (r.decompose ok, d.decompose null).
    ef.vm = VarianceMax();
  }
  return ef;
}

inline void push_event(tramo_seats_model_result& r, int depth, const std::string& kind,
                const std::string& msg, const std::string& detail) {
  tramo_seats_model_result::event e;
  e.seq = static_cast<int>(r.events.size()) + 1;
  e.kind = kind;
  e.depth = depth;
  e.msg = msg;
  e.detail = detail;
  r.events.push_back(e);
}


// The 24-A driver, verbatim (moved out of validate_tramo_seats_model; the
// public entry point is now a one-line wrapper).
inline tramo_seats_model_result engine_run(
    const tramo_seats_model_request& request) {
  if (request.period != 12 && request.period != 4)
    refuse("frequency " + std::to_string(request.period) +
           " outside the captured corpus shapes {4, 12}");
  if (request.frequency != 12 && request.frequency != 4)
    refuse("series frequency " + std::to_string(request.frequency) +
           " outside the captured corpus shapes {4, 12}");
  if (request.d < 0 || request.bd < 0) refuse("negative differencing order");
  for (const vd* v : {&request.phi, &request.bphi, &request.theta, &request.btheta})
    for (double x : *v)
      if (is_nan(x) || std::fabs(x) == std::numeric_limits<double>::infinity())
        refuse("non-finite model parameter");
  if (is_nan(request.xl) || is_nan(request.epsphi) || is_nan(request.rmod) ||
      is_nan(request.smod) || is_nan(request.smod_pi) ||
      is_nan(request.innovation_variance))
    refuse("non-finite setting");

  tramo_seats_model_result out;
  out.xl = request.xl;

  // ---- SeatsKernel.validate seam (DefaultModelValidator) ------------------
  push_event(out, 0, "push", "seats", "-");
  push_event(out, 1, "push", "validation", "-");
  Sarima m0;
  m0.period = request.period;
  m0.d = request.d;
  m0.bd = request.bd;
  m0.phi = request.phi;
  m0.bphi = request.bphi;
  m0.theta = request.theta;
  m0.btheta = request.btheta;
  Validator val;
  val.xl = request.xl;
  val.eps = 1e-4;  // DefaultModelValidator.Builder DEF_EPS (the toolkit sets xl only)
  bool untouched = val.validate(m0);
  out.parameters_cut_off = !untouched;
  if (!untouched) push_event(out, 2, "remark", "arima parameters cut off", "-");
  push_event(out, 1, "pop", "-", "-");
  const Sarima& m1 = val.model;
  out.orders_p = m1.p();
  out.orders_d = m1.d;
  out.orders_q = m1.q();
  out.orders_bp = m1.bp();
  out.orders_bd = m1.bd;
  out.orders_bq = m1.bq();
  out.phi = m1.phi;
  out.bphi = m1.bphi;
  out.theta = m1.theta;
  out.btheta = m1.btheta;
  out.innovation_variance = request.innovation_variance;
  out.period = m1.period;

  // ---- SeatsKernel.decomposeModel seam -------------------------------------
  const bool noisy = request.mode == seats_approx_mode::noisy;
  const bool approximator_present = request.mode != seats_approx_mode::none;
  push_event(out, 1, "push", "decomposition", "-");
  EngineFacts ef;
  try {
    ef = decompose(m1, request.rmod, request.epsphi,
                   request.smod, request.smod_pi, noisy);
  } catch (const std::invalid_argument&) {
    // Envelope violations propagate to the caller (never swallowed): the live
    // path would have used machinery this unit does not model.
    throw;
  }
  out.decomposition_threw = ef.chain_threw;
  out.variance_floor = ef.vm.sum;
  out.variance_max_ok = ef.vm.ok;
  for (const Arima& c : ef.raw) {
    seats_raw_component rc;
    rc.sar = c.stationary_ar().c;
    rc.ur = c.non_stationary_ar().c;
    rc.sfd = c.symmetric_ma().a;
    out.raw_components.push_back(rc);
  }

  // The kernel loop's first iteration decides everything this unit models.
  const bool um_null = !ef.vm.ok;
  if (um_null) push_event(out, 2, "warning", "non decomposable model", "-");
  if (um_null && !approximator_present) {
    push_event(out, 2, "error", "canonical decomposition failed", "-");
    push_event(out, 1, "pop", "-", "-");
    push_event(out, 0, "pop", "-", "-");
    out.seats_exception = true;
    return out;
  }
  if (um_null) {
    out.approximation_required = true;  // D6: the kernel calls the approximator here
    return out;
  }
  if (noisy && lt_g(ef.vm.sum, 0.0)) {
    // setVarianceMax replaced the model with the noisy ArimaModel; the kernel
    // compares um.getModel() against the SarimaModel, sees them differ and
    // marks the model changed.
    out.noisy_model_used = true;
    out.model_changed = true;
    push_event(out, 2, "warning", "noisy model used", "-");
  }
  push_event(out, 1, "pop", "-", "-");
  push_event(out, 0, "pop", "-", "-");
  return out;
}

}  // namespace detail
}  // namespace kalos
