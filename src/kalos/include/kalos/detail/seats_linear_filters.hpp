#pragma once

// kalos M5 slice 24-C.1, cluster 1: the linearfilters mirror backing the WK
// extraction machine. Transcription of the pinned 3.9.0 bytecode of
//   jdplus.toolkit.base.core.math.polynomials.RationalFunction
//   jdplus.toolkit.base.core.math.polynomials.Polynomial$SimplifyingTool
//   jdplus.toolkit.base.core.math.linearfilters.{BackFilter, ForeFilter,
//     SymmetricFilter (core algebra), FiniteFilter, RationalBackFilter,
//     RationalForeFilter, RationalFilter (+ $Builder, $Decomposition)}
//   jdplus.toolkit.base.core.math.linearfilters.FilterUtility.frequencyResponse
//   jdplus.toolkit.base.core.math.ComplexComputer (accumulator algebra)
//   jdplus.toolkit.base.api.math.{Complex.div/plus, DoubleSeq.equals(DDD)}
//   jdplus.toolkit.base.core.math.linearsystem.LULinearSystemSolver2(1.0E-13)
//     (LinearSystemSolver.robustSolver() per the LS_Factory BootstrapMethods
//     table, entry 1 -> lambda$static$1)
// javap dumps: /tmp/opencode/ssf/{RF_cluster.txt, jdplus...RationalFilter.txt,
// Rational{Back,Fore}Filter.txt, BackFilter.txt, FilterUtility.txt,
// Polynomial.txt, Polynomial_SimplifyingTool.txt, ComplexComputer.txt,
// Complex.txt, DoubleSeq.txt, LULinearSystemSolver2.txt, SymmetricFilter.txt,
// RFBuilder.txt, Polynomials_shaded.txt}.
//
// Fidelity notes (pinned against the dumps, 2026-11-22):
//  * RfSimp below is the WK mirror's OWN faithful copy of
//    Polynomial$SimplifyingTool. It differs from the 24-B seam's
//    poly_simplify_private (seats_ucarima_model.hpp) in two points that the
//    24-B corpus never exercises; the seam stays untouched this slice and
//    the divergences are recorded in the pinning doc SECTION 15 and the WK
//    contract doc:
//      (a) 3.9.0 bytecode takes the roots of ARG2 (the SMALLER polynomial)
//          in the private simplify (the roots parameter is null from the
//          2-arg caller); the seam takes them of `big`
//          (seats_ucarima_model.hpp line ~1206).
//      (b) For a NaN imaginary root part the pinned dcmpl-ifle refreshes
//          f2 (to NaNs); the seam keeps the stale f2 there. Both agree for
//          im <= 0 (stale f2) and for the real-root branch, where the f2
//          shared tail DOES run after the f1 division (149: goto 184,
//          re-verified fresh from the pinned jar 2026-11-22 - an earlier
//          decode note claiming "goto 221" was wrong).
//  * Polynomial.equals(Polynomial, eps): the JVM reference-equality
//    shortcut cannot be mirrored under value semantics; it only diverges
//    for two aliased objects carrying NaN coefficients, which this corpus
//    never builds.
//  * Polynomial.get(int) is the raw daload (ArrayIndexOutOfBounds in the
//    JVM); pg() refuses instead.
//  * Math.cos/sin inside frequencyResponse dispatch through the live C2
//    stubs (host-libm lineage, CI-FDLIBM ledger 2026-11-21) - raw std
//    calls, same precedent as rfr() in ucarima_engine.hpp.
//  * Decomposition.of assembles the dense system with the bytecode's own
//    indexing (column-major FastMatrix, set(i,j) = storage[i + j*n]) and
//    solves it in place with lu2_solve below. MatrixException surfaces as
//    linear_filter_exception("Invalid decomposition of rational filter").
//  * SymmetricFilter.decompose (SymmetricFilterAlgorithms.decomposer) is
//    NOT mirrored yet - only reachable through
//    RationalFilter.rationalSymmetricFilter, which stays unported until the
//    WK cluster proves whether the builder path needs it. apply/applyTo
//    (IFiniteFilter defaults) likewise port with the consuming cluster.
//
// NOTE (linkage): this header includes seats_ucarima_model.hpp; a handful
// of its free functions lack `inline'. Linking two translation units that
// include the seam together needs a separate behavior-neutral commit to
// add the keywords (tracked for cluster 5+).

#include "kalos/detail/hostlineage.hpp"
#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include "kalos/detail/seats_ucarima_model.hpp"

namespace kalos {
namespace detail {
namespace seats_lfilters {

using seats_model::eval_complex;  // Polynomial.evaluateAt(Complex)
using seats_model::ls_divide;     // LeastSquaresDivision.divide
using seats_model::ls_exact;      // LeastSquaresDivision.isExact
using seats_model::idx;           // narrowing-free index cast
using seats_model::LsDiv;
using seats_model::c_add;         // Complex.plus accumulation

// Exception mirrors (jdplus exception surface).
struct polynomial_exception : std::runtime_error {
  explicit polynomial_exception(const std::string& m) : std::runtime_error(m) {}
};
struct matrix_exception : std::runtime_error {
  explicit matrix_exception(const std::string& m) : std::runtime_error(m) {}
};
struct linear_filter_exception : std::runtime_error {
  explicit linear_filter_exception(const std::string& m) : std::runtime_error(m) {}
};

// Polynomial.get(int), raw daload; JVM AIOOBE paths refuse instead.
inline double pg(const Poly& p, int i) {
  if (i < 0 || static_cast<std::size_t>(i) >= p.c.size())
    throw std::out_of_range("seats_linear_filters: Polynomial.get index (AIOOBE parity)");
  return p.c[idx(i)];
}

// DoubleSeq.equals(double a, double b, double eps) (pinned static @1098).
inline bool dseq_eq(double a, double b, double eps) {
  if (!gt_l(a, b)) return !gt_g(b - eps, a);  // dcmpl ifle branch (NaN goes here)
  return !gt_g(a - eps, b);
}

// Polynomial.equals(Polynomial, double) without the reference shortcut.
inline bool poly_equals_eps(const Poly& p, const Poly& q, double eps) {
  if (p.c.size() != q.c.size()) return false;
  for (std::size_t i = 0; i < p.c.size(); ++i)
    if (!dseq_eq(p.c[i], q.c[i], eps)) return false;
  return true;
}

// Polynomials.plus: swap to the longer base (clone), add the shorter
// elementwise, tail copied verbatim (no +0.0 on the tail).
inline vd pl_plus(const vd& a, const vd& b) {
  const vd* base = &a;
  const vd* add = &b;
  if (a.size() < b.size()) std::swap(base, add);
  vd r = *base;
  for (std::size_t i = 0; i < add->size(); ++i) r[i] = r[i] + (*add)[i];
  return r;
}

// Polynomial.plus(Polynomial): arg.isZero() shortcut, then trim(plus).
inline Poly p_plus_p(const Poly& a, const Poly& b) {
  if (b.is_zero()) return a;
  return p_of(pl_plus(a.c, b.c));
}

// Polynomial.plus(double) / minus(double): c == 0.0 returns this, else the
// raw copy with coeff[0] += c.
inline Poly p_plus_d(const Poly& a, double c) {
  if (c == 0.0) return a;
  vd r = a.c;
  r[0] = r[0] + c;
  return p_raw(std::move(r));
}
inline Poly p_minus_d(const Poly& a, double c) { return p_plus_d(a, -c); }

// ---------------------------------------------------------------------------
// ComplexComputer (mutable accumulator) + Complex.div(Complex).
// ---------------------------------------------------------------------------
struct CC {
  double re = 0.0, im = 0.0;
  void mul_d(double s) {  // mul(double): dcmpl ifne; s == 0 (incl. -0.0) zeroes
    if (!(s != 0.0)) {
      re = 0.0;
      im = 0.0;
    } else {
      re = re * s;
      im = im * s;
    }
  }
  void add_ac(double a, cplx z) {  // addAC(a, z)
    re = re + a * z.re;
    im = im + a * z.im;
  }
  void sub_cc(cplx z) {  // sub(ComplexType)
    re = re - z.re;
    im = im - z.im;
  }
};

// Complex.div(Complex) == ComplexComputer.div(re, im) (pinned @420-507):
// scaled division branching on |br| >= |bi| (dcmpg - NaN takes the ge side).
inline cplx cq_div(cplx a, cplx b) {
  double br = b.re, bi = b.im;
  if (!lt_g(std::fabs(br), std::fabs(bi))) {
    double rat = 1.0 / (br + bi * (bi / br));
    double re = rat * (a.re + a.im * (bi / br));
    double im = rat * (a.im - a.re * (bi / br));
    return {re, im};
  }
  double rat = 1.0 / (bi + br * (br / bi));
  double re = rat * (a.re * (br / bi) + a.im);
  double im = rat * (a.im * (br / bi) - a.re);
  return {re, im};
}

// ---------------------------------------------------------------------------
// FilterUtility.frequencyResponse(g, lo, hi, w) (pinned @576-662): first
// factor from trig, the rest by the three-term recurrence
// z_{n+1} = (2 cos w) z_n - z_{n-1} through the ComplexComputer. Trig calls
// are host-libm (C2 stub lineage).
// ---------------------------------------------------------------------------
inline cplx fu_freq(const std::function<double(int)>& g, int lo, int hi, double w) {
  double cs = kalos::detail::glibc_cos(w);
  int i = lo;
  cplx z = (lo == 0) ? cplx{1.0, 0.0}
                     : c_cart(kalos::detail::glibc_cos(w * static_cast<double>(lo)),
                              kalos::detail::glibc_sin(w * static_cast<double>(lo)));
  CC acc{z.re, z.im};
  acc.mul_d(g(i++));
  if (i > hi) return c_cart(acc.re, acc.im);
  cplx zi = c_cart(kalos::detail::glibc_cos(w * static_cast<double>(i)),
                   kalos::detail::glibc_sin(w * static_cast<double>(i)));
  acc.add_ac(g(i++), zi);
  while (i <= hi) {
    CC t{zi.re, zi.im};
    t.mul_d(2.0 * cs);
    t.sub_cc(z);
    cplx zn = c_cart(t.re, t.im);
    acc.add_ac(g(i++), zn);
    z = zi;
    zi = zn;
  }
  return c_cart(acc.re, acc.im);
}

// ---------------------------------------------------------------------------
// LULinearSystemSolver2 (robustSolver(), eps = 1.0E-13) for a square
// n x n column-major matrix A and one right-hand side y (solve(FastMatrix,
// DataBlock) == columnOf(y) round trip). Flat-array indexing transcribed
// verbatim, including the mixed views of the augmented storage.
// ---------------------------------------------------------------------------
inline void lu2_solve(const vd& a, int n, vd& y) {
  if (static_cast<int>(y.size()) != n) throw matrix_exception("m_err_dim");
  int ncols = n + 1;
  vd data = a;
  data.insert(data.end(), y.begin(), y.end());
  std::vector<int> perm(static_cast<std::size_t>(ncols), -1);
  constexpr double eps = 1.0e-13;
  for (int k = 0; k < n; ++k) {
    double max = std::numeric_limits<double>::denorm_min();  // 4.9E-324
    int i_max = -1;
    for (int i = 0, ii = k; i < n; ++i, ii += n) {
      if (perm[idx(i)] == -1) {
        double v = std::fabs(data[idx(ii)]);
        if (gt_l(v, max)) {
          i_max = i;
          max = v;
        }
      }
    }
    if (i_max < 0) throw matrix_exception("");  // bare MatrixException()
    perm[idx(i_max)] = k;
    double piv = data[idx(k + i_max * n)];
    if (lt_g(std::fabs(piv), eps)) throw matrix_exception("m_err_sing");
    double inv = 1.0 / piv;
    for (int j = 0, rb = i_max * n; j < n; ++j, ++rb) {
      if (j == k) continue;
      double h = data[idx(rb)];
      if (h == 0.0) continue;
      double hh = inv * h;
      for (int c = 0, src = k, dst = j; c < ncols; ++c, src += n, dst += n)
        if (perm[idx(c)] == -1) data[idx(dst)] = data[idx(dst)] - hh * data[idx(src)];
    }
    for (int c = 0, ik = k; c < ncols; ++c, ik += n)
      if (perm[idx(c)] == -1) data[idx(ik)] = data[idx(ik)] * inv;
  }
  for (int i = 0; i < n; ++i)
    if (perm[idx(i)] != -1) y[idx(i)] = data[idx(perm[idx(i)] + n * n)];
}

// ---------------------------------------------------------------------------
// IFiniteFilter core: the three primitives + the weightsToArray /
// frequencyResponse defaults (pinned IFiniteFilter defaults).
// ---------------------------------------------------------------------------
struct FinFilter {
  virtual ~FinFilter() = default;
  virtual int get_lower_bound() const = 0;
  virtual int get_upper_bound() const = 0;
  virtual double weight(int t) const = 0;
  virtual cplx frequency_response(double w) const {
    return fu_freq([this](int t) { return weight(t); }, get_lower_bound(),
                   get_upper_bound(), w);
  }
  virtual vd weights_to_array() const {
    std::size_t len = static_cast<std::size_t>(get_upper_bound() - get_lower_bound() + 1);
    vd w(len, 0.0);
    int t = get_lower_bound();
    for (std::size_t k = 0; k < len; ++k, ++t) w[k] = weight(t);
    return w;
  }
};

// ---------------------------------------------------------------------------
// Polynomial$SimplifyingTool, faithful 3.9.0 copy (see header notes).
// ---------------------------------------------------------------------------
struct RfSimp {
  bool has_lr = false;
  Poly left, right;  // placeholders until set (Java's null fields)
  Poly common{vd{1.0}};

  bool simplify_priv(const Poly& big, const Poly& small) {
    if (small.c.size() == 1) return false;
    LsDiv dv = ls_divide(big, small);
    if (dv.have_quotient && ls_exact(dv.err)) {  // simplifyExact
      has_lr = true;
      left = dv.quotient;
      right = Poly(vd{1.0});
      common = small;
      return true;
    }
    std::vector<cplx> zroots = p_roots(small, "rf.simplify");  // roots of ARG2
    vd f1{0.0, 1.0};       // Coefficients.fromDegree(1)
    vd f2{0.0, 0.0, 1.0};  // Coefficients.fromDegree(2)
    Poly cm{vd{1.0}};
    for (const cplx& z : zroots) {
      if (lt_g(c_abs(eval_complex(big, z)), 1.0e-12)) {  // dcmpg ifge: NaN skips
        double re = z.re, im = z.im;
        if (im == 0.0) {  // dcmpl ifne (-0.0 counts as real)
          f1[0] = -re;
          cm = p_raw(pl_times(cm.c, pl_divide(big.c, f1)));
        }  // 149: goto 184 - the f2 shared tail ALSO runs after real roots
        if (!(im <= 0.0)) {  // dcmpl ifle: im <= 0 keeps the stale f2; NaN im
                             // (dcmpl -> 1) DOES refresh it (with NaNs)
          f2[0] = re * re + im * im;
          f2[1] = -2.0 * re;
        }
        cm = p_raw(pl_times(cm.c, pl_divide(big.c, f2)));
      }
    }
    if (cm.degree() > 0) {
      has_lr = true;
      left = big;
      right = p_quotient(small, cm);
      common = cm;
      return true;
    }
    common = cm;
    return false;
  }

  bool simplify(const Poly& p1, const Poly& p2) {
    has_lr = false;
    left = Poly();
    right = Poly();
    common = Poly(vd{1.0});
    if (p1.c.size() >= p2.c.size()) {
      simplify_priv(p1, p2);
    } else {
      simplify_priv(p2, p1);
      std::swap(left, right);  // the JVM swaps the (possibly unset) fields
    }
    return common.c.size() > 1;
  }
};

// ---------------------------------------------------------------------------
// RationalFunction (lazy power-series cache; the JVM `synchronized' /
// volatile pair is behavior-free in the single-threaded replay).
// ---------------------------------------------------------------------------
struct RatFunc {
  Poly num{vd{0.0}};   // private ctor: Polynomial.ZERO
  Poly den{vd{1.0}};   // private ctor: Polynomial.ONE
  std::optional<vd> coeff;

  RatFunc() = default;

  static RatFunc zero() { return RatFunc(); }
  static RatFunc one() {
    RatFunc r;
    r.num = Poly(vd{1.0});
    r.den = Poly(vd{1.0});
    return r;  // of(ONE, ONE), simplify = false
  }
  static RatFunc of(Poly n, Poly d) { return RatFunc(std::move(n), std::move(d), false); }

  RatFunc(Poly n, Poly d, bool simplify) {
    if (simplify) {
      RfSimp t;
      if (t.simplify(n, d)) {
        num = std::move(t.left);
        den = std::move(t.right);
        return;
      }
    }
    num = std::move(n);
    den = std::move(d);
  }

  bool is_finite() const { return den.degree() == 0; }

  void prepare(int n) {
    if (coeff.has_value() && static_cast<int>(coeff->size()) > n) return;
    int dd = den.degree(), dn = num.degree();
    double d0 = den.c[0];
    if (dd == 0) {
      vd c = num.c;  // num.coefficients().toArray()
      if (!(d0 == 1.0))  // dcmpl ifeq: only d0 == 1.0 skips the division
        for (std::size_t k = 0; k < c.size(); ++k) c[k] = c[k] / d0;
      coeff = std::move(c);
      return;
    }
    int k = std::max(dd, dn);
    if (n < k)
      n = k;
    else
      n = ((n - 1) / 32 + 1) * 32;  // Java int division; C++ matches (trunc)
    int start = 1;
    if (!coeff.has_value()) {
      vd c(static_cast<std::size_t>(n) + 1, 0.0);
      c[0] = num.c[0] / d0;
      for (int i = 1; i <= dd; ++i) {
        double s = 0.0;
        for (int j = 1; j <= i; ++j) s += den.c[idx(j)] * c[idx(i - j)];
        if (i <= dn)
          c[idx(i)] = (num.c[idx(i)] - s) / d0;
        else
          c[idx(i)] = (-s) / d0;
      }
      if (dn > dd) {  // fresh numerator continuation
        for (int i = dd + 1; i <= dn; ++i) {
          double s = 0.0;
          for (int j = 1; j <= dd; ++j) s += den.c[idx(j)] * c[idx(i - j)];
          c[idx(i)] = (num.c[idx(i)] - s) / den.c[0];
        }
      }
      start = k + 1;
      coeff = std::move(c);
    } else {
      vd nc(static_cast<std::size_t>(n) + 1, 0.0);
      start = static_cast<int>(coeff->size());
      std::copy(coeff->begin(), coeff->end(), nc.begin());
      coeff = std::move(nc);
    }
    for (int i = start; i <= n; ++i) {
      double s = 0.0;
      for (int j = 1; j <= dd; ++j) s += den.c[idx(j)] * (*coeff)[idx(i - j)];
      (*coeff)[idx(i)] = (-s) / d0;
    }
  }

  double get(int i) {
    prepare(i);
    if (i < 0)  // the JVM daloads coeff[i] -> AIOOBE (after prepare ran)
      throw std::out_of_range("seats_linear_filters: get(i < 0) (AIOOBE parity)");
    if (i >= static_cast<int>(coeff->size())) return 0.0;  // dconst_0
    return (*coeff)[idx(i)];
  }

  vd coefficients(int n) {
    if (n < 0) throw std::out_of_range("seats_linear_filters: negative coefficients(n)");
    prepare(n - 1);
    vd out(static_cast<std::size_t>(n), 0.0);
    int k = static_cast<int>(coeff->size()) < n ? static_cast<int>(coeff->size()) : n;
    for (int i = 0; i < k; ++i) out[idx(i)] = (*coeff)[idx(i)];
    return out;
  }

  RatFunc plus(const RatFunc& o) const {
    Poly n = p_plus_p(p_times(num, o.den), p_times(den, o.num));
    Poly d = p_times(den, o.den);
    return RatFunc(std::move(n), std::move(d), true);
  }
  RatFunc minus(const RatFunc& o) const {
    Poly n = p_minus_p(p_times(num, o.den), p_times(den, o.num));
    Poly d = p_times(den, o.den);
    return RatFunc(std::move(n), std::move(d), true);
  }
  RatFunc times(const RatFunc& o) const {
    Poly n = p_times(num, o.num);
    Poly d = p_times(den, o.den);
    return RatFunc(std::move(n), std::move(d), true);
  }
  RatFunc divide(const RatFunc& o) const {
    Poly n = p_times(num, o.den);
    Poly d = p_times(den, o.num);
    return RatFunc(std::move(n), std::move(d), true);
  }

  RatFunc drop(int d) {
    int i4 = den.degree() + 1;
    int i5 = i4;
    if (num.degree() + 1 - d > i5) i5 = num.degree() + 1 - d;
    if (i5 == 0) return zero();
    vd a(static_cast<std::size_t>(i5), 0.0), q(static_cast<std::size_t>(i5), 0.0);
    for (int k = 0; k < i5; ++k) a[idx(k)] = get(k + d);  // the SERIES
    q[0] = a[0] * den.c[0];
    for (int i = 1; i < i5; ++i) {
      double s = 0.0;
      int i11 = i < i5 ? i : i5 - 1;  // always i; kept for transcription parity
      for (int k = 0; k <= i11; ++k) {
        int i13 = i - k;
        if (i13 < i4) s += pg(den, i13) * a[idx(k)];
      }
      if (lt_g(std::fabs(s), kFeEps)) s = 0.0;  // dcmpg keeps NaN
      q[idx(i)] = s;
    }
    int m = i5;
    while (m > 0 && q[idx(m - 1)] == 0.0) --m;  // dcmpl: NaN stops the trim
    if (m == 0) return zero();
    if (m != i5) q = vd(q.begin(), q.begin() + m);
    RatFunc r;
    r.num = p_of(q);  // Polynomial.ofInternal (trim)
    r.den = den;
    return r;  // ctor(nq, den, false)
  }

  double evaluate_at(double x) const {
    double n = p_value_at(num, x);
    double d = p_value_at(den, x);
    if (!gt_g(std::fabs(d - 0.0), kFeEps))
      throw polynomial_exception("rf_err_pole");  // dcmpg: NaN denom passes
    return n / d;
  }
  cplx evaluate_at(cplx z) const {
    cplx n = eval_complex(num, z);
    cplx d = eval_complex(den, z);
    if (!gt_g(c_abs(d), kFeEps)) throw polynomial_exception("rf_err_pole");
    return cq_div(n, d);
  }

  std::vector<cplx> roots() const { return p_roots(num, "rf.roots"); }
  std::vector<cplx> poles() const { return p_roots(den, "rf.poles"); }
};

// ---------------------------------------------------------------------------
// BackFilter / ForeFilter (one-sided filters == their Polynomial, raw).
// ---------------------------------------------------------------------------
class ForeFilter;
class BackFilter : public FinFilter {
 public:
  Poly p;
  explicit BackFilter(Poly q) : p(std::move(q)) {}
  static BackFilter zero() { return BackFilter(Poly(vd{0.0})); }
  static BackFilter one() { return BackFilter(Poly(vd{1.0})); }
  static BackFilter of_internal(const vd& w) {  // ofInternal(double...)
    if (w.size() == 1) {
      if (w[0] == 1.0) return one();  // dcmpl
      if (w[0] == 0.0) return zero();
    }
    return BackFilter(p_of(w));  // Polynomial.ofInternal (trim)
  }
  int degree() const { return p.degree(); }
  int length() const { return p.degree() + 1; }
  int get_lower_bound() const override { return -p.degree(); }
  int get_upper_bound() const override { return 0; }
  double weight(int t) const override { return pg(p, -t); }
  double get(int i) const { return pg(p, i); }
  cplx frequency_response(double w) const override {
    return fu_freq([this](int i) { return pg(p, i); }, 0, p.degree(), -w);
  }
  vd weights_to_array() const override {  // weights(-deg..0) = reversed coeffs
    vd w = p.c;
    std::reverse(w.begin(), w.end());
    return w;
  }
  const Poly& as_polynomial() const { return p; }
  BackFilter times(const BackFilter& o) const { return BackFilter(p_times(p, o.p)); }
  BackFilter times(double d) const { return BackFilter(p_times_d(p, d)); }
  BackFilter plus(double c) const { return BackFilter(p_plus_d(p, c)); }
  BackFilter minus(double c) const { return BackFilter(p_minus_d(p, c)); }
  BackFilter normalize() const {  // get(0) != 0 && != 1 (dcmpl; NaN divides)
    double p0 = p.c[0];
    if (p0 == 0.0 || p0 == 1.0) return *this;
    return BackFilter(p_times_d(p, 1.0 / p0));
  }
  // BackFilter.mirror (pinned BackFilter.txt 230): new ForeFilter(the SAME
  // Polynomial) - no reversal; the ForeFilter convention indexes it forward.
  ForeFilter mirror() const;  // defined after ForeFilter
  // BackFilter.divide (pinned BackFilter.txt 82): Polynomial.divide(a, b)
  // classic long division (Polynomial.txt 1412 - NOT the least-squares
  // divide), Division(m_r, m_q) = (remainder, quotient), quotient kept.
  // The quotient array (length deg(a)-deg(b)+1) is Arrays.copyOf'd and fed
  // to the RAW Polynomial ctor - trailing zeros survive the invariant;
  // deg(b) > deg(a) yields the ZERO quotient (Division(a, ZERO) at 15-26);
  // the remainder (copyOf(c, deg(b)) through Coefficients.ofInternal, ZERO
  // when deg(b) == 0) is discarded here. Inner update is c[j] - q*b.get().
  BackFilter divide(const BackFilter& o) const {
    int da = p.degree(), db = o.p.degree();
    if (db > da) return BackFilter(Poly(vd{0.0}));
    vd c(static_cast<std::size_t>(da) + 1), qq(static_cast<std::size_t>(da) + 1, 0.0);
    for (int i = 0; i <= da; ++i) c[idx(i)] = pg(p, i);  // qq[i] = 0.0 by ctor
    for (int i = da - db; i >= 0; --i) {
      qq[idx(i)] = c[idx(db + i)] / pg(o.p, db);
      for (int j = db + i - 1; j >= i; --j)
        c[idx(j)] = c[idx(j)] - (qq[idx(i)] * pg(o.p, j - i));
    }
    vd qc(qq.begin(), qq.begin() + static_cast<std::ptrdiff_t>(da - db + 1));
    return BackFilter(p_raw(std::move(qc)));
  }
};

class ForeFilter : public FinFilter {
 public:
  Poly p;
  explicit ForeFilter(Poly q) : p(std::move(q)) {}
  static ForeFilter zero() { return ForeFilter(Poly(vd{0.0})); }
  static ForeFilter one() { return ForeFilter(Poly(vd{1.0})); }
  static ForeFilter of_internal(const vd& w) {
    if (w.size() == 1) {
      if (w[0] == 1.0) return one();
      if (w[0] == 0.0) return zero();
    }
    return ForeFilter(p_of(w));
  }
  int degree() const { return p.degree(); }
  int length() const { return p.degree() + 1; }
  int get_lower_bound() const override { return 0; }
  int get_upper_bound() const override { return p.degree(); }
  double weight(int t) const override { return pg(p, t); }
  double get(int i) const { return pg(p, i); }
  cplx frequency_response(double w) const override {
    return fu_freq([this](int i) { return pg(p, i); }, 0, p.degree(), w);
  }
  const Poly& get_polynomial() const { return p; }
  ForeFilter times(const ForeFilter& o) const { return ForeFilter(p_times(p, o.p)); }
  ForeFilter times(double d) const { return ForeFilter(p_times_d(p, d)); }
  ForeFilter plus(double c) const { return ForeFilter(p_plus_d(p, c)); }
  ForeFilter minus(double c) const { return ForeFilter(p_minus_d(p, c)); }
  ForeFilter normalize() const {
    double p0 = p.c[0];
    if (p0 == 0.0 || p0 == 1.0) return *this;
    return ForeFilter(p_times_d(p, 1.0 / p0));
  }
  BackFilter mirror() const { return BackFilter(p); }  // new BackFilter(polynomial)
};

inline ForeFilter BackFilter::mirror() const { return ForeFilter(p); }

// ---------------------------------------------------------------------------
// SymmetricFilter core: constructor, ofInternal, convolutionOf, bounds,
// weight(t) = c|t|, frequencyResponse = cart(realFrequencyResponse(w)).
// (decompose / factorize chains stay in the 24-B seam; not re-exported here.)
// ---------------------------------------------------------------------------
class SymmetricFilter : public FinFilter {
 public:
  Poly p;
  explicit SymmetricFilter(Poly q) : p(std::move(q)) {}
  static SymmetricFilter zero() { return SymmetricFilter(Poly(vd{0.0})); }
  static SymmetricFilter one() { return SymmetricFilter(Poly(vd{1.0})); }
  static SymmetricFilter of_internal(const vd& w) {
    if (w.size() == 1) {
      if (w[0] == 1.0) return one();
      if (w[0] == 0.0) return zero();
    }
    return SymmetricFilter(p_of(w));
  }
  static SymmetricFilter conv_of(const FinFilter& f) { return conv_of_d(f, 1.0); }
  static SymmetricFilter conv_of_d(const FinFilter& f, double scale) {
    vd w = f.weights_to_array();
    return conv_of_arr(w, scale);
  }
  static SymmetricFilter conv_of_poly(const Poly& q, double scale) {
    return conv_of_arr(q.c, scale);  // Polynomial.toArray()
  }
  static SymmetricFilter conv_of_arr(const vd& w, double scale) {
    vd r(w.size(), 0.0);
    for (std::size_t i = 0; i < w.size(); ++i)
      for (std::size_t j = i; j < w.size(); ++j) r[j - i] = r[j - i] + w[i] * w[j];
    if (!(scale == 1.0))  // dcmpl ifeq: only scale == 1.0 skips
      for (std::size_t i = 0; i < r.size(); ++i) r[i] = r[i] * scale;
    return of_internal(r);
  }
  int degree() const { return p.degree(); }
  // SymmetricFilter.times(double) (pinned SymmetricFilter.txt 661): wrap of
  // polynomial.times(d).
  SymmetricFilter times(double d) const { return SymmetricFilter(p_times_d(p, d)); }
  // SymmetricFilter.times(SymmetricFilter) (pinned SymmetricFilter.txt 672).
  // Product of two symmetric filters stored as one-sided polynomials; the
  // r[0] term of an (i, i) pair is added as 2.0 * m (bytecode 259-271), the
  // loop-1/loop-2 zero tests are dcmpl against 0.0 (-0.0 skips, NaN
  // participates) and the loop-3 inner tests skip only exact +0.0.
  SymmetricFilter times(const SymmetricFilter& f) const {
    int d1 = p.degree(), d2 = f.p.degree();
    vd r(static_cast<std::size_t>(d1 + d2 + 1), 0.0);
    if (!(pg(p, 0) == 0.0))
      for (int i = 0; i <= d1; ++i)
        r[idx(i)] = r[idx(i)] + (pg(p, i) * pg(f.p, 0));
    if (!(pg(f.p, 0) == 0.0))
      for (int i = 1; i <= d2; ++i)
        r[idx(i)] = r[idx(i)] + (pg(f.p, i) * pg(p, 0));
    for (int i = 1; i <= d1; ++i) {
      if (pg(p, i) == 0.0) continue;
      for (int j = 1; j <= d2; ++j) {
        if (pg(f.p, j) == 0.0) continue;
        double m = pg(p, i) * pg(f.p, j);
        r[idx(i + j)] = r[idx(i + j)] + m;
        if (i > j) r[idx(i - j)] = r[idx(i - j)] + m;
        else if (i < j) r[idx(j - i)] = r[idx(j - i)] + m;
        else r[0] = r[0] + (2.0 * m);
      }
    }
    return of_internal(r);
  }
  int get_lower_bound() const override { return -p.degree(); }
  int get_upper_bound() const override { return p.degree(); }
  double weight(int t) const override { return pg(p, t < 0 ? -t : t); }
  vd weights_to_array() const override {
    int deg = p.degree();
    vd w(static_cast<std::size_t>(2 * deg + 1), 0.0);
    for (int k = 0; k <= 2 * deg; ++k) w[idx(k)] = pg(p, k - deg < 0 ? deg - k : k - deg);
    return w;
  }
  cplx frequency_response(double w) const override {
    return c_cart(rfr(SymF{p.c}, w));
  }
  const Poly& coefficients_as_polynomial() const { return p; }
};

// ---------------------------------------------------------------------------
// FiniteFilter + the sum/product algebra (statics of FiniteFilter).
// ---------------------------------------------------------------------------
class FiniteFilter : public FinFilter {
 public:
  vd w;
  int lbv = 0;
  FiniteFilter() : w{0.0} {}
  FiniteFilter(vd ww, int lb) : w(std::move(ww)), lbv(lb) {}
  explicit FiniteFilter(const FinFilter& f)  // ctor(IFiniteFilter)
      : w(f.weights_to_array()), lbv(f.get_lower_bound()) {}
  FiniteFilter(const Poly& q, int lb) : w(q.c), lbv(lb) {}  // ctor(Polynomial, int)
  static FiniteFilter of_internal(vd ww, int lb) { return FiniteFilter(std::move(ww), lb); }
  static FiniteFilter of(const vd& ww, int lb) { return FiniteFilter(ww, lb); }  // clone
  int length() const { return static_cast<int>(w.size()); }
  int get_lower_bound() const override { return lbv; }
  int get_upper_bound() const override { return lbv + static_cast<int>(w.size()) - 1; }
  double weight(int t) const override { return w[idx(t - lbv)]; }  // raw daload
  FiniteFilter mirror() const {  // reverse + lb' = -lb - len + 1
    vd r = w;
    std::reverse(r.begin(), r.end());
    return FiniteFilter(std::move(r), -lbv - static_cast<int>(w.size()) + 1);
  }
  bool is_symmetric() const {  // pinned: (len-1) % 2 == 0 && len-1 == -lb,
    int n = static_cast<int>(w.size()) - 1;  // pairwise |w[i]-w[n-i]| <= 1e-4
    if ((n % 2) != 0) return false;          // (dcmpl ifle: NaN -> false)
    if (n != -lbv) return false;
    for (int i = 0; i < n / 2; ++i)
      if (!(std::fabs(w[idx(i)] - w[idx(n - i)]) <= 1.0e-4)) return false;
    return true;
  }
};

// add(f, c): zero-padded span [min(lb,0), max(ub,0)] with c += at lag 0.
inline FiniteFilter ff_add_d(const FinFilter& f, double c) {
  int lo = f.get_lower_bound(), hi = f.get_upper_bound();
  int lo5 = lo < 0 ? lo : 0;
  int up6 = hi < 0 ? 0 : hi;
  vd w(static_cast<std::size_t>(up6 - lo5 + 1), 0.0);
  for (int t = lo; t <= hi; ++t) w[idx(t - lo5)] = f.weight(t);
  w[idx(-lo5)] = w[idx(-lo5)] + c;
  return FiniteFilter::of_internal(std::move(w), lo5);
}

inline FiniteFilter ff_add(const FinFilter& a, const FinFilter& b) {
  int lo1 = a.get_lower_bound(), hi1 = a.get_upper_bound();
  int lo2 = b.get_lower_bound(), hi2 = b.get_upper_bound();
  int lo = lo1 < lo2 ? lo1 : lo2;
  int hi = hi1 < hi2 ? hi2 : hi1;
  vd w(static_cast<std::size_t>(hi - lo + 1), 0.0);
  for (int t = lo1; t <= hi1; ++t) w[idx(t - lo)] = a.weight(t);
  for (int t = lo2; t <= hi2; ++t) w[idx(t - lo)] = w[idx(t - lo)] + b.weight(t);
  return FiniteFilter::of_internal(std::move(w), lo);
}

inline FiniteFilter ff_multiply_d(const FinFilter& f, double c) {
  vd ws = f.weights_to_array();
  for (std::size_t k = 0; k < ws.size(); ++k) ws[k] = ws[k] * c;
  return FiniteFilter::of_internal(std::move(ws), f.get_lower_bound());
}

// multiply(a, b): trim(ofInternal(weights))s through Polynomial.times;
// lb = lb(a) + lb(b); the ctor copies poly.toArray().
inline FiniteFilter ff_multiply(const FinFilter& a, const FinFilter& b) {
  int lb = a.get_lower_bound() + b.get_lower_bound();
  Poly pa = p_of(a.weights_to_array());
  Poly pb = p_of(b.weights_to_array());
  Poly t = p_times(pa, pb);
  return FiniteFilter(t.c, lb);
}

inline FiniteFilter ff_negate(const FinFilter& f) {
  vd ws = f.weights_to_array();
  for (std::size_t k = 0; k < ws.size(); ++k) ws[k] = -ws[k];
  return FiniteFilter::of_internal(std::move(ws), f.get_lower_bound());
}

inline FiniteFilter ff_subtract_d(const FinFilter& f, double c) {
  return ff_add_d(f, -c);
}

inline FiniteFilter ff_subtract(const FinFilter& a, const FinFilter& b) {
  int lo1 = a.get_lower_bound(), hi1 = a.get_upper_bound();
  int lo2 = b.get_lower_bound(), hi2 = b.get_upper_bound();
  int lo = lo1 < lo2 ? lo1 : lo2;
  int hi = hi1 < hi2 ? hi2 : hi1;
  vd w(static_cast<std::size_t>(hi - lo + 1), 0.0);
  for (int t = lo1; t <= hi1; ++t) w[idx(t - lo)] = a.weight(t);
  for (int t = lo2; t <= hi2; ++t) w[idx(t - lo)] = w[idx(t - lo)] - b.weight(t);
  return FiniteFilter::of_internal(std::move(w), lo);
}

// ---------------------------------------------------------------------------
// RationalBackFilter / RationalForeFilter.
// ---------------------------------------------------------------------------
struct RatForeFilter;

struct RatBackFilter {
  RatFunc rf;
  int bshift = 0;

  RatBackFilter() = default;  // ZERO: RationalFunction.zero(), shift 0
  RatBackFilter(RatFunc r, int s) : rf(std::move(r)), bshift(s) {}
  RatBackFilter(const BackFilter& n, const BackFilter& d, int s)
      : rf(RatFunc::of(n.p, d.p)), bshift(s) {}
  static RatBackFilter zero() { return RatBackFilter(); }

  BackFilter get_numerator() const { return BackFilter(rf.num); }
  BackFilter get_denominator() const { return BackFilter(rf.den); }
  int get_l_bound() const { return rf.is_finite() ? -rf.num.degree() : INT_MIN; }
  int get_u_bound() const { return -bshift; }
  bool has_lower_bound() const { return rf.is_finite(); }
  bool has_upper_bound() const { return true; }
  double weight(int t) { return rf.get(-bshift - t); }
  void prepare(int n) { rf.prepare(n); }
  vd get_weights(int n) { return rf.coefficients(n); }
  RatBackFilter drop(int d) { return RatBackFilter(rf.drop(d), bshift + d); }

  cplx frequency_response(double w) {
    Poly n = rf.num, d = rf.den;
    int sh = bshift;
    cplx num = fu_freq([n, sh](int i) { return pg(n, i - sh); }, sh,
                       sh + n.degree(), -w);
    cplx den = fu_freq([d](int i) { return pg(d, i); }, 0, d.degree(), -w);
    return cq_div(num, den);
  }

  RatBackFilter times(const RatBackFilter& o) {
    Poly n1 = rf.num, n2 = o.rf.num, d1 = rf.den, d2 = o.rf.den;
    RfSimp t;
    if (t.simplify(n1, d2)) {
      n1 = t.left;
      d2 = t.right;
    }
    if (t.simplify(n2, d1)) {
      n2 = t.left;
      d1 = t.right;
    }
    Poly n = p_times(n1, n2), d = p_times(d1, d2);
    return RatBackFilter(RatFunc::of(n, d), bshift + o.bshift);
  }

  RatForeFilter mirror() const;  // after RatForeFilter
};

struct RatForeFilter {
  RatFunc rf;
  int fshift = 0;

  RatForeFilter() = default;
  RatForeFilter(RatFunc r, int s) : rf(std::move(r)), fshift(s) {}
  RatForeFilter(const ForeFilter& n, const ForeFilter& d, int s)
      : rf(RatFunc::of(n.p, d.p)), fshift(s) {}
  static RatForeFilter zero() { return RatForeFilter(); }

  ForeFilter get_numerator() const { return ForeFilter(rf.num); }
  ForeFilter get_denominator() const { return ForeFilter(rf.den); }
  int get_l_bound() const { return 0; }
  int get_u_bound() const { return rf.is_finite() ? rf.num.degree() : INT_MAX; }
  bool has_lower_bound() const { return true; }
  bool has_upper_bound() const { return rf.is_finite(); }
  double weight(int t) { return rf.get(t); }  // no fshift in the bytecode
  void prepare(int n) { rf.prepare(n); }
  vd get_weights(int n) { return rf.coefficients(n); }
  RatForeFilter drop(int d) { return RatForeFilter(rf.drop(d), fshift + d); }

  cplx frequency_response(double w) {
    Poly n = rf.num, d = rf.den;
    int sh = fshift;
    cplx num = fu_freq([n, sh](int i) { return pg(n, i - sh); }, sh,
                       sh + n.degree(), w);
    cplx den = fu_freq([d](int i) { return pg(d, i); }, 0, d.degree(), w);
    return cq_div(num, den);
  }

  RatForeFilter times(const RatForeFilter& o) {
    Poly n1 = rf.num, n2 = o.rf.num, d1 = rf.den, d2 = o.rf.den;
    RfSimp t;
    if (t.simplify(n1, d2)) {
      n1 = t.left;
      d2 = t.right;
    }
    if (t.simplify(n2, d1)) {
      n2 = t.left;
      d1 = t.right;
    }
    Poly n = p_times(n1, n2), d = p_times(d1, d2);
    double k = d.c[0];
    if (!(k == 1.0)) {  // dcmpl ifeq; NaN divides through p_divide_d
      n = p_divide_d(n, k);
      d = p_divide_d(d, k);
    }
    return RatForeFilter(RatFunc::of(n, d), fshift + o.fshift);
  }

  RatBackFilter get_mirror() const { return RatBackFilter(rf, fshift); }
};

inline RatForeFilter RatBackFilter::mirror() const { return RatForeFilter(rf, bshift); }

// ---------------------------------------------------------------------------
// RationalFilter$Decomposition: the Toeplitz-like system with the bytecode's
// own indexing (set(i, j) -> A[i + j*n]), solved by lu2_solve.
// ---------------------------------------------------------------------------
struct Decomposition {
  RatBackFilter back;
  RatForeFilter fore;

  static Decomposition of(const FiniteFilter& w0, const BackFilter& b, const ForeFilter& f) {
    int p = -w0.get_lower_bound();
    int q = w0.get_upper_bound();
    vd ws = w0.weights_to_array();
    int v6 = p < 0 ? 0 : p;
    int v7 = q < 0 ? 0 : q;
    int b8 = -b.get_lower_bound();
    int f9 = f.get_upper_bound();
    int nfore = std::max(v7, f9);
    int nback = std::max(v6, b8);
    int n = nback + nfore + 1;
    if (static_cast<int>(ws.size()) != n) {  // center the filter weights
      // System.arraycopy(ws, 0, nw, nback - p, ws.length): the dst offset is
      // UNCLAMPED (AIOOBE in the JVM when p > nback or the tail overruns).
      int off = nback - p;
      if (off < 0 || off + static_cast<int>(ws.size()) > n)
        throw std::out_of_range("seats_linear_filters: arraycopy parity (AIOOBE)");
      vd nw(static_cast<std::size_t>(n), 0.0);
      std::copy(ws.begin(), ws.end(), nw.begin() + off);
      ws = std::move(nw);
    }
    vd back_w(static_cast<std::size_t>(nback) + 1, 0.0);
    vd fore_w(static_cast<std::size_t>(nfore) + 1, 0.0);
    fore_w[0] = 0.0;  // explicit dastore in the bytecode
    vd b_arr = b.weights_to_array(), f_arr = f.weights_to_array();
    vd A(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
    // f loop (bytecode 163-206): set(i + k, k, fArr[i]) for i = 0..fDeg,
    // k = 0..nback  ->  storage[(i + k) + k * n].
    for (int i = 0; i <= f9; ++i)
      for (int k = 0; k <= nback; ++k) A[idx(i + k + k * n)] = f_arr[idx(i)];
    // b loop (bytecode 209-273): row starts at 1 - bDeg + nback, advanced
    // with i; set(row + k, nback + 1 + k, bArr[i]), k = 0..nfore-1 strict.
    int row = 1 - b8 + nback;
    for (int i = 0; i <= b8; ++i, ++row)
      for (int k = 0; k < nfore; ++k)
        A[idx((row + k) + (k + nback + 1) * n)] = b_arr[idx(i)];
    try {
      lu2_solve(A, n, ws);
    } catch (const matrix_exception&) {
      throw linear_filter_exception("Invalid decomposition of rational filter");
    }
    for (int i = 0; i <= nback; ++i) back_w[idx(i)] = ws[idx(i)];
    for (int i = 1; i <= nfore; ++i) fore_w[idx(i)] = ws[idx(i + nback)];
    std::reverse(back_w.begin(), back_w.end());
    Decomposition d;
    d.back = RatBackFilter(BackFilter::of_internal(back_w), b, 0);
    d.fore = RatForeFilter(ForeFilter::of_internal(fore_w), f, 0);
    return d;
  }
};

// ---------------------------------------------------------------------------
// RationalFilter (Builder == direct ctor per $Builder bytecode). The lazy
// numerator/denominator hold the IFiniteFilter results (FiniteFilter or
// SymmetricFilter). rationalSymmetricFilter is not ported yet (needs the
// SymmetricFilterAlgorithms decomposer; see the header notes).
// ---------------------------------------------------------------------------
class RationalFilter {
 public:
  RatBackFilter back;
  RatForeFilter fore;
  std::shared_ptr<FinFilter> numerator, denominator;  // nullable fields

  RationalFilter() = default;
  RationalFilter(RatBackFilter b, RatForeFilter f, std::shared_ptr<FinFilter> n,
                 std::shared_ptr<FinFilter> d)
      : back(std::move(b)), fore(std::move(f)), numerator(std::move(n)),
        denominator(std::move(d)) {}

  static std::shared_ptr<FinFilter> sp(FiniteFilter f) {
    return std::make_shared<FiniteFilter>(std::move(f));
  }

  // of(BackFilter, BackFilter, ForeFilter, ForeFilter)
  static RationalFilter of_bbff(const BackFilter& n, const BackFilter& d,
                                const ForeFilter& fn, const ForeFilter& fd) {
    FiniteFilter w = ff_multiply(n, fn);
    Decomposition dec = Decomposition::of(w, d, fd);
    return RationalFilter(dec.back, dec.fore, sp(w),
                          sp(ff_multiply(d, fd)));
  }
  // of(RationalBackFilter)
  static RationalFilter of_back(const RatBackFilter& fb) {
    return RationalFilter(fb, RatForeFilter::zero(), sp(FiniteFilter(fb.get_numerator())),
                          sp(FiniteFilter(fb.get_denominator())));
  }
  // of(RationalForeFilter)
  static RationalFilter of_fore(const RatForeFilter& ff) {
    return RationalFilter(RatBackFilter::zero(), ff, sp(FiniteFilter(ff.get_numerator())),
                          sp(FiniteFilter(ff.get_denominator())));
  }
  // of(IFiniteFilter, BackFilter, ForeFilter)
  static RationalFilter of_wbf(const FinFilter& w0, const BackFilter& b, const ForeFilter& f) {
    FiniteFilter fw(w0);
    FiniteFilter den = ff_multiply(FiniteFilter(b), f);
    Decomposition dec = Decomposition::of(fw, b, f);
    return RationalFilter(dec.back, dec.fore, sp(fw), sp(std::move(den)));
  }

  cplx frequency_response(double w) {
    if (numerator && denominator)
      return cq_div(numerator->frequency_response(w), denominator->frequency_response(w));
    return c_add(back.frequency_response(w), fore.frequency_response(w));
  }

  std::shared_ptr<FinFilter> get_denominator() {
    if (!denominator) {
      BackFilter d1 = back.get_denominator();
      ForeFilter d2 = fore.get_denominator();
      if (poly_equals_eps(d1.p, d2.p, 0.0))
        denominator = std::make_shared<SymmetricFilter>(SymmetricFilter::conv_of(d1));
      else
        denominator = sp(ff_multiply(FiniteFilter(d1), FiniteFilter(d2)));
    }
    return denominator;
  }

  std::shared_ptr<FinFilter> get_numerator() {
    if (!numerator) {
      FiniteFilter n1(back.get_numerator()), n2(fore.get_numerator());
      FiniteFilter d1(back.get_denominator()), d2(fore.get_denominator());
      numerator = sp(ff_add(ff_multiply(n1, d2), ff_multiply(n2, d1)));
    }
    return numerator;
  }

  int get_l_bound() const { return back.get_u_bound(); }
  int get_u_bound() const { return fore.get_u_bound(); }
  bool has_lower_bound() const { return back.has_lower_bound(); }
  bool has_upper_bound() const { return fore.has_upper_bound(); }

  double weight(int t) {
    double w = 0.0;
    if (t <= back.get_u_bound()) w = back.weight(t);
    if (t >= fore.get_l_bound()) w = w + fore.weight(t);
    return w;
  }

  void prepare(int a, int b) {
    back.prepare(a);
    fore.prepare(b);
  }

  // RationalFilter.weights() (pinned RationalFilter.txt 357 + lambda
  // lambda$weights$0 at 426): the IntToDoubleFunction bound to weight(t) -
  // i.e. exactly weight(t). Mirrored for surface completeness; the WK
  // builder calls weight(t) directly.
  std::function<double(int)> weights() { return [this](int t) { return weight(t); }; }
};

}  // namespace seats_lfilters
}  // namespace detail
}  // namespace kalos
