#pragma once

// kalos M5 slice 24-C.1, cluster 2: the acf/solver mirror backing the WK
// extraction machine. Transcription of the pinned 3.9.0 bytecode of
//   jdplus.toolkit.base.core.arima.AutoCovarianceFunction
//   internal.toolkit.base.core.arima.AutoCovarianceComputers
//     (defaultComputer / defaultSymmetricComputer lambda bodies)
//   jdplus.toolkit.base.core.arima.LinearProcess
//   jdplus.toolkit.base.core.arima.AbstractArimaModel.getAutoCovarianceFunction
//     + isStationary
//   internal.toolkit.base.core.math.linearfilters.SymmetricFilterAlgorithms
//     .decomposer(null) lambda  (SymmetricFilter.decompose(BackFilter) -
//     the need for this was PROVEN in cluster 2: WienerKolmogorovEstimators
//     calls SymmetricFilter.decompose directly, pinning note "prove the
//     need in cluster 3/5" from the seats_linear_filters.hpp header closed)
//   LinearSystemSolver.robustSolver() == LULinearSystemSolver2(1.0E-13)
//     (LS_Factory BootstrapMethods entry 1 -> lambda$static$1 - the SAME
//     solver instance the cluster-1 lu2_solve mirrors; fastSolver() is the
//     LULinearSystemSolver builder path and is NOT used by any call site
//     ported here).
// javap dumps: /tmp/opencode/ssf/{jdplus.toolkit.base.core.arima.
// AutoCovarianceFunction.txt, AcfComputers.txt, jdplus.toolkit.base.core.
// arima.LinearProcess.txt, AbstractArimaModel.txt, SymFilterAlgos.txt,
// LS_Factory.txt, LUL2.txt, ArimaModel.txt (regenerated 2026-11-22 from the
// pinned jar - the staged RationalFunction.txt was a broken 1-line artifact;
// `of(Polynomial, Polynomial)' re-verified FRESH: iconst_0 -> NO simplify,
// so the cluster-1 RatFunc::of(simplify = false) is pinned correct and the
// SECTION 14 memo claiming of() simplifies was wrong).
//
// Fidelity notes (pinned against the dumps, 2026-11-22):
//  * AutoCovarianceFunction caches in 36-quantized blocks (BLOCK = 36):
//    prepare(0) -> 36, else round n UP to the next multiple of 36 (Java %-
//    sign semantics kept: C++ % truncates identically). get(i) returns 0.0
//    past the cache; get(i < 0) passes the cache compare and the JVM
//    daload AIOOBEs - refused here AFTER prepare ran (AIOOBE parity).
//  * The default computer scales by ivar only when ivar != 1.0 (dcmpl ifeq:
//    NaN ivar multiplies); the extension recurrence ac[i] =
//    -sum(ar[j] * ac[i-j]) runs on the COPIED array so new entries feed
//    later ones, matching the JVM's use of the new buffer.
//  * The Toeplitz rhs is a DataBlock VIEW of the acf array: the solver
//    writes ac[0..m) only. Mirrored through a temporary of exactly m
//    entries (lu2_solve's y.size() == n contract).
//  * FastMatrix.add(row, col, d): storage[row + col*ncols] += d (cluster-1
//    column-major indexing). Zero-skip is `c != 0.0' with dcmpl ifeq - a
//    NaN coefficient DOES enter the matrix.
//  * Computer failure path: MatrixException -> ArimaException("acgf of
//    non stationary model") (exception table 215-244). LinearProcess never
//    wraps; SymmetricFilter.decompose lets MatrixException escape.
//  * LinearProcess.doStationary()/getSpectrum() are NOT on the WK path
//    (WKToolkit call-site inventory, 2026-11-22) and stay unported; the
//    consuming cluster fails closed by their absence.
//  * AutoCovarianceFunction(SymmetricFilter, Polynomial) (the 2-arg ctor)
//    is referenced by NO class in the pinned toolkit/sa/tramoseats jars
//    (constant-pool scan, 2026-11-22) - ported for mirror completeness
//    because the SymmetricComputer's q > 0 branch is the only consumer of
//    sf_decompose inside the acf surface. It fails closed where the JVM
//    would hit null fields (values()/getBound() with ma == null).
//  * RationalFunction.coefficients(n) length is ALWAYS n (cluster-1 mirror,
//    golden-verified), so the computers' `w.size() == p + 1' and the q > 0
//    symmetric branch's `buf.length == n + 1 > n -> return buf' hold.

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "kalos/detail/seats_linear_filters.hpp"

namespace kalos {
namespace detail {
namespace seats_acf {

using seats_lfilters::BackFilter;
using seats_lfilters::FinFilter;
using seats_lfilters::lu2_solve;
using seats_lfilters::matrix_exception;
using seats_lfilters::pg;
using seats_lfilters::RatFunc;
using seats_lfilters::RationalFilter;
using seats_lfilters::SymmetricFilter;
using seats_model::AModel;
using seats_model::BFilter;
using seats_model::idx;
using seats_model::JvmNpe;
using seats_model::p_of;
using seats_model::p_times;

// ArimaException mirror (jdplus ...core.arima.ArimaException).
struct arima_exception : std::runtime_error {
  explicit arima_exception(const std::string& m) : std::runtime_error(m) {}
};

// ---------------------------------------------------------------------------
// SymmetricFilterAlgorithms.decomposer(null) lambda: spectral decomposition
// of a SymmetricFilter against a BackFilter (lambda$decomposer$0). Matrix-
// Exception (matrix_exception) escapes, as in SymmetricFilter.decompose.
// ---------------------------------------------------------------------------
inline BackFilter sf_decompose(const SymmetricFilter& sma, const BackFilter& fb) {
  if (fb.length() == 1) {
    vd c = sma.coefficients_as_polynomial().c;  // coefficientsAsPolynomial().toArray()
    c[0] = c[0] / 2.0;
    double b = std::fabs(fb.get(0));  // Math.abs(fb.get(0))
    if (!(b == 1.0))                  // dcmpl ifeq: NaN divides through
      for (std::size_t i = 0; i < c.size(); ++i) c[i] = c[i] / b;
    return BackFilter::of_internal(c);
  }
  Poly pfb = fb.as_polynomial();
  Poly smp = sma.coefficients_as_polynomial();
  int p = pfb.degree(), q = smp.degree();
  int n = p > q ? p : q;
  int dim = n + 1;
  vd A(static_cast<std::size_t>(dim) * static_cast<std::size_t>(dim), 0.0);
  vd rhs(static_cast<std::size_t>(dim), 0.0);
  for (int i = 0; i <= n; ++i) {
    rhs[idx(n - i)] = (i <= q) ? pg(smp, i) : 0.0;
    for (int j = 0; j <= i; ++j) {
      if (i - j <= p) {
        double c = pg(pfb, i - j);
        A[idx(i + j * dim)] = A[idx(i + j * dim)] + c;  // FastMatrix.add(i, j, c)
      }
      int u = (n - i) + j;
      if (u <= p) {
        double c = pg(pfb, u);
        A[idx(i + (n - j) * dim)] = A[idx(i + (n - j) * dim)] + c;  // add(i, n-j, c)
      }
    }
  }
  lu2_solve(A, dim, rhs);  // robustSolver(); MatrixException escapes
  std::reverse(rhs.begin(), rhs.end());  // DataBlock.reverse().toArray()
  return BackFilter::of_internal(rhs);
}

// ---------------------------------------------------------------------------
// AutoCovarianceComputers.defaultComputer(null): Computer.ac(num, den, n).
// Returns the fresh (pre-ivar-scale) acf block of length n + 1 with the AR
// extension. Throws arima_exception("acgf of non stationary model") when
// the Toeplitz system is singular (the pinned MatrixException translation).
// ---------------------------------------------------------------------------
inline vd ac_computer(const Poly& ma, const Poly& ar, int n) {
  int q = ar.degree(), p = ma.degree();
  int m = (q > p ? q : p) + 1;
  if (n < m) n = m;
  vd out(static_cast<std::size_t>(n) + 1, 0.0);
  RatFunc rf = RatFunc::of(ma, ar);  // of(): NO simplify (fresh-verified)
  vd w = rf.coefficients(p + 1);     // length p + 1 exactly
  for (int i = 0; i <= p; ++i) {     // rhs via the DataBlock view (out[0..m))
    double s = 0.0;
    for (int j = i; j <= p; ++j) s += pg(ma, j) * w[idx(j - i)];
    out[idx(i)] = s;
  }
  vd A(static_cast<std::size_t>(m) * static_cast<std::size_t>(m), 0.0);
  for (int i = 0; i < m; ++i)
    for (int k = 0; k <= q; ++k) {
      double c = pg(ar, k);
      if (c != 0.0) {  // dcmpl ifeq skip; NaN enters
        int col = (k > i) ? k - i : i - k;
        A[idx(i + col * m)] = A[idx(i + col * m)] + c;
      }
    }
  vd y(out.begin(), out.begin() + m);  // the DataBlock view
  try {
    lu2_solve(A, m, y);
  } catch (const matrix_exception&) {
    throw arima_exception("acgf of non stationary model");
  }
  for (int i = 0; i < m; ++i) out[idx(i)] = y[idx(i)];
  for (int i = m; i <= n; ++i) {
    double s = 0.0;
    for (int j = 1; j <= q; ++j) s += pg(ar, j) * out[idx(i - j)];
    out[idx(i)] = -s;
  }
  return out;
}

// ---------------------------------------------------------------------------
// AutoCovarianceComputers.defaultSymmetricComputer(null): SymmetricComputer
// .ac(sma, den, n). See the header notes for the q > 0 return-early shape.
// ---------------------------------------------------------------------------
inline vd ac_sym_computer(const SymmetricFilter& sma, const Poly& den, int n) {
  int q = den.degree();
  int ub = sma.get_upper_bound();
  int m = (q > ub ? q : ub) + 1;
  if (n < m) n = m;
  vd buf;
  if (q == 0) {
    buf = sma.coefficients_as_polynomial().c;  // toArray()
  } else {
    BackFilter ar_bf{den};  // new BackFilter(Polynomial): raw ctor
    BackFilter bf = sf_decompose(sma, ar_bf);
    RatFunc rf = RatFunc::of(bf.p, den);
    buf = rf.coefficients(n + 1);  // length n + 1
    buf[0] = buf[0] * 2.0;
  }
  if (n < static_cast<int>(buf.size())) return buf;
  int len = static_cast<int>(buf.size());
  vd out(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < len; ++i) out[idx(i)] = buf[idx(i)];
  if (q > 0)
    for (int i = len; i < n; ++i) {
      double s = 0.0;
      for (int j = 1; j <= q; ++j) s += pg(den, j) * out[idx(i - j)];
      out[idx(i)] = -s;
    }
  return out;
}

// ---------------------------------------------------------------------------
// AutoCovarianceFunction: lazy acf with the 36-quantized cache, the ivar
// scale, and the AR-recurrence extension. The JVM `synchronized'/volatile
// pair is behavior-free in the single-threaded replay.
// ---------------------------------------------------------------------------
struct Acf {
  bool has_ma = true;
  Poly ma;           // valid when has_ma (null in the 2-arg ctor)
  Poly ar;           // denominator polynomial
  bool has_sma = false;
  SymmetricFilter sma{Poly(vd{0.0})};  // valid when has_sma
  double ivar = 1.0;
  std::optional<vd> ac;

  // AutoCovarianceFunction(Polynomial ma, Polynomial ar, double ivar).
  static Acf of_polys(const Poly& ma_, const Poly& ar_, double ivar) {
    Acf f;
    f.ma = ma_;
    f.ar = ar_;
    if (ar_.degree() == 0) {
      f.has_sma = true;
      BackFilter b{ma_};  // raw BackFilter(Polynomial) ctor
      f.sma = SymmetricFilter::conv_of_d(b, ivar);
      f.ivar = 1.0;
    } else {
      f.ivar = ivar;
    }
    return f;
  }
  // AutoCovarianceFunction(SymmetricFilter, Polynomial): no live caller
  // (constant-pool scan); mirrored for surface completeness.
  static Acf of_sym(const SymmetricFilter& sma_, const Poly& ar_) {
    Acf f;
    f.has_ma = false;
    f.ar = ar_;
    f.has_sma = true;
    f.sma = sma_;
    f.ivar = 1.0;
    return f;
  }

  void prepare(int n) {
    if (n == 0) {
      n = 36;
    } else {
      int r = n % 36;  // Java % (truncating) matches C++11
      if (r != 0) n = n + (36 - r);
    }
    if (ac.has_value() && static_cast<int>(ac->size()) > n) return;
    ac = ac_extend(*this, ac, n);
  }

  double get(int i) {
    prepare(i + 1);
    if (i < 0)  // the JVM passes the cache compare, then daload AIOOBEs
      throw std::out_of_range("seats_acf_solver: get(i < 0) (AIOOBE parity)");
    if (i >= static_cast<int>(ac->size())) return 0.0;  // dconst_0
    return (*ac)[idx(i)];
  }

  vd values(int n) {
    if (n < 0) throw std::out_of_range("seats_acf_solver: negative values(n)");
    prepare(n);
    vd out(static_cast<std::size_t>(n), 0.0);
    int k = static_cast<int>(ac->size()) < n ? static_cast<int>(ac->size()) : n;
    for (int i = 0; i < k; ++i) out[idx(i)] = (*ac)[idx(i)];
    return out;
  }

  bool has_bound() const { return ar.degree() + 1 == 1; }
  int get_bound() const {
    if (!has_ma) throw JvmNpe{};  // ma.degree() on the null field
    return has_bound() ? ma.degree() + 1 : -1;
  }

 private:
  static vd ac_extend(Acf& f, const std::optional<vd>& old, int n);
};

inline vd Acf::ac_extend(Acf& self, const std::optional<vd>& old, int n) {
  vd r;
  if (!old.has_value()) {
    if (!self.has_sma) {
      r = ac_computer(self.ma, self.ar, n);
      if (!(self.ivar == 1.0))  // dcmpl ifeq: NaN ivar multiplies
        for (std::size_t i = 0; i < r.size(); ++i) r[i] = r[i] * self.ivar;
    } else {
      r = ac_sym_computer(self.sma, self.ar, n);
    }
  } else {
    r = *old;
  }
  if (static_cast<int>(r.size()) > n) return r;
  vd out(static_cast<std::size_t>(n) + 1, 0.0);
  for (std::size_t i = 0; i < r.size(); ++i) out[idx(i)] = r[idx(i)];
  int d = self.ar.degree();
  for (int i = static_cast<int>(r.size()); i <= n; ++i) {
    double s = 0.0;
    for (int j = 1; j <= d; ++j) s += pg(self.ar, j) * out[idx(i - j)];
    out[idx(i)] = -s;
  }
  return out;
}

// ---------------------------------------------------------------------------
// LinearProcess (RationalFilter + innovation variance): the lazily built
// acf is Acf(numerator.weightsToArray, back.den * fore.den, var). doStatio-
// nary()/getSpectrum() are not on the WK path (header notes).
// ---------------------------------------------------------------------------
struct LinearProc {
  RationalFilter rf;
  double var = 0.0;
  std::optional<Acf> acf;

  LinearProc() = default;
  LinearProc(RationalFilter r, double v) : rf(std::move(r)), var(v) {}

  RationalFilter& get_filter() { return rf; }
  double get_innovation_variance() const { return var; }
  bool is_invertible() const { return true; }
  bool is_stationary() const { return true; }
  bool is_null() const { return var == 0.0; }  // dcmpl ifne: NaN -> false

  Acf& get_auto_covariance_function() {
    if (!acf.has_value()) {
      Poly num = p_of(rf.get_numerator()->weights_to_array());  // ofInternal
      Poly d1 = rf.back.get_denominator().as_polynomial();
      Poly d2 = rf.fore.get_denominator().get_polynomial();
      Poly den = p_times(d1, d2);
      acf = Acf::of_polys(num, den, var);
    }
    return *acf;
  }
};

// ---------------------------------------------------------------------------
// AbstractArimaModel.getAutoCovarianceFunction on the seam's AModel (24-B
// types): isStationary() = getNonStationaryArOrder() == 0, then the lazy
// Acf(getMa().asPolynomial(), getStationaryAr().asPolynomial(), getIV()).
// get_ma()/get_iv() carry the pinned spectral-factorization caches (MaError
// = ArimaException("invalid ARIMA"), JvmNpe on the null Factorization).
// ---------------------------------------------------------------------------
inline Acf am_acf(const AModel& a) {
  if (a.non_stationary_ar().degree() != 0)
    throw arima_exception("acgf of non stationary model");
  BFilter ma = a.get_ma();
  return Acf::of_polys(ma.p, a.stationary_ar(), a.get_iv());
}

}  // namespace seats_acf
}  // namespace detail
}  // namespace kalos
