#ifndef KALOS_DETAIL_SEATS_EXACT_FORECASTS_HPP
#define KALOS_DETAIL_SEATS_EXACT_FORECASTS_HPP
// ---------------------------------------------------------------------------
// Item 24-C.1 cluster 4: jdplus.toolkit.base.core.ssf.arima.ExactArimaForecasts
// (the SSF-based forecaster BurmanEstimates delegates to) as a bit-exact C++
// mirror. Every loop, guard and float-op ORDER below was re-read from fresh
// javap -p -c dumps of the pinned 3.9.0 jars (/tmp/opencode/ssf, 2026-11-22):
// ExactArimaForecasts.txt, SsfArima.txt (stateComponent/ofStationary/
// ofNonStationary/fastInitializer+diffuseInitializer lambdas/stInitialize/
// dInitialize), jdplus...SsfArima_{ArimaData,ArimaInitialization,
// ArmaInitialization,ArimaDynamics}.txt, CkmsFilter.txt, CkmsState.txt,
// jdplus...ckms.CkmsDiffuseInitializer.txt, jdplus...ssf.UpdateInformation.txt,
// Loading1.txt, jdplus...univariate.SsfData.txt, ArimaForecasts.txt,
// ArimaModel.txt (stationaryTransformation), IFinF.txt (IFiniteFilter.apply),
// jdplus...data.DataBlock.txt, FastMatrix.txt, SymMat.txt, GenMat.txt (AB =
// plain product via the pure-java GEMM kernel - same accumulation order the
// 24-B fe_gemm mirror pins; this header carries its literal twin xf_gemm).
//
// Machine (verified against the 1160-tag xf_goldens.inc corpus):
//  * prepare(model, bmean): bmean -> bar = stat x nonstat and the model is
//    rebuilt as (stat, D1 x nonstat, D1 x ma, ivar); else bar = getAr()
//    (= stat x nonstat, same product). prepare ALWAYS returns true.
//  * fcasts(data, n) rebuilds a fresh CkmsFilter(fastInitializer(model)) on
//    every call and runs process() forward (PED(false) is a read-only sink).
//    The filter keeps a = state mean and l = the CKMS sqrt-factor of the
//    prediction covariance (P = l l^T is exact in this companion form); the
//    steady-state shortcut needs steadypos == 50 (kSteadyPos) which is
//    unreachable for our corpora (steadypos starts at -1 and each non-
//    missing step advances it by one only AFTER the steady check... it CAN
//    reach 50 for long series - mirrored exactly, including the norm2 test
//    that sets steadypos, after which only TX(a) runs).
//  * Stationary dispatch is model.isStationary() == (nonStationaryArOrder
//    == 0); the arma init is M = l = acf.values(dim), TX(l), f = acf[0],
//    t0 = 0, a = 0; the arima init pre-runs the diffuse "fast" initializer
//    (CkmsDiffuseInitializer over the SsfArima diffuseInitializer lambda:
//    a from the first d obs through the (n+d) x d unit-root filter matrix,
//    P = Sigma X Sigma^T with X = stVar(psi, acgf, var), Sigma = the
//    Toeplitz lift of 1/dif) at t0 = d, then builds l from the largest-
//    diagonal column of -f(TP T^T + V - P - l_0 l_0^T/f) (pivot > 0 path),
//    with the isZero(1e-6) -> l = 0 degradation and the SsfException
//    "Invalid fast filter" refusal mirrored as xf_ssf_exception.
//  * forecasts(n): if n >= dim the FIRST n >= dim outputs are the final
//    state a copy-then-TX-loop (the copy aliases y[0..dim-1] BEFORE the TX
//    loop mutates a - pinned); for n < dim y = a[0..n-1] raw. bmean mean is
//    the LAST output of bar.filter(final a) - recomputed inside every
//    fcasts call (the mean-overwrite quirk the {z,m}.meang goldens pin).
//  * backcasts = interface default: reverse(data) -> forecasts -> reverse.
// NaN in the working data propagates exactly like the JVM (the CkmsFilter
// here has no missing-value handling at all - error(i) subtracts blindly).
// ---------------------------------------------------------------------------
#include "kalos/detail/seats_linear_filters.hpp"
#include "kalos/detail/seats_acf_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace kalos {
namespace detail {
namespace seats_xf {

using seats_lfilters::BackFilter;
using seats_lfilters::RatFunc;
using seats_acf::Acf;

inline std::size_t idx(int i) { return static_cast<std::size_t>(i); }

// SsfException from the diffuse fast initializer (never expected on the
// corpus; the replay reports an unexpected throw as a failure).
// JVM abort classes raised by the PredictionErrorDecomposition sink that
// ExactArimaForecasts.fcasts attaches to every CkmsFilter.process run
// (ResidualsCumulator.add / DeterminantalTerm.update). The 24-C extraction
// unit surfaces them as `status error:<class>[:<msg>]' rows; no xf/bm
// golden case ever triggers them (the corpus pins that).
struct eco_abort {
  std::string message;  // "Invalid variance" | "Negative variance"
};
struct demetra_abort {};  // DemetraException(), null message

// jdplus...ssf.univariate.PredictionErrorDecomposition(false) in its sink
// role: per-OBSERVATION (eps, var) cumulation. The Burman path never
// consumes the accumulated values (the innovation variance is the spec
// echo) - only the throw behaviour is live.
struct PedSink {
  double detman = 1.0;  // DeterminantalTerm fields
  double detcar = 0.0;
  double ssqerr = 0.0;
  int n = 0;

  void det_update(double x) {
    if (!std::isfinite(x)) throw eco_abort{"Invalid variance"};
    if (!(x > 0.0)) throw eco_abort{"Negative variance"};  // dcmpg ifgt: >0
                                                           // or NaN passes
    detman *= x;
    while (detman >= 1.0) {  // dcmpl(detman,1) iflt: NaN exits
      detman *= 0.0625;
      detcar += 4.0;
    }
    if (detman != 0.0)
      while (detman <= 0.0625) {  // dcmpg ifgt: NaN exits
        detman *= 16.0;
        detcar -= 4.0;
      }
  }

  void add(double eps, double var) {  // ResidualsCumulator.add
    if (std::fabs(var) < 1e-9) {  // Math.abs dcmpg ifge: NaN takes the big
                                  // branch
      if (std::fabs(eps) < 1e-9) return;
      throw demetra_abort{};
    }
    det_update(var);
    ssqerr += eps * eps / var;
    ++n;
  }
};

struct xf_ssf_exception {
  std::string msg;
};

// ---------------------------------------------------------------------------
// Dense square matrix with the JVM column-major storage (FastMatrix square
// layout, lda == n) plus the exact-op-order kernels the filter needs.
// ---------------------------------------------------------------------------
struct SMat {
  int n = 0;
  vd a;  // a[i + j*n]
  SMat() = default;
  explicit SMat(int n_) : n(n_), a(static_cast<std::size_t>(n_) * idx(n_), 0.0) {}
  double get(int i, int j) const { return a[idx(i) + idx(j) * idx(n)]; }
  void set(int i, int j, double v) { a[idx(i) + idx(j) * idx(n)] = v; }
};

// GeneralMatrix.AB via the pure-java GEMM.addAB kernel (twin of the 24-B
// seats_model::fe_gemm pin: per column fill 0, ascending inner index
// ii: s = alpha * B[ii, j], s == 0.0 skips, col += s * A(:, ii)).
inline void xf_gemm(SMat& c, const SMat& A, const SMat& B, double alpha) {
  for (int j = 0; j < c.n; ++j) {
    for (int i = 0; i < c.n; ++i) c.set(i, j, 0.0);  // beta == 0
    for (int ii = 0; ii < c.n; ++ii) {
      double s = alpha * B.get(ii, j);
      if (s == 0.0) continue;  // dcmpl ifeq: NaN scale does NOT skip
      for (int i = 0; i < c.n; ++i) c.set(i, j, c.get(i, j) + s * A.get(i, ii));
    }
  }
}

// DataBlock.dot fast paths: ascending s += x_i * y_i (all inc variants share
// the order).
inline double db_dot(const double* x, int ix, const double* y, int iy, int len) {
  double s = 0.0;
  for (int k = 0; k < len; ++k) s += x[idx(k) * idx(ix)] * y[idx(k) * idx(iy)];
  return s;
}

// DataBlock.norm2: the BLAS2-scaled loop. Branch test is `scale dcmpg |v|,
// ifge -> av/scale branch': a NaN |v| takes the av/scale side (dcmpg pushes
// +1), so only scale < |v| (NaN-free) takes the rescale branch where
// acc = 1 + ((acc*q)*q) (bytecode groups as (acc*q)*q).
inline double db_norm2(const double* x, int ix, int len) {
  if (len == 0) return 0.0;
  if (len == 1) return std::fabs(x[0]);
  double m = 0.0, acc = 1.0;
  for (int k = 0; k < len; ++k) {
    double v = x[idx(k) * idx(ix)];
    if (v != 0.0) {  // dcmpl ifeq: NaN enters
      double av = std::fabs(v);
      double q;
      if (m < av) {  // dcmpg ifge: NaN av -> else branch (q = av/m = NaN)
        q = m / av;
        acc = 1.0 + (acc * q) * q;
        m = av;
      } else {
        q = av / m;
        acc = acc + (q * q);
      }
    }
  }
  return m * std::sqrt(acc);
}

// DataBlock.addAY(a, y): a == 0.0 returns (NaN proceeds); a == 1 add; a ==
// -1 sub; else x += a * y (a first operand).
inline void db_adday(vd& x, double a, const double* y, int iy) {
  if (a == 0.0) return;
  const int n = static_cast<int>(x.size());
  if (a == 1.0) {
    for (int i = 0; i < n; ++i) x[idx(i)] = x[idx(i)] + y[idx(i) * idx(iy)];
    return;
  }
  if (a == -1.0) {
    for (int i = 0; i < n; ++i) x[idx(i)] = x[idx(i)] - y[idx(i) * idx(iy)];
    return;
  }
  for (int i = 0; i < n; ++i) x[idx(i)] = x[idx(i)] + a * y[idx(i) * idx(iy)];
}

// DataBlock.mul(d): the ONLY skip is d == 1.0 (dcmpl ifne).
inline void db_mul(vd& x, double d) {
  if (d == 1.0) return;
  for (std::size_t i = 0; i < x.size(); ++i) x[i] = x[i] * d;
}

// SymmetricMatrix.fromLower (pure copy).
inline void from_lower(SMat& m) {
  for (int j = 0; j < m.n; ++j)
    for (int i = 0; i < j; ++i) m.set(i, j, m.get(j, i));
}

// SymmetricMatrix.xxt(x) * var (ArmaInitialization.v): for each column j,
// rows i >= j: val = x_j * x_i; mirror into (j, i).
inline SMat v_matrix(const vd& psi, double var) {
  const int n = static_cast<int>(psi.size());
  SMat v(n);
  for (int j = 0; j < n; ++j)
    for (int i = j; i < n; ++i) {
      double val = psi[idx(j)] * psi[idx(i)];
      v.set(i, j, val);
      if (i != j) v.set(j, i, val);
    }
  // FastMatrix.mul(var): d == 0 fills storage with d, else per-entry * d.
  if (var == 0.0) {
    std::fill(v.a.begin(), v.a.end(), 0.0 * var);  // preserves +-0 semantics
  } else {
    for (std::size_t i = 0; i < v.a.size(); ++i) v.a[i] = v.a[i] * var;
  }
  return v;
}

// ArmaInitialization.p0(var, acgf, psi) == ArimaInitialization.stVar:
// col0 = acgf; X(i+1,i+1) = X(i,i) - psi_i * (psi_i * var) [bytecode:
// (psi_i * psi_i) * var]; X(i+1,k+1) = X(i,k) - (psi_i * psi_k) * var;
// fromLower.
inline SMat p0_matrix(double var, const vd& acgf, const vd& psi) {
  const int n = static_cast<int>(acgf.size());
  SMat x(n);
  for (int i = 0; i < n; ++i) x.set(i, 0, acgf[idx(i)]);
  for (int i = 0; i + 1 < n; ++i) {
    double a = psi[idx(i)];
    x.set(i + 1, i + 1, x.get(i, i) - (a * a) * var);
    for (int k = 0; k < i; ++k)
      x.set(i + 1, k + 1, x.get(i, k) - (a * psi[idx(k)]) * var);
  }
  from_lower(x);
  return x;
}

// ArimaInitialization.sigma(S, dif): w = RationalFunction.of(ONE,
// Polynomial.of(dif)).coefficients(n); lower-triangular Toeplitz.
inline SMat sigma_matrix(const vd& dif, int n) {
  RatFunc rf = RatFunc::of(Poly(vd{1.0}), Poly(dif));
  vd w = rf.coefficients(n);
  SMat s(n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) s.set(i, j, w[idx(i - j)]);
  return s;
}

// SymmetricMatrix.XSXt(A, B, dst) = B * A * B^T: tmp = AB(B, A) via the GEMM
// kernel; for each dst column k the rows i >= k hold dot(row_i(tmp),
// row_k(B)); fromLower.
inline SMat xsxt(const SMat& A, const SMat& B) {
  SMat tmp(B.n);
  xf_gemm(tmp, B, A, 1.0);  // GeneralMatrix.AB(B, A) = B . A
  SMat dst(B.n);
  for (int k = 0; k < B.n; ++k)
    for (int i = k; i < B.n; ++i) {
      double d = db_dot(&tmp.a[idx(i)], B.n, &B.a[idx(k)], B.n, B.n);
      dst.set(i, k, d);
    }
  from_lower(dst);
  return dst;
}

// FastMatrix.addXaXt(a, x) -> addXaYt(a, x, x): FULL rank-1 update. Column
// k (y-storage ascending): w = x_k * a; x_k == 0.0 skips the column (dcmpl
// ifeq, NaN proceeds); every row i: M(i,k) = M(i,k) + w * x_i. NOT a
// triangular mirror write - the JVM updates both sides independently with
// different association ((x_k*a)*x_i vs (x_i*a)*x_k), so the result is
// asymmetric at the ulp level (pinned by the quarterly-m l vector). a ==
// 0.0 returns early (dcmpl ifne; NaN proceeds).
inline void add_xaxt(SMat& m, double a, const vd& x) {
  if (a == 0.0) return;
  const int n = m.n;
  for (int k = 0; k < n; ++k) {
    double u = x[idx(k)];
    if (u == 0.0) continue;
    double w = u * a;
    for (int i = 0; i < n; ++i) m.set(i, k, m.get(i, k) + w * x[idx(i)]);
  }
}

// FastMatrix.upLeftShift(1) on a square matrix.
inline void up_left_shift(SMat& m) {
  for (int j = 0; j + 1 < m.n; ++j)
    for (int i = 0; i + 1 < m.n; ++i) m.set(i, j, m.get(i + 1, j + 1));
}

// FastMatrix.isZero(tol) = every entry passes |v| <= tol; NaN fails
// (Math.abs(NaN) <= tol is false, like the dcmpg-guarded predicate).
inline bool mat_is_zero(const SMat& m, double tol) {
  for (double v : m.a)
    if (!(std::fabs(v) <= tol)) return false;
  return true;
}

// ---------------------------------------------------------------------------
// ArimaDynamics (shared by the stationary and diffuse paths).
// ---------------------------------------------------------------------------
struct ArimaDyn {
  vd phi;   // ArimaData.phi = getAr().asPolynomial().toArray()
  int dim;  // ArimaData.dim = max(arPoly.degree, maPoly.degree + 1)

  // ISsfDynamics.TX: s = 0; for i in 1..len(phi)-1: s -= phi[i] * x[len-i]
  // (reverseReader cursor starts at the last entry); bshift(1) (elements
  // slide left); x[dim-1] = s. NOTE phi.length may exceed dim: the JVM
  // indexes x storage-relative, and dim >= phi.length - 1 always holds.
  void tx(double* x) const {
    double s = 0.0;
    for (int i = 1; i < static_cast<int>(phi.size()); ++i) s -= phi[idx(i)] * x[idx(dim - i)];
    for (int k = 0; k + 1 < dim; ++k) x[idx(k)] = x[idx(k + 1)];
    x[idx(dim - 1)] = s;
  }

  // ISsfDynamics.TVT(M): phi.length == 1 -> upLeftShift + zero the last row
  // and column; else build z = sum_{i=1..p} (-phi_i) * col(dim - i), advance
  // z through TX, upLeftShift, then write z into the last column and row.
  void tvt(SMat& m) const {
    if (phi.size() == 1) {
      up_left_shift(m);
      for (int i = 0; i < dim; ++i) {
        m.set(i, dim - 1, 0.0);
        m.set(dim - 1, i, 0.0);
      }
      return;
    }
    vd z(static_cast<std::size_t>(dim), 0.0);
    int p = static_cast<int>(phi.size());
    for (int i = 1; i < p; ++i) {
      // reverseColumnsIterator.next() walks columns dim-1, dim-2, ...
      int col = dim - i;
      db_adday(z, -phi[idx(i)], &m.a[idx(col)], m.n);
    }
    tx(z.data());
    up_left_shift(m);
    for (int r = 0; r < dim; ++r) {
      m.set(r, dim - 1, z[idx(r)]);
      m.set(dim - 1, r, z[idx(r)]);
    }
  }
};

// ---------------------------------------------------------------------------
// ExactArimaForecasts.
// ---------------------------------------------------------------------------
struct XfModel {
  Poly stat, nonstat, ma;  // getStationaryAr/getNonStationaryAr/getMa
  double var = 0.0;        // getInnovationVariance
  // SarimaModel.getAr multiplies the differencing filter first.
  // A mean-augmented model is rebuilt as ArimaModel and uses the inverse order.
  bool sarima = false;
  Poly full_ar() const { return sarima ? p_times(nonstat, stat) : p_times(stat, nonstat); }
};

class ExactFc {
 public:
  // prepare(model, bmean) - always true (pinned).
  bool prepare(const XfModel& model, bool bmean) {
    bmean_ = bmean;
    XfModel m = model;
    bar_ = bmean ? p_times(m.stat, m.nonstat) : m.full_ar();
    if (bmean) {
      Poly d1(vd{1.0, -1.0});  // BackFilter.D1.asPolynomial()
      // Bytecode: nonstat.times(D1), ma.times(D1) (D1 is the RIGHT operand).
      m.sarima = false;
      m.nonstat = p_times(m.nonstat, d1);
      m.ma = p_times(m.ma, d1);
    }
    m_ = m;
    return true;
  }

  double get_mean() const { return mean_; }

  // forecasts(DoubleSeq, n) (fresh filter each call).
  vd forecasts(const vd& data, int n) { return fcasts(data, n); }

  // ArimaForecasts.backcasts default: reverse -> forecasts -> reverse.
  vd backcasts(const vd& data, int n) {
    vd rev(data.rbegin(), data.rend());
    vd f = fcasts(rev, n);
    return vd(f.rbegin(), f.rend());
  }

  int state_dim() const {
    Poly full = m_.full_ar();
    int d = std::max(full.degree(), m_.ma.degree() + 1);
    return d < 1 ? 1 : d;
  }

 private:
  Poly bar_;
  bool bmean_ = false;
  XfModel m_;
  double mean_ = 0.0;

  vd fcasts(const vd& data, int n) {
    Poly full = m_.full_ar();
    ArimaDyn dyn;
    dyn.phi = full.c;
    dyn.dim = std::max(full.degree(), m_.ma.degree() + 1);
    const int dim = dyn.dim;
    const int ndata = static_cast<int>(data.size());

    vd a(static_cast<std::size_t>(dim), 0.0);  // CkmsState.a
    vd l(static_cast<std::size_t>(dim), 0.0);  // CkmsState.l
    vd M(static_cast<std::size_t>(dim), 0.0);  // UpdateInformation.M
    double f = 0.0;
    int t0 = 0;

    if (m_.nonstat.degree() == 0) {
      // stInitialize: acf.values(dim); M = l = acf; TX(l); f = acf[0].
      Acf acf = Acf::of_polys(m_.ma, m_.stat, m_.var);
      vd ac = acf.values(dim);
      M = ac;
      l = ac;
      dyn.tx(l.data());
      f = ac[idx(0)];
    } else {
      diffuse_init(data, dyn, a, l, M, f, t0);
    }

    // CkmsFilter.process loop body (initialize(): steadypos = -1, eps =
    // 1e-15, neps = eps * pe.variance). ExactArimaForecasts.fcasts passes
    // a FRESH PredictionErrorDecomposition(false) per call; process saves
    // (error(i), pre-update variance) with OBSERVATION status before
    // update(i) (no missing values on the ported domain, pinned 24-C).
    const double eps = 1.0e-15;
    const double neps = eps * f;
    int steadypos = -1;
    PedSink ped;
    for (int i = t0; i < ndata; ++i) {
      double e = data[idx(i)] - a[idx(0)];  // error(i) (no missing handling)
      ped.add(e, f);                        // PED.save(OBS) before update
      db_adday(a, e / f, M.data(), 1);      // update(i)
      if (steadypos < 0) {
        double zx = l[idx(0)];
        double azx = std::fabs(zx);
        if (azx > neps) {  // dcmpl ifle: NaN zx takes the else branch
          double v = f;
          double h = zx / v;
          v = v - zx * h;
          for (int k = 0; k < dim; ++k) {
            double lk = l[idx(k)];
            l[idx(k)] = l[idx(k)] - M[idx(k)] * h;
            M[idx(k)] = M[idx(k)] - lk * h;
          }
          f = v;
        } else {
          if (!(db_norm2(l.data(), 1, dim) > eps)) steadypos = i;
        }
        dyn.tx(l.data());
      }
      dyn.tx(a.data());
    }

    // Output extraction (fcasts bytecode).
    vd y(static_cast<std::size_t>(n < 0 ? 0 : n), 0.0);
    if (n >= dim) {
      for (int i = 0; i < dim; ++i) y[idx(i)] = a[idx(i)];
      for (int i = dim; i < n; ++i) {
        dyn.tx(a.data());
        y[idx(i)] = a[idx(dim - 1)];
      }
    } else {
      for (int i = 0; i < n; ++i) y[idx(i)] = a[idx(i)];
    }

    if (bmean_) {
      // mean = last entry of bar applied to the FINAL state a. apply() works
      // through weightsToArray() = REVERSED polynomial coefficients
      // (BackFilter.getLowerBound = -deg, weights k -> p.get(-k)), i.e.
      // tmp[m] = sum_j p_{deg-j} * a[m+j]. Only palindromic bars (airline,
      // quarterly) read the same forwards - the seasonal/trend_ma/ar1 goldens
      // pin the reversal (a plain c[j]-order dot flips their mean's sign).
      const int deg = bar_.degree();
      const int mlen = dim - deg;
      if (mlen <= 0)
        throw std::out_of_range("seats_exact_forecasts: bar longer than state");
      double last = 0.0;
      for (int k = 0; k < mlen; ++k) {
        double s = 0.0;
        for (int j = 0; j <= deg; ++j) s += bar_.c[idx(deg - j)] * a[idx(k + j)];
        last = s;
      }
      mean_ = last;
    }
    return y;
  }

  // dInitialize: CkmsDiffuseInitializer(SsfArima.diffuseInitializer(model))
  // .initializeFilter(...), both bodies decoded in SsfArima.txt.
  void diffuse_init(const vd& data, ArimaDyn& dyn, vd& a, vd& l, vd& M,
                    double& f, int& t0) {
    const int dim = dyn.dim;
    const vd& dif = m_.nonstat.c;  // getUnitRoots().asPolynomial().toArray()
    const int d = m_.nonstat.degree();
    if (d == 0) throw std::logic_error("seats_exact_forecasts: diffuse d=0");

    // ArimaInitialization fields (stationaryTransformation(): the stationary
    // model is ArimaModel(stat, ONE, ma, var)).
    Acf stacf = Acf::of_polys(m_.ma, m_.stat, m_.var);
    vd stacgf = stacf.values(dim);
    RatFunc rfpsi = RatFunc::of(m_.ma, m_.stat);
    vd stpsi = rfpsi.coefficients(dim);

    // diffuseInitializer lambda: the (d+n) x d unit-root filter matrix.
    vd Mf(static_cast<std::size_t>(dim + d) * idx(d), 0.0);
    auto mf = [&](int r, int j) -> double& {
      return Mf[idx(r) + idx(j) * idx(dim + d)];  // column-major, stride = nrows
    };
    for (int j = 0; j < d; ++j) mf(j, j) = 1.0;
    for (int j = 0; j < d; ++j) {
      for (int rl = d; rl < d + dim; ++rl) {
        double s = 0.0;
        for (int i = 1; i <= d; ++i) s -= dif[idx(i)] * mf(rl - i, j);
        mf(rl, j) = s;
      }
    }
    for (int i = 0; i < dim; ++i) {
      double s = 0.0;
      for (int j = 0; j < d; ++j) {
        if (j >= static_cast<int>(data.size()))
          throw std::out_of_range("seats_exact_forecasts: diffuse seed (IOOBE parity)");
        s += mf(i + d, j) * data[idx(j)];
      }
      a[idx(i)] = s;
    }
    SMat P = xsxt(p0_matrix(m_.var, stacgf, stpsi), sigma_matrix(dif, dim));
    t0 = d;

    // CkmsDiffuseInitializer.initializeFilter tail.
    double var = P.get(0, 0);  // ZVZ(0, P) (measurement error is null)
    f = var;
    for (int i = 0; i < dim; ++i) M[idx(i)] = P.get(0, i);  // ZM -> row 0
    l = M;
    dyn.tx(l.data());  // dynamics.TX(0, l) (pe.M storage is independent)

    SMat Pc = P;  // deepClone
    dyn.tvt(Pc);  // TVT(0, Pc)
    {
      RatFunc rf = RatFunc::of(m_.ma, Poly(dyn.phi));  // ArimaData.psi: ma/ar
      vd psi = rf.coefficients(dim);
      SMat V = v_matrix(psi, m_.var);
      for (std::size_t q = 0; q < Pc.a.size(); ++q) Pc.a[q] = Pc.a[q] + V.a[q];
    }
    for (std::size_t q = 0; q < Pc.a.size(); ++q) Pc.a[q] = Pc.a[q] - P.a[q];
    add_xaxt(Pc, -1.0 / var, l);
    const double negvar = -var;
    if (negvar == 0.0) {
      std::fill(Pc.a.begin(), Pc.a.end(), 0.0 * negvar);
    } else {
      for (std::size_t q = 0; q < Pc.a.size(); ++q) Pc.a[q] = Pc.a[q] * negvar;
    }

    int pIdx = 0;
    double pivot = Pc.get(0, 0);
    for (int j = 1; j < dim; ++j) {
      double v = Pc.get(j, j);
      if (v > pivot) {  // dcmpl ifle: strict >, NaN skips
        pIdx = j;
        pivot = v;
      }
    }
    if (pivot > 0.0) {
      for (int i = 0; i < dim; ++i) l[idx(i)] = Pc.get(i, pIdx);
      db_mul(l, std::sqrt(1.0 / pivot));
    } else {
      if (!mat_is_zero(Pc, 1.0e-6))
        throw xf_ssf_exception{"Invalid fast filter"};
      std::fill(l.begin(), l.end(), 0.0);
    }
  }
};

}  // namespace seats_xf
}  // namespace detail
}  // namespace kalos

#endif  // KALOS_DETAIL_SEATS_EXACT_FORECASTS_HPP
