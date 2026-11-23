// Kalos - TRAMO-SEATS wk extraction (roadmap item 24-C.1, cluster 5).
//
// Bit-exact mirror of jdplus.tramoseats.base.core.seats.BurmanEstimates
// (BurmanEstimates.txt, freshly dumped via javap -p -c from the pinned
// jdplus-tramoseats-base-core-3.9.0.jar, 2026-11-22; the Builder defaults
// per BurmanBuilder.txt: ser = 1.0, mcmp = 0, bmean = false, nf = nb = 0).
// Auxiliary kernels re-read from fresh dumps the same day:
//   Gauss.txt                  Gauss.decompose(m) / decompose(m, scale)
//   LUDecomposition_fresh.txt  constructor / solve(DataBlock) / swap
//   LowerTri.txt               solveL1x(FastMatrix, double[], int, int)
//   UpperTri.txt               solveUx(FastMatrix, double[], int, int, d)
//   FastMatrix.txt             toArray / top / square / subDiagonal
//   MatrixWindow.txt           constructor / vnext (MUTATES the shared cur)
//   DataBlock.txt              set(double) fill / of(double[])
//   Polynomial.txt             isIdentity / evaluateAt / times(double)
//   BackFilter.txt             isIdentity / divide / times / SimplifyingTool
//
// Numerics doctrine: every loop order, guard and float-op pairing below is
// transcribed from the pinned bytecode; nothing is "improved":
//  * Gauss.decompose(m) passes Constants.getEpsilon() as scale but the
//    3.9.0 body NEVER READS the scale local (verified: no dload_2 in the
//    method) - the singularity test is the exact `pivot == 0.0' (dcmpl
//    ifne; NaN pivots pass through and divide). Row pivoting picks the
//    FIRST strict max via dcmpl ifle (NaN candidates never promoted).
//  * UpperTriangularMatrix.solveUx with diagEps 0.0: |u| > 0.0 test is
//    dcmpg (NaN u takes the DIVISION branch); the singular branch throws
//    MatrixException m_err_sing unless |x_i| dcmpl-lt 0 (only NaN takes
//    the store-0 escape).
//  * calc() mirrors both Burman solves, the g-filter applications (one
//    forward over z, one backward), the in-place MA recursions (divide by
//    ma[0], never multiply), the AR recursions around the ±mean/2 seeds
//    and the DoubleSeq.of VIEW semantics (the meanc tail-addition at
//    bytecode 1114 mutates the buffer the estimates view was created on -
//    the mirror mutates before materializing views, same final values).
//  * JVM exceptions modeled as data: null estimates/forecasts/backcasts
//    slots (kernel-refused components) return null on the *T getters and
//    NPE (JvmNpe) on the *F fastOp path; stdev getters catch exactly
//    {ArimaException, MatrixException} -> EMPTY (estimates) / null
//    (forecasts/backcasts) per their Exception tables.

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "seats_exact_forecasts.hpp"
#include "seats_wk_estimators.hpp"

namespace kalos {
namespace detail {
namespace seats_burman {

using seats_acf::arima_exception;
using seats_lfilters::BackFilter;
using seats_lfilters::matrix_exception;
using seats_lfilters::SymmetricFilter;
using seats_model::AModel;
using seats_model::BackSimp;
using seats_model::back_simplify;
using seats_model::BFilter;
using seats_model::idx;
using seats_model::JvmNpe;
using seats_model::p_times_d;
using seats_model::sf_is_null;
using seats_model::UcarimaModel;
using seats_wk::am_ma_bf;
using seats_wk::am_non_stationary_ar;
using seats_wk::WkEstimators;
using seats_xf::ExactFc;
using seats_xf::XfModel;

// BackFilter.isIdentity -> Polynomial.isIdentity: coeff.length == 1 &&
// coeff[0] == 1.0 (dcmpl ifne: NaN is not the identity).
inline bool bf_is_identity(const BackFilter& f) {
  return f.p.c.size() == 1 && f.p.c[0] == 1.0;
}

// ---------------------------------------------------------------------------
// Gauss.decompose(m) -> LUDecomposition over the clone of m's column-major
// storage (FastMatrix.toArray(); the result wraps the SAME array, builder
// nrows == ncolumns == n). The scale argument of the 2-arg overload is
// provably unused by the 3.9.0 bytecode; decompose(m) passes
// Constants.getEpsilon() - irrelevant here.
// ---------------------------------------------------------------------------
struct BurmanLu {
  int n = 0;
  vd lu;                // column-major, element (i, j) = lu[i + j*n]
  std::vector<int> piv; // EMPTY models the JVM null pivot array

  // Gauss.decompose(FastMatrix.square(p+q)) over `a` (taken by value: the
  // JVM runs on toArray()'s clone). Throws matrix_exception("m_err_sing").
  void decompose(vd a, int n_) {
    n = n_;
    const int ncols = n;
    std::vector<int> pv;  // null until the first row swap
    int k = 0;
    int start = 0;
    while (k < ncols) {
      int pi = k;
      double d = std::fabs(a[idx(start + k)]);
      for (int i = k + 1; i < ncols; ++i) {
        double d2 = std::fabs(a[idx(start + i)]);
        if (d2 > d) {  // dcmpl ifle: NaN candidate never promoted
          pi = i;
          d = d2;
        }
      }
      if (pi != k) {
        if (pv.empty()) {
          pv.resize(idx(ncols));
          for (int i = 0; i < ncols; ++i) pv[idx(i)] = i;
        }
        int ia = pi, ib = k;
        for (int j = 0; j < ncols; ++j) {  // full rows (all columns)
          double t = a[idx(ia)];
          a[idx(ia)] = a[idx(ib)];
          a[idx(ib)] = t;
          ia += ncols;
          ib += ncols;
        }
        int t = pv[idx(pi)];
        pv[idx(pi)] = pv[idx(k)];
        pv[idx(k)] = t;
      }
      double pivval = a[idx(start + k)];
      if (pivval == 0.0)  // dcmpl ifne: NaN passes and divides
        throw matrix_exception{"m_err_sing"};
      for (int i = k + 1; i < ncols; ++i) {
        a[idx(start + i)] = a[idx(start + i)] / pivval;  // (i, col k) / pivot
        const double f = a[idx(start + i)];
        int off = (k + 1) * ncols;  // bytecodes read a[13 + 17] WITHOUT the
                                    // start offset: rows k+1.., cols k+1..
        for (int j = k + 1; j < ncols; ++j) {
          a[idx(i + off)] = a[idx(i + off)] - f * a[idx(k + off)];
          off += ncols;
        }
      }
      ++k;
      start += ncols;
    }
    lu = std::move(a);
    piv = std::move(pv);
  }

  // LUDecomposition.solve(DataBlock.of(x)) with DataBlock (storage, 0, 1):
  // swap (only when a pivot array exists), then solveL1x + solveUx(0.0),
  // both on the inc == 1 fast paths.
  void solve(vd& x) const {
    if (!piv.empty()) {  // swap: x_i = copy(pivot_i)
      vd cp = x;
      for (int i = 0; i < n; ++i) x[idx(i)] = cp[idx(piv[idx(i)])];
    }
    for (int r = 0; r < n; ++r) {  // solveL1x (unit lower, ascending rows)
      double xi = x[idx(r)];
      for (int j = r + 1; j < n; ++j)
        x[idx(j)] = x[idx(j)] - xi * lu[idx(j + idx(r) * idx(n))];
    }
    for (int r = n - 1; r >= 0; --r) {  // solveUx, diagEps = 0.0
      double xi = x[idx(r)];
      double u = lu[idx(r + idx(r) * idx(n))];
      if (!(std::fabs(u) > 0.0)) {  // dcmpg ifgt: NaN u takes else? NO -
                                    // dcmpg(NaN, 0) = +1 -> ifgt TAKEN;
                                    // this branch is u == 0.0 only.
        if (!(std::fabs(xi) >= 0.0)) {  // |xi| dcmpl-lt 0 (NaN only): store 0
          x[idx(r)] = 0.0;
        } else {
          throw matrix_exception{"m_err_sing"};
        }
        continue;
      }
      xi = xi / u;
      x[idx(r)] = xi;
      for (int j = r - 1; j >= 0; --j)
        x[idx(j)] = x[idx(j)] - xi * lu[idx(j + idx(r) * idx(n))];
    }
  }
};

// ---------------------------------------------------------------------------
// jdplus.tramoseats.base.core.seats.BurmanEstimates.
// Construction = the private constructor: field copy, wk, initModel(),
// extendSeries(), then calc(c) for every component (bytecode order).
// ---------------------------------------------------------------------------
class BurmanEstimates {
 public:
  BurmanEstimates(const vd& data, WkEstimators& wk, bool bmean, int mcmp,
                  double ser, int nfcasts, int nbcasts)
      : data_(data), wk_(wk), bmean_(bmean), mcmp_(mcmp), ser_(ser),
        nfcasts_(nfcasts), nbcasts_(nbcasts) {
    init_model();
    extend_series();
    const int ncomp = n_components();
    for (int c = 0; c < ncomp; ++c) calc(c);
  }

  // ---- captured privates (BurmanProbe .priv row) --------------------------
  double mean = 0.0;
  double meanc = 0.0;
  int nf = 0;
  int nbc = 0;
  int nfc = 0;

  bool is_mean_correction() const { return bmean_; }

  // ---- public getters (null slots = nullopt; *F on null throws JvmNpe) ----
  std::optional<vd> estimates(int c, bool transformed) const {
    if (transformed) return est_[idx(c)];
    const vd& e = req(est_[idx(c)]);
    vd out(idx(e.size()));
    for (std::size_t k = 0; k < e.size(); ++k)
      out[k] = data_[k] - e[k];  // DoubleSeq.fastOp(data, est, a - b)
    return out;
  }
  std::optional<vd> forecasts(int c, bool transformed) const {
    if (transformed) return fc_[idx(c)];
    const vd& f = req(fc_[idx(c)]);
    vd out(idx(f.size()));
    for (std::size_t k = 0; k < f.size(); ++k)  // receiver - xfcasts[0..nfc)
      out[k] = f[k] - xfcasts_[idx(k)];
    return out;
  }
  std::optional<vd> backcasts(int c, bool transformed) const {
    if (transformed) return bc_[idx(c)];
    const vd& b = req(bc_[idx(c)]);
    const std::size_t s0 = xbcasts_.size() - idx(nbcasts_);
    vd out(idx(b.size()));
    for (std::size_t k = 0; k < b.size(); ++k)  // receiver - xbcasts tail
      out[k] = b[k] - xbcasts_[s0 + k];
    return out;
  }
  vd series_forecasts() const {  // xbcasts-free: xfcasts.range(0, nfcasts)
    return vd(xfcasts_.begin(), xfcasts_.begin() + idx(nfcasts_));
  }
  vd series_backcasts() const {  // xbcasts.drop(0, len - nbcasts)
    return vd(xbcasts_.end() - idx(nbcasts_), xbcasts_.end());
  }

  // ---- stdev getters (D11: conditional-error variance only) ---------------
  // Doubles.EMPTY on null component / caught {Arima,Matrix}Exception.
  vd stdev_estimates(int c) const {
    if (wk_.component(c).is_null()) return vd{};
    const int n = static_cast<int>(data_.size());
    try {
      std::optional<vd> var = wk_.total_error_variance(c, true, 0, (n + 1) / 2);
      if (!var.has_value()) throw JvmNpe{};  // var.length on null
      vd out(idx(n));
      for (std::size_t k = 0; k < var->size(); ++k) {
        double s = ser_ * std::sqrt((*var)[k]);
        out[k] = s;
        out[idx(n) - 1 - k] = s;
      }
      return out;
    } catch (const arima_exception&) {
      return vd{};
    } catch (const matrix_exception&) {
      return vd{};
    }
  }
  std::optional<vd> stdev_forecasts(int c, bool b) const {
    if (wk_.component(c).is_null() || nfcasts_ == 0) return std::nullopt;
    try {
      std::optional<vd> var =
          wk_.total_error_variance(c, b, -nfcasts_, nfcasts_);
      if (!var.has_value()) throw JvmNpe{};
      vd out(idx(nfcasts_));
      for (int k = 0; k < nfcasts_; ++k)
        out[idx(k)] = ser_ * std::sqrt((*var)[idx(nfcasts_ - 1 - k)]);
      return out;
    } catch (const arima_exception&) {
      return std::nullopt;
    } catch (const matrix_exception&) {
      return std::nullopt;
    }
  }
  std::optional<vd> stdev_backcasts(int c, bool b) const {
    if (wk_.component(c).is_null() || nbcasts_ == 0) return std::nullopt;
    try {
      std::optional<vd> var =
          wk_.total_error_variance(c, b, -nbcasts_, nbcasts_);
      if (!var.has_value()) throw JvmNpe{};
      vd out(idx(nbcasts_));
      for (int i = nbcasts_ - 1; i >= 0; --i)
        out[idx(i)] = ser_ * std::sqrt((*var)[idx(i)]);
      return out;
    } catch (const arima_exception&) {
      return std::nullopt;
    } catch (const matrix_exception&) {
      return std::nullopt;
    }
  }

  // Burman-order surfaces captured for cross-checks.
  const vd& xfcasts() const { return xfcasts_; }
  const vd& xbcasts() const { return xbcasts_; }

  // Test-visible internals (debug-twin parity; JVM reflection sees these).
  const std::optional<vd>& g_at(int c) const { return g_[idx(c)]; }
  const vd& z() const { return z_; }
  const vd& ar_arr() const { return ar_; }
  const vd& ma_arr() const { return ma_; }
  const BurmanLu& lu() const { return lu_; }

 private:
  static const vd& req(const std::optional<vd>& o) {
    if (!o.has_value()) throw JvmNpe{};  // fastOp on the null slot
    return *o;
  }

  int n_components() const { return static_cast<int>(wk_.um->comps.size()); }

  // ArimaModel.getStationaryAr of a component (seam: BackFilter(a.ar)).
  static BackFilter comp_stat(const AModel& a) { return BackFilter(a.ar); }

  // ---- initModel ----------------------------------------------------------
  void init_model() {
    const AModel& outer = wk_.outer;
    const int ncomp = n_components();
    est_.assign(idx(ncomp), std::nullopt);
    fc_.assign(idx(ncomp), std::nullopt);
    bc_.assign(idx(ncomp), std::nullopt);
    g_.assign(idx(ncomp), std::nullopt);

    Poly mapoly = am_ma_bf(outer).p;  // getMa().asPolynomial()
    double iv = outer.get_iv();       // getInnovationVariance()
    if (iv != 1.0)                    // dcmpl ifeq: NaN scales
      mapoly = p_times_d(mapoly, std::sqrt(iv));
    Poly arpoly = BackFilter(p_times_(outer)).p;  // getAr() = ar.times(delta)

    for (int c = 0; c < ncomp; ++c) {
      const AModel& a = wk_.component(c);
      if (a.is_null()) continue;  // g[c] stays null
      const seats_model::SymF smf = a.symmetric_ma();
      if (sf_is_null(smf)) {
        g_[idx(c)] = vd{1.0};  // ONE
        continue;
      }
      SymmetricFilter sma{Poly{smf.a}};
      BackFilter arn = am_non_stationary_ar(outer);       // outer nonstat
      BackFilter arnc = am_non_stationary_ar(a);          // comp nonstat
      BackFilter num12 = arn.divide(arnc);
      BackFilter st1 = comp_stat(outer);
      BackFilter st2 = comp_stat(a);
      BackSimp tool = back_simplify(BFilter{st1.p}, BFilter{st2.p});
      if (tool.ret) {
        st1 = BackFilter{tool.left.p};
        st2 = BackFilter{tool.right.p};
      }
      BackFilter num = num12.times(st1);
      BackFilter den = BackFilter(mapoly).times(st2);  // raw maPoly
      // Bytecode 245: RECEIVER is sma - sma.times(convolutionOf(num)). The
      // convolution accumulation order of SymmetricFilter.times is receiver-
      // dependent, so the argument order is load-bearing for bit parity.
      SymmetricFilter sf = sma.times(SymmetricFilter::conv_of(num));
      double iv2 = outer.get_iv();
      if (iv2 != 1.0) sf = sf.times(1.0 / iv2);
      BackFilter gbf = seats_acf::sf_decompose(sf, den);
      g_[idx(c)] = gbf.p.c;  // asPolynomial().toArray()
    }
    ma_ = mapoly.c;  // the SQRT-SCALED ma
    ar_ = arpoly.c;
    init_solver();
  }

  static BackFilter p_times_(const AModel& a) {
    return BackFilter(seats_model::p_times(a.ar, a.delta));
  }

  // ---- initSolver ---------------------------------------------------------
  // MatrixWindow quirk: top(0) yields cur (start 0, nrows 0); each vnext(k)
  // advances start by the PREVIOUS nrows then sets nrows = k, so the ma
  // window is (start 0, nrows p) and the ar window (start p, nrows q) on
  // the SAME n x n storage (lda = n). subDiagonal(j) fills the j-th lower
  // diagonal of the window, len = min(nrows, ncols - j), EMPTY past bounds.
  void init_solver() {
    const int q = static_cast<int>(ma_.size()) - 1;
    const int p = static_cast<int>(ar_.size()) - 1;
    const int n = p + q;
    vd s(static_cast<std::size_t>(n) * idx(n), 0.0);
    for (int i = 0; i <= q; ++i) {  // window (0, nrows p, ncols n)
      if (i >= n) continue;         // subDiagonal EMPTY guard (i >= ncols)
      int len = p < n - i ? p : n - i;
      double v = ma_[idx(i)];
      for (int r = 0; r < len; ++r) s[idx(i * n + r + idx(r) * idx(n))] = v;
    }
    for (int j = 0; j <= p; ++j) {  // window (start p, nrows q, ncols n)
      if (j >= n) continue;
      int len = q < n - j ? q : n - j;
      double v = ar_[idx(p - j)];
      for (int r = 0; r < len; ++r)
        s[idx(p + j * n + r + idx(r) * idx(n))] = v;
    }
    lu_.decompose(std::move(s), n);
  }

  // ---- extendSeries -------------------------------------------------------
  void extend_series() {
    const AModel& outer = wk_.outer;
    const int q = static_cast<int>(ma_.size()) - 1;
    const int p = static_cast<int>(ar_.size()) - 1;
    const auto padding = q > p ? 2LL * q : static_cast<long long>(p) + q;
    if (padding > std::numeric_limits<int>::max())
      throw std::invalid_argument("seats_extract_wk: filter padding exceeds the index range");
    nf = static_cast<int>(padding);
    nbc = nf > nbcasts_ ? nf : nbcasts_;
    nfc = nf > nfcasts_ ? nf : nfcasts_;
    const auto span = static_cast<std::uint64_t>(nbc) + nfc + data_.size();
    if (span > std::numeric_limits<int>::max())
      throw std::invalid_argument("seats_extract_wk: extended series exceeds the monthly index range");

    XfModel xm;
    xm.stat = outer.ar;
    xm.nonstat = outer.delta;
    xm.ma = am_ma_bf(outer).p;
    xm.var = outer.get_iv();
    ExactFc fx;
    fx.prepare(xm, bmean_);  // return value popped by the JVM
    xfcasts_ = fx.forecasts(data_, nfc);
    xbcasts_ = fx.backcasts(data_, nbc);
    mean = bmean_ ? fx.get_mean() : 0.0;

    const int n = static_cast<int>(data_.size());
    z_.assign(idx(n) + idx(nbc) + idx(nfc), 0.0);
    for (int k = 0; k < n; ++k) z_[idx(nbc + k)] = data_[idx(k)];
    for (int k = 0; k < nfc; ++k) z_[idx(nbc) + idx(n) + idx(k)] = xfcasts_[idx(k)];
    for (int k = 0; k < nbc; ++k) z_[idx(k)] = xbcasts_[idx(k)];

    meanc = mean / p_value_at(comp_stat(outer).p, 1.0);
    if (use_mean())  // bmean && outer nonstat order 0
      for (std::size_t k = 0; k < z_.size(); ++k) z_[k] = z_[k] - meanc;
  }

  bool use_mean() const {
    return bmean_ && wk_.outer.delta.degree() == 0;
  }

  bool is_trend_constant() const { return wk_.component(mcmp_).is_null(); }

  // ---- calc(int) ----------------------------------------------------------
  void calc(int comp) {
    const int n = static_cast<int>(data_.size());
    double d3 = 0.0, d5 = 0.0, d7 = 0.0;
    if (mcmp_ == comp) {
      if (is_trend_constant()) {
        est_[idx(comp)] = vd(idx(n), meanc);
        fc_[idx(comp)] = vd(idx(nfcasts_), meanc);
        bc_[idx(comp)] = vd(idx(nbcasts_), meanc);
        return;
      }
      BackFilter arn = am_non_stationary_ar(wk_.component(comp));
      if (mean != 0.0) {  // dcmpl ifeq; NaN proceeds
        if (bf_is_identity(arn)) {
          d7 = meanc;
        } else {
          d5 = mean;
          double prod = arn.get(0) * arn.get(arn.degree());
          // dcmpg ifge at bytecode 158-161: dcmpg(NaN, 0) = +1 takes the
          // ifge branch -> NaN KEEPS +mean; only a strict negative product
          // negates.
          d3 = prod < 0.0 ? -mean : mean;
        }
      }
    }
    if (!g_[idx(comp)].has_value()) return;
    const vd& gv = *g_[idx(comp)];
    const int q = static_cast<int>(ma_.size()) - 1;
    const int p = static_cast<int>(ar_.size()) - 1;
    const int pq = p + q;
    const int gmax = static_cast<int>(gv.size()) - 1;

    const int wlen = nbc + n + q;
    vd w(idx(wlen));
    for (int j = 0; j < wlen; ++j) {  // FORWARD g application
      double s = gv[0] * z_[idx(j)];
      for (int k = 1; k <= gmax; ++k) s += gv[idx(k)] * z_[idx(j + k)];
      w[idx(j)] = s;
    }
    vd x(idx(pq));
    int i19 = nbc + n + q - p;
    for (int j = 0; j < p; ++j) x[idx(j)] = w[idx(i19 + j)];
    for (int j = p; j < pq; ++j) x[idx(j)] = d3 / 2.0;
    lu_.solve(x);
    vd y(idx(nbc) + idx(n) + idx(nfc), 0.0);
    const int base = nbc + n + q - p;
    for (int j = 0; j < pq; ++j) y[idx(base + j)] = x[idx(j)];
    for (int i = base - 1; i >= 0; --i) {
      double s = w[idx(i)];
      for (int k = 1; k <= q; ++k) s -= y[idx(i + k)] * ma_[idx(k)];
      y[idx(i)] = s / ma_[0];
    }
    // Bytecode 556: the AR-forward loop starts at nbc + n + 2q (== base+pq),
    // never touching the data region or the solved block.
    for (int i = nbc + n + 2 * q; i < nbc + n + nfc; ++i) {
      double s = d3 / 2.0;
      for (int k = 1; k <= p; ++k) s -= ar_[idx(k)] * y[idx(i - k)];
      y[idx(i)] = s;
    }

    vd u(idx(n) + idx(nbc) + idx(nfc));
    const int ustart = nbc - q;
    for (int i = ustart; i < n + nbc + nfc; ++i) {  // BACKWARD g application
      double s = gv[0] * z_[idx(i)];
      for (int k = 1; k <= gmax; ++k) s += gv[idx(k)] * z_[idx(i - k)];
      u[idx(i)] = s;
    }
    vd x2(idx(pq));
    int i21 = nbc + p - q;
    for (int j = 0; j < p; ++j) x2[idx(j)] = u[idx(--i21)];
    for (int j = p; j < pq; ++j) x2[idx(j)] = d5 / 2.0;
    lu_.solve(x2);
    vd v(idx(n) + idx(nbc) + idx(nfc), 0.0);
    const int base2 = nbc + p - q;
    {
      int w2 = base2;
      for (int j = 0; j < pq; ++j) v[idx(--w2)] = x2[idx(j)];
    }
    for (int i = base2; i < n + nbc + nfc; ++i) {
      double s = u[idx(i)];
      for (int k = 1; k <= q; ++k) s -= v[idx(i - k)] * ma_[idx(k)];
      v[idx(i)] = s / ma_[0];
    }
    for (int i = nbc - 2 * q - 1; i >= 0; --i) {
      double s = d5 / 2.0;
      for (int k = 1; k <= p; ++k) s -= ar_[idx(k)] * v[idx(i + k)];
      v[idx(i)] = s;
    }

    vd est(idx(n) + idx(nfc) + idx(nbc));
    for (int k = 0; k < nbc + n + nfc; ++k)
      est[idx(k)] = y[idx(k)] + v[idx(k)];
    if (d7 != 0.0)  // dcmpl ifeq; NaN adds
      for (std::size_t k = 0; k < est.size(); ++k) est[k] = est[k] + d7;
    // Materialize the DoubleSeq.of views AFTER the tail-addition (view
    // semantics of the JVM: reads observe the mutated buffer).
    est_[idx(comp)] =
        vd(est.begin() + idx(nbc), est.begin() + idx(nbc) + idx(n));
    if (nfcasts_ > 0)
      fc_[idx(comp)] = vd(est.begin() + idx(n) + idx(nbc),
                          est.begin() + idx(n) + idx(nbc) + idx(nfcasts_));
    if (nbcasts_ > 0)
      bc_[idx(comp)] =
          vd(est.begin() + idx(nbc - nbcasts_),
             est.begin() + idx(nbc - nbcasts_) + idx(nbcasts_));
  }

  vd data_;
  WkEstimators& wk_;  // JVM mutates the lazy WK grid from the stdev getters
  bool bmean_ = false;
  int mcmp_ = 0;
  double ser_ = 1.0;
  int nfcasts_ = 0;
  int nbcasts_ = 0;

  vd ma_, ar_;                        // sqrt-scaled ma / ar (initModel)
  std::vector<std::optional<vd>> g_;  // null = component refused
  BurmanLu lu_;
  vd xfcasts_, xbcasts_, z_;
  std::vector<std::optional<vd>> est_, fc_, bc_;
};

}  // namespace seats_burman
}  // namespace detail
}  // namespace kalos
