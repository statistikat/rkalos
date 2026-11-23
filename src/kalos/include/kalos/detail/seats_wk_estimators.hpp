// Kalos - seasonal adjustment toolkit.
// 24-C.1 cluster 3: bit-faithful C++ mirror of
// jdplus.toolkit.base.core.ucarima.WienerKolmogorovEstimators (the WK
// estimator builder) plus its result object ucarima.WienerKolmogorovEstimator
// and the StationaryTransformation holder produced by
// finalStationaryEstimator. Source dumps: /tmp/opencode/ssf/{WKToolkit,
// WKIface}.txt (javap, jdplus 3.9.0, staged 2026-11-22); the builder
// composes purely on the cluster-1/2 mirrors (seats_linear_filters.hpp,
// seats_acf_solver.hpp) and the 24-B seam (seats_ucarima_model.hpp).
//
// Surface pinned here:
//  - calcestimator (WKToolkit.txt 16): the isSpecial split, the BackFilter
//    SimplifyingTool / divide / times chain, SymmetricFilter.convolutionOf +
//    times(iv-ratio) + decompose(BackFilter), the RationalBackFilter ctor
//    with delay 0 and RationalFilter ctor (back, fore, num, den) whose lazy
//    fields are pre-set from the decomposition, and the estimator Linear-
//    Process through RationalFilter.of(numBack, denBack, numFore, denFore);
//  - isSpecial (208): outer getArOrder minus the sum over ALL components
//    (raw aaload surface), true when the difference is non-zero;
//  - finalErrorModel (243): the 3-arg ArimaModel ctor with the outer MA in
//    the AR slot, delta = BackFilter.ONE and the sma product;
//  - finalEstimator (290) / finalStationaryEstimator (323) - the latter
//    always runs the simplify tool, re-reads cm's UNsimplified stationary
//    AR as the denominator slot, and wraps the process in the
//    StationaryTransformation holder;
//  - revisionModel (444) - null for h < 0, drop(h+1) on the estimator fore
//    filter, rebuilt through RationalFilter.of(RationalForeFilter);
//  - totalErrorAcf (488) - the 4-arg ArimaModel ctor (mirrored fore
//    denominator, ONE, normalized mirrored numerator, iv*d*d) folded by
//    ArimaModel.plus(ArimaModel) = pm(add=true, trim=true); the trim flag
//    upgrades the pm tail simplifyAr() -> simplify() =
//    simplifyUr(simplifyAr()), mirrored by am_plus_trim;
//  - revisionVariance (544) - including the dead-looking prelude loop that
//    walks i = -1 DOWN to pos, stores ar[i - pos] and then has acc and i
//    RESET at bytecode 150 (the stores survive); the dcmpg-ifgt tests make
//    a NaN running sum CONTINUE the loop;
//  - relativeRevisionVariance (689), totalErrorVariance (747),
//    variationPrecision (789), variationRevisionVariance (863).
//
// The finals[comp][b ? 1 : 0] cache is mirrored by value; the JVM shares
// the RationalFilter/LinearProcess instances, the mirror shares their lazy
// caches through shared_ptr (RationalFilter) and recomputes the pure
// LinearProcess acf - behavior-free differences only.
// StationaryTransformation.transform() and LinearProcess.doStationary()/
// getSpectrum() have no ported caller and stay unported by design
// (fail-closed outside the ported domain).

#pragma once

#include "kalos/detail/seats_acf_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace kalos {
namespace detail {
namespace seats_wk {

using seats_acf::Acf;
using seats_acf::am_acf;
using seats_acf::arima_exception;
using seats_acf::LinearProc;
using seats_acf::sf_decompose;
using seats_lfilters::BackFilter;
using seats_lfilters::FinFilter;
using seats_lfilters::FiniteFilter;
using seats_lfilters::ForeFilter;
using seats_lfilters::RatBackFilter;
using seats_lfilters::RatForeFilter;
using seats_lfilters::RationalFilter;
using seats_lfilters::RfSimp;
using seats_lfilters::SymmetricFilter;
using seats_model::AModel;
using seats_model::BackSimp;
using seats_model::back_simplify;
using seats_model::BFilter;
using seats_model::idx;
using seats_model::JvmNpe;
using seats_model::MaError;
using seats_model::p_of;
using seats_model::Poly;
using seats_model::p_times;
using seats_model::sf_times;
using seats_model::SymF;
using seats_model::um_get_complement;
using seats_model::UcarimaModel;
using seats_model::ur_solver_run;
using seats_model::UrSolverResult;

// ---------------------------------------------------------------------------
// ArimaModel glue on the seam's AModel in linearfilters types (the builder
// works with BackFilter/SymmetricFilter; the seam stores Polynomials).
// ---------------------------------------------------------------------------
inline BackFilter am_bf(const BFilter& b) { return BackFilter(b.p); }
inline BackFilter am_stationary_ar(const AModel& a) { return BackFilter(a.ar); }
inline BackFilter am_non_stationary_ar(const AModel& a) { return BackFilter(a.delta); }
inline BackFilter am_ma_bf(const AModel& a) { return am_bf(a.get_ma()); }
inline BackFilter am_ar_full(const AModel& a) {  // getAr() = ar.times(delta)
  return BackFilter(p_times(a.ar, a.delta));
}
// getArOrder() = getStationaryArOrder() + getNonStationaryArOrder(), each
// BackFilter.length() - 1 == the stored polynomial degree.
inline int am_ar_order(const AModel& a) { return a.ar.degree() + a.delta.degree(); }

// ArimaModel.simplifyUr() (pinned ArimaModel.txt 944): UnitRootsSolver(0)
// over getMa().asPolynomial(); no roots -> this. Polynomial$SimplifyingTool
// over (unitRoots.asPolynomial(), delta.asPolynomial()); false -> this.
// New model: (getStationaryAr(), new BackFilter(right),
// new BackFilter(left.times(solver.remainder())), getInnovationVariance()).
// Catches ArimaException ONLY (the getMa factorization failure -> this);
// the getIV NPE and the tool's PolynomialException propagate.
inline AModel am_simplify_ur(const AModel& a) {
  try {
    BFilter ma = a.get_ma();
    UrSolverResult solver = ur_solver_run(ma.p);
    if (solver.roots.roots_count() == 0) return a;  // factorize() false
    RfSimp tool;
    if (!tool.simplify(solver.roots.as_polynomial(), a.delta)) return a;
    AModel r;
    r.ar = a.ar;  // getStationaryAr() (the ar field)
    r.delta = tool.right;
    r.symmetric = false;
    r.has_ma = true;
    r.ma1 = p_times(tool.left, solver.remainder).c;
    r.var = a.get_iv();  // NPE not caught
    return r;
  } catch (const MaError&) {
    return a;
  }
}

// ArimaModel.plus(ArimaModel) (pinned ArimaModel.txt 789): pm(this, b,
// add = true, trim = true); the trim tail is model.simplify() =
// simplifyUr(simplifyAr()) - the seam's am_pm carries the simplifyAr tail,
// so trim=true chains simplifyUr on its result (same two calls, same
// receivers - no re-run of simplifyAr).
inline AModel am_plus_trim(const AModel& a, const AModel& b, bool trim) {
  AModel r = seats_model::am_pm(a, b, true, false);
  if (trim) r = am_simplify_ur(r);
  return r;
}

// ---------------------------------------------------------------------------
// jdplus ...core.ucarima.WienerKolmogorovEstimator (WKIface.txt): the
// per-(component, side) result pair (wienerKolmogorovFilter,
// estimatorModel). equals/hashCode/toString are not numerics - skipped.
// ---------------------------------------------------------------------------
struct WkEstimator {
  RationalFilter wkf;
  LinearProc est;
};

// jdplus ...core.arima.StationaryTransformation<ILinearProcess> holder:
// (process, delta) as built by finalStationaryEstimator. The transform()
// machinery has no ported caller.
struct StationaryTransformation {
  LinearProc proc;
  BackFilter delta;
};

// ---------------------------------------------------------------------------
// jdplus ...core.ucarima.WienerKolmogorovEstimators (WKToolkit.txt).
// The JVM holds `ucm' (a UcarimaModel) and the lazy finals grid; this
// mirror additionally carries the outer model as a seam AModel - on the
// corpus path the JVM's ucm.getModel() is the SarimaModel and its IArimaModel
// getters are the 21-era detail::Sarima/Arima mapping, converted to this
// AModel by the consuming cluster (the probe cases use plain ArimaModel
// outers built through UcarimaModel.Builder.model).
// ---------------------------------------------------------------------------
struct WkEstimators {
  const UcarimaModel* um = nullptr;  // null ucm -> every entry null (pinned)
  AModel outer;
  std::vector<std::array<std::optional<WkEstimator>, 2>> finals;  // lazy [n][2]

  WkEstimators() = default;
  WkEstimators(const UcarimaModel* u, AModel outer_) : um(u), outer(std::move(outer_)) {}

  static double jvm_nan() { return std::numeric_limits<double>::quiet_NaN(); }

  // UcarimaModel.getComponent(i): raw aaload - out-of-range mirrors AIOOBE.
  const AModel& component(int i) const {
    int n = static_cast<int>(um->comps.size());
    if (i < 0 || i >= n)
      throw std::out_of_range("seats_wk_estimators: getComponent (AIOOBE parity)");
    return um->comps[idx(i)];
  }

  // isSpecial(ucm) (pinned 208): outer.getArOrder() minus every component's
  // getArOrder(); true iff the remainder is non-zero.
  bool is_special() const {
    int ord = am_ar_order(outer);
    for (const AModel& c : um->comps) ord -= am_ar_order(c);
    return ord != 0;
  }

  // private calcestimator(int, boolean) (pinned 16).
  void calcestimator(int comp, bool b) {
    if (um == nullptr) return;
    int k = b ? 1 : 0;
    if (finals.empty()) {
      finals.resize(um->comps.size());
    } else {
      if (comp < 0 || comp >= static_cast<int>(finals.size()))
        throw std::out_of_range("seats_wk_estimators: finals[comp] (AIOOBE parity)");
      if (finals[idx(comp)][idx(k)].has_value()) return;
    }
    const AModel& m = component(comp);
    AModel cm = b ? m : um_get_complement(*um, comp);  // slot 5
    BackFilter ar = am_stationary_ar(outer);   // slot 6
    BackFilter car = am_stationary_ar(cm);     // slot 7
    BackFilter delta = am_non_stationary_ar(outer);  // slot 8
    BackFilter cdelta = am_non_stationary_ar(cm);    // slot 9
    BackFilter ma{Poly(vd{0.0})};  // slot 10 (placeholder; both branches assign)
    if (is_special()) {
      BackSimp t = back_simplify(BFilter{ar.p}, BFilter{car.p});
      if (t.ret) {
        ar = BackFilter{t.left.p};
        car = BackFilter{t.right.p};
      }
      delta = delta.divide(cdelta);
      ar = ar.times(delta);
      ma = am_ma_bf(outer).times(car);  // special branch folds car in
    } else {
      delta = delta.divide(cdelta);
      ar = ar.divide(car);
      ar = ar.times(delta);
      ma = am_ma_bf(outer);
    }
    BackFilter num = am_ma_bf(cm).times(ar);  // slot 11
    SymmetricFilter sma = SymmetricFilter::conv_of(num);  // slot 12
    double cmvar = cm.get_iv();               // slot 13 (dstore 13)
    double mvar = outer.get_iv();             // slot 15 (dstore 15)
    sma = sma.times(cmvar / mvar);            // times(div - raw IEEE divide)
    BackFilter den = sf_decompose(sma, ma);   // SymmetricFilter.decompose(10)
    RatBackFilter rbf(den, ma, 0);            // slot 18
    RationalFilter est(rbf, rbf.mirror(), RationalFilter::sp(FiniteFilter(sma)),
                       RationalFilter::sp(FiniteFilter(SymmetricFilter::conv_of(ma))));  // slot 19
    RationalFilter rf2 = RationalFilter::of_bbff(  // slot 20
        am_ma_bf(cm).times(cmvar / mvar), am_ar_full(cm), num.mirror(), ma.mirror());
    finals[idx(comp)][idx(k)] = WkEstimator{est, LinearProc(rf2, mvar)};
  }

  // public finalEstimator(int, boolean) (pinned 290).
  std::optional<WkEstimator> final_estimator(int comp, bool b) {
    if (um == nullptr) return std::nullopt;
    if (component(comp).is_null()) return std::nullopt;
    calcestimator(comp, b);
    return finals[idx(comp)][idx(b ? 1 : 0)];
  }

  // public finalErrorModel(int) (pinned 243): the 3-arg ArimaModel ctor
  // (ar = outer MA!, delta = ONE, sma = symmetricMa * symmetricMa; ma null,
  // var 0 on that ctor).
  std::optional<AModel> final_error_model(int comp) {
    if (um == nullptr) return std::nullopt;
    const AModel& m = component(comp);
    if (m.is_null()) return std::nullopt;
    AModel cm = um_get_complement(*um, comp);
    SymF sma1 = m.symmetric_ma();
    SymF sma2 = cm.symmetric_ma();
    SymF prod = sf_times(sma1, sma2);
    AModel r;
    r.ar = outer.get_ma().p;  // BackFilter arg 1 -> the ar field
    r.delta = p_of(vd{1.0});  // BackFilter.ONE
    r.symmetric = true;
    r.has_ma = false;
    r.sma = std::move(prod);
    r.var = 0.0;
    return r;
  }

  // public finalStationaryEstimator(int, boolean) (pinned 323): always
  // simplifies, denBack is cm's ORIGINAL stationary AR (re-read at 206),
  // delta3 = cm's non-stationary AR.
  std::optional<StationaryTransformation> final_stationary_estimator(int comp, bool b) {
    if (um == nullptr) return std::nullopt;
    AModel cm = b ? component(comp) : um_get_complement(*um, comp);  // slot 5
    if (cm.is_null()) return std::nullopt;
    BackFilter ar = am_stationary_ar(outer);   // slot 6
    BackFilter car = am_stationary_ar(cm);     // slot 7
    BackFilter delta = am_non_stationary_ar(outer);  // slot 8
    BackFilter cdelta = am_non_stationary_ar(cm);    // slot 9
    BackSimp t = back_simplify(BFilter{ar.p}, BFilter{car.p});
    if (t.ret) {
      ar = BackFilter{t.left.p};
      car = BackFilter{t.right.p};
    }
    delta = delta.divide(cdelta);
    ar = ar.times(delta);
    BackFilter ma = am_ma_bf(outer).times(car);   // slot 11
    BackFilter num = am_ma_bf(cm).times(ar);      // slot 12
    double cmvar = cm.get_iv();                   // slot 13
    double mvar = outer.get_iv();                 // slot 15
    BackFilter ar17 = am_stationary_ar(cm);       // slot 17 - un-simplified
    BackFilter delta3 = am_non_stationary_ar(cm);  // slot 3
    RationalFilter rf = RationalFilter::of_bbff(
        am_ma_bf(cm).times(cmvar / mvar), ar17, num.mirror(), ma.mirror());
    return StationaryTransformation{LinearProc(rf, mvar), delta3};
  }

  // public revisionModel(int, int) (pinned 444): null for h < 0; otherwise
  // drop(h + 1) of the FINAL estimator's fore filter, rebuilt as a fresh
  // RationalForeFilter(delay 0) through RationalFilter.of(RationalForeFilter).
  // A null component's estimator NPEs on .getEstimatorModel() in the JVM.
  std::optional<LinearProc> revision_model(int comp, int h) {
    if (h < 0) return std::nullopt;
    std::optional<WkEstimator> est = final_estimator(comp, true);
    if (!est.has_value()) throw JvmNpe{};
    LinearProc lp = est->est;
    RatForeFilter rff = lp.get_filter().fore.drop(h + 1);
    RationalFilter g = RationalFilter::of_fore(
        RatForeFilter(rff.get_numerator(), rff.get_denominator(), 0));
    return LinearProc(g, lp.get_innovation_variance());
  }

  // public totalErrorAcf(int, int) (pinned 488): the revision model's fore
  // denominator/numerator mirrored into the 4-arg ArimaModel (var scaled by
  // d*d with d the mirrored numerator's get(0)), folded onto the error
  // model with the TRIMMED plus, then AbstractArimaModel.
  // getAutoCovarianceFunction (am_acf).
  Acf total_error_acf(int comp, int h) {
    std::optional<LinearProc> lpo = revision_model(comp, h);
    std::optional<AModel> erro = final_error_model(comp);
    if (!lpo.has_value()) throw JvmNpe{};    // null.getFilter()
    if (!erro.has_value()) throw JvmNpe{};   // null.plus(...)
    LinearProc& lp = *lpo;
    RatForeFilter& rff = lp.get_filter().fore;
    BackFilter num = rff.get_numerator().mirror();  // ForeFilter.mirror()
    double d0 = num.get(0);
    AModel m2;  // 4-arg ctor: (ar, delta, ma, var), sma null, one-sided ma
    m2.ar = rff.get_denominator().mirror().p;
    m2.delta = p_of(vd{1.0});  // BackFilter.ONE
    m2.symmetric = false;
    m2.has_ma = true;
    m2.ma1 = num.normalize().p.c;
    double iv = lp.get_innovation_variance();
    m2.var = (iv * d0) * d0;
    AModel sum = am_plus_trim(*erro, m2, true);
    return am_acf(sum);
  }

  // public revisionAcf(int, int) (pinned 533).
  Acf revision_acf(int comp, int h) {
    std::optional<LinearProc> lp = revision_model(comp, h);
    if (!lp.has_value()) throw JvmNpe{};  // null.getAutoCovarianceFunction()
    return lp->get_auto_covariance_function();
  }

  // public revisionVariance(int, boolean, int, int) (pinned 544). The
  // prelude loop walks i from -1 DOWN while i >= pos storing ar[i - pos];
  // at bytecode 150 acc REVERTS to d6 and i JUMPS to 1 - only the prelude
  // STORES survive (measured on pos < -1 probes). The 0-comparisons are
  // dcmpg: a NaN running sum is > 0 for the loop and CONTINUES.
  std::optional<vd> revision_variance(int comp, bool b, int pos, int n) {
    if (um == nullptr) throw JvmNpe{};  // ucm.getComponent on the null field
    if (component(comp).is_null()) return std::nullopt;
    vd ar(static_cast<std::size_t>(n), 0.0);
    LinearProc lp0 = *revision_model(comp, 0);
    double d6 = lp0.get_auto_covariance_function().get(0);
    std::optional<WkEstimator> est = final_estimator(comp, b);
    if (!est.has_value()) throw JvmNpe{};
    RationalFilter f = est->est.get_filter();
    double iv = est->est.get_innovation_variance();
    if (pos <= 0 && -pos < n) ar[idx(-pos)] = d6;
    double acc = d6;
    int i = -1;
    for (;;) {  // prelude (94): loop while i >= pos
      if (i < pos) break;
      double wi = f.weight(i + 1);
      acc = acc + (wi * wi * iv);
      if (i - pos < n) ar[idx(i - pos)] = acc;
      --i;
    }
    acc = d6;  // 150: reset (the prelude stores survive)
    i = 1;
    for (; i < pos; ++i) {
      double wi = f.weight(i);
      acc = acc - (wi * wi * iv);
      if (!(std::isnan(acc) || acc > 0.0)) return ar;  // dcmpg ifgt: NaN stays
    }
    for (i = std::max(pos, 1); i < n + pos; ++i) {
      double wi = f.weight(i);
      acc = acc - (wi * wi * iv);
      if (!(std::isnan(acc) || acc > 0.0)) break;  // 243-250: return ar
      ar[idx(i - pos)] = acc;
    }
    return ar;
  }

  // public relativeRevisionVariance(int, boolean, int, int) (pinned 689):
  // ar[0] stays 0.0; ar[i] = cumsum_{j=1..i} w(j + pos)^2 * iv.
  std::optional<vd> relative_revision_variance(int comp, bool b, int pos, int n) {
    if (um == nullptr) throw JvmNpe{};
    if (component(comp).is_null()) return std::nullopt;
    vd ar(static_cast<std::size_t>(n), 0.0);
    std::optional<WkEstimator> est = final_estimator(comp, b);
    if (!est.has_value()) throw JvmNpe{};
    RationalFilter f = est->est.get_filter();
    double iv = est->est.get_innovation_variance();
    double acc = 0.0;
    for (int i = 1; i < n; ++i) {
      double wi = f.weight(i + pos);
      acc = acc + (wi * wi * iv);
      ar[idx(i)] = acc;
    }
    return ar;
  }

  // public totalErrorVariance(int, boolean, int, int) (pinned 747):
  // revisionVariance(...) + finalErrorModel(...).getAutoCovarianceFunction()
  // .get(0) added to EVERY entry (the component-null guard is doubled:
  // this method's own and revisionVariance's).
  std::optional<vd> total_error_variance(int comp, bool b, int pos, int n) {
    if (um == nullptr) throw JvmNpe{};
    if (component(comp).is_null()) return std::nullopt;
    std::optional<vd> ar = revision_variance(comp, b, pos, n);
    if (!ar.has_value()) return ar;
    std::optional<AModel> erro = final_error_model(comp);
    if (!erro.has_value()) throw JvmNpe{};
    Acf acf = am_acf(*erro);
    double d = acf.get(0);
    for (std::size_t k = 0; k < ar->size(); ++k) (*ar)[k] = (*ar)[k] + d;
    return ar;
  }

  // public variationPrecision(int, int, int, boolean) (pinned 789):
  // 2.0 * (acf(0) - acf(k)) minus sum_{i=0..k-1} w(i + h + 1)^2 * iv with
  // the FINAL estimator's filter and iv; acf = totalErrorAcf or revisionAcf.
  double variation_precision(int comp, int h, int k, bool total) {
    if (um == nullptr) throw JvmNpe{};
    if (component(comp).is_null()) return jvm_nan();
    std::optional<WkEstimator> est = final_estimator(comp, true);
    if (!est.has_value()) throw JvmNpe{};
    RationalFilter f = est->est.get_filter();
    double iv = est->est.get_innovation_variance();
    Acf acf = total ? total_error_acf(comp, h) : revision_acf(comp, h);
    double v = 2.0 * (acf.get(0) - acf.get(k));
    for (int i = 0; i < k; ++i) {
      double wi = f.weight(i + h + 1);
      v = v - (wi * wi * iv);
    }
    return v;
  }

  // public variationRevisionVariance(int, int, int, int) (pinned 863):
  // iv * sum_{i=1..n} (w(i+h)^2 + w(i+h+k)^2 - (2.0*w(i+h))*w(i+h+k)).
  double variation_revision_variance(int comp, int h, int k, int n) {
    if (um == nullptr) throw JvmNpe{};
    if (component(comp).is_null()) return jvm_nan();
    std::optional<WkEstimator> est = final_estimator(comp, true);
    if (!est.has_value()) throw JvmNpe{};
    RationalFilter f = est->est.get_filter();
    double iv = est->est.get_innovation_variance();
    double acc = 0.0;
    for (int i = 1; i <= n; ++i) {
      double w1 = f.weight(i + h);
      double w2 = f.weight(i + h + k);
      acc = acc + ((w1 * w1 + w2 * w2) - ((2.0 * w1) * w2));
    }
    return acc * iv;
  }
};

}  // namespace seats_wk
}  // namespace detail
}  // namespace kalos
