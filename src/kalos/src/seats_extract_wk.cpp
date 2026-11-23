// Item 24-C.1: the SEATS Wiener-Kolmogorov extraction machine as a
// bit-exact native mirror - the SeatsKernel path from the validated model
// over the decomposition seam into the WienerKolmogorovEstimator initial
// surface and the DefaultBiasCorrector final surface.
//
//   validate (24-A)  ->  selector chain + setVarianceMax (24-B)
//                    ->  UcarimaModel.compact(2, 2)
//                    ->  BurmanEstimates (cluster 5 mirror)
//                    ->  WienerKolmogorovEstimator.decompose (init map)
//                    ->  DefaultBiasCorrector.correctBias (final map)
//                    ->  status / dom / init / final fact rows
//
// Every routine transcribes the pinned bytecode recorded in /tmp/opencode/
// ssf (fresh dumps of WienerKolmogorovEstimator, DefaultBiasCorrector,
// SeatsKernel, SeatsToolkit, SeatsModel, SeatsResults, UcarimaModel,
// TsData, TsDomain/TsPeriod, DoublesMath, DoubleSeq, dstats.LogNormal and
// the PED chain); the mirrored classes are named in the comments. Operation
// order, comparison flavour (dcmpg/dcmpl) and rounding sequence are
// preserved; the build uses -ffp-contract=off so no FMA is introduced.
// Math.log/Math.exp map to the fdlibm port (the HotSpot implementations).

#include <exception>
#include "kalos/seats_extract_wk.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/seats_burman_estimates.hpp"
#include "kalos/detail/seats_exact_forecasts.hpp"
#include "kalos/detail/seats_ucarima_model.hpp"
#include "kalos/detail/seats_wk_estimators.hpp"
#include "kalos/detail/ucarima_engine.hpp"

namespace kalos {
namespace {

namespace sb = detail::seats_burman;
namespace slf = detail::seats_lfilters;
namespace sm = detail::seats_model;
namespace swk = detail::seats_wk;
using detail::seats_xf::demetra_abort;
using detail::seats_xf::eco_abort;
using sm::AModel;
using sm::UcarimaModel;
using vd = std::vector<double>;

// ---------------------------------------------------------------------------
// TsData, reduced to the monthly-indexed window the extraction path uses
// (start position + dense values). Every op reproduces the pinned domain
// intersection; on the corpus path operands always share their domain.
// ---------------------------------------------------------------------------
struct Ts {
  std::int64_t m0 = 0;  // absolute month index (year * 12 + period - 1)
  vd v;
  std::int64_t end() const { return m0 + static_cast<std::int64_t>(v.size()); }
};
using Opt = std::optional<Ts>;
using OptV = std::optional<vd>;

[[noreturn]] void jvm_npe() { throw sm::JvmNpe{}; }

struct wrong_dimensions_abort {};  // IllegalArgumentException("wrong ...")

template <class F>
Ts ts_fn(const Ts& a, F op) {
  Ts r{a.m0, vd(a.v.size(), 0.0)};
  for (std::size_t k = 0; k < a.v.size(); ++k) r.v[k] = op(a.v[k]);
  return r;
}

template <class F>
Ts ts_fn2(const Ts& a, const Ts& b, F op) {
  // TsData.fn/fastFn: intersection domain (lambda$fastFn$1 offset mapping);
  // a domain miss (hi < lo) is the JVM null return.
  const std::int64_t lo = a.m0 > b.m0 ? a.m0 : b.m0;
  const std::int64_t hi = a.end() < b.end() ? a.end() : b.end();
  if (hi < lo) throw sm::JvmNpe{};  // caller-side NPE parity (unreachable:
                                    // the extraction path aligns)
  Ts r{lo, vd(hi - lo, 0.0)};
  for (std::int64_t k = 0; k < hi - lo; ++k) r.v[k] = op(a.v[k + (lo - a.m0)], b.v[k + (lo - b.m0)]);
  return r;
}

template <class F>
Opt opt_fn(const Opt& a, F op) {
  if (!a) return std::nullopt;
  return ts_fn(*a, op);
}

double plus(double x, double y) { return x + y; }

Opt ts_add(const Opt& a, const Opt& b) {  // TsData.add static (null quirks)
  if (!a) return b;
  if (!b) return a;
  return ts_fn2(*a, *b, plus);
}

// TsData.add(a, {b, c}) (varargs, bytecode 861-958): k = first non-null arg;
// none -> return a; a null && k last -> return args[k] (alias); acc = a or
// args[k], then left fold acc = acc.fastFn(args[j], +) over j > k (k for a
// non-null) - a null argument at j >= k dereferences null in fastFn.
Opt ts_add_v(const Opt& a, const Opt& b, const Opt& c) {
  const Opt* args[2] = {&b, &c};
  int k = args[0] ? 0 : (args[1] ? 1 : -1);
  if (k < 0) return a;
  Opt acc = a;
  int j = k;
  if (!acc) {
    acc = *args[k];
    j = k + 1;
    if (j == 2) return acc;  // alias return
  }
  for (; j < 2; ++j) {
    if (!*args[j]) jvm_npe();
    acc = ts_fn2(*acc, **args[j], plus);
  }
  return acc;
}

Opt ts_sub(const Opt& a, const Opt& b) {  // TsData.subtract static
  if (!a) return b ? opt_fn(b, [](double x) { return -x; }) : Opt{};
  if (!b) return a;
  return ts_fn2(*a, *b, [](double x, double y) { return x - y; });
}

Opt ts_div(const Opt& a, const Opt& b) {  // TsData.divide static
  if (!a) return opt_fn(b, [](double x) { return 1.0 / x; });  // inv()
  if (!b) return a;
  return ts_fn2(*a, *b, [](double x, double y) { return x / y; });
}

// TsData.multiply(double) (shortcut + zero- and generic lambdas in order).
Opt ts_mult(const Opt& a, double c) {
  if (!a) return std::nullopt;
  if (c == 1.0) return a;
  if (c == 0.0) return ts_fn(*a, [](double x) { return 0.0 * x; });
  return ts_fn(*a, [c](double x) { return x * c; });
}

Ts ts_exp(const Ts& a) { return ts_fn(a, detail::fdlibm_exp); }
Opt ts_exp(const Opt& a) { return opt_fn(a, detail::fdlibm_exp); }

// DoubleSeq.reduce(0.0, +).average-then-divide (DoublesMath.average).
double seq_average(const vd& v) {
  double s = 0.0;
  for (std::size_t k = 0; k < v.size(); ++k) s += v[k];
  return s / static_cast<double>(v.size());
}

// LogNormal.stdev2(level, stdev) applied elementwise over the fastFn
// intersection: correctStdevForLog(stdev, level) (lambda$correctStdev-
// ForLog$0 passes (d2=level, d1=stdev)).
Ts correct_stdev_for_log(const Ts& stdev, const Ts& level) {
  const std::int64_t lo = stdev.m0 > level.m0 ? stdev.m0 : level.m0;
  const std::int64_t hi = stdev.end() < level.end() ? stdev.end() : level.end();
  if (hi < lo) throw sm::JvmNpe{};
  Ts r{lo, vd(hi - lo, 0.0)};
  for (std::int64_t k = 0; k < hi - lo; ++k) {
    const double m = level.v[k + (lo - level.m0)];
    const double s = stdev.v[k + (lo - stdev.m0)];
    const double s2 = s * s;
    r.v[k] = (m * detail::fdlibm_exp(0.5 * s2)) *
             std::sqrt(detail::fdlibm_exp(s2) - 1.0);
  }
  return r;
}

// DoublesMath.add(acc, [v]) with the pinned null handling (null acc -> v
// alias; length mismatch -> IllegalArgumentException("wrong dimensions")).
OptV dm_add(const OptV& acc, const vd& v) {
  if (!acc) return v;
  if (v.size() != acc->size()) throw wrong_dimensions_abort{};
  OptV r{vd(*acc)};
  for (std::size_t k = 0; k < v.size(); ++k) (*r)[k] += v[k];
  return r;
}

// DoublesMath.subtract(a, b): b null -> a; a null -> -b; length mismatch ->
// IllegalArgumentException("wrong dimensions").
vd dm_sub(const vd& a, const OptV& b) {
  if (!b) return a;
  if (a.size() != b->size()) throw wrong_dimensions_abort{};
  vd r(a.size());
  for (std::size_t k = 0; k < a.size(); ++k) r[k] = a[k] - (*b)[k];
  return r;
}

// lambda$decompose$2..5: dcmpg ifgt - x > 0 or NaN take sqrt, else 0.
double sqrt_plus(double x) {
  if (x > 0.0 || std::isnan(x)) return std::sqrt(x);
  return 0.0;
}
double square(double x) { return x * x; }

// ---------------------------------------------------------------------------
// SeriesDecomposition: the six-by-six EnumMap of raw TsData items. Builder.
// add is EnumMap.put (a repeated key OVERWRITES); getSeries is the raw
// Map.get - no derivation, no domain alignment (pinned SeriesDecomposition
// fresh dump).
// ---------------------------------------------------------------------------
enum { kSeries, kTrend, kSeasonal, kAdjusted, kIrregular, kCalendar, kNT };
enum {
  kValue, kStdev, kForecast, kStdevForecast, kBackcast, kStdevBackcast, kNI
};
const char* const kTypeNames[kNT] = {"Series",        "Trend",
                                     "Seasonal",      "SeasonallyAdjusted",
                                     "Irregular",     "CalendarEffect"};
const char* const kInfNames[kNI] = {"Value",
                                    "Stdev",
                                    "Forecast",
                                    "StdevForecast",
                                    "Backcast",
                                    "StdevBackcast"};

struct Dec {
  Opt items[kNT][kNI];
  void put(int t, int i, const Opt& x) { items[t][i] = x; }
  const Opt& get(int t, int i) const { return items[t][i]; }
};

struct MonthStart {
  std::int64_t year;
  int position;  // Java annualPosition(), January = 0
};

MonthStart month_start(std::int64_t index) {
  // TsPeriod uses floor division for a backcast crossing year zero.
  const int remainder = static_cast<int>(index % 12);
  return {index / 12 - (remainder < 0), remainder < 0 ? remainder + 12 : remainder};
}

std::vector<seats_component_series> collect_decomposition(
    std::vector<std::string>& out, const char* key, Dec& d) {
  constexpr seats_component types[] = {
      seats_component::series, seats_component::trend, seats_component::seasonal,
      seats_component::seasonally_adjusted, seats_component::irregular,
      seats_component::calendar_effect};
  constexpr seats_information kinds[] = {
      seats_information::value, seats_information::stdev, seats_information::forecast,
      seats_information::stdev_forecast, seats_information::backcast,
      seats_information::stdev_backcast};
  std::vector<seats_component_series> rows;
  // probe dumpDecomposition: fixed enum enumeration, absent rows skipped,
  // year + 0-based annualPosition, %.17g doubles (seq).
  for (int t = 0; t < kNT; ++t)
    for (int i = 0; i < kNI; ++i) {
      Opt& x = d.items[t][i];
      if (!x) continue;
      const auto start = month_start(x->m0);
      std::string line = std::string(key) + " " + kTypeNames[t] + " " +
                         kInfNames[i] + " n " + std::to_string(x->v.size()) +
                         " start " + std::to_string(start.year) + " " +
                         std::to_string(start.position) + " vals ";
      // probe seq(): elements space-separated, NO trailing space; the
      // " vals " prefix carries one (empty rows end with a space).
      for (std::size_t k = 0; k < x->v.size(); ++k) {
        if (k != 0) line += " ";
        line += sm::fmt_g(x->v[k]);
      }
      out.push_back(std::move(line));
      rows.push_back({types[t], kinds[i], start.year, start.position + 1, std::move(x->v)});
    }
  return rows;
}

void dump_domain(std::vector<std::string>& out, const char* key,
                 const Opt& x) {
  const auto start = x ? month_start(x->m0) : MonthStart{1970, 0};
  out.push_back(std::string(key) + " n " + std::to_string(x ? x->v.size() : 0) +
                " start " + std::to_string(start.year) + " " + std::to_string(start.position));
}

// ---------------------------------------------------------------------------
// WienerKolmogorovEstimator.decompose(SeatsModel) - the initial surface
// (fresh dump WKE_full.txt). data/start/end/bstart are the transformed
// series frame; be carries the Burman arrays.
// ---------------------------------------------------------------------------
Dec wke_decompose(sb::BurmanEstimates& be, const UcarimaModel& um3,
                  const Ts& data, int nbc) {
  const int ncomp = static_cast<int>(um3.comps.size());
  OptV est[kNT], fc[3], se[3], sf[3], bc[3], sbc[3];
  for (int c = 0; c < ncomp; ++c) {
    if (c != 0 && um3.comps[sm::idx(c)].is_null()) continue;  // WKE 175-189
    est[c] = be.estimates(c, true);
    se[c] = be.stdev_estimates(c);  // never null on a fetched component
    fc[c] = be.forecasts(c, true);
    sf[c] = be.stdev_forecasts(c, true);
    bc[c] = be.backcasts(c, true);
    sbc[c] = be.stdev_backcasts(c, true);
  }
  const vd seriesFc = be.series_forecasts();
  const vd seriesBc = be.series_backcasts();

  // Variance accumulators (335-450): v = sf[c].fn(x*x); accSeries over all
  // c, accSA iff c != 1; same for the sb[] side.
  OptV accSf, accSaFc, accSb, accSaBc;
  for (int c = 0; c < ncomp; ++c) {
    if (!sf[c]) continue;
    OptV v{vd(sf[c]->size())};
    for (std::size_t k = 0; k < sf[c]->size(); ++k) (*v)[k] = square((*sf[c])[k]);
    accSf = dm_add(accSf, *v);
    if (c != 1) accSaFc = dm_add(accSaFc, *v);
  }
  for (int c = 0; c < ncomp; ++c) {
    if (!sbc[c]) continue;
    OptV v{vd(sbc[c]->size())};
    for (std::size_t k = 0; k < sbc[c]->size(); ++k)
      (*v)[k] = square((*sbc[c])[k]);
    accSb = dm_add(accSb, *v);
    if (c != 1) accSaBc = dm_add(accSaBc, *v);
  }

  const std::int64_t start = data.m0;
  const std::int64_t end = data.end();
  const std::int64_t bstart = start - nbc;
  auto at = [](std::int64_t m, const OptV& x) -> Opt {
    if (!x) return std::nullopt;
    return Ts{m, *x};
  };
  auto at_sqrt = [](std::int64_t m, const OptV& x) -> Opt {
    if (!x) return std::nullopt;
    Ts r{m, vd(x->size())};
    for (std::size_t k = 0; k < x->size(); ++k) r.v[k] = sqrt_plus((*x)[k]);
    return r;
  };
  auto val = [](std::int64_t m, const OptV& x) -> Ts {
    if (!x) jvm_npe();  // TsData.of(period, null)
    return Ts{m, *x};
  };

  Dec d;
  d.put(kSeries, kValue, data);
  d.put(kSeries, kForecast, at(end, OptV{seriesFc}));
  d.put(kSeries, kStdevForecast, at_sqrt(end, accSf));
  d.put(kSeries, kBackcast, at(bstart, OptV{seriesBc}));
  d.put(kSeries, kStdevBackcast, at_sqrt(bstart, accSb));
  d.put(kTrend, kValue, at(start, est[0]));
  d.put(kSeasonal, kValue, at(start, est[1]));
  d.put(kTrend, kForecast, at(end, fc[0]));
  d.put(kSeasonal, kForecast, at(end, fc[1]));
  d.put(kAdjusted, kForecast, Ts{end, dm_sub(seriesFc, fc[1])});
  d.put(kTrend, kStdev, at(start, se[0]));
  d.put(kTrend, kStdevForecast, at(end, sf[0]));
  d.put(kAdjusted, kValue, Ts{start, dm_sub(data.v, est[1])});
  if (se[1]) {  // ONE guard feeds Seasonal/Stdev AND SA/Stdev (789-835)
    d.put(kSeasonal, kStdev, at(start, se[1]));
    d.put(kAdjusted, kStdev, at(start, se[1]));
  }
  d.put(kSeasonal, kStdevForecast, at(end, sf[1]));
  d.put(kIrregular, kValue, val(start, est[2]));  // UNCONDITIONAL add
  d.put(kIrregular, kForecast, at(end, fc[2]));
  d.put(kIrregular, kStdev, at(start, se[2]));
  d.put(kIrregular, kStdevForecast, at(end, sf[2]));
  d.put(kTrend, kBackcast, at(bstart, bc[0]));
  d.put(kSeasonal, kBackcast, at(bstart, bc[1]));
  d.put(kAdjusted, kStdevForecast, at_sqrt(end, accSaFc));
  d.put(kAdjusted, kBackcast, Ts{bstart, dm_sub(seriesBc, bc[1])});
  d.put(kTrend, kStdevBackcast, at(bstart, sbc[0]));
  d.put(kSeasonal, kStdevBackcast, at(bstart, sbc[1]));
  d.put(kIrregular, kBackcast, at(bstart, bc[2]));
  d.put(kIrregular, kStdevBackcast, at(bstart, sbc[2]));
  d.put(kAdjusted, kStdevBackcast, at_sqrt(bstart, accSaBc));
  return d;
}

// ---------------------------------------------------------------------------
// DefaultBiasCorrector (fresh dump DefaultBiasCorrector_fresh.txt). The
// final map is built by the same put-semantics builder; insertion order is
// invisible through the probe enumeration but kept as in the JVM.
// ---------------------------------------------------------------------------
struct BiasContext {
  Ts original;   // SeatsModel.getOriginalSeries
  Dec init;      // getInitialComponents
  bool bias;     // DefaultBiasCorrector(bias)
  int period;    // SeatsModel.getPeriod (original model period)
};

void fill_fc_bc(Dec& out, const BiasContext& cx, bool forecasts) {
  // fill{Fore,Back}casts(init, builder) - the 2-arg ADDITIVE passes (bytecode
  // 158-279 / 280-436). Slot order: Series, Seasonal, Trend, Irregular.
  const int I = forecasts ? kForecast : kBackcast;
  const int IS = forecasts ? kStdevForecast : kStdevBackcast;
  const Opt fcS = cx.init.get(kSeries, I);
  const Opt fcSe = cx.init.get(kSeasonal, I);
  const Opt fcT = cx.init.get(kTrend, I);
  const Opt fcI = cx.init.get(kIrregular, I);
  Opt series = fcS, saFc;
  if (!fcS) {
    series = ts_add_v(fcSe, fcT, fcI);
    saFc = ts_add(fcT, fcI);
  } else {
    saFc = ts_sub(fcS, fcSe);
  }
  if (!series || series->v.empty()) return;  // the whole pass aborts (91-102)
  out.put(kSeries, I, series);
  if (cx.init.get(kSeries, IS)) out.put(kSeries, IS, cx.init.get(kSeries, IS));
  if (fcSe) {
    out.put(kSeasonal, I, fcSe);
    if (cx.init.get(kSeasonal, IS))
      out.put(kSeasonal, IS, cx.init.get(kSeasonal, IS));
  }
  if (fcT) {
    out.put(kTrend, I, fcT);
    if (cx.init.get(kTrend, IS)) out.put(kTrend, IS, cx.init.get(kTrend, IS));
  }
  if (saFc) {
    out.put(kAdjusted, I, saFc);
    if (cx.init.get(kAdjusted, IS))
      out.put(kAdjusted, IS, cx.init.get(kAdjusted, IS));
  }
  if (fcI) {
    out.put(kIrregular, I, fcI);
    if (cx.init.get(kIrregular, IS))
      out.put(kIrregular, IS, cx.init.get(kIrregular, IS));
  }
}

Dec fill_levels(const BiasContext& cx) {
  // fillLevels(SeatsModel) (bytecode 36-157): values first (SA/I derived by
  // TsData.subtract with the null quirks), then the two additive passes into
  // the SAME builder.
  Dec out;
  const Ts& orig = cx.original;
  out.put(kSeries, kValue, orig);
  const Opt seasV = cx.init.get(kSeasonal, kValue);
  if (seasV) {
    out.put(kSeasonal, kValue, seasV);
    if (cx.init.get(kSeasonal, kStdev))
      out.put(kSeasonal, kStdev, cx.init.get(kSeasonal, kStdev));
  }
  const Opt trendV = cx.init.get(kTrend, kValue);
  if (trendV) {
    out.put(kTrend, kValue, trendV);
    if (cx.init.get(kTrend, kStdev))
      out.put(kTrend, kStdev, cx.init.get(kTrend, kStdev));
  }
  const Opt sa = ts_sub(Opt{orig}, seasV);
  out.put(kAdjusted, kValue, sa);
  if (cx.init.get(kAdjusted, kStdev))
    out.put(kAdjusted, kStdev, cx.init.get(kAdjusted, kStdev));
  const Opt irr = ts_sub(sa, trendV);
  out.put(kIrregular, kValue, irr);
  if (cx.init.get(kIrregular, kStdev))
    out.put(kIrregular, kStdev, cx.init.get(kIrregular, kStdev));
  fill_fc_bc(out, cx, true);
  fill_fc_bc(out, cx, false);
  return out;
}

void fill_fc3(Dec& out, const BiasContext& cx, double saMean, double irMean) {
  // fillForecasts(init, builder, saMean, irMean) - the LOG path (fresh dump
  // bytecode 687-888); fillBackcasts is the exact same shape with the
  // Backcast information (889-1090, lambdas $3 = exp/saMean, $4 = exp*prod).
  const int I = kForecast;
  const int IS = kStdevForecast;
  Opt fcS = cx.init.get(kSeries, I);
  const Opt& fcSe = cx.init.get(kSeasonal, I);
  const Opt& fcT = cx.init.get(kTrend, I);
  const Opt& fcI = cx.init.get(kIrregular, I);
  if (!fcS) {
    fcS = ts_add_v(fcSe, fcT, fcI);
  }
  if (!fcS || fcS->v.empty()) return;
  fcS = ts_exp(fcS);
  out.put(kSeries, I, fcS);
  if (cx.init.get(kSeries, IS))
    out.put(kSeries, IS, correct_stdev_for_log(*cx.init.get(kSeries, IS), *fcS));
  Opt seScaled;
  if (fcSe) {
    seScaled = cx.bias
                   ? opt_fn(fcSe, [saMean](double x) {
                       return detail::fdlibm_exp(x) / saMean;
                     })
                   : ts_exp(fcSe);
    out.put(kSeasonal, I, seScaled);
    if (cx.init.get(kSeasonal, IS))
      out.put(kSeasonal, IS,
              correct_stdev_for_log(*cx.init.get(kSeasonal, IS), *seScaled));
  }
  Opt tScaled;
  if (fcT) {
    const double m12 = saMean * irMean;
    tScaled = cx.bias
                  ? opt_fn(fcT, [m12](double x) {
                      return detail::fdlibm_exp(x) * m12;
                    })
                  : ts_exp(fcT);
    out.put(kTrend, I, tScaled);
    if (cx.init.get(kTrend, IS))
      out.put(kTrend, IS,
              correct_stdev_for_log(*cx.init.get(kTrend, IS), *tScaled));
  }
  Opt saFc = ts_div(fcS, seScaled);
  if (saFc) {
    out.put(kAdjusted, I, saFc);
    if (cx.init.get(kAdjusted, IS))
      out.put(kAdjusted, IS,
              correct_stdev_for_log(*cx.init.get(kAdjusted, IS), *saFc));
  }
  Opt irrFc = ts_div(saFc, tScaled);
  if (irrFc) {
    out.put(kIrregular, I, irrFc);
    if (cx.init.get(kIrregular, IS))
      out.put(kIrregular, IS,
              correct_stdev_for_log(*cx.init.get(kIrregular, IS), *irrFc));
  }
}

void fill_bc3(Dec& out, const BiasContext& cx, double saMean, double irMean) {
  const int I = kBackcast;
  const int IS = kStdevBackcast;
  Opt bcS = cx.init.get(kSeries, I);
  const Opt& bcSe = cx.init.get(kSeasonal, I);
  const Opt& bcT = cx.init.get(kTrend, I);
  const Opt& bcI = cx.init.get(kIrregular, I);
  if (!bcS) {
    bcS = ts_add_v(bcSe, bcT, bcI);
  }
  if (!bcS || bcS->v.empty()) return;
  bcS = ts_exp(bcS);
  out.put(kSeries, I, bcS);
  if (cx.init.get(kSeries, IS))
    out.put(kSeries, IS, correct_stdev_for_log(*cx.init.get(kSeries, IS), *bcS));
  Opt seScaled;
  if (bcSe) {
    seScaled = cx.bias
                   ? opt_fn(bcSe, [saMean](double x) {
                       return detail::fdlibm_exp(x) / saMean;
                     })
                   : ts_exp(bcSe);
    out.put(kSeasonal, I, seScaled);
    if (cx.init.get(kSeasonal, IS))
      out.put(kSeasonal, IS,
              correct_stdev_for_log(*cx.init.get(kSeasonal, IS), *seScaled));
  }
  Opt tScaled;
  if (bcT) {
    const double m12 = saMean * irMean;
    tScaled = cx.bias
                  ? opt_fn(bcT, [m12](double x) {
                      return detail::fdlibm_exp(x) * m12;
                    })
                  : ts_exp(bcT);
    out.put(kTrend, I, tScaled);
    if (cx.init.get(kTrend, IS))
      out.put(kTrend, IS,
              correct_stdev_for_log(*cx.init.get(kTrend, IS), *tScaled));
  }
  Opt saBc = ts_div(bcS, seScaled);
  if (saBc) {
    out.put(kAdjusted, I, saBc);
    if (cx.init.get(kAdjusted, IS))
      out.put(kAdjusted, IS,
              correct_stdev_for_log(*cx.init.get(kAdjusted, IS), *saBc));
  }
  Opt irrBc = ts_div(saBc, tScaled);
  if (irrBc) {
    out.put(kIrregular, I, irrBc);
    if (cx.init.get(kIrregular, IS))
      out.put(kIrregular, IS,
              correct_stdev_for_log(*cx.init.get(kIrregular, IS), *irrBc));
  }
}

Dec correct_logs(const BiasContext& cx) {
  Dec out;
  const Ts& orig = cx.original;
  out.put(kSeries, kValue, orig);
  const int n = static_cast<int>(orig.v.size());
  const int trunc = n - n % cx.period;
  double saMean = 1.0, irMean = 1.0;
  Opt seas;
  if (cx.init.get(kSeasonal, kValue)) {
    seas = ts_exp(*cx.init.get(kSeasonal, kValue));
    if (cx.bias) {
      saMean = seq_average(vd(seas->v.begin(), seas->v.begin() + trunc));
      seas = ts_mult(seas, 1.0 / saMean);
    }
    out.put(kSeasonal, kValue, seas);
    if (cx.init.get(kSeasonal, kStdev))
      out.put(kSeasonal, kStdev,
              correct_stdev_for_log(*cx.init.get(kSeasonal, kStdev), *seas));
  }
  Opt irr;
  if (cx.init.get(kIrregular, kValue)) {
    irr = ts_exp(*cx.init.get(kIrregular, kValue));
    if (cx.bias) {
      irMean = seq_average(irr->v);
      irr = ts_mult(irr, 1.0 / irMean);
    }
    out.put(kIrregular, kValue, irr);
    if (cx.init.get(kIrregular, kStdev))
      out.put(kIrregular, kStdev,
              correct_stdev_for_log(*cx.init.get(kIrregular, kStdev), *irr));
  }
  Opt trend;
  if (cx.init.get(kTrend, kValue)) {
    trend = ts_exp(*cx.init.get(kTrend, kValue));
    if (cx.bias) trend = ts_mult(trend, saMean * irMean);
    out.put(kTrend, kValue, trend);
    if (cx.init.get(kTrend, kStdev))
      out.put(kTrend, kStdev,
              correct_stdev_for_log(*cx.init.get(kTrend, kStdev), *trend));
  }
  Opt sa = ts_div(Opt{orig}, seas);
  out.put(kAdjusted, kValue, sa);
  if (cx.init.get(kAdjusted, kStdev))
    out.put(kAdjusted, kStdev,
            correct_stdev_for_log(*cx.init.get(kAdjusted, kStdev), *sa));
  out.put(kIrregular, kValue, ts_div(sa, trend));  // overwrites the exp key
  fill_fc3(out, cx, saMean, irMean);
  fill_bc3(out, cx, saMean, irMean);
  return out;
}

// ---------------------------------------------------------------------------
// The extraction driver.
// ---------------------------------------------------------------------------
int extrap(int n, int period) { return n < 0 ? -period * n : n; }

seats_extract_wk_result run_extraction(const tramo_seats_model_request& request,
                                       const tramo_seats_model_result& dec) {
  // The 24-B decomposition seam, live (the same engine the seam unit and
  // the validator run; dec guarantees a non-null d inside the envelope).
  detail::Sarima m1;
  m1.period = dec.period;
  m1.d = dec.orders_d;
  m1.bd = dec.orders_bd;
  m1.phi = dec.phi;
  m1.bphi = dec.bphi;
  m1.theta = dec.theta;
  m1.btheta = dec.btheta;
  detail::EngineFacts ef;
  try {
    ef = detail::decompose(m1, request.rmod, request.epsphi, request.smod,
                           request.smod_pi, false);
  } catch (const detail::jvm_decompose_abort&) {
    ef.chain_threw = true;
  }
  UcarimaModel raw;
  detail::Arima ma1 = m1.to_arima();
  raw.m_sar = ma1.stationary_ar().c;
  raw.m_ur = ma1.non_stationary_ar().c;
  raw.m_ma = ma1.ma1;
  raw.m_var = 1.0;  // SarimaModel.getInnovationVariance() on every corpus
  for (const detail::Arima& c : ef.raw) {
    AModel a;
    a.ar = c.ar;
    a.delta = c.delta;
    a.symmetric = c.symmetric;
    if (c.symmetric) {
      a.sma = c.sma;
    } else {
      a.has_ma = true;
      a.ma1 = c.ma1;
      a.var = c.var;
    }
    raw.comps.push_back(a);
  }
  std::optional<UcarimaModel> d;
  if (!ef.chain_threw) {
    try {
      d = sm::um_set_variance_max(raw, -1, false);
    } catch (const detail::jvm_decompose_abort&) {
      d = std::nullopt;
    }
  }
  if (!d) {  // unreachable inside the envelope (both null branches refused)
    return seats_extract_wk_result({"status error:jdplus.tramoseats.base.api.seats."
                                    "SeatsException:Invalid decomposition"});
  }
  // WKE.decompose 13-15: UcarimaModel.compact(2, 2) is UNCONDITIONAL; the
  // JVM indexes comps[2] (plus the s-minus-1 merges) without a guard - a
  // short model is an ArrayIndexOutOfBoundsException (unreachable inside
  // the envelope: every corpus d has four components).
  if (static_cast<int>(d->comps.size()) < 3) throw sm::Aioobe{};
  UcarimaModel um3 = sm::um_compact(*d, 2, 2);

  AModel outer;
  outer.ar = sm::Poly{raw.m_sar};
  outer.delta = sm::Poly{raw.m_ur};
  outer.symmetric = false;
  outer.has_ma = true;
  outer.ma1 = raw.m_ma;
  outer.var = raw.m_var;  // BurmanEstimates.initModel reads ucm.getModel():
  // the UcarimaModel's OWN model (var 1 on the corpus path), NOT the spec
  // echo. The spec variance enters the extraction only through the builder's
  // innovationStdev = sqrt(model.getInnovationVariance()) below.

  // The SeatsModel frame: transformed series, extrapolation counts, ser.
  const std::int64_t m0 =
      static_cast<std::int64_t>(request.start_year) * 12 + (request.start_month - 1);
  Ts data{m0, request.series};
  if (request.log) data = ts_fn(data, detail::fdlibm_log);  // TsData.log
  const int nfc = extrap(request.forecast_count, dec.period);
  const int nbc = extrap(request.backcast_count, dec.period);
  const double ser =
      std::sqrt(dec.innovation_variance);  // Math.sqrt(getInnovationVariance)

  swk::WkEstimators wk(&um3, outer);
  sb::BurmanEstimates be(data.v, wk, request.mean_correction, 0, ser, nfc,
                         nbc);
  Dec init = wke_decompose(be, um3, data, nbc);

  BiasContext cx;
  cx.original = Ts{m0, request.series};
  cx.init = std::move(init);
  cx.bias = request.bias != seats_bias_correction::none;
  cx.period = dec.period;
  Dec final = request.log ? correct_logs(cx) : fill_levels(cx);

  std::vector<std::string> out{"status ok"};
  dump_domain(out, "dom actual", cx.init.get(kAdjusted, kValue));
  dump_domain(out, "dom backcast", cx.init.get(kAdjusted, kBackcast));
  dump_domain(out, "dom forecast", cx.init.get(kAdjusted, kForecast));
  auto initial_rows = collect_decomposition(out, "init", cx.init);
  auto final_rows = collect_decomposition(out, "final", final);
  return {std::move(out), std::move(initial_rows), std::move(final_rows)};
}

std::string abort_status(const std::exception_ptr& e) {
  try {
    std::rethrow_exception(e);
  } catch (const eco_abort& a) {
    return "status error:jdplus.toolkit.base.api.eco.EcoException:" +
           a.message;
  } catch (const demetra_abort&) {
    return "status error:jdplus.toolkit.base.api.DemetraException";
  } catch (const wrong_dimensions_abort&) {
    return "status error:java.lang.IllegalArgumentException:wrong "
           "dimensions";
  } catch (const slf::matrix_exception&) {
    return "status error:jdplus.toolkit.base.core.math.matrices."
           "MatrixException";
  } catch (...) {
    return "status error:" + std::string(sm::exception_class(e));
  }
}

}  // namespace

std::string seats_extract_wk_result::fact(const std::string& key) const {
  for (const std::string& line : facts_) {
    if (line.size() >= key.size() && line.compare(0, key.size(), key) == 0 &&
        (line.size() == key.size() || line[key.size()] == ' ')) {
      const std::size_t v = line.find_first_not_of(' ', key.size());
      return v == std::string::npos ? std::string() : line.substr(v);
    }
  }
  return std::string();
}

seats_extract_wk_result seats_extract_wk(
    const tramo_seats_model_request& request) {
  tramo_seats_model_result dec = validate_tramo_seats_model(request);
  if (request.frequency != 12)
    throw std::invalid_argument(
        "seats_extract_wk: annual frequency " +
        std::to_string(request.frequency) +
        " is outside the D11 monthly-first envelope");
  if (request.components_method != seats_components_method::burman)
    throw std::invalid_argument(
        "seats_extract_wk: components_method KalmanSmoother is outside the "
        "D1 envelope (the Kalman estimator is not ported)");
  if (request.innovation_variance == 0.0)
    throw std::invalid_argument(
        "seats_extract_wk: innovation variance 0 selects the JVM "
        "ConcentratedLikelihood variance path (outside the ported domain)");
  if (request.mode == seats_approx_mode::noisy)
    throw std::invalid_argument(
        "seats_extract_wk: Noisy approximation mode is outside the D6 "
        "envelope");
  if (dec.approximation_required)
    throw std::invalid_argument(
        "seats_extract_wk: approximation_required decision is outside the "
        "D6 envelope (the DefaultModelApproximator retry loop is not "
        "ported)");

  std::vector<std::string> facts;
  if (dec.seats_exception) {
    // The JVM kernel throws SeatsException before estimation whenever the
    // decomposition came back null under mode None (the approximation loop
    // is not available); no dom/init/final rows follow.
    facts.push_back("status error:jdplus.tramoseats.base.api.seats."
                    "SeatsException:Invalid decomposition");
    return seats_extract_wk_result(std::move(facts));
  }
  try {
    return run_extraction(request, dec);
  } catch (const std::invalid_argument&) {
    throw;  // envelope violation, never swallowed
  } catch (...) {
    facts.push_back(abort_status(std::current_exception()));
  }
  return seats_extract_wk_result(std::move(facts));
}

}  // namespace kalos
