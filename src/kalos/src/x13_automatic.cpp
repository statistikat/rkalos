// kalos :: x13_automatic.cpp
//
// Item 21-H.1: the automatic-spec X-13 pipeline driver (stage Z). This file
// is the pinned RegArimaKernel.calc transcription
// (docs/contracts/x13-automatic-pipeline-contract.md, bytecode re-derived
// from the current-jar dump cur/RegArimaKernel_full.txt per the pinning
// notes SECTION 12) composed on the qualified stage-Z units: the 21-B
// decide_transformation, the 21-G run_mean_controller /
// run_regression_variable_tests / final_mean_rule_adds_mean, the 21-C
// detect_outliers with the sweep-qualified callback recipe, the 21-F
// refine_arma_model, the slice-1 model_controller_accept accept machine, the
// seven-filter pass2 gate with checkUnitRoots/checkMA/redoEstimation/
// lastSolution, and the FinalEstimator cascade on final_estimator_test over
// the 21-H.0 final_information producer. The downstream composition reuses
// the fixed fit's qualified seams in the pinned X13Kernel order
// (x13_preadjust -> x11_selection_pass -> results-builder -> x13_restore);
// the fixed path is untouched by construction (x13_fit.cpp carries no
// change; the two file-local helpers below are verbatim copies of its
// pinned stages).
//
// Implementation-actual deltas from the contract sketch are recorded in the
// contract's "Implementation-actual deltas at landing" section: the G3
// LS -> Trend resolution, the row-3 sweep recipe (Multiple + Bfgs, not the
// sketch's PinnedLM + HannanRissanen), the fresh-detection ODM adapter (the
// 21-C unit carries no addOutlier seeding), the design-row cascade deviance
// mapping (the concentrated RegressionLikelihood deviance under the 21-H.0
// ARMA-parameter stencil - the bit-pinned deviances machine covers only the
// design-free domain), the accept-refresh chain tags, and the JVM
// differenced-model estimation surface: all driver-side estimations and
// likelihood profiles run on the RegArmaModel.of surface (the filtered
// response, the Griffies mean column of raw ones first, the filtered outlier
// columns, spec_d/spec_bd carrying the differencing - probe
// jdprobe/item21h1_meantrial.jsh against the current jar), while the 21-C
// outlier machine keeps its sweep-qualified raw-domain surface with a
// mean-free base (the JVM Griffies mean is inexpressible there, recorded per
// entry with the "mean on: excluded from the ODM base (delta)" note).

#include "kalos/x13_automatic.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/hostlineage.hpp"

#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/final_information.hpp"
#include "kalos/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace kalos {
namespace {

// ---------------------------------------------------------------------------
// Transcription constants (JVM call-site literals; the decoder-derived
// scalars are request data, these are the kernel's own ldc operands).
// ---------------------------------------------------------------------------

constexpr double kCheckMuCalc = 1.96;   // calc's checkMu(m, 1.96)
constexpr double kCheckMuRedo = 1.6;    // redoEstimation's checkMu(m, 1.6)
constexpr double kPcrLoop1 = 0.025;     // pass2's pcr increment, loop == 1
constexpr double kPcrOther = 0.015;     // pass2's pcr increment, else
constexpr double kPlboxMid = 0.75;      // pass2 filters 1/3
constexpr double kFct = 1.0126582278481011;  // ModelController.FCT (filter 3)
constexpr double kFilter5Phi = -0.82;   // filter 5 phi(1) shape bound
constexpr double kFilter6Bphi = -0.65;  // filter 6 bphi(1) shape bound
constexpr double kMinCurva = 2.8;       // MINCV, the curva reduction floor
constexpr double kImWindow = 0.05;      // searchur |im| window
constexpr double kMaSingular = 1e-3;    // checkMA |MA(1)| rejection bound
constexpr int kMaxRounds = 5;           // the round < 5 loop bound
constexpr int kMaxCascadeRounds = 5;    // FinalEstimator.estimate's loop
constexpr int kAnnualFrequency = 12;
constexpr double kFeps = 1e-4;  // RegSarimaComputer DEF_INTERNAL_EPS

// ---------------------------------------------------------------------------
// JVM comparison flavours. dcmpg pushes +1 for a NaN operand and dcmpl -1
// (JVM20 §6.5), so every guarded comparison of the pinned kernel - "dcmpg
// (a, b) ifge skip" (the a < b guard), "dcmpl(a, b) iflt skip" (the a >= b
// guard), "dcmpg(a, b) ifgt take" (the a > b take) - reduces to the plain
// IEEE spelling with NaN on the FALSE (skip/take-the-other-branch) side:
// jvm_lt/jvm_ge/jvm_gt below are those readings. The only NaN-true readings
// are the bare dcmpg/dcmpl stores (none here) and Math.max, where a NaN
// operand propagates (java_max below; std::max would silently drop it).
// ---------------------------------------------------------------------------

bool jvm_lt(double a, double b) { return a < b; }   // false for NaN
bool jvm_ge(double a, double b) { return a >= b; }  // false for NaN
bool jvm_gt(double a, double b) { return a > b; }   // false for NaN

// Java Math.max: NaN on either side returns NaN.
double java_max(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::max(a, b);
}

// ---------------------------------------------------------------------------
// Polynomials. The Kalos convention: AR block phi with polynomial 1 - phi(B),
// MA block theta with 1 + theta(B). The JVM SarimaModel stores the AR with
// the negated signs ("1 + phi(B)" storage with phi = -phiK, per the
// ssq_deviances header), so the raw JVM coordinate vector
// [phi | bphi | theta | btheta] is [-arK | -sarK | maK | smaK].
// ---------------------------------------------------------------------------

// SarimaModel.getRegularAR().asPolynomial() (ascending): [1, -phiK ...].
std::vector<double> jvm_ar_ascending(const std::vector<double>& arK) {
  std::vector<double> out(arK.size() + 1, 0.);
  out[0] = 1.;
  for (std::size_t i = 0; i < arK.size(); ++i) {
    out[i + 1] = -arK[i];
  }
  return out;
}

// SarimaModel.getSeasonalAR().asPolynomial() over B (ascending):
// prod_j (1 - sarK_j B^{(j+1) P}).
std::vector<double> jvm_seasonal_ar_ascending(const std::vector<double>& sarK,
                                              int period) {
  std::vector<double> poly{1.};
  for (std::size_t j = 0; j < sarK.size(); ++j) {
    std::vector<double> factor((j + 1) * static_cast<std::size_t>(period) + 1,
                               0.);
    factor[0] = 1.;
    factor[(j + 1) * static_cast<std::size_t>(period)] = -sarK[j];
    std::vector<double> next(poly.size() + factor.size() - 1, 0.);
    for (std::size_t i = 0; i < poly.size(); ++i) {
      for (std::size_t k = 0; k < factor.size(); ++k) {
        next[i + k] += poly[i] * factor[k];
      }
    }
    poly = std::move(next);
  }
  return poly;
}

// Polynomial.mirror(): the reversed coefficient sequence (ascending input,
// ascending output of the mirrored polynomial).
std::vector<double> mirror_ascending(const std::vector<double>& ascending) {
  return std::vector<double>(ascending.rbegin(), ascending.rend());
}

// SarimaModel.getRegularMA().evaluateAt(1.0) = 1 + sum theta (the Kalos MA
// block matches the JVM MA storage: 1 + theta(B)).
double ma_evaluate_at_1(const std::vector<double>& maK) {
  double value = 1.;
  for (const double theta : maK) {
    value += theta;
  }
  return value;
}

int arma_dim_of(const std::array<std::vector<double>, 4>& blocks) {
  return static_cast<int>(blocks[0].size() + blocks[1].size() +
                          blocks[2].size() + blocks[3].size());
}

// ---------------------------------------------------------------------------
// Period labels / outlier names. IOutlier.defaultName(code, period) probed
// through jshell: "AO (1951-05-01)" for TsPeriod.monthly(1951, 5); the
// monthly day is always 01.
// ---------------------------------------------------------------------------

std::string outlier_name(OutlierType type, std::size_t position,
                         int first_year, int first_period) {
  const char* code = "AO";
  switch (type) {
    case OutlierType::Additive: code = "AO"; break;
    case OutlierType::LevelShift: code = "LS"; break;
    case OutlierType::TransitoryChange: code = "TC"; break;
    case OutlierType::SeasonalOutlier: code = "SO"; break;
  }
  const int index = first_period - 1 + static_cast<int>(position);
  const int year = first_year + index / 12;
  const int month = index % 12 + 1;
  std::ostringstream stream;
  if (year < 0) {
    stream << '-';
  }
  const unsigned magnitude = year < 0 ? static_cast<unsigned>(-year)
                                      : static_cast<unsigned>(year);
  if (magnitude < 10) stream << "000";
  else if (magnitude < 100) stream << "00";
  else if (magnitude < 1000) stream << '0';
  stream << magnitude << '-';
  if (month < 10) stream << '0';
  stream << month << "-01";
  return std::string(code) + " (" + stream.str() + ")";
}

// SaVariable.defaultComponentTypeOf(IOutlier), re-derived from the bytecode
// (G3 probe, 21-H.1: AO/TransitoryChange/SwitchOutlier -> Irregular,
// LevelShift -> Trend, PeriodicOutlier -> Seasonal; the contract's
// "LS -> LevelShift" expectation was superseded by the probe - LevelShift is
// not a ComponentType).
X13Component outlier_effect(OutlierType type) {
  switch (type) {
    case OutlierType::Additive:
    case OutlierType::TransitoryChange:
      return X13Component::irregular;
    case OutlierType::LevelShift:
      return X13Component::trend;
    case OutlierType::SeasonalOutlier:
      return X13Component::seasonal;
  }
  return X13Component::undefined;
}

const char* component_name(X13Component effect) {
  switch (effect) {
    case X13Component::trend: return "Trend";
    case X13Component::seasonal: return "Seasonal";
    case X13Component::irregular: return "Irregular";
    default: return "Undefined";
  }
}

// The pinned X13Kernel.updateSpec forcing (verbatim copy of the x13_fit.cpp
// file-local stage; the fixed path keeps its own byte-identical copy).
X11DecompositionMode forceMode(bool log, X11DecompositionMode mode) {
    if(!log) return X11DecompositionMode::additive;
    if(mode==X11DecompositionMode::log_additive
        ||mode==X11DecompositionMode::pseudo_additive) return mode;
    return X11DecompositionMode::multiplicative;
}

// The pinned X11Kernel results-builder stage (verbatim x13_fit.cpp copy).
void logAdditiveResults(const std::vector<double>& step10,
    const std::vector<double>& step12, const std::vector<double>& step13,
    const std::vector<double>& input, std::vector<double>& d10,
    std::vector<double>& d11, std::vector<double>& d12,
    std::vector<double>& d13) {
    const std::size_t n=step10.size();
    std::vector<double> A(n),B(n),C(n);
    double sa=0.,sc=0.;
    for(std::size_t i=0;i<n;++i){
        // Raw libm, intentionally (CI-FDLIBM ledger): third copy of the
        // X13Kernel log-additive results builder (A/B/C). Its bit contract is
        // the x13-allocation/x13-fit fixture corpora, captured WITHOUT the
        // DisableIntrinsic recipe, so those fixtures froze the HotSpot _dexp
        // stub bits (IBM-APML lineage: glibc and the macOS libm derive from
        // it, which is why host-libm bits hold the pinned 1e-15 exp-boundary
        // class on every qualified platform, while fdlibm_exp moves 1-2 ulp
        // at the q_log export arguments). The bytecode-recipe consumers
        // (oracle/sweep) reach these tables through tolerance bands only.
        // See src/x13_allocation.cpp for the recorded measurement.
        A[i]=kalos::detail::glibc_exp(step10[i]); B[i]=kalos::detail::glibc_exp(step12[i]); C[i]=kalos::detail::glibc_exp(step13[i]);
        sa+=A[i]; sc+=C[i];
    }
    const double v9=sa/static_cast<double>(n),
        v13=v9*(sc/static_cast<double>(n));
    d10.resize(n); d11.resize(n); d12.resize(n); d13.resize(n);
    for(std::size_t i=0;i<n;++i){
        d10[i]=A[i]/v9;
        d12[i]=B[i]*v13;
        d11[i]=input[i]/d10[i];
        d13[i]=(input[i]/d10[i])/d12[i];
    }
}

// ---------------------------------------------------------------------------
// The modelling state (the pinned RegSarimaModelling + ModelDescription
// fields this kernel reads and writes: log, mean, specification, ami
// outliers, the retained candidate).
// ---------------------------------------------------------------------------

struct OutlierVariable {
  OutlierType type = OutlierType::Additive;
  std::size_t position = 0;
  std::string name;
  X13Component effect = X13Component::irregular;
};

struct Estimation {
  std::array<std::vector<double>, 4> coefficients;  // [ar, ma, sar, sma]
  std::string status;
  double log_likelihood = 0.;
};

struct Description {
  bool log = false;
  bool mean = false;
  // The specification (SarimaOrders): stationary part + differencing. The
  // X13ModelBuilder start is the AIRLINE (0,1,1)(0,1,1) - q and bq are 1
  // from the start (odmscan spec echoes: "0 1 1 0 1 1"); the ODM/refine
  // units size their whitening and search from these orders.
  int p = 0, q = 1, bp = 0, bq = 1;
  int d = 1, bd = 1;  // the airline differencing (freq != 1)
  std::vector<OutlierVariable> variables;  // ami outliers (no user
                                           // variables inside the envelope)
};

struct Modelling {
  Description desc;
  std::optional<Estimation> estimation;  // nullopt == needEstimation()
  std::array<std::vector<double>, 4> last_coefficients;  // desc.arima()
  bool has_coefficients = false;
};

// ---------------------------------------------------------------------------
// The driver.
// ---------------------------------------------------------------------------

class AutomaticPipeline {
 public:
  AutomaticPipeline(const X13AutomaticRequest& request,
                    X13AutomaticFitResult& result)
      : request_(request), result_(result) {}

  // The failed: stage/round of the last calc position (entry-point use).
  const std::string& stage() const { return stage_; }
  int round() const { return round_; }

  void run() {
    transitory_rate_ = transitory_change_rate(kAnnualFrequency,
                                             request_.tcrate);
    y_ = request_.input;  // the level domain; log below if chosen

    // (init) the X13ModelBuilder airline start + the critical-value rule
    // (process(): va0 = options.va > 0 ? va : X13Utility.calcCv(n)).
    model_.desc = Description{};  // log=false, mean=false, (0,1,1)(0,1,1)
    pcr_ = request_.ljung_box_limit;  // clear(): pcr = options.ljungBoxLimit
    va0_ = request_.va > 0. ? request_.va
                            : default_critical_value(y_.size());
    curva_ = va0_;
    need_outliers_ = request_.outliers_used;
    need_auto_modelling_ = true;  // autoModel != null for this kind
    {
      AutoDecisionRecord entry;
      entry.stage = "init";
      entry.scalars = {{"curva", curva_},
                       {"pcr", pcr_},
                       {"n", static_cast<double>(y_.size())}};
      entry.notes = {request_.va > 0. ? "curva=va" : "curva=calcCv(n)",
                     "airline start (0,1,1)(0,1,1), mean off"};
      result_.chain.push_back(std::move(entry));
    }

    // calc: desc.isLogTransformation() && any(v <= 0) -> error. The log flag
    // only turns on through the transformation decision below, which
    // short-circuits to levels on any non-positive value; the transcription
    // keeps the guard for the log-chosen state.
    if (model_.desc.log) {
      for (const double value : request_.input) {
        if (!(value > 0.)) {
          throw std::runtime_error(
              "can't apply log transformation 'some obs. are <= 0");
        }
      }
    }

    stage_ = "transformation";
    transformation_prepass();

    stage_ = "meantest";
    // regAIC: calendarTest/easterTest fail-closed (21-I) -> Unchanged.
    check_mu(kCheckMuCalc, "meantest0", 0, nullptr);

    if (need_outliers_) {
      stage_ = "outliers";
      const bool changed =
          process_outliers(curva_, 0, "outliers0", nullptr);
      if (changed) {
        stage_ = "estimation";
        estimate_if_needed();
        stage_ = "regtest";
        run_regtest(request_.mean_test_initial, "regtest0", 0);
      }
    }

    // isAutoModelling() is always true for this kind: the non-auto branch
    // (finalEstimator.estimate + build) is unreachable here.
    stage_ = "estimation";
    estimate_if_needed();
    reference_ = std::make_unique<Modelling>(model_);
    stage_ = "accept";
    accept_refresh("accept0", 0, nullptr);
    const bool accepted_airline = last_accept_ && last_accept_->accepted;

    if (!(request_.accept_airline && accepted_airline)) {
      loop_ = 1;
      for (int round = 1; round <= kMaxRounds; ++round) {
        round_ = round;
        stage_ = "round";
        if (calc_round(round)) {
          break;  // the cascade ended the loop
        }
      }
    }

    stage_ = "final";
    build_result();
  }

 private:
  // ---- chain plumbing -------------------------------------------------------

  // Chain entries are appended by value (references into result_.chain must
  // never be held across appends: the vector reallocates).
  void push(AutoDecisionRecord entry) {
    result_.chain.push_back(std::move(entry));
  }

  // ---- the working data ------------------------------------------------------

  std::vector<double> y_;  // the transformed (level or log) working series
  double transitory_rate_ = 0.;

  // Whether the retained estimation was fitted with the mean column in the
  // design. The transformation prepass retains the CHOSEN-branch estimation,
  // which is always a mean-free fit; a mean-on description paired with it is
  // design-mismatched (JVM estimate-if-needed semantics re-fit before the
  // mean controller reads the cl statistics).
  bool estimation_has_mean_ = false;

  // The JVM differenced-model surface (RegArmaModel.of, current jar, probed
  // 21-H.1): the response is the unit-root-filtered series and the mean is a
  // GRIFFIES mean - the differenced design's FIRST column carries raw 1.0
  // (bytecode 328-342/421-435: isMean() writes iterator.next().set(1.0)
  // ahead of the BackFilter.apply of every X column). The Kalos estimation
  // and profile units take this surface as the pre-differenced mode (d = 0,
  // spec_d/spec_bd carry the spec differencing); the raw-domain mean column
  // would difference to zero and the G2 whitening reports it rank-deficient.
  std::vector<double> yd_;   // the differenced response (n - d - 12*bd rows)
  std::vector<double> yd_filter_;  // (1 - B)^d (1 - B^12)^bd, ascending

  void ensure_domain() {
    yd_filter_ = detail::pinned_differencing_polynomial(
        kAnnualFrequency, model_.desc.d, model_.desc.bd);
    yd_ = (model_.desc.d == 0 && model_.desc.bd == 0)
              ? y_
              : detail::pinned_backfilter_apply(yd_filter_, y_);
  }

  std::vector<double> diff_column(std::vector<double> column) const {
    if (model_.desc.d == 0 && model_.desc.bd == 0) {
      return column;
    }
    return detail::pinned_backfilter_apply(yd_filter_, column);
  }

  std::vector<double> outlier_column(const OutlierVariable& variable) const {
    std::vector<double> column(y_.size(), 0.);
    fill_outlier_column(variable.type, variable.position, kAnnualFrequency,
                        transitory_rate_, column);
    return column;
  }

  // The design (mean column first, then one column per ami outlier in model
  // order). All coefficients free.
  RegressionDesign build_design(
      const std::vector<OutlierVariable>* variables = nullptr) {
    ensure_domain();
    const std::vector<OutlierVariable>* set =
        variables != nullptr ? variables : &model_.desc.variables;
    const std::size_t n = yd_.size();
    RegressionDesign design;
    std::vector<std::vector<double>> columns;
    if (model_.desc.mean) {
      columns.emplace_back(n, 1.);  // the Griffies mean column (raw ones)
    }
    for (const OutlierVariable& variable : *set) {
      columns.push_back(diff_column(outlier_column(variable)));
    }
    design.columns = columns.size();
    design.values.assign(design.columns * n, 0.);
    design.coefficients.assign(design.columns, 0.);
    design.fixed.assign(design.columns, 0);
    for (std::size_t j = 0; j < columns.size(); ++j) {
      for (std::size_t i = 0; i < n; ++i) {
        design.values[j * n + i] = columns[j][i];
      }
    }
    return design;
  }

  // The raw-domain design of the 21-C sweep-qualified surface (contract
  // delta): the detection base NEVER carries the mean column - the JVM mean
  // lives as raw ones in the differenced domain and cannot be represented on
  // this surface (a raw constant annihilates; 21-D ruled a drift column out).
  RegressionDesign build_raw_design(
      const std::vector<OutlierVariable>* variables = nullptr) const {
    const std::vector<OutlierVariable>* set =
        variables != nullptr ? variables : &model_.desc.variables;
    const std::size_t n = y_.size();
    RegressionDesign design;
    std::vector<std::vector<double>> columns;
    for (const OutlierVariable& variable : *set) {
      columns.push_back(outlier_column(variable));
    }
    design.columns = columns.size();
    design.values.assign(design.columns * n, 0.);
    design.coefficients.assign(design.columns, 0.);
    design.fixed.assign(design.columns, 0);
    for (std::size_t j = 0; j < columns.size(); ++j) {
      for (std::size_t i = 0; i < n; ++i) {
        design.values[j * n + i] = columns[j][i];
      }
    }
    return design;
  }

  // The base estimation request over the current design (the JVM
  // differenced-model surface): the series carries the filtered response and
  // the design the differenced-domain columns, the differencing travels in
  // spec_d/spec_bd (the estimation units' pre-differenced mode).
  SarimaEstimationRequest base_request(
      double eps, EstimationStartingPoint start, bool fast_derivatives,
      const RegressionDesign* design = nullptr) {
    ensure_domain();
    SarimaEstimationRequest out;
    out.series.values = yd_;
    out.series.missing.assign(yd_.size(), 0);
    out.design = design != nullptr ? *design : build_design();
    out.blocks[0].initial.assign(static_cast<std::size_t>(model_.desc.p), 0.);
    out.blocks[1].initial.assign(static_cast<std::size_t>(model_.desc.q), 0.);
    out.blocks[2].initial.assign(static_cast<std::size_t>(model_.desc.bp), 0.);
    out.blocks[3].initial.assign(static_cast<std::size_t>(model_.desc.bq), 0.);
    for (auto& block : out.blocks) {
      block.fixed.assign(block.initial.size(), 0);
    }
    out.period = kAnnualFrequency;
    out.d = 0;
    out.seasonal_d = 0;
    out.spec_d = model_.desc.d;
    out.spec_bd = model_.desc.bd;
    out.control.eps = eps;
    out.control.feps = kFeps;
    out.control.maximum_likelihood = true;
    out.control.unit_root_search = false;
    out.control.corrected_degrees_of_freedom = true;
    out.control.fast_derivatives = fast_derivatives;
    out.control.starting_point = start;
    // CI-FDLIBM: the landing estimate runs on the transcribed JVM LM machine
    // (RegArmaSsqFunction + Levenberg-Marquardt, detail/ssq_levmar.hpp — the
    // same machine the differencing and initial-ARMA paths already replay
    // bit-exactly), superseding the Bfgs product default of maintainer
    // ruling #2. Rationale (measured 2026-11-21; contract "CI-FDLIBM routing"
    // section): with the transcendental surface pinned to the JVM's bytecode
    // semantics, Bfgs landing deterministically lands in a different basin
    // than the JVM's LM landing on the four knife-edge cases (8/36 out of the
    // 1e-2 coefficient band, worst-rel 0.0622 — the CI run-#17 signature,
    // machine-independent); on the LM machine the family collapses to 0 bad,
    // 15/36 bit-exact, worst-rel 5.75e-4, legacy suite 140/140, shim/hwcaps
    // invariant. The solver remains explicitly set — the "never implicit"
    // half of ruling #2 stands.
    out.control.solver = EstimationSolver::PinnedLevenbergMarquardt;
    return out;
  }

  // The 21-C sweep-qualified raw-domain request (the outlier unit generates
  // its candidate columns on its series' domain; the pinned differencing is
  // carried on d/seasonal_d exactly like the sweep's base requests).
  SarimaEstimationRequest raw_request(
      double eps, EstimationStartingPoint start, bool fast_derivatives,
      const RegressionDesign* design = nullptr) const {
    SarimaEstimationRequest out;
    out.series.values = y_;
    out.series.missing.assign(y_.size(), 0);
    out.design = design != nullptr ? *design : build_raw_design();
    out.blocks[0].initial.assign(static_cast<std::size_t>(model_.desc.p), 0.);
    out.blocks[1].initial.assign(static_cast<std::size_t>(model_.desc.q), 0.);
    out.blocks[2].initial.assign(static_cast<std::size_t>(model_.desc.bp), 0.);
    out.blocks[3].initial.assign(static_cast<std::size_t>(model_.desc.bq), 0.);
    for (auto& block : out.blocks) {
      block.fixed.assign(block.initial.size(), 0);
    }
    out.period = kAnnualFrequency;
    out.d = model_.desc.d;
    out.seasonal_d = model_.desc.bd;
    out.control.eps = eps;
    out.control.feps = kFeps;
    out.control.maximum_likelihood = true;
    out.control.unit_root_search = false;
    out.control.corrected_degrees_of_freedom = true;
    out.control.fast_derivatives = fast_derivatives;
    out.control.starting_point = start;
    out.control.solver = EstimationSolver::Bfgs;
    return out;
  }

  MaskedSeries make_series() const {
    MaskedSeries series;
    series.values = y_;
    series.missing.assign(y_.size(), 0);
    return series;
  }

  MaskedSeries make_differenced_series() {
    ensure_domain();
    MaskedSeries series;
    series.values = yd_;
    series.missing.assign(yd_.size(), 0);
    return series;
  }

  // ---- the retained candidate -------------------------------------------------

  // The row-1 kernel-round estimation (m.estimate(options.precision):
  // RegSarimaComputer HR start, PinnedLevenbergMarquardt per the contract
  // adapter row 1 — CI-FDLIBM ruling; the Bfgs default it replaced is
  // recorded in the contract's CI-FDLIBM section).
  void estimate_row1() {
    auto request =
        base_request(request_.precision, EstimationStartingPoint::HannanRissanen,
                     true);
    SarimaEstimationResult fit = estimate_sarima_fixed(request);
    Estimation estimation;
    estimation.coefficients = fit.coefficients;
    estimation.status = std::move(fit.status);
    estimation.log_likelihood = fit.log_likelihood;
    model_.last_coefficients = estimation.coefficients;
    model_.has_coefficients = true;
    model_.estimation = std::move(estimation);
    estimation_has_mean_ = model_.desc.mean;  // build_design carries the mean
  }

  void estimate_if_needed() {
    if (!model_.estimation) {
      estimate_row1();
    }
  }

  void clear_estimation() {
    model_.estimation.reset();
    estimation_has_mean_ = false;
  }

  // The retained candidate must exist wherever the JVM dereferences the
  // estimation (its NPE maps onto our failed: status).
  const Estimation& estimation_required() const {
    if (!model_.estimation) {
      throw std::runtime_error("estimation required but cleared");
    }
    return *model_.estimation;
  }

  // The G2 landing composition: the regression likelihood at the retained
  // ARMA candidate on the differenced-model surface (the 21-C sweep recipe,
  // basin-tier qualified; the Griffies mean column of ones is the JVM one).
  RegressionLikelihood landing_fit() {
    const Estimation& estimation = estimation_required();
    return sarima_regression_likelihood(
        make_differenced_series(), build_design(), estimation.coefficients[0],
        estimation.coefficients[1], estimation.coefficients[2],
        estimation.coefficients[3], kAnnualFrequency, 0, 0);
  }

  // The concentrated-likelihood statistics view (RegressionTestStatistics)
  // at a likelihood: the dimension recovery follows the sweep recipe,
  // round(ssq / profiled variance), the finite count otherwise.
  static RegressionTestStatistics statistics_of(
      const RegressionLikelihood& fit) {
    RegressionTestStatistics statistics;
    statistics.coefficients = fit.coefficients;
    statistics.unscaled_covariance = fit.unscaled_covariance;  // symmetric
    statistics.ssq = fit.sum_squared_standardized;
    int dim = 0;
    if (std::isfinite(fit.profiled_variance) && fit.profiled_variance > 0. &&
        std::isfinite(fit.sum_squared_standardized)) {
      const double rounded =
          std::round(fit.sum_squared_standardized / fit.profiled_variance);
      if (rounded >= 0. && rounded <= static_cast<double>(fit.finite_count)) {
        dim = static_cast<int>(rounded);
      }
    }
    if (dim <= 0) {
      dim = static_cast<int>(fit.finite_count);
    }
    statistics.dim = dim;
    return statistics;
  }

  RegressionTestStatistics landing_statistics() {
    RegressionLikelihood fit = landing_fit();
    if (fit.status != "available") {
      throw std::runtime_error("landing likelihood unavailable: " +
                               fit.status);
    }
    return statistics_of(fit);
  }

  // ---- the pre-round machines --------------------------------------------------

  // LogLevelModule.process: the 21-B decision; the CHOSEN branch estimation
  // is retained in the modelling (pinning note 12.8).
  void transformation_prepass() {
    SarimaEstimationRequest branch;
    branch.series.values = request_.input;  // the LEVEL data; the unit
                                            // builds the log twin
    branch.series.missing.assign(request_.input.size(), 0);
    branch.design.columns = 0;  // mean off, no variables
    branch.blocks[0].initial.clear();
    branch.blocks[1].initial.assign(1, 0.);  // airline q = 1
    branch.blocks[2].initial.clear();
    branch.blocks[3].initial.assign(1, 0.);  // airline bq = 1
    for (auto& block : branch.blocks) {
      block.fixed.assign(block.initial.size(), 0);
    }
    branch.period = kAnnualFrequency;
    branch.d = 1;
    branch.seasonal_d = 1;
    branch.control.eps = request_.precision;  // readTransformation 1e-7
    branch.control.feps = kFeps;
    branch.control.maximum_likelihood = true;
    branch.control.unit_root_search = false;
    branch.control.corrected_degrees_of_freedom = true;
    branch.control.fast_derivatives = true;
    branch.control.starting_point = EstimationStartingPoint::Multiple;
    branch.control.solver = EstimationSolver::Bfgs;

    TransformationDecisionRequest decision_request;
    decision_request.branch = std::move(branch);
    decision_request.aicc_log_correction = request_.aicc_log_correction;
    TransformationDecisionResult decision =
        decide_transformation(decision_request);

    AutoDecisionRecord entry;
    entry.stage = "transformation";
    entry.scalars = {{"log_chosen", decision.log_chosen ? 1. : 0.},
                     {"aicc_level", decision.aicc_level},
                     {"aicc_log", decision.aicc_log}};
    entry.notes = {decision.outcome};
    entry.transformation = decision;
    push(std::move(entry));

    if (decision.outcome == "decided") {
      const TransformationBranchReport& chosen =
          decision.log_chosen ? decision.log : decision.level;
      if (decision.log_chosen) {
        model_.desc.log = true;
        y_.clear();
        y_.reserve(request_.input.size());
        for (const double value : request_.input) {
          y_.push_back(detail::fdlibm_log(value));  // decided-log implies all > 0
        }
      }
      if (chosen.available) {
        Estimation estimation;
        estimation.coefficients = chosen.estimation.coefficients;
        estimation.status = chosen.estimation.status;
        estimation.log_likelihood = chosen.log_likelihood;
        model_.last_coefficients = estimation.coefficients;
        model_.has_coefficients = true;
        model_.estimation = std::move(estimation);
        estimation_has_mean_ = false;  // the branch fit is always mean-free
      }
    }
  }

  // checkMu(m, cval): MeanController(cval).test(m) through the 21-G unit.
  // With the mean off the machine estimates the mean-added trial at
  // mean_precision (the row-3 recipe); with the mean on it reads the
  // CURRENT estimation's cl.
  bool check_mu(double cval, const char* tag, int round, const char* note) {
    if (!request_.check_mu) {
      return false;  // Unchanged, the machine never runs
    }
    MeanControllerRequest controller;
    if (model_.desc.mean && model_.estimation && estimation_has_mean_) {
      // JVM checkMu reads the CURRENT estimation's cl. The retained
      // estimation exists exactly where the landing composition is
      // 21-G-qualified to reproduce it (redoEstimation's checkMu(1.6),
      // the final mean rule); elsewhere — the round-0 test precedes the
      // first estimate — the JVM's "current estimation" is the mean-included
      // differenced-design fit, which the mean-on landing surface cannot
      // express (the Griffies column is inexpressible there).
      controller.statistics = landing_statistics();
    } else {
      // RegArimaUtility.processor(true, 1e-5) on the mean-added
      // description: a trial estimation with the mean column forced on.
      const bool mean_before = model_.desc.mean;
      model_.desc.mean = true;
      RegressionDesign trial_design = build_design();
      model_.desc.mean = mean_before;
      auto trial = base_request(request_.mean_precision,
                                EstimationStartingPoint::Multiple, true,
                                &trial_design);
      SarimaEstimationResult fit = estimate_sarima_fixed(trial);
      RegressionLikelihood likelihood = sarima_regression_likelihood(
          make_differenced_series(), trial_design, fit.coefficients[0],
          fit.coefficients[1], fit.coefficients[2], fit.coefficients[3],
          kAnnualFrequency, 0, 0);
      if (likelihood.status != "available") {
        throw std::runtime_error("mean-added trial likelihood unavailable: " +
                                 likelihood.status);
      }
      controller.statistics = statistics_of(likelihood);
    }
    controller.mean = model_.desc.mean;
    controller.cval = cval;
    MeanControllerResult controller_result = run_mean_controller(controller);
    if (controller_result.changed) {
      model_.desc.mean = controller_result.mean;
      clear_estimation();
    } else if (model_.desc.mean && !estimation_has_mean_) {
      // The mean stayed on and the retained estimation is design-mismatched
      // (mean-free branch fit): the JVM's estimate-if-needed re-fit replaced
      // it before the machine read the cl; drop the stale carrier so the
      // next estimate_if_needed rebuilds the current-design estimation.
      clear_estimation();
    }
    AutoDecisionRecord entry;
    entry.stage = tag;
    entry.round = round;
    entry.scalars = {{"cval", cval},
                     {"tstat", controller_result.tstat},
                     {"mean_after", controller_result.mean ? 1. : 0.},
                     {"changed", controller_result.changed ? 1. : 0.}};
    if (note != nullptr) {
      entry.notes = {note};
    }
    entry.mean = controller_result;
    push(std::move(entry));
    return controller_result.changed;
  }

  // ---- the outlier machine (fresh-detection adapter, contract delta) ----------

  bool process_outliers(double cv, int round, const char* tag,
                        const char* note) {
    if (!request_.outliers_used) {
      return false;
    }
    // The base design carries no outlier columns (the 21-C unit appends the
    // candidates itself) and never the mean (contract delta: the JVM mean is
    // a differenced-domain Griffies column, inexpressible on this surface;
    // the sweep-qualified base is the mean-free raw design).
    RegressionDesign base_design = build_raw_design(&empty_variables_);
    OutlierScoringRequest scoring;
    scoring.base = raw_request(request_.precision,
                               EstimationStartingPoint::Multiple, true,
                               &base_design);
    // The pinned ODM chain estimates first when the model holds no
    // estimation (mod.estimate before the detector reads estimation.getArima()),
    // so the whitening parameters are the RETAINED estimation - never the
    // spec-change-surviving desc.arima() carrier (which may carry stale
    // block shapes after a refinement changed the orders).
    estimate_if_needed();
    scoring.initial_coefficients = estimation_required().coefficients;
    scoring.types = request_.outlier_types;
    scoring.critical_value = cv;
    scoring.transitory_rate = transitory_rate_;
    scoring.lower_bound = 0;
    scoring.upper_bound = y_.size();
    scoring.max_round = request_.outlier_max_round;
    scoring.max_outliers = request_.outlier_max_outliers;

    // The sweep-qualified callback recipe (contract adapter row 3).
    const auto estimate = [](const SarimaEstimationRequest& candidate,
                             const Interrupt&) {
      return callback_estimation(estimate_sarima_fixed(candidate), candidate);
    };
    const auto optimize = [](const SarimaEstimationRequest& candidate,
                             const std::array<std::vector<double>, 4>& start,
                             const Interrupt&) {
      return callback_estimation(
          optimize_sarima_from(candidate, start), candidate);
    };
    OutlierScoringResult outcome =
        detect_outliers(scoring, estimate, optimize);

    const bool changed = outcome.changed;
    if (changed) {
      // ODM.process's Changed path: remove all ami outliers, re-add the
      // selected set with the pinned attributes, clearEstimation.
      model_.desc.variables.clear();
      for (const OutlierDecision& decision : outcome.outliers) {
        OutlierVariable variable;
        variable.type = decision.type;
        variable.position = decision.position;
        variable.effect = outlier_effect(decision.type);
        variable.name =
            outlier_name(decision.type, decision.position,
                         request_.first_year, request_.first_period);
        model_.desc.variables.push_back(std::move(variable));
      }
      clear_estimation();
    }
    AutoDecisionRecord entry;
    entry.stage = tag;
    entry.round = round;
    entry.scalars = {
        {"cv", cv},
        {"changed", changed ? 1. : 0.},
        {"n_outliers", static_cast<double>(outcome.outliers.size())},
        {"rounds", static_cast<double>(outcome.rounds)}};
    entry.notes = {outcome.stop_reason};
    if (model_.desc.mean) {
      entry.notes.push_back("mean on: excluded from the ODM base (delta)");
    }
    if (note != nullptr) {
      entry.notes.push_back(note);
    }
    entry.outliers = std::move(outcome);
    push(std::move(entry));
    return changed;
  }

  // DriverCallbacks::from from the 21-C sweep (verbatim recipe).
  static OutlierEstimation callback_estimation(
      const SarimaEstimationResult& fit_result,
      const SarimaEstimationRequest& candidate) {
    OutlierEstimation out;
    out.available = false;
    out.status = fit_result.status;
    out.coefficients = fit_result.coefficients;
    const RegressionLikelihood fit = [&] {
      try {
        return sarima_regression_likelihood(
            candidate.series, candidate.design, fit_result.coefficients[0],
            fit_result.coefficients[1], fit_result.coefficients[2],
            fit_result.coefficients[3], candidate.period, candidate.d,
            candidate.seasonal_d);
      } catch (const std::exception&) {
        RegressionLikelihood bad;  // mirrors the pinned null estimation
        return bad;
      }
    }();
    out.status += "/" + fit.status;
    if (fit.status != "available" ||
        fit.free_columns.size() != candidate.design.columns) {
      return out;
    }
    double dim = 0.0;
    if (std::isfinite(fit.profiled_variance) && fit.profiled_variance > 0.0) {
      dim = std::llround(fit.sum_squared_standardized / fit.profiled_variance);
    }
    if (!(dim > 0.0)) {
      dim = static_cast<double>(fit.finite_count);
    }
    out.t_statistics =
        regression_t_statistics(fit, static_cast<std::size_t>(dim));
    out.available = true;
    return out;
  }

  // ---- the accept machine -------------------------------------------------------

  static void fold_column(std::vector<double>& residuals,
                          const std::vector<double>& filter,
                          const std::vector<double>& column,
                          double coefficient) {
    if (coefficient == 0.) {
      return;  // the pinned addAY zero skip
    }
    const bool identity = filter.size() == 1 && filter[0] == 1.;
    const std::vector<double> filtered =
        identity ? column : detail::pinned_backfilter_apply(filter, column);
    if (coefficient == 1.) {
      for (std::size_t i = 0; i < residuals.size(); ++i) {
        residuals[i] -= filtered[i];  // the exact x - y fold
      }
      return;
    }
    for (std::size_t i = 0; i < residuals.size(); ++i) {
      residuals[i] -= coefficient * filtered[i];  // pinned mul then sub
    }
  }

  // ModelController(ljung_box_limit).accept: the differenced linear-model
  // residuals (diff(y) - X b, the pinned differencing filter applied to
  // every column) drive the qualified slice-1 unit.
  ModelControllerResult accept_model() {
    const Estimation& estimation = estimation_required();
    RegressionLikelihood fit = landing_fit();
    if (fit.status != "available") {
      throw std::runtime_error("accept likelihood unavailable: " +
                               fit.status);
    }
    ensure_domain();
    // The LinearModel residuals on the differenced-model surface: yd minus
    // the Griffies design fold (identity filter - the columns enter already
    // filtered; the mean term is b * 1, never annihilated).
    const std::vector<double> identity{1.};
    std::vector<double> residuals = yd_;
    const std::vector<double>& coefficients = fit.coefficients;  // allCoefficients
    std::size_t col = 0;
    if (model_.desc.mean) {
      fold_column(residuals, identity, std::vector<double>(yd_.size(), 1.),
                  coefficients[col++]);
    }
    for (const OutlierVariable& variable : model_.desc.variables) {
      fold_column(residuals, identity, diff_column(outlier_column(variable)),
                  coefficients[col++]);
    }

    ModelControllerRequest controller;
    controller.residuals = std::move(residuals);
    controller.phi = estimation.coefficients[0];
    controller.theta = estimation.coefficients[1];
    controller.bphi = estimation.coefficients[2];
    controller.btheta = estimation.coefficients[3];
    controller.period = kAnnualFrequency;
    controller.effective_observations = static_cast<int>(yd_.size());
    controller.annual_frequency = kAnnualFrequency;
    controller.free_parameters = arma_dim_of(estimation.coefficients);
    controller.lb_limit = request_.ljung_box_limit;
    return model_controller_accept(controller);
  }

  void accept_refresh(const char* tag, int round, const char* note) {
    ModelControllerResult controller = accept_model();
    plbox0_ = 1.0 - controller.lb_p;
    rvr0_ = controller.rvr;
    rtval0_ = controller.rtval;
    last_accept_ = controller;
    AutoDecisionRecord entry;
    entry.stage = tag;
    entry.round = round;
    entry.scalars = {{"plbox0", plbox0_},
                     {"rvr0", rvr0_},
                     {"rtval0", rtval0_},
                     {"accepted", controller.accepted ? 1. : 0.}};
    if (note != nullptr) {
      entry.notes = {note};
    }
    entry.accept = std::move(controller);
    push(std::move(entry));
  }

  // ---- the RegressionVariablesTest machine ---------------------------------------

  bool run_regtest(double mean_cutoff, const char* tag, int round) {
    RegressionVariableTestRequest test;
    test.statistics = landing_statistics();
    int position = 0;
    for (std::size_t i = 0; i < model_.desc.variables.size(); ++i) {
      // ami outliers are never tested by the live machine (kind other); the
      // layout keeps the description order for the audit fields.
      test.variables.push_back(
          RegressionVariableLayout{RegressionVariableKind::other, position++, 1});
    }
    test.mean_test = mean_cutoff;
    test.td_test = 0.;   // decoder: no tdTest inside the envelope
    test.mh_test = 0.;   // moving holidays only with a live easter (21-I)
    test.derived_test = false;
    test.td_joint = 0.;  // the dead joint branch, unit-rejects nonzero
    test.mean = model_.desc.mean;
    test.adjusted = false;  // log && lp != None: lp is None in scope
    test.preadjustment = 0;
    RegressionVariableTestResult outcome = run_regression_variable_tests(test);
    const bool changed = outcome.changed;
    if (changed) {
      model_.desc.mean = outcome.mean;
      clear_estimation();
    }
    AutoDecisionRecord entry;
    entry.stage = tag;
    entry.round = round;
    entry.scalars = {{"mean_cutoff", mean_cutoff},
                     {"changed", changed ? 1. : 0.},
                     {"mean_after", outcome.mean ? 1. : 0.}};
    entry.regtest = std::move(outcome);
    push(std::move(entry));
    return changed;
  }

  // ---- searchur / checkUnitRoots / checkMA / redoEstimation -----------------------

  struct SearchurResult {
    int count = 0;
    double margin = std::numeric_limits<double>::quiet_NaN();
  };

  // The pinned searchur count over the AR mirror roots: the JVM's
  // Random-seeded Polynomial.roots replaced by the deterministic native
  // port; the margin feeds the slice-3 boundary-removal rule.
  SearchurResult searchur(const std::vector<double>& mirror) {
    SearchurResult result;
    if (mirror.size() < 2) {
      return result;  // a constant polynomial has no zeros
    }
    for (const ComplexRoot& root : polynomial_roots(mirror)) {
      const double im = std::abs(root.imag);
      const double modulus =
          std::sqrt(root.real * root.real + root.imag * root.imag);
      const double distance =
          std::min(std::abs(im - kImWindow),
                   std::abs(modulus - request_.ur_limit));
      if (std::isnan(result.margin) || distance < result.margin) {
        result.margin = distance;
      }
      if (jvm_gt(modulus, request_.ur_limit) && im <= kImWindow &&
          root.real > 0.) {
        ++result.count;
      }
    }
    return result;
  }

  bool check_unit_roots(int round) {
    const Estimation& estimation = estimation_required();
    AutoUnitRootsRecord record;
    Description& spec = model_.desc;
    bool fired = false;
    if (spec.p > 0 && spec.d < 2) {
      const auto result = searchur(
          mirror_ascending(jvm_ar_ascending(estimation.coefficients[0])));
      record.regular_count = result.count;
      record.regular_margin = result.margin;
      if (result.count != 0) {
        record.regular_fired = true;
        --spec.p;
        ++spec.d;
        fired = true;
      }
    }
    if (spec.bp > 0 && spec.bd < 1) {
      const auto result = searchur(mirror_ascending(
          jvm_seasonal_ar_ascending(estimation.coefficients[2],
                                    kAnnualFrequency)));
      record.seasonal_count = result.count;
      record.seasonal_margin = result.margin;
      if (result.count != 0) {
        record.seasonal_fired = true;
        --spec.bp;
        ++spec.bd;
        fired = true;
      }
    }
    record.p_after = spec.p;
    record.d_after = spec.d;
    record.bp_after = spec.bp;
    record.bd_after = spec.bd;
    if (fired) {
      record.redo_estimation = true;
      redo_estimation(round);
    }
    AutoDecisionRecord entry;
    entry.stage = "round" + std::to_string(round) + ".unitroots";
    entry.round = round;
    entry.scalars = {{"regular_count", record.regular_count},
                     {"seasonal_count", record.seasonal_count},
                     {"regular_fired", record.regular_fired ? 1. : 0.},
                     {"seasonal_fired", record.seasonal_fired ? 1. : 0.},
                     {"regular_margin", record.regular_margin},
                     {"seasonal_margin", record.seasonal_margin},
                     {"p_after", record.p_after},
                     {"d_after", record.d_after},
                     {"bp_after", record.bp_after},
                     {"bd_after", record.bd_after}};
    entry.unitroots = record;
    push(std::move(entry));
    return !fired;
  }

  bool check_ma(int round) {
    Description& spec = model_.desc;
    AutoMaTestRecord record;
    if (spec.q == 0 || spec.d == 0) {
      record.q_after = spec.q;
      record.d_after = spec.d;
      AutoDecisionRecord entry;
      entry.stage = "round" + std::to_string(round) + ".matest";
      entry.round = round;
      entry.scalars = {{"scanned", 0.},
                       {"q_after", record.q_after},
                       {"d_after", record.d_after}};
      entry.matest = record;
      push(std::move(entry));
      return true;
    }
    record.scanned = true;
    const double value =
        ma_evaluate_at_1(estimation_required().coefficients[1]);
    record.abs_theta_at_1 = std::abs(value);
    if (std::abs(value) < kMaSingular) {
      record.fired = true;
      --spec.q;
      --spec.d;
      model_.desc.mean = true;
      record.mean_set = true;
      record.redo_estimation = true;
      redo_estimation(round);
    }
    record.q_after = spec.q;
    record.d_after = spec.d;
    AutoDecisionRecord entry;
    entry.stage = "round" + std::to_string(round) + ".matest";
    entry.round = round;
    entry.scalars = {{"scanned", 1.},
                     {"abs_theta_at_1", record.abs_theta_at_1},
                     {"fired", record.fired ? 1. : 0.},
                     {"q_after", record.q_after},
                     {"d_after", record.d_after},
                     {"mean_set", record.mean_set ? 1. : 0.}};
    entry.matest = record;
    push(std::move(entry));
    return !record.fired;
  }

  void redo_estimation(int round) {
    stage_ = "estimation";
    estimate_row1();  // the pinned unconditional m.estimate(precision)
    if (model_.desc.mean) {
      stage_ = "meantest";
      check_mu(kCheckMuRedo,
               ("round" + std::to_string(round) + ".meantest").c_str(), round,
               "redo");
    }
    if (!model_.desc.variables.empty()) {
      model_.desc.variables.clear();
      clear_estimation();
    }
    stage_ = "estimation";
    estimate_if_needed();
    if (request_.outliers_used) {
      stage_ = "outliers";
      process_outliers(curva_, round,
                       ("round" + std::to_string(round) + ".outliers").c_str(),
                       "redo");
      stage_ = "estimation";
      estimate_if_needed();
    }
    stage_ = "accept";
    accept_refresh(("round" + std::to_string(round) + ".accept").c_str(),
                   round, "redo");
    reference_ = std::make_unique<Modelling>(model_);
  }

  // ---- pass2 (the seven-filter reference-restore gate) --------------------------

  bool pass2(bool stable_in, int round) {
    if (!reference_) {
      throw std::runtime_error("pass2 without a reference");
    }
    stage_ = "accept";
    ModelControllerResult controller = accept_model();
    const double rtval = controller.rtval;
    const double rvr = controller.rvr;
    double plbox = 1.0 - controller.lb_p;
    stage_ = "pass2";

    const Description& ref_desc = reference_->desc;
    const int n_out = static_cast<int>(model_.desc.variables.size());
    const int n0 = static_cast<int>(ref_desc.variables.size());
    const bool mean0 = ref_desc.mean;
    const bool mean = model_.desc.mean;
    const Description spec = model_.desc;

    int filter = 0;
    const bool same_model =
        ref_desc.p == spec.p && ref_desc.q == spec.q &&
        ref_desc.bp == spec.bp && ref_desc.bq == spec.bq &&
        ref_desc.d == spec.d && ref_desc.bd == spec.bd && mean0 == mean;
    if (!(n0 > n_out) && !same_model) {
      if (jvm_lt(plbox, request_.ljung_box_limit) &&
          jvm_lt(plbox, kPlboxMid) && jvm_lt(rvr0_, rvr)) {
        filter = 1;
      } else if (loop_ == 1 &&
                 jvm_ge(plbox, request_.ljung_box_limit) &&
                 jvm_lt(plbox0_, request_.ljung_box_limit)) {
        filter = 2;
      } else if (jvm_lt(plbox, request_.ljung_box_limit) &&
                 jvm_lt(plbox0_, kPlboxMid) && jvm_lt(plbox0_, plbox) &&
                 jvm_lt(rvr0_, kFct * rvr)) {
        filter = 3;
      } else if (jvm_ge(plbox, request_.ljung_box_limit) &&
                 jvm_lt(plbox0_, request_.ljung_box_limit) &&
                 jvm_lt(rvr0_, 1.0 * rvr)) {
        filter = 4;
      } else if (spec.d == 0 && spec.bd == 1 && spec.p == 1 &&
                 spec.bp == 0 && spec.q == 1 && spec.bq == 1 &&
                 model_.has_coefficients &&
                 jvm_lt(-model_.last_coefficients[0][0], kFilter5Phi)) {
        filter = 5;
      } else if (spec.d == 1 && spec.bd == 0 && spec.p == 0 &&
                 spec.bp == 1 && spec.q == 1 && spec.bq == 1 &&
                 model_.has_coefficients &&
                 jvm_lt(-model_.last_coefficients[2][0], kFilter6Bphi)) {
        filter = 6;
      }
    }

    AutoPass2Record record;
    record.filter = filter;
    bool stable = stable_in;
    if (filter > 0) {
      record.restored = true;
      model_ = *reference_;
      plbox = plbox0_;
      stable = true;
    } else {
      rtval0_ = rtval;
      rvr0_ = rvr;
      plbox0_ = plbox;
      reference_ = std::make_unique<Modelling>(model_);
    }
    pcr_ += loop_ == 1 ? kPcrLoop1 : kPcrOther;
    ++loop_;
    record.plbox = plbox;
    record.plbox0 = plbox0_;
    record.rvr = rvr;
    record.rtval = rtval;
    record.pcr_after = pcr_;
    record.loop_after = loop_;

    if (plbox <= pcr_) {
      // The dcmpg/ifgt reading (bytecode 566: dcmpg, ifgt 572): plbox > pcr
      // takes the false path and a NaN plbox/pcr takes it too (dcmpg -> +1);
      // only a non-NaN plbox <= pcr returns true (calc keeps the round body).
      record.continue_round = true;
      record.need_outliers_after = need_outliers_;
      record.need_auto_modelling_after = need_auto_modelling_;
      record.curva_after = curva_;
      attach_pass2(record, round);
      return true;
    }
    bool var21 = false;
    if (loop_ == 2 && request_.outliers_used) {
      if (jvm_gt(curva_, kMinCurva)) {
        // Bytecode 590: dcmpl(curva, 2.8) ifle skips, so a NaN curva skips
        // the reduction; Math.max then propagates a NaN reduceVa (java_max).
        curva_ = java_max(kMinCurva, curva_ * (1. - request_.reduce_va));
        var21 = true;
        record.curva_reduced = true;
      }
      need_outliers_ = true;  // unconditional in this branch (bytecode 640)
    }
    need_auto_modelling_ = !stable;
    if (loop_ > 2 || !var21) {
      record.last_solution = true;
      last_solution();
      if (loop_ == 2) {
        loop_ = 3;
      }
    }
    record.need_outliers_after = need_outliers_;
    record.need_auto_modelling_after = need_auto_modelling_;
    record.curva_after = curva_;
    attach_pass2(record, round);
    return false;
  }

  void attach_pass2(const AutoPass2Record& record, int round) {
    AutoDecisionRecord entry;
    entry.stage = "round" + std::to_string(round) + ".pass2";
    entry.round = round;
    entry.scalars = {{"filter", record.filter},
                     {"plbox", record.plbox},
                     {"plbox0", record.plbox0},
                     {"rvr", record.rvr},
                     {"rtval", record.rtval},
                     {"pcr_after", record.pcr_after},
                     {"loop_after", record.loop_after},
                     {"curva_after", record.curva_after},
                     {"continue_round", record.continue_round ? 1. : 0.}};
    if (record.last_solution) {
      entry.notes = {"last chance model"};
    }
    if (record.curva_reduced) {
      entry.notes.push_back("curva reduction");
    }
    entry.pass2 = record;
    push(std::move(entry));
  }

  void last_solution() {
    Description& spec = model_.desc;
    spec.p = kAnnualFrequency == 2 ? 1 : (kAnnualFrequency == 3 ? 2 : 3);
    if (spec.bd > 0 || kAnnualFrequency == 1) {
      spec.bp = 0;  // the pinned bytecode condition
    }
    spec.q = request_.mixed ? 1 : 0;
    if (kAnnualFrequency > 1) {
      spec.bq = 1;
    }
    stage_ = "estimation";
    clear_estimation();  // the spec changed
    estimate_row1();     // the pinned m.estimate(options.precision)
    if (request_.outliers_used) {
      curva_ = va0_;
      need_outliers_ = true;
    }
    need_auto_modelling_ = false;
  }

  // ---- the round body --------------------------------------------------------------

  // One calc loop round; returns true iff the cascade ended the loop.
  bool calc_round(int round) {
    const std::string prefix = "round" + std::to_string(round) + ".";
    bool stable = false;
    bool refine_failed_no_estimation = false;

    if (need_auto_modelling_) {
      stage_ = "auto modelling";
      if (!model_.estimation) {
        // The pinned module entry guard (needEstimation -> Failed): the
        // machine never runs, stable stays false and the removal branch
        // below mirrors the JVM's Failed handling.
        refine_failed_no_estimation = true;
        AutoDecisionRecord entry;
        entry.stage = prefix + "refine";
        entry.round = round;
        entry.notes = {"failed:no estimation"};
        push(std::move(entry));
      } else {
        ArmaRefinementRequest refinement;
        refinement.series = y_;
        refinement.mean = model_.desc.mean;
        // linearizedData removes the non-mean design columns with the
        // concentrated-likelihood coefficients; the 21-F unit consumes
        // exactly those inputs (regressors = the outlier columns).
        for (const OutlierVariable& variable : model_.desc.variables) {
          refinement.regressors.push_back(outlier_column(variable));
        }
        {
          RegressionLikelihood fit = landing_fit();
          if (fit.status != "available") {
            throw std::runtime_error("refine likelihood unavailable: " +
                                     fit.status);
          }
          refinement.regression_coefficients = fit.coefficients;
        }
        refinement.p = model_.desc.p;
        refinement.q = model_.desc.q;
        refinement.bp = model_.desc.bp;
        refinement.bq = model_.desc.bq;
        refinement.d = model_.desc.d;
        refinement.seasonal_d = model_.desc.bd;
        refinement.period = kAnnualFrequency;
        refinement.max_d = request_.max_d;
        refinement.max_seasonal_d = request_.max_bd;
        refinement.ub1 = request_.ub1;
        refinement.ub2 = request_.ub2;
        refinement.cancel = request_.cancel;
        refinement.precision = request_.differencing_precision;
        refinement.models_count = request_.models_count;
        refinement.max_p = request_.max_p;
        refinement.max_q = request_.max_q;
        refinement.max_bp = request_.max_bp;
        refinement.max_bq = request_.max_bq;
        refinement.balanced = request_.balanced;
        refinement.mixed = request_.mixed;
        refinement.accept_white_noise = request_.arma_accept_white_noise;
        refinement.eps = request_.arma_eps;

        ArmaRefinement outcome = refine_arma_model(refinement);
        stable = !outcome.changed;
        AutoDecisionRecord entry;
        entry.stage = prefix + "refine";
        entry.round = round;
        entry.scalars = {{"p", outcome.p},
                         {"q", outcome.q},
                         {"bp", outcome.bp},
                         {"bq", outcome.bq},
                         {"d", outcome.d},
                         {"bd", outcome.seasonal_d},
                         {"mean", outcome.mean ? 1. : 0.},
                         {"changed", outcome.changed ? 1. : 0.},
                         {"fallback", outcome.fallback ? 1. : 0.}};
        if (outcome.fallback) {
          entry.notes = {"fallback:" + outcome.fallback_stage};
        }
        entry.refine = std::move(outcome);
        push(std::move(entry));
        const ArmaRefinement& applied = *result_.chain.back().refine;
        model_.desc.p = applied.p;
        model_.desc.q = applied.q;
        model_.desc.bp = applied.bp;
        model_.desc.bq = applied.bq;
        model_.desc.d = applied.d;
        model_.desc.bd = applied.seasonal_d;
        model_.desc.mean = applied.mean;
        if (!stable) {
          // AutoModellingModule.process's Changed path: the
          // specification/mean change drives clearEstimation (contract
          // 21-F: "specification_changed || mean_changed drives
          // clearEstimation + Changed") BEFORE calc's removal branch - the
          // stale candidate (old orders) must never survive into the
          // estimate_if_needed below.
          clear_estimation();
        }
      }
      if (!stable) {
        // removeVariable(isOutlier(v, true)): all ami outliers.
        if (!model_.desc.variables.empty()) {
          model_.desc.variables.clear();
          clear_estimation();
          if (!result_.chain.empty() && refine_failed_no_estimation == false) {
            result_.chain.back().notes.push_back("outliers removed");
          }
        }
        need_outliers_ = request_.outliers_used;
      }
      stage_ = "estimation";
      estimate_if_needed();
      // regAIC (stable && loop <= 1): fail-closed no-op (21-I).
    }

    if (need_outliers_) {
      stage_ = "outliers";
      process_outliers(curva_, round,
                       (prefix + "outliers").c_str(), nullptr);
      stage_ = "estimation";
      estimate_if_needed();
    }

    if (request_.outliers_used && loop_ <= 2) {
      stage_ = "pass2";
      if (!pass2(stable, round)) {
        return false;  // calc TAIL (round++)
      }
    } else {
      stage_ = "accept";
      accept_refresh((prefix + "accept").c_str(), round, nullptr);
    }

    stage_ = "regtest";
    {
      const bool regtest_changed =
          run_regtest(request_.mean_test_final,
                      (prefix + "regtest").c_str(), round);
      if (regtest_changed) {
        if (loop_ < 3) {
          loop_ = 3;
        }
        stage_ = "estimation";
        estimate_if_needed();
        stage_ = "accept";
        accept_refresh((prefix + "accept").c_str(), round, "regtest");
      }
    }

    bool ok = false;
    while (!ok) {
      stage_ = "unit roots";
      const bool unit_ok = check_unit_roots(round);
      stage_ = "ma test";
      const bool ma_ok = check_ma(round);
      ok = unit_ok && ma_ok;
    }

    stage_ = "mean rule";
    {
      const bool adds = final_mean_rule_adds_mean(true, model_.desc.mean,
                                                  request_.check_mu, rtval0_);
      AutoMeanRuleRecord record;
      record.rtval0 = rtval0_;
      record.mean_added = adds;
      if (adds) {
        model_.desc.mean = true;
        clear_estimation();
      }
      AutoDecisionRecord entry;
      entry.stage = prefix + "meanrule";
      entry.round = round;
      entry.scalars = {{"rtval0", rtval0_}, {"mean_added", adds ? 1. : 0.}};
      entry.meanrule = record;
      push(std::move(entry));
    }

    stage_ = "final estimation";
    {
      std::vector<AutoCascadeRoundRecord> records;
      const bool cascade_ok = run_cascade(records);
      AutoDecisionRecord entry;
      entry.stage = prefix + "cascade";
      entry.round = round;
      for (const AutoCascadeRoundRecord& record : records) {
        entry.scalars.push_back({"status", record.status});
      }
      entry.cascade = std::move(records);
      push(std::move(entry));
      if (cascade_ok) {
        return true;
      }
    }
    if (loop_ <= 2 && request_.outliers_used) {
      curva_ = java_max(kMinCurva, curva_ * (1. - request_.reduce_va));
      need_outliers_ = true;
      need_auto_modelling_ = true;
    }
    return false;
  }

  // ---- the FinalEstimator cascade ---------------------------------------------------

  SarimaOrders final_orders() const {
    SarimaOrders orders;
    orders.p = model_.desc.p;
    orders.bp = model_.desc.bp;
    orders.q = model_.desc.q;
    orders.bq = model_.desc.bq;
    return orders;
  }

  // The raw JVM coordinate vector [phi | bphi | theta | btheta] of the
  // retained candidate.
  std::vector<double> raw_point() const {
    const Estimation& estimation = estimation_required();
    std::vector<double> point;
    for (const double phi : estimation.coefficients[0]) {
      point.push_back(-phi);
    }
    for (const double bphi : estimation.coefficients[2]) {
      point.push_back(-bphi);
    }
    for (const double theta : estimation.coefficients[1]) {
      point.push_back(theta);
    }
    for (const double btheta : estimation.coefficients[3]) {
      point.push_back(btheta);
    }
    return point;
  }

  // Raw [phi | bphi | theta | btheta] -> Kalos blocks [ar, ma, sar, sma]
  // under the current orders (the AR sign flip of the storage conventions).
  std::array<std::vector<double>, 4> raw_blocks(
      const std::vector<double>& raw) const {
    std::array<std::vector<double>, 4> blocks;
    std::size_t at = 0;
    auto take = [&](int n) {
      std::vector<double> out;
      for (int i = 0; i < n && at < raw.size(); ++i, ++at) {
        out.push_back(raw[at]);
      }
      return out;
    };
    std::vector<double> phi = take(model_.desc.p);
    std::vector<double> bphi = take(model_.desc.bp);
    std::vector<double> theta = take(model_.desc.q);
    std::vector<double> btheta = take(model_.desc.bq);
    for (double& v : phi) v = -v;
    for (double& v : bphi) v = -v;
    blocks[0] = std::move(phi);
    blocks[1] = std::move(theta);
    blocks[2] = std::move(bphi);
    blocks[3] = std::move(btheta);
    return blocks;
  }

  bool run_cascade(std::vector<AutoCascadeRoundRecord>& records) {
    for (int cround = 0; cround < kMaxCascadeRounds; ++cround) {
      AutoCascadeRoundRecord record;
      record.round = cround;
      record.fi_status = "none";
      record.p_before = model_.desc.p;
      record.bp_before = model_.desc.bp;
      record.q_before = model_.desc.q;
      record.bq_before = model_.desc.bq;
      record.p_after = record.p_before;
      record.bp_after = record.bp_before;
      record.q_after = record.q_before;
      record.bq_after = record.bq_before;
      try {
        // Row 4: the FinalEstimator computer (RegSarimaComputer + LM
        // minimizer + Multiple start + exact final derivatives). Kalos:
        // since CI-FDLIBM the same PinnedLevenbergMarquardt machine as row 1
        // (base_request), fast_derivatives = false; the recorded solver
        // delta against the JVM's LM minimizer is retired by the ruling.
        stage_ = "final estimation";
        auto request = base_request(request_.precision,
                                    EstimationStartingPoint::Multiple, false);
        SarimaEstimationResult fit_result = estimate_sarima_fixed(request);
        Estimation estimation;
        estimation.coefficients = fit_result.coefficients;
        estimation.status = fit_result.status;
        estimation.log_likelihood = fit_result.log_likelihood;
        model_.last_coefficients = estimation.coefficients;
        model_.has_coefficients = true;
        model_.estimation = std::move(estimation);
        estimation_has_mean_ = model_.desc.mean;
        record.estimation_status = std::move(fit_result.status);
        record.estimation_ok = true;

        // FinalEstimator.estimate's dim is desc.mapping().getDim() and
        // ModelDescription.mapping() (no fixed parameters, freeArimaParameters
        // already ran) is SarimaMapping2.of(orders): the dim is the ARMA
        // parameter count (the design never enters it - it is the X matrix),
        // so the dim == 0 early-true fires on a parameter-less ARMA, not on
        // an empty regression design.
        const int arma_parameters = model_.desc.p + model_.desc.bp +
                                    model_.desc.q + model_.desc.bq;
        if (arma_parameters == 0) {
          record.status = -2;  // the pinned dim == 0 early-true
          record.cascade_ok = true;
          records.push_back(std::move(record));
          return true;
        }
        // ami == (autoModel != null) == true for this kind: the test runs.

        stage_ = "final information";
        final_information_request information;
        information.orders = final_orders();
        information.point = raw_point();
        ensure_domain();
        const std::vector<double>& y_d = yd_;
        // RegSarimaComputer's cdf correction (RegSarimaComputer bytecode
        // 58-96, builder default cdf = true): the nobs fed to the
        // information/gradient scalings is y.length() - X.getColumnsCount()
        // - mapping.getDim() on the differenced-model surface.
        const int xcols = static_cast<int>(build_design().columns);
        const int arma_dim = arma_parameters;
        information.nobs = static_cast<int>(y_d.size()) - xcols - arma_dim;
        record.nobs = information.nobs;
        if (xcols == 0) {
          // The bit-pinned deviances domain (the design-free machine).
          const SarimaOrders orders = information.orders;
          const std::vector<double> differenced = y_d;
          information.deviance =
              [orders,
               differenced](const std::vector<double>& x) -> std::optional<double> {
            const auto point = detail::pinned_deviances_evaluate(
                differenced, static_cast<std::size_t>(orders.p),
                static_cast<std::size_t>(orders.bp),
                static_cast<std::size_t>(orders.q),
                static_cast<std::size_t>(orders.bq), kAnnualFrequency, {}, x);
            if (!point) {
              return std::nullopt;
            }
            return point->ssq_e;
          };
        } else {
          // The concentrated deviance mapping at the candidate ARMA on the
          // differenced-model surface (contract delta: the CompoundMapping
          // jacobian block lands as measured-tier data, the G2 recipe, never
          // bit-claimed).
          const RegressionDesign design = build_design();
          const MaskedSeries series = make_differenced_series();
          const Description desc = model_.desc;
          information.deviance =
              [this, design, series,
               desc](const std::vector<double>& raw) -> std::optional<double> {
            const Description restore_desc = model_.desc;
            model_.desc = desc;  // raw_blocks reads the order sizes
            const auto blocks = raw_blocks(raw);
            model_.desc = restore_desc;
            RegressionLikelihood fit = sarima_regression_likelihood(
                series, design, blocks[0], blocks[1], blocks[2], blocks[3],
                kAnnualFrequency, 0, 0);
            if (fit.status != "available" ||
                !std::isfinite(fit.sum_squared_standardized) ||
                fit.finite_count == 0 ||
                !std::isfinite(fit.sum_log_finite_variances)) {
              return std::nullopt;
            }
            return fit.sum_squared_standardized *
                   detail::fdlibm_exp(fit.sum_log_finite_variances /
                            static_cast<double>(fit.finite_count));
          };
        }
        final_information_result information_result =
            final_information(information);
        record.fx = information_result.fx;
        switch (information_result.status) {
          case final_information_status::available:
            record.fi_status = "available";
            break;
          case final_information_status::inverse_failed:
            record.fi_status = "inverse_failed";
            break;
          case final_information_status::final_point_rejected:
            record.fi_status = "final_point_rejected";
            break;
        }
        if (information_result.status != final_information_status::available) {
          // The JVM throw analogue (no LLPoint): the pinned catch maps the
          // cascade round onto failure.
          record.status = -1;
          record.cascade_ok = false;
          records.push_back(std::move(record));
          break;
        }
        record.asy_diagonal = information_result.asymptotic_covariance_diagonal;

        stage_ = "cascade test";
        FinalEstimatorRound test = final_estimator_test(
            information.orders, information.point, record.asy_diagonal,
            static_cast<int>(y_.size()), request_.tsig);
        record.status = test.status;
        record.p_after = test.orders.p;
        record.bp_after = test.orders.bp;
        record.q_after = test.orders.q;
        record.bq_after = test.orders.bq;
        records.push_back(std::move(record));
        if (test.status == 0) {
          records.back().cascade_ok = true;
          return true;
        }
        if (test.status > 1) {
          break;  // the pinned return false (round failure)
        }
        // status == 1: the reduced specification re-estimates next round
        // (setSpecification; the retained candidate stays, the next row-4
        // estimation overwrites it exactly as the JVM loop re-estimates).
        model_.desc.p = test.orders.p;
        model_.desc.q = test.orders.q;
        model_.desc.bp = test.orders.bp;
        model_.desc.bq = test.orders.bq;
      } catch (const std::exception&) {
        // The pinned catch(RuntimeException): "estimation faied", the round
        // fails; the modelling keeps whatever the partial run retained.
        record.status = -1;
        record.cascade_ok = false;
        records.push_back(std::move(record));
        break;
      }
    }
    return false;  // the five inconclusive rounds
  }

  // ---- the final description + downstream composition ------------------------------

  void build_result() {
    X13AutomaticFitResult& out = result_;
    out.log = model_.desc.log;
    out.mean = model_.desc.mean;
    out.p = model_.desc.p;
    out.q = model_.desc.q;
    out.bp = model_.desc.bp;
    out.bq = model_.desc.bq;
    out.d = model_.desc.d;
    out.bd = model_.desc.bd;

    const Estimation& estimation = estimation_required();
    // Export convention (pinned by the oracle's arma facts, 21-H.3): the
    // JVM surface is SarimaSpec.parameters(), where the AR and SAR blocks
    // carry the all-plus polynomial signs (a = -phi); MA/SMA already match
    // the driver's theta convention. The internal estimation keeps phi.
    out.coefficients = estimation.coefficients;
    for (const int block : {0, 2}) {
      for (double& v : out.coefficients[static_cast<std::size_t>(block)]) {
        v = -v;
      }
    }
    out.estimation_status = estimation.status;
    out.log_likelihood = estimation.log_likelihood;
    out.effective_observations = static_cast<std::size_t>(
        static_cast<int>(y_.size()) -
        (model_.desc.d + model_.desc.bd * kAnnualFrequency));

    // The landing composition at the final point (the G2 recipe; never
    // bit-claimed): the outlier coefficients and the AICC echo.
    RegressionLikelihood fit = landing_fit();
    if (fit.status != "available") {
      throw std::runtime_error("final likelihood unavailable: " + fit.status);
    }
    std::size_t col = 0;
    if (model_.desc.mean) {
      ++col;
    }
    for (const OutlierVariable& variable : model_.desc.variables) {
      AutoVariableRecord record;
      record.name = variable.name;
      record.type = variable.type;
      record.regeffect = component_name(variable.effect);
      record.position = variable.position;
      record.coefficient =
          col < fit.coefficients.size() ? fit.coefficients[col] : 0.;
      record.column = outlier_column(variable);
      out.variables.push_back(std::move(record));
      ++col;
    }
    {
      const double neff = static_cast<double>(out.effective_observations);
      const double parameters =
          static_cast<double>(arma_dim_of(estimation.coefficients) +
                              fit.coefficients.size() + 1);
      // AICC rides the ADJUSTED log-likelihood (LikelihoodStatistics.calc
      // reads adjustedLogLikelihood; the Builder adds the transformation
      // Jacobian): in log mode the JVM subtracts sum(log y_t) over the
      // estimation domain - the last n - d - 12*bd original values (pinned
      // by the oracle stats facts: ao -333.0415626 with d=0, bd=1; rw
      // -450.5425654 with d=1, bd=0; level cases carry no adjustment).
      double adjusted_log_likelihood = out.log_likelihood;
      if (out.log) {
        const std::size_t skip = static_cast<std::size_t>(
            model_.desc.d + model_.desc.bd * kAnnualFrequency);
        for (std::size_t i = skip; i < request_.input.size(); ++i) {
          adjusted_log_likelihood -= detail::fdlibm_log(request_.input[i]);
        }
      }
      out.aicc = -2. *
                 (adjusted_log_likelihood -
                  neff * parameters / (neff - parameters - 1.));
    }

    // The downstream composition in the pinned X13Kernel order (the fixed
    // fit's seam, decoder-shaped X-11 block: Msr seasonal requests and the
    // automatic Henderson length pass straight to the qualified unit).
    X13FitResult& fit_result = out.fit;
    fit_result.frequency = kAnnualFrequency;
    fit_result.first_period = request_.first_period;
    fit_result.log = out.log;
    fit_result.seasonal = request_.seasonal;
    fit_result.trend_length = request_.trend_length;
    fit_result.sigma_set = request_.sigma.has_value();
    if (request_.sigma) {
      fit_result.sigma = *request_.sigma;
    }
    fit_result.lower_sigma = request_.lower_sigma;
    fit_result.upper_sigma = request_.upper_sigma;
    fit_result.msr_policy = request_.msr_policy;
    fit_result.mode_requested = request_.mode;
    fit_result.mode_effective = forceMode(out.log, request_.mode);
    // The pinned process semantics (item 21-H.4): the spec horizons are
    // month counts; a negative value means |h| annual frequencies.
    const int nb = request_.backcast < 0
                       ? -request_.backcast * kAnnualFrequency
                       : request_.backcast;
    const int nf = request_.forecast < 0
                       ? -request_.forecast * kAnnualFrequency
                       : request_.forecast;
    fit_result.backcast = nb;
    fit_result.forecast = nf;

    X13PreadjustInput preadjust;
    preadjust.input = request_.input;
    preadjust.log = out.log;
    preadjust.backcast = nb;
    preadjust.forecast = nf;
    if (nb != 0 || nf != 0) {
      // The landing stochastic component on the JVM SarimaModel surface
      // (the AR blocks carry the all-plus polynomial signs, the same export
      // convention as the oracle's arma facts) plus the mean-correction
      // flag, the exact inputs of the pinned ExactArimaForecasts route.
      ArimaForecastModel forecasts_model;
      forecasts_model.period = kAnnualFrequency;
      forecasts_model.p = out.p;
      forecasts_model.q = out.q;
      forecasts_model.bp = out.bp;
      forecasts_model.bq = out.bq;
      forecasts_model.d = out.d;
      forecasts_model.bd = out.bd;
      forecasts_model.mean_correction = out.mean;
      // The estimation block layout is [ar, ma, sar, sma] (the export flips
      // the two AR blocks {0,2}); the JVM SarimaModel parameter surface is
      // [ar, sar, ma, sma] with the all-plus polynomial signs on the AR
      // blocks (SarimaModel.parameters(), probe-verified in the pinning
      // notes section E).
      const std::array<std::vector<double>, 4>& blocks = estimation.coefficients;
      forecasts_model.phi.assign(blocks[0].size(), 0.);
      forecasts_model.bphi.assign(blocks[2].size(), 0.);
      for (std::size_t i = 0; i < blocks[0].size(); ++i)
        forecasts_model.phi[i] = -blocks[0][i];
      for (std::size_t i = 0; i < blocks[2].size(); ++i)
        forecasts_model.bphi[i] = -blocks[2][i];
      forecasts_model.theta = blocks[1];
      forecasts_model.btheta = blocks[3];
      preadjust.arima = forecasts_model;
    }
    for (const AutoVariableRecord& variable : out.variables) {
      X13AllocationVariable allocation;
      allocation.cls = X13VariableClass::outlier;
      allocation.effect = outlier_effect(variable.type);
      allocation.coefficient = {variable.coefficient};
      // The pinned factory fill over the extended calendar: the zero-ended
      // LS/SO patterns reach back across the position into the backcast
      // window, the TC ramp runs into the forecast one (the JVM reg-column
      // probe, pinning notes section J).
      std::vector<double> column(
          static_cast<std::size_t>(nb) + request_.input.size() +
              static_cast<std::size_t>(nf),
          0.);
      fill_outlier_column(variable.type,
          static_cast<std::size_t>(nb) + variable.position, kAnnualFrequency,
          transitory_rate_, column);
      allocation.columns = {std::move(column)};
      preadjust.variables.push_back(std::move(allocation));
    }
    fit_result.preadjustment = x13_preadjust(preadjust);
    if (fit_result.preadjustment.status != "available") {
      fit_result.status = fit_result.preadjustment.status;  // verbatim
      out.status = fit_result.status;
      AutoDecisionRecord entry;
      entry.stage = "final";
      entry.notes = {"preadjust refused", fit_result.status};
      push(std::move(entry));
      return;
    }

    const int extended_first_period =
        ((request_.first_period - 1 - nb) % kAnnualFrequency +
            kAnnualFrequency) %
            kAnnualFrequency +
        1;
    MaskedSeries series;
    series.values = fit_result.preadjustment.out;
    series.missing.assign(fit_result.preadjustment.out.size(), 0);
    fit_result.selection = x11_selection_pass(
        series, kAnnualFrequency, extended_first_period, request_.seasonal,
        request_.trend_length, *request_.sigma, request_.msr_policy,
        request_.lower_sigma, request_.upper_sigma, nb, nf,
        fit_result.mode_effective);
    if (fit_result.selection.status != "available") {
      fit_result.status = fit_result.selection.status;  // verbatim
      out.status = fit_result.status;
      AutoDecisionRecord entry;
      entry.stage = "final";
      entry.notes = {"downstream refused", fit_result.status};
      push(std::move(entry));
      return;
    }
    if (fit_result.mode_effective == X11DecompositionMode::log_additive) {
      logAdditiveResults(fit_result.selection.d.d10.values,
                         fit_result.selection.d.d12.values,
                         fit_result.selection.d.d13.values,
                         fit_result.preadjustment.out, fit_result.d10,
                         fit_result.d11, fit_result.d12, fit_result.d13);
    } else {
      fit_result.d10 = fit_result.selection.d.d10.values;
      fit_result.d11 = fit_result.selection.d.d11.values;
      fit_result.d12 = fit_result.selection.d.d12.values;
      fit_result.d13 = fit_result.selection.d.d13.values;
    }
    fit_result.c17 = fit_result.selection.c.c17.values;
    X13RestoreInput restore;
    restore.pre = fit_result.preadjustment;
    restore.d10 = fit_result.d10;
    restore.d11 = fit_result.d11;
    restore.d12 = fit_result.d12;
    restore.d13 = fit_result.d13;
    restore.c17 = fit_result.c17;
    restore.mode = fit_result.mode_effective;
    restore.backcast = nb;
    restore.forecast = nf;
    fit_result.finals = x13_restore(restore);
    fit_result.finalSA = fit_result.finals.d11final;
    fit_result.finalTrend = fit_result.finals.d12final;
    fit_result.finalIrregular = fit_result.finals.d13final;
    fit_result.status = "available";

    AutoDecisionRecord entry;
    entry.stage = "final";
    entry.scalars = {{"p", out.p},
                     {"q", out.q},
                     {"bp", out.bp},
                     {"bq", out.bq},
                     {"d", out.d},
                     {"bd", out.bd},
                     {"mean", out.mean ? 1. : 0.},
                     {"log", out.log ? 1. : 0.},
                     {"n_outliers", static_cast<double>(out.variables.size())},
                     {"log_likelihood", out.log_likelihood},
                     {"aicc", out.aicc}};
    entry.notes = {"downstream composed", "available"};
    push(std::move(entry));
    out.status = "available";
  }

  // ---- state ------------------------------------------------------------------------

  const X13AutomaticRequest& request_;
  X13AutomaticFitResult& result_;
  std::vector<OutlierVariable> empty_variables_;

  Modelling model_;
  std::unique_ptr<Modelling> reference_;

  double va0_ = 0.;
  double curva_ = 0.;
  double pcr_ = 0.;
  int loop_ = 0;
  int round_ = 0;
  double plbox0_ = 0.;
  double rvr0_ = 0.;
  double rtval0_ = 0.;
  bool need_outliers_ = false;
  bool need_auto_modelling_ = false;
  std::optional<ModelControllerResult> last_accept_;
  std::string stage_ = "init";
};

}  // namespace

// ---------------------------------------------------------------------------
// The public entry point.
// ---------------------------------------------------------------------------

X13AutomaticFitResult x13_automatic_fit(const X13AutomaticRequest& request) {
  X13AutomaticFitResult out;

  // Shape validation (throws std::invalid_argument, as the seams below do).
  // The frequency refusal of the fail-closed matrix precedes the shape
  // checks that interpret the request as monthly (the contract matrix reads
  // "frequency != 12" as a status answer for EVERY non-monthly request,
  // before anything inspects its shape).
  if (request.input.empty()) {
    throw std::invalid_argument("x13 automatic fit requires a non-empty series");
  }
  if (request.frequency != 12) {
    out.status = "unsupported:x13 automatic fit monthly series only";
    return out;
  }
  if (request.first_period < 1 || request.first_period > 12) {
    throw std::invalid_argument(
        "x13 automatic fit first_period is a 1-based month index");
  }
  if (request.seasonal.size() != static_cast<std::size_t>(request.frequency)) {
    throw std::invalid_argument(
        "x13 automatic fit seasonal requests have one entry per calendar period");
  }
  if (request.input.size() <= 13) {
    throw std::invalid_argument(
        "x13 automatic fit series too short for the airline start");
  }
  if (request.outliers_used && request.outlier_types.empty()) {
    throw std::invalid_argument(
        "x13 automatic fit outlier detection needs at least one type");
  }

  // The fail-closed matrix, in the contract order (the monthly refusal
  // already ran above, as the first shape-independent domain gate).
  for (const double value : request.input) {
    if (!std::isfinite(value)) {
      out.status =
          "unsupported:x13 automatic fit requires a complete finite series";
      return out;
    }
  }
  if (!request.sigma) {
    out.status = "unsupported:x13 fit requires an explicit sigma policy";
    return out;
  }
  if (request.use_calendar_detection || request.use_easter_detection ||
      request.automatic_trading_days || request.use_user_variables) {
    out.status =
        "unsupported:x13 automatic fit defers calendar/easter/TD/user "
        "detection (21-I)";
    return out;
  }
  if (request.outlier_span &&
      (request.outlier_span->first != 0 ||
       request.outlier_span->second !=
           static_cast<int>(request.input.size()) - 1)) {
    out.status = "unsupported:x13 automatic fit spans the full series only";
    return out;
  }

  AutomaticPipeline pipeline(request, out);
  try {
    pipeline.run();
  } catch (const std::runtime_error& error) {
    // The pinned catch(Exception) -> return null: the failed: status with
    // the chain kept so far (the invalid_argument shape violations of the
    // seams below propagate to the caller, as documented).
    out.status = "failed:" + pipeline.stage();
    AutoDecisionRecord entry;
    entry.stage = "failed";
    entry.round = pipeline.round();
    entry.notes = {error.what()};
    out.chain.push_back(std::move(entry));
  }
  return out;
}

}  // namespace kalos
