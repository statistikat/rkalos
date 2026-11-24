#ifndef KALOS_X13_AUTOMATIC_HPP
#define KALOS_X13_AUTOMATIC_HPP
#include "x13_estimated_calendar.hpp"
// Item 21-H.1: the automatic-spec X-13 pipeline driver (stage Z). This is
// the pinned RegArimaKernel.calc transcription (docs/contracts/
// x13-automatic-pipeline-contract.md) composed on the qualified stage-Z
// units: the 21-B transformation decision, the 21-G mean controller /
// regression-variable tests / final mean rule, the 21-C outlier rounds with
// the isOutlier(v, true) forward-add vocabulary, the 21-F refinement
// wrapping the 21-D/21-E machines, the slice-1 ModelController accept
// machine, the seven-filter pass2 gate, checkUnitRoots/checkMA/
// redoEstimation/lastSolution, and the FinalEstimator cascade on top of
// final_estimator_test and the 21-H.0 final_information producer. The
// result reuses the restricted fit's X13FitResult verbatim as its
// downstream composition block (preadjust -> selection pass -> finals over
// the final description); the fixed path is untouched by construction -
// this kind is a distinct request/result pair behind its own entry point
// and x13_fit.hpp/cpp carry no change.
//
// Delivered scope (contract envelope): frequency 12, checked backcast and forecast horizons (21-H.4), the
// four decomposition modes through the pinned forcing, a complete finite
// input, outlier span = the full estimation domain, transform.function =
// Auto (the decoder builds the LogLevelModule only in that case).
// Calendar/easter-detection/automatic-TD/user-variable detection are
// deferred to 21-I and answered with the explicit unsupported status;
// SeasonalityDetector and the FRegressionTest joint branch are dead in the
// current jar and not wired; the estimation adapter uses the qualified
// call-site recipes (row 1: HannanRissanen + PinnedLevenbergMarquardt —
// the landing estimate runs on the transcribed JVM LM machine per the
// CI-FDLIBM ruling recorded in the pipeline contract; rows 3/4: the 21-C
// sweep recipe with the Multiple start; the PinnedLevenbergMarquardt solver
// stays inside the qualified 21-D/E/F units). Every module scalar arrives
// as request data - the core promotes nothing.
//
// Non-determinism envelope (slice-3 policy): Polynomial.roots is
// Random-seeded in the JVM for degree >= 3; the driver's root checks run
// the deterministic native ports and record, per checkUnitRoots fire, the
// distance of the closest root to the |im| = 0.05 / |z| = ur_limit
// boundaries so corpus cases landing within 0.01 of a boundary can be
// removed at the oracle tier; the final cascade stays within checkRoots
// degrees 0-2 (the final_estimator unit's envelope).
#include "arma_diagnostics.hpp"
#include "arma_refinement.hpp"
#include "estimation_control.hpp"
#include "final_estimator.hpp"
#include "outlier_scoring.hpp"
#include "regression_variable_tests.hpp"
#include "transformation_decision.hpp"
#include "x11_selection_pass.hpp"
#include "x13_fit.hpp"
#include <array>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>
namespace kalos {
// The fit-spec kind of docs/contracts/x13-automatic-pipeline-contract.md.
// The fixed kind IS kalos::x13_fit and keeps its exact request/result;
// the automatic kind is the X13AutomaticRequest / X13AutomaticFitResult
// pair below. The enum is informational (it lives on this side so the
// fixed header stays byte-identical - the schema-change note of the fit
// contract).
enum class X13FitSpecKind { fixed_spec, automatic_spec };
// One automatic-spec request: the decoder-equivalent module state, ALL
// required (the core applies no promoted default to a module scalar).
// The scalars mirror the X13SpecDecoder call sites of pinning note 12.1
// exactly: DifferencingModule(cancel, ub1, ub2, precision), ArmaModule
// (models_count, grid bounds, balanced, mixed, eps), the x13
// OutliersDetectionModule (type set, max round/outliers, tcrate), the
// AmiOptions (precision, va, reduceVa, ljungBoxLimit, checkMu, mixedModel,
// urLimit, acceptAirline), the RegressionVariablesTest builder cutoffs
// (mean_test_initial 1.96 / mean_test_final 1.0) and the FinalEstimator
// tsig (builder default 1.0). The X-11 block mirrors the decoder X11Spec
// (filters Msr allowed, Henderson 0 = automatic, REQUIRED sigma policy,
// sigma cuttings and horizons in months or negative annual counts).
struct X13AutomaticRequest {
    std::vector<double> input; // complete, finite, monthly
    int frequency=12;
    int first_period=1;        // 1-based month index of input[0]
    int first_year=1970;       // calendar year of input[0] (the IOutlier
                               // defaultName dates; the JVM period carries it)
    // ---- LogLevelModule (transformation.function == Auto) ----
    double aicc_log_correction;       // TransformSpec aicDiff (-2.0 in corpus)
    // ---- DifferencingModule (21-D state) ----
    int max_d;
    int max_bd;
    double ub1, ub2, cancel, differencing_precision;
    // ---- ArmaModule (21-E state) ----
    int models_count, max_p, max_q, max_bp, max_bq;
    bool balanced, mixed, arma_accept_white_noise;
    double arma_eps;
    // ---- OutliersDetectionModule (21-C state) ----
    bool outliers_used;
    std::vector<OutlierType> outlier_types;    // pinned factory order subset
    std::optional<std::pair<int,int>> outlier_span; // [first,last] positions;
                                                    // unset or the full series
    std::size_t outlier_max_round, outlier_max_outliers;
    double tcrate;                              // monthly TC rate
    // ---- AmiOptions ----
    double precision;               // kernel-round eps (1e-7)
    double mean_precision;          // MeanController processor eps (1e-5)
    double va;                      // 0 -> X13Utility.calcCv(n) rule
    double reduce_va, ljung_box_limit;
    bool check_mu, accept_airline;
    double ur_limit;
    // ---- RegressionVariablesTest / FinalEstimator cutoffs ----
    double mean_test_initial, mean_test_final;  // 1.96 / 1.0 (0 = no test)
    double tsig;                                // cascade ratioLimit (1.0)
    // ---- deferred-module flags (21-I): any true answers the explicit
    // unsupported status; the pinned drivers are not transcribed ----
    bool use_calendar_detection=false, use_easter_detection=false;
    bool automatic_trading_days=false, use_user_variables=false;
    // ---- the decoder X-11 block (G4 surface) ----
    std::vector<X11SeasonalRequest> seasonal;   // Msr allowed here
    int trend_length=0;                         // 0 = automatic Henderson
    std::optional<X11SigmaPolicy> sigma;        // required; no default
    double lower_sigma=1.5, upper_sigma=2.5;
    X11MsrIterationPolicy msr_policy=X11MsrIterationPolicy::jdplus;
    X11DecompositionMode mode=X11DecompositionMode::additive; // Undefined
    // The X-11 spec horizons (item 21-H.4): month counts, a negative value
    // meaning |h| annual frequencies, exactly as the pinned process reads
    // them; the initialStep extension route then drives a1b/a1a, the
    // extended out calendar and the split finals.
    int backcast=0, forecast=0;
};
// ---- decision-record vocabulary (the H.2 diff surface) ----
// One pass2 invocation: the fired filter (0 = keep), the state after the
// gate and the loop bookkeeping, all pinned by RegArimaKernel.pass2.
struct AutoPass2Record {
    int filter=0;                 // 0-6 (0 = no restore)
    bool restored=false;          // filter > 0 -> reference restore
    bool continue_round=false;    // return value: plbox <= pcr
    double plbox=0., plbox0=0., rvr=0., rtval=0.;
    double pcr_after=0.; int loop_after=0;
    double curva_after=0.; bool curva_reduced=false;
    bool need_outliers_after=false; bool need_auto_modelling_after=false;
    bool last_solution=false;     // the "last chance model" branch fired
};
// One checkUnitRoots invocation: searchur counts per block and the
// boundary-distance margins of the slice-3 non-determinism policy.
struct AutoUnitRootsRecord {
    int regular_count=0, seasonal_count=0;
    bool regular_fired=false, seasonal_fired=false;
    int p_after=0,d_after=0,bp_after=0,bd_after=0;
    // min over the inspected roots of the distance to the |im| = 0.05 and
    // |z| = ur_limit boundaries; NaN when the block was not scanned.
    double regular_margin=std::numeric_limits<double>::quiet_NaN();
    double seasonal_margin=std::numeric_limits<double>::quiet_NaN();
    bool redo_estimation=false;
};
// One checkMA invocation.
struct AutoMaTestRecord {
    bool scanned=false;           // q > 0 && d > 0 (the guard)
    double abs_theta_at_1=0.;     // |regular MA evaluated at 1|
    bool fired=false;             // < 1e-3: q--, d--, mean set
    int q_after=0, d_after=0; bool mean_set=false; bool redo_estimation=false;
};
// The final mean rule (rtval0 > 2.5 with the pinned dcmpl flavour).
struct AutoMeanRuleRecord {
    double rtval0=0.; bool mean_added=false;
};
// One FinalEstimator round: the row-4 estimation, the final-information
// producer inputs and the cascade status (one entry per cascade.rounds).
struct AutoCascadeRoundRecord {
    int round=0;                  // 0-based cascade round
    std::string estimation_status;// row-4 estimation status
    bool estimation_ok=false;     // the pinned null-estimation guard
    std::string fi_status;        // final_information_status name, or "none"
    double fx=0.; int nobs=0;
    int p_before=0,bp_before=0,q_before=0,bq_before=0;
    int p_after=0,bp_after=0,q_after=0,bq_after=0;
    int status=0;                 // final_estimator_test status (or <0:
                                  // -1 estimation throw, -2 dim == 0,
                                  // -3 non-ami, never for this kind)
    std::vector<double> asy_diagonal; // the consumed asy cov diagonal
    bool cascade_ok=false;        // estimate() returned true in this round
};
// One entry of the decision chain, in execution order. Every machine
// execution appends exactly one entry; entries are never reordered or
// merged (contract chain invariant). Exactly one payload mirrors the
// machine's qualified record for its stage tag; scalars/notes echo the
// pinned state transitions.
struct AutoDecisionRecord {
    // "init" | "transformation" | "meantest0" | "outliers0" | "regtest0" |
    // "accept0" | "round<k>.refine" | "round<k>.outliers" |
    // "round<k>.pass2" | "round<k>.regtest" | "round<k>.unitroots" |
    // "round<k>.matest" | "round<k>.meantest" | "round<k>.accept" |
    // "round<k>.meanrule" | "round<k>.cascade" | "final" (k is 1-based).
    // The redoEstimation sub-machines re-enter the mean controller, the
    // outlier machine and the accept refresh: those entries repeat the
    // base tag with round context and carry the note "redo" (vocabulary
    // extension recorded in the pipeline contract - the accept-refreshes
    // of the pass2-skip and regtest-Changed branches are plain
    // "round<k>.accept" entries).
    std::string stage;
    int round=0;                  // loop round (0 for the pre-round stages)
    // Pinned scalar echoes (curva/pcr/loop/plbox0/rvr0/rtval0/... per
    // stage; insertion order preserved).
    std::vector<std::pair<std::string,double>> scalars;
    std::vector<std::string> notes;
    // The consumed unit's qualified record, verbatim.
    std::optional<TransformationDecisionResult> transformation;
    std::optional<MeanControllerResult> mean;
    std::optional<OutlierScoringResult> outliers;
    std::optional<RegressionVariableTestResult> regtest;
    std::optional<ModelControllerResult> accept;
    std::optional<ArmaRefinement> refine;
    std::optional<AutoPass2Record> pass2;
    std::optional<AutoUnitRootsRecord> unitroots;
    std::optional<AutoMaTestRecord> matest;
    std::optional<AutoMeanRuleRecord> meanrule;
    std::vector<AutoCascadeRoundRecord> cascade; // one per cascade round
};
// One selected regression variable of the final description: the pinned
// name, class attributes and the estimated coefficient on its column.
struct AutoVariableRecord {
    std::string name;             // IOutlier.defaultName(code, period)
    OutlierType type;             // factory code (AO/LS/TC/SO)
    std::string regeffect;        // the pinned defaultComponentTypeOf name
    std::size_t position=0;       // estimation-domain position of the outlier
    double coefficient=0.;
    std::vector<double> column;   // the factory column over the model span
};
// The public RegArimaModel likelihood and innovations of the final fitted
// model: populated by the stock estimated templates and, since the M7
// cutover entry, by every available automatic fit.
struct X13StockLikelihood {
    double residual_ssq = 0.;
    std::size_t estimated_parameters = 0;
    std::vector<double> full_residuals; // starts after d + 12*bd observations
    int start_year = 0, start_month = 0; // input domain
    std::vector<double> regression_covariance; // row-major, in fitted model order
};

// Stock estimated RSA2c calendar variables precede the outlier variables in model order.
struct X13StockCalendarFit {
    X13Rsa2cCalendar selection;
    std::vector<std::string> names; // td, lp, easter, in fitted order
    std::vector<double> coefficients;
};

// The result of one automatic-spec or stock estimated fit.
struct X13AutomaticFitResult {
    // "available" | "unsupported:<...>" (fail-closed matrix) |
    // "failed:<stage>" (the pinned return-null chain outcomes, chain kept).
    std::string status;
    // The final description: orders (stationary part + differencing), the
    // log/mean flags and the landed ARMA coefficients (Kalos convention,
    // blocks [regular AR, regular MA, seasonal AR, seasonal MA]).
    bool log=false, mean=false;
    int p=0, q=0, bp=0, bq=0, d=0, bd=0;
    std::array<std::vector<double>, 4> coefficients;
    // Selected regression variables (ami outliers) with their columns and
    // final coefficients, in model order.
    std::vector<AutoVariableRecord> variables;
    std::optional<X13StockCalendarFit> stock_calendar;
    // Landing statistics of the final estimation (the G2 recipe at the
    // final point; never bit-claimed).
    std::string estimation_status;
    double log_likelihood=0., aicc=0.;
    std::size_t effective_observations=0;
    std::optional<X13StockLikelihood> stock_likelihood;
    // The decision chain, in pinned execution order (provenance of the
    // whole run: "the final description plus the chain is the full
    // regeneration recipe of the fit").
    std::vector<AutoDecisionRecord> chain;
    // The downstream composition: preadjust -> x11_selection_pass over the
    // decoder X-11 surface (forced mode) -> results-builder -> finals.
    X13FitResult fit;
};
// The automatic-spec fit. Unsupported-domain requests answer the explicit
// status of the contract matrix before any state is built; a throwing
// machine inside the pinned calc body maps onto the failed: status with
// the chain kept so far. Shape violations throw std::invalid_argument
// exactly as the seams below do.
X13AutomaticFitResult x13_automatic_fit(const X13AutomaticRequest& request);

// Stock monthly X-13 RSA0 with estimated airline parameters and a one-year
// forecast extension. The observed series starts at the given year and month.
X13AutomaticFitResult x13_estimated_rsa0(const std::vector<double>& observations, int year, int month);
// Stock RSA1 selects the transformation and detects AO, LS and TC outliers.
X13AutomaticFitResult x13_estimated_rsa1(const std::vector<double>& observations, int year, int month);
// Stock RSA2c estimates the working-day and Easter candidates before the outlier pass.
X13AutomaticFitResult x13_estimated_rsa2c(const std::vector<double>& observations, int year, int month);

struct X13StockResidualDiagnostics {
    double standard_error;
    LjungBoxTest ljung_box;
    std::optional<ResidualNormalityTest> normality;
};
X13StockResidualDiagnostics x13_stock_residual_diagnostics(const X13AutomaticFitResult& fit);

// The public automatic-route residual diagnostics (M7, item 27 cutover):
// the residuals.ser / residuals.lb / residuals.doornikhansen getters over
// the final-model likelihood of ANY complete automatic fit (the surface the
// stock consumer above pins for the airline-only estimated templates).
// The fit populates the likelihood surface on the success path; the
// Ljung-Box degrees of freedom mirror the fitted regular+seasonal ARMA
// parameter count, and the Doornik-Hansen arm mirrors the bridge capture
// domain of this surface (tramo_residual_normality), where the stock
// consumer pins the estimated-template capture domain
// (portable_residual_normality). Guards: status "available", the surface
// present with its full effective-observation residual length.
X13StockResidualDiagnostics x13_automatic_residual_diagnostics(const X13AutomaticFitResult& fit);

struct X13StockForecastErrors {
    std::string status;
    std::vector<double> standard_errors;
};
// Public original-scale preprocessing forecast uncertainty (y_ef). A nonnegative
// month count begins immediately after the input. Zero requests an empty result.
X13StockForecastErrors x13_stock_forecast_errors(const X13AutomaticFitResult& fit, int months);
} // namespace kalos
#endif
