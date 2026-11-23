#ifndef KALOS_TRAMO_FIT_HPP
#define KALOS_TRAMO_FIT_HPP

// Item 23-F F.4: the TRAMO automatic-fit assembly - the x13_automatic_fit
// analogue for the pinned jdplus 3.9.0 tramo engine (docs/contracts/
// tramo-fit-contract.md). This is the TramoKernel.process composition of the
// landed 23-series units, NOTHING ELSE: the stage order is the one the F.3
// pass-execution witnesses proved (tmp_opencode/kalos-m5-pinning.md SECTIONS
// 29.4/30.8/31), every decision is owned by its landed unit and this file
// re-implements no stage logic:
//
//   stage                     landed unit (witness tier)
//   builder registration      build_tramo_model        23-A (tramo-model-builder)
//   seasonality + log/level   decide_tramo_transform   23-B (tramo-transform-decision)
//   differencing selection    decide_tramo_differencing 23-C (tramo-differencing-decision)
//   arma selection            select_tramo_arma        23-D (tramo-arma-selection)
//   outliers selection        tramo_outlier_rounds     23-F F.2 (tramo-outlier-rounds)
//   ami() round loop          tramo_run_kernel_rounds  23-F F.3 (tramo-kernel-rounds)
//   refStats gate / mean rule tramo_ref_stats_gate et al 23-E (inside F.3)
//   final information         final_information        21-H.0 (final_information.hpp)
//   coefficient covariance    the 20c.1b.2 tramo_residual covariance algebra (D9)
//
// Like the F.3 driver, the heavy stages ENTER THE ROUND LOOP AS ORDERED
// SEAM FACTS (the request embeds the very TramoKernelRoundsRequest of the
// replay corpus - the differencing/arma/scan/estimation-module outcomes and
// module log blocks the JVM produced). What F.4 adds on top of F.3 is the
// process() assembly around that loop:
//  * the LIVE pre-stage: the 23-A builder mirror and the 23-B transform
//    decision run natively on the request series and are cross-checked
//    bit-exactly against the captured pre-stage state (fail-closed throw on
//    any disagreement - the same echo doctrine as F.3);
//  * the LIVE per-seam zero-drift cross-checks (contract "witness tiers"):
//    for every differencing/arma/scan seam whose design stays inside the
//    landed stage unit's qualified envelope, that unit is RUN on facts
//    derived from the seam stream and its result compared bit-exactly with
//    the seam fact the kernel machine already validated (result,
//    specification, mean, picks, module events). Envelope-excluded seams
//    (outlier-column warm-ML differencing states, the estimateModel
//    attempts outside the F.2 family domain) are recorded fact-only,
//    never approximated;
//  * the final surface with the D9 coefficient covariance: the regression
//    block covariance is recomputed natively as unscaledCovariance *
//    ssq/(dim - nx - k) - the ConcentratedLikelihood.covariance(k, true)
//    chain the JVM runs (bytecode-verified, pinning SECTION 32), the same
//    algebra as the 20c.1b.2 tramo_residual accurate covariance - against
//    the double-captured `cov` bit rows, kept DISTINCT from the SEATS
//    extraction variances (register D9) and from the Point-tail
//    asymptoticCovariance surface the 21-H.0 producer replays (fe_est).
//
// Refusals are fail-closed with pinned reasons BEFORE any state is built:
// non-12 annual frequency (D11), every calendar-family request - trading
// days automatic/tested, Easter tests, user variables, the TRfull/rsa4/
// rsa5c-style TD/Easter specs and the automatic AIC/BIC/WALD regression
// methods (register D10, calendar-roadmap I-10 - the 23-E calendar seam is
// the joint 21-I surface and the TRAMO auto-method family NPEs on strong
// TD inputs), the transform outliersCorrection pre-correction (+oc: the
// robust AO+LS module has never been witnessed THROUGH the fit seam),
// decomposition/approximation semantics (register D6 - the fit never
// decomposes; decomposition policy belongs to decompose_ucarima and the
// 24-D/24-E consumers), missing values in the series (the AverageInterpolator
// plan stays qualified in the 23-A/23-B unit tiers), and struct-violating
// requests. Every refusal row carries its corpus witness
// (validation/evidence/native/tramo-fit/, ctest #79).
//
// Dated 2026-11-23 (F.4 debug): the live seam re-runs consume the EVOLVING
// seasonality context - the jar module factories read ctx.seasonal at module
// creation, so stage 4 tracks ctx_seasonal from the pre-stage value through
// the seasonality_controller seams (TramoKernelSeam.seasonal_after); the
// pre-stage flag alone is the pass-1 reading only. The ami() branch is the
// decoded isFullySpecified (transform decided && !autoModel && !outliers &&
// !meanTest); the empty seam stream is a cross-checked consequence, never
// the authority. See docs/contracts/tramo-fit-contract.md ("Debug
// correctives").
//
// All doubles are raw IEEE-754 values (the probe's %.17g / bit columns;
// every replay comparison is exact bit equality). Statuses: "available",
// "refused:<...>" (the fail-closed matrix), plus the fail-closed throws
// (std::invalid_argument) for structural violations - the x13_automatic_fit
// grammar adapted to the TRAMO seam doctrine.

#include "kalos/tramo_kernel_rounds.hpp"
#include "kalos/tramo_model_builder.hpp"
#include "kalos/tramo_outlier_scan.hpp"
#include "kalos/tramo_transform_decision.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace kalos {

// One final-model outlier (ami) variable: the pinned defaultName(code,
// period), factory type, estimation-domain position, the estimated
// coefficient (the design-column slot of the cl coefficient vector) and the
// factory column over the model span (request data; the JVM fit surface
// exposes names and coefficients, the columns ride for the consumer-facing
// variable record of the x13 AutoVariableRecord analogue).
struct TramoFitVariable {
    std::string name;
    TramoOutlierKind type = TramoOutlierKind::Ao;
    int position = 0;
    double coefficient = 0.0;
    std::vector<double> column;
};

// The resolved TramoSpec scalars of the fit request (one-for-one with the
// probe's "KER <id> resolved" line, which pins each default from the JVM
// spec getters - calendar fields exist here ONLY to carry refusal state).
struct TramoFitSpec {
    std::string template_name;      // "TR0".."TR5", "TRfull", "K"(=TR3) echo
    int annual_frequency = 12;      // D11: 12 only
    int observations = 0;           // actual (non-missing) count: the calcCv input
    bool using_auto_model = true;
    // TransformSpec.
    tramo_transform_function transform_function = tramo_transform_function::automatic;
    double fct = 0.0;               // getFct(): 0 -> Math.log(0.95) default
    double log_preference = 0.0;    // resolved fct (the JVM g17 echo)
    bool outliers_correction = false;   // TransformSpec.isOutliersCorrection (+oc)
    // OutlierSpec.
    bool outliers_used = true, ao = true, ls = true, tc = true, so = true;
    double delta_tc = 0.7, critical_value = 0.0;  // va (0 -> calcCv init)
    bool outliers_ml = false;
    // AutoModelSpec.
    bool accept_airline = false;
    // (Header amendment 2026-11-23: resolved captures show every TR template
    // getCancel() = 0.05 and getTol() = 1e-7 — the b512ec2 defaults 0.044 /
    // 1e-4 mis-spoke; request data arrives from the resolved line anyway.
    // The outliers module's estimation precision is options.
    // intermediatePrecision, javap-pinned at the outliersModule() factory.)
    double pcr = 0.95, tsig = 1.0, reduce_va = 0.12, cancel = 0.05;
    double ub1 = 0.97, ub2 = 0.91;
    // EstimateSpec.
    double precision = 1e-7;            // EstimateSpec.getTol()
    double intermediate_precision = 1e-7;
    // RegressionSpec.
    bool mean_is_test = false;          // RegressionSpec.getMean().isTest()
    bool mean_used = true;              // RegressionSpec.getMean().isUsed()
    // (+nomean rows register MeanSpec.none(); header amendment 2026-11-23,
    // post-corpus: the fitbuild row of fit_nomean/fit_pure witnesses the
    // mean flag the builder consumed - the b512ec2 shape could not express
    // it.)
    // Refusal-state fields (all default false / unused = the admitted
    // envelope; any true value lands a pinned fail-closed refusal).
    bool td_automatic = false, td_test = false, easter_test = false;
    bool user_variables = false;
    tramo_auto_method auto_method = tramo_auto_method::unused; // AIC/BIC/WALD/FTEST
    // Decomposition/approximation request knob (D6): -1 = unset (admitted);
    // any other value refuses - the fit never decomposes.
    int seats_decomposition_mode = -1;
    // (Header amendment 2026-11-23, post-corpus: the builder mirror needs
    // the fixed-template SarimaSpec when using_auto_model is false - the TR
    // templates register SarimaSpec.AIRLINE, which is request data the
    // b512ec2 shape did not carry. Unused (period 0) under auto modelling.)
    tramo_sarima_spec fixed_arima;
};

// One differencing/arma/scan seam of the cross-check record: which landed
// unit the fit re-ran live, and whether the live unit AGREED bit-exactly
// with the JVM seam fact the kernel machine validated. tier:
//   "live"  - the landed unit ran on seam-derived facts and matched;
//   "fact"  - envelope-excluded seam (pinned note why) - fact-only.
struct TramoFitSeamCheck {
    std::size_t seam = 0;               // index into kernel.seams
    std::string unit;                   // "tramo_differencing" | "tramo_arma_selection"
                                        // | "tramo_outlier_rounds" | "final_estimator_run"
    std::string tier;                   // "live" | "fact"
    std::string note;                   // match echo / pinned envelope note
};

// The final-model record of the request: the rfinal-analogue capture lines
// (the JVM bits the replay pins) plus the covariance inputs. Every field is
// request data captured through the live public path (Scratch23f probe); the
// unit validates, cross-checks and exposes them - it never invents a value.
struct TramoFitFinalRecord {
    // rfinal spec/log/mean/outliers (orders as the JVM `spec` line,
    // p,d,q,bp,bd,bq).
    TramoOrders orders;
    bool log = false;
    bool mean = false;
    std::vector<std::string> outlier_names;
    std::vector<TramoFitVariable> variables;   // names/positions/coefficients
    // Model bits (the probe bit columns): desc.arima().parameters(),
    // the estimation design coefficients ge.getCoefficients (the
    // regression block, design order: [mean], then the outlier columns in
    // sortVariables order; the cl allCoefficients row additionally carries
    // the missing-value estimates - the fit_missing witness - and is not
    // the fit record carrier), the differenced
    // response the landing machine consumed (log-applied, interpolated -
    // the fit recomputes it and asserts equality), and fullResiduals.
    std::vector<double> parameters;            // JVM [phi|bphi|theta|btheta]
    std::vector<double> coefficients;          // nx regression entries
    std::vector<double> differenced_response;
    std::vector<double> full_residuals;
    // ConcentratedLikelihoodWithMissing views (the covariance algebra
    // inputs; Section 32): dim, nx, nmissing, ssq(), logDeterminant(),
    // the unscaledCovariance() (nx x nx, row-major) as the record carries
    // it, the isScalingFactor() flag, and the free-ARMA-parameter count k
    // the RegSarimaModel.of call passes to covariance(k, true)
    // (= desc.getArimaSpec().freeParametersCount()).
    int dim = 0, nx = 0, nmissing = 0;
    bool scaling_factor = false;
    double ssq = 0.0, log_determinant = 0.0;
    std::vector<double> unscaled_covariance;   // nx*nx row-major
    int free_arma_parameters = 0;
    // The JVM covariance of the regression block (the `cov` upper triangle
    // as captured, bit-pinned; replay target of the D9 recomputation).
    std::vector<double> captured_covariance_upper; // i<=j order, nx*(nx+1)/2
    // LikelihoodStatistics (the `stats` row).
    int ls_obs = 0, ls_eff = 0, ls_npar = 0;
    double ls_ll = 0.0, ls_tadj = 0.0, ls_all = 0.0, ls_ssq = 0.0;
    double ls_aic = 0.0, ls_aicc = 0.0, ls_bic = 0.0, ls_bicc = 0.0;
    double ls_bic2 = 0.0, ls_hq = 0.0;
};

// Request data for the live per-seam cross-checks the stage envelopes
// cannot recompute themselves (the probe's fitseam/fitresid lines; header
// amendment 2026-11-23, pre-capture: the landed b512ec2 shape carried only
// the arma residual hand-over, but the differencing/scan live checks read
// the module ENTRY state - what the module saw on entry, as opposed to the
// kr seam fact's post-state - so one record per module-entry seam is the
// faithful carrier).
struct TramoFitSeamInputs {
    std::size_t seam = 0;               // the module-entry seam index
    // Description state at the module call (fitseam line):
    TramoOrders entry_orders;
    bool entry_mean = false;
    std::vector<std::string> entry_outliers;   // names ("-"/empty = none)
    int pass = 0, round = 0;            // kernel fields at the call
    int entry_nobs = 0;                 // regarima().getObservationsCount()
    // arma seams (the 23-D contract owns them as the JVM-pinned hand-over):
    std::vector<double> ols_residuals;  // optional verification when the entry has only a mean or no regressors
    int nobs = 0;                       // regarima().getObservationsCount()
    int spc = -1;                       // the comespa gate readout (-1: none)
};

struct TramoFitRequest {
    std::string id;                     // corpus id echo
    std::vector<double> series;         // complete, finite, monthly
    int start_year = 1990, start_month = 1;
    TramoFitSpec spec;
    // The builder-request extras the fit mirror feeds the 23-A unit
    // (preset outlier/ramp/intervention variables; empty in the corpus).
    std::vector<tramo_outlier_request> builder_outliers;
    // The kernel-loop request (spec facts, pre-stage facts, seam facts) -
    // exactly the F.3 request the kr corpus replays.
    TramoKernelRoundsRequest kernel;
    // Per-module-entry seam records (diff/arma/scan), sorted by seam index.
    std::vector<TramoFitSeamInputs> seam_inputs;
    // The final-model record.
    TramoFitFinalRecord final;
};

// The coefficient-covariance block of the result (D9): the native
// recomputation of the JVM regression-block covariance
// V = unscaledCovariance * (ssq / (dim - nx - k)) - the Concentrated
// Likelihood.covariance(k, true) chain, byte-identical algebra to the
// 20c.1b.2 tramo_residual accurate covariance - alongside the captured JVM
// upper triangle. `matches_capture` is computed by the unit itself: true
// only after the raw-bit comparison of every upper-triangle entry.
struct TramoFitCovariance {
    std::size_t n = 0;                          // nx
    std::vector<double> covariance;             // n x n row-major (native)
    std::vector<double> covariance_upper;       // i<=j (native)
    std::vector<double> captured_upper;         // JVM bits
    bool matches_capture = false;
};

// The fit result: the full regeneration surface of one automatic TRAMO run.
struct TramoFitResult {
    // "available" | "refused:<reason>" (fail-closed matrix; chain empty).
    std::string status;
    std::string refusal_note;           // the pinned reason detail, if any
    // The final description.
    bool log = false, mean = false;
    TramoOrders orders;
    std::vector<double> parameters;     // JVM [phi|bphi|theta|btheta] bits
    std::vector<TramoFitVariable> variables;
    TramoFitCovariance covariance;
    // sigma2 bookkeeping (the cl record views): ssq, the effective count
    // (cl dim), nx, nmissing, logDeterminant and the derived
    // scale = ssq/(dim - nx - k) of the covariance algebra.
    double ssq = 0.0, covariance_scale = 0.0;
    int cl_dim = 0, cl_nx = 0, cl_nmissing = 0;
    double log_determinant = 0.0;
    // LikelihoodStatistics echo (estimation metadata).
    int ls_obs = 0, ls_eff = 0, ls_npar = 0;
    double ls_ll = 0.0, ls_tadj = 0.0, ls_all = 0.0, ls_ssq = 0.0;
    double ls_aic = 0.0, ls_aicc = 0.0, ls_bic = 0.0, ls_bicc = 0.0;
    double ls_bic2 = 0.0, ls_hq = 0.0;
    std::vector<double> full_residuals;
    // The live pre-stage records (the transformation + differencing
    // decisions the x13 analogue's chain carries).
    tramo_build_result builder;
    TramoTransformDecisionResult transformation;
    bool transform_matches_capture = false;  // log/seasonality/spec vs pre-stage facts
    // The ami() loop: the F.3 run verbatim (events, snapshots, decisions,
    // init bits, final state) plus the per-seam witness tiers.
    TramoKernelRoundsRun kernel;
    std::vector<TramoFitSeamCheck> seam_checks;
    std::size_t seams_live = 0, seams_fact = 0;
};

// The automatic TRAMO fit. Order of operations (contract "stage order"):
// validate/refuse -> build_tramo_model (cross-checked) ->
// decide_tramo_transform (cross-checked) -> tramo_run_kernel_rounds ->
// per-seam live cross-checks -> final surface + D9 covariance recomputation.
// Structural violations and every desync between a live unit and the
// captured stream throw std::invalid_argument naming the site; out-of-
// envelope REQUESTS answer the refused: status instead where the matrix
// says so (the refusals are decisions, not shape errors).
TramoFitResult tramo_fit(const TramoFitRequest& request);

}  // namespace kalos

#endif  // KALOS_TRAMO_FIT_HPP
