#ifndef KALOS_TRAMO_MODEL_BUILDER_HPP
#define KALOS_TRAMO_MODEL_BUILDER_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace kalos {

// Item 23-A: the TRAMO model-builder decision mirror. This unit reproduces,
// operation for operation, what the LIVE pinned jdplus 3.9.0 builder
// (jdplus.tramoseats.base.core.tramo.TramoModelBuilder.build, reached via
// TramoSpec + ModellingContext) decides when it turns an expanded TramoSpec
// and a series into a ModelDescription: the start rule (airline or the
// re-periodized fixed spec), the initial SARIMA spec, the log flag and the
// length-of-period preadjustment, the estimated-mean flag with its
// fixed-coefficient "const" quirk, the AverageInterpolator missing plan, the
// registered regression-variable list with the regeffect attribute
// vocabulary and the ModelDescription.sortVariables order, and the builder's
// own ProcessingLog warnings. It mirrors the DECISIONS and the registered
// variable metadata only: no calendar design algebra (the day-clustering
// contrasts, LoP and Easter column values - deferred to the 21-I joint
// surface per register D10) and no collinearity check (the pinned check() is
// a Householder QR over the design matrices; the qualified corpus is
// full-rank and the unit makes no collinearity claim).
//
// Envelope (register D11, monthly-first, faithful-to-capture): the unit
// boundary accepts monthly frequency only; the annual/quarterly corpus
// cases are JVM-qualified record rows the unit REFUSES. Requests the pinned
// JVM path never demonstrates in the double-captured corpus are refused
// fail-closed (std::invalid_argument), never silently reinterpreted:
// non-finite values other than the NaN missing marker, an all-missing
// series, non-default transform/estimate spans (fixed as "all" here), a
// resolvable user variable (its design data lives in 21-I), stock/holidays/
// user-trading-day calendar requests, trading-day types other than TD2 with
// a LeapYear length-of-period rule, Easter variables other than the standard
// (non-Julian, non-include) form with an Undefined-type registration,
// non-default declared SARIMA periods (only 0/absent and 12 are accepted on
// the fixed-order path). A fully unused regression block
// (tramo_regression_used false) is NOT a refusal: the pinned builder skips
// initializeVariables entirely and the mean stays off (mean_off row).

// Pinned enum mirrors (constant names printed by the oracle getters).
enum class tramo_transform_function { none, automatic, log };        // TransformationType
enum class tramo_lop { none, leap_year, length_of_period };          // LengthOfPeriodType
enum class tramo_trading_days_type { none, td2, td3, td4, td6, td7, working_days };
enum class tramo_auto_method { unused, ftest, wald, bic, aic };      // TradingDaysSpec.AutoMethod
enum class tramo_regression_test { none, separate_t, joint_f };      // RegressionTestType
enum class tramo_easter_type { unused, standard, include_easter, include_easter_monday };
// Java ParameterType: isFree() is "type != Fixed", so the pinned getters
// render Undefined/Initial/Estimated identically as free values.
enum class tramo_parameter_kind { undefined, initial, fixed, estimated };

struct tramo_parameter {
    tramo_parameter_kind kind = tramo_parameter_kind::undefined;
    double value = 0.0;
};

// Pinned SarimaSpec shape: the orders p,q,bp,bq are the PARAMETER ARRAY
// LENGTHS (the Java constructor takes only period, d, bd and the four
// arrays); p = phi.size(), q = theta.size(), bp = bphi.size(),
// bq = btheta.size().
struct tramo_sarima_spec {
    int period = 0;  // 0 = unspecified (template data); withPeriod(12) applies it
    int d = 0;
    int bd = 0;
    std::vector<tramo_parameter> phi, theta, bphi, btheta;

    // Pinned SarimaSpec.isAirline(): p == 0 && d == 1 && q == 1 && bp == 0
    // && bd == 1 && bq == 1 (period and parameter kinds irrelevant).
    bool is_airline() const;
    // Pinned SarimaSpec.isUndefined(): every parameter of all four arrays
    // is the Undefined singleton (kind == undefined here; empty arrays are
    // default-true).
    bool is_undefined() const;
};

struct tramo_mean_spec {
    bool used = false;
    bool trend_constant = false;
    bool test = false;
    // MeanSpec.getCoefficient() may be null (DEFAULT_UNUSED); null and a
    // value are distinct request states (isFixed(null) is false).
    bool coefficient_present = false;
    tramo_parameter coefficient;
};

struct tramo_trading_days_spec {
    bool used = false;
    bool test = false;
    bool automatic = false;
    bool stock = false;
    tramo_trading_days_type type = tramo_trading_days_type::none;
    tramo_auto_method auto_method = tramo_auto_method::unused;
    int stock_days = 0;
    tramo_lop lop = tramo_lop::none;
};

struct tramo_easter_spec {
    bool used = false;
    bool test = false;
    tramo_easter_type type = tramo_easter_type::unused;
    int duration = 6;
    bool julian = false;
};

// The pinned template-expansion table (TramoSpec.fromString / the TR*
// statics). Accepted names, exactly the fromString spellings (case-
// sensitive): tr0 tr1 tr2 tr3 tr4 tr5 trfull and TR0..TR5, TRfull. Any
// other name - including "Trfull" and every rsa* name - yields nullptr
// (fail-closed). The fields are the raw getter values captured from the
// JVM specs (double-captured, pinned evidence); the calendar fields arrive
// as data and are builder-gated (test or automatic paths never register).
// Coefficient arrays are NOT part of the template row: the JVM spec
// defaults them to a single Undefined parameter, which is what the request
// defaults below carry as well.
struct tramo_template_fields {
    const char* name;  // canonical uppercase spelling
    bool using_auto_model;
    tramo_transform_function transform_function;
    tramo_lop transform_adjust;
    tramo_sarima_spec arima;
    tramo_mean_spec mean;
    // OutlierSpec flags (detection policy for later slices; the builder
    // itself only consumes delta_tc).
    bool outliers_ao, outliers_ls, outliers_tc, outliers_so;
    double delta_tc;
    tramo_trading_days_spec trading_days;
    tramo_easter_spec easter;
};
const tramo_template_fields* tramo_template(const std::string& name);

// ---- build request -----------------------------------------------------

struct tramo_outlier_request {
    std::string name;
    std::string code;  // exact case: AO LS TC SO register; any other code
                       // is silently skipped exactly as the builder's
                       // case-sensitive switch does
    int position_year = 1990;
    int position_month = 1;  // 1..12
    // Preset "regeffect" attribute (Variable.setAttribute before the
    // builder would set its factory default); when present the builder
    // keeps it verbatim.
    bool has_regeffect_override = false;
    std::string regeffect_override;
    // Variable coefficient list; the JVM Variable.builder() default is a
    // single Undefined-rank parameter at value 0.
    std::vector<tramo_parameter> coefficients;
};

struct tramo_user_variable_request {
    std::string name;
    // ModellingContext.instantiateFrom result. The first envelope cannot
    // carry resolved user data (21-I), so resolved == true is refused;
    // resolved == false takes the pinned warning path.
    bool resolved = false;
};

struct tramo_ramp_request {
    std::string name;
    int start_year = 1990, start_month = 1;
    int end_year = 1990, end_month = 1;
    std::vector<tramo_parameter> coefficients;
};

struct tramo_intervention_request {
    std::string name;
    double delta = 0.0;
    double delta_seasonal = 0.0;
    // InterventionVariable.builder().sequence(Range.of(start, end))...
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> sequences;
    std::vector<tramo_parameter> coefficients;
};

struct tramo_build_request {
    // Series. NaN is the missing marker; +-Inf is refused. The frequency
    // is part of the request so the fail-closed matrix can pin it: only
    // annual_frequency == 12 builds.
    std::vector<double> values;
    int annual_frequency = 12;
    int start_year = 1990;
    int start_month = 1;  // 1..12
    // Spans are always the full domain in this envelope (the TimeSelector
    // knobs are refused by construction: no fields exist for them).
    bool using_auto_model = false;
    tramo_transform_function transform_function = tramo_transform_function::none;
    tramo_lop transform_adjust = tramo_lop::none;
    tramo_sarima_spec arima;  // read only when !using_auto_model
    tramo_mean_spec mean;
    tramo_trading_days_spec trading_days;
    tramo_easter_spec easter;
    double delta_tc = 0.7;
    std::vector<tramo_outlier_request> outliers;
    std::vector<tramo_user_variable_request> user_variables;
    std::vector<tramo_ramp_request> ramps;
    std::vector<tramo_intervention_request> interventions;
};

// RegressionSpec.isUsed(): the mean spec, the calendar spec (the used flags
// carry the pinned equals-vs-DEFAULT derivation of the expanded spec) or any
// of the four variable lists being non-empty. When false, the pinned
// initializeVariables returns BEFORE the mean block, so nothing registers
// and the mean flag keeps its default off (the mean_off corpus row; the
// autoModel overwrite in initializeArima still applies - mean_auto_off).
bool tramo_regression_used(const tramo_build_request& request);

// ---- built model description mirror -------------------------------------

enum class tramo_variable_kind {
    trend_constant,
    additive_outlier,
    level_shift,
    transitory_change,
    periodic_outlier,
    trading_days,
    stock_trading_days,
    length_of_period,
    easter_variable,
    julian_easter_variable,
    user_variable,
    intervention,
    ramp,
};

struct tramo_model_variable {
    std::string name;
    tramo_variable_kind kind;
    std::size_t dim = 1;
    bool has_regeffect = false;
    std::string regeffect;
    bool has_ami = false;  // "ami" attribute (set by 23-B detection, never here)
    std::vector<tramo_parameter> coefficients;

    // Core metadata by kind (mirrors the pinned getters the oracle read):
    // outliers: code + position (+ SO period);
    // ramp: start/end; intervention: deltas + sequence ranges;
    // trading_days: clustering type + group count;
    // length_of_period: rule; easter_variable: duration + end position
    // (Simple mean correction); julian_easter_variable: duration +
    // gregorian flag; trend_constant: d + bd + reference period.
    std::string outlier_code;
    int position_year = 0, position_month = 0;
    int period = 0;  // PeriodicOutlier
    int start_year = 0, start_month = 0, end_year = 0, end_month = 0;
    double delta = 0.0, delta_seasonal = 0.0;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> sequences;
    tramo_trading_days_type td_type = tramo_trading_days_type::none;
    int td_groups = 0;
    tramo_lop lop = tramo_lop::none;
    int easter_duration = 0, easter_end_position = 0;
    bool easter_gregorian = false;
    int constant_d = 0, constant_bd = 0, constant_ref_year = 0, constant_ref_month = 0;
};

struct tramo_build_result {
    int annual_frequency = 12;
    // TsDomain echoes: year : annualPosition (0-based) : frequency : length.
    int domain_year = 0, domain_position = 0, domain_length = 0;
    int estimation_year = 0, estimation_position = 0, estimation_length = 0;
    // initializeMissing + AverageInterpolator plan: every non-finite
    // position (ascending, pinned Arrays.sort) and the FULL interpolated
    // array (missing slots replaced per the pinned fill rules: interior
    // run (left + right) / 2 with the run's outer neighbours, leading run
    // = first finite value, trailing run = last finite value to its left).
    std::vector<int> missing;
    std::vector<double> interpolated;
    bool log_transformation = false;  // function == Log (Auto is decided later)
    tramo_lop preadjustment = tramo_lop::none;
    bool mean = false;
    tramo_sarima_spec arima;
    // ModelDescription.sortVariables order: TrendConstant, user category
    // (users, interventions, ramps in registration order), days-related
    // (td, lp), moving holidays (easter), manual outliers sorted by
    // position, auto-identified outliers sorted by position.
    std::vector<tramo_model_variable> variables;
    // ProcessingLog warnings in emission order ("u1: not found (removed
    // from the model)").
    std::vector<std::string> warnings;
};

// Builds the mirror. Throws std::invalid_argument outside the envelope
// above; otherwise reproduces the pinned builder facts bit for bit on the
// monthly corpus (decision fields exact, doubles through IEEE-exact
// operations only: the copied data and the (a + b) / 2.0 average).
tramo_build_result build_tramo_model(const tramo_build_request& request);

}  // namespace kalos

#endif
