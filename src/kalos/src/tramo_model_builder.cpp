#include "kalos/tramo_model_builder.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace kalos {
namespace {

// ---- template-expansion table (pinned data) -----------------------------
// Rows captured from the pinned JVM specs through the double-captured req
// echo (validation/evidence/native/tramo-model-builder/inputs.csv). All
// seven templates share the undefined-airline SARIMA shape with period 0
// (d 1, bd 1, theta/btheta Undefined), a full transform/estimate span and
// deltaTC 0.7; TR0 disables all outlier detection flags, TR1..TRfull enable
// AO/LS/TC (SO stays off); TR0..TR3 carry the trend-constant mean spec with
// an Undefined coefficient, TR4/TR5/TRfull the DEFAULT_UNUSED mean (null
// coefficient). Calendar data: TR2/TR4 register nothing (TD test), TR5 is
// TD7-tested, TRfull is automatic-FTEST plus a tested IncludeEaster - the
// builder gates them all away from registration.
const tramo_sarima_spec kTemplateArima{
    0, 1, 1, {}, {tramo_parameter{}}, {}, {tramo_parameter{}}};

const tramo_mean_spec kMeanTrendConstant{true, true, false, true, tramo_parameter{}};
const tramo_mean_spec kMeanDefaultUnused{false, false, false, false, tramo_parameter{}};

const tramo_trading_days_spec kTdUnused{};
const tramo_trading_days_spec kTd2Tested{true,  true, false, false,
                                             tramo_trading_days_type::td2,
                                             tramo_auto_method::unused, 0,
                                             tramo_lop::leap_year};
const tramo_trading_days_spec kTd7Tested{true,  true, false, false,
                                             tramo_trading_days_type::td7,
                                             tramo_auto_method::unused, 0,
                                             tramo_lop::leap_year};
const tramo_trading_days_spec kTd7AutoFtest{true,  false, true, false,
                                                tramo_trading_days_type::td7,
                                                tramo_auto_method::ftest, 0,
                                                tramo_lop::leap_year};

const tramo_easter_spec kEasterUnused{};
const tramo_easter_spec kEasterStandardTested{true, true,
                                                  tramo_easter_type::standard, 6,
                                                  false};
const tramo_easter_spec kEasterIncludeTested{true, true,
                                                 tramo_easter_type::include_easter, 6,
                                                 false};

const tramo_template_fields kTemplates[7] = {
    {"TR0", false, tramo_transform_function::none, tramo_lop::none, kTemplateArima,
     kMeanTrendConstant, false, false, false, false, 0.7, kTdUnused, kEasterUnused},
    {"TR1", false, tramo_transform_function::automatic, tramo_lop::none, kTemplateArima,
     kMeanTrendConstant, true, true, true, false, 0.7, kTdUnused, kEasterUnused},
    {"TR2", false, tramo_transform_function::automatic, tramo_lop::none, kTemplateArima,
     kMeanTrendConstant, true, true, true, false, 0.7, kTd2Tested, kEasterStandardTested},
    {"TR3", true, tramo_transform_function::automatic, tramo_lop::none, kTemplateArima,
     kMeanTrendConstant, true, true, true, false, 0.7, kTdUnused, kEasterUnused},
    {"TR4", true, tramo_transform_function::automatic, tramo_lop::none, kTemplateArima,
     kMeanDefaultUnused, true, true, true, false, 0.7, kTd2Tested, kEasterStandardTested},
    {"TR5", true, tramo_transform_function::automatic, tramo_lop::none, kTemplateArima,
     kMeanDefaultUnused, true, true, true, false, 0.7, kTd7Tested, kEasterStandardTested},
    {"TRfull", true, tramo_transform_function::automatic, tramo_lop::none, kTemplateArima,
     kMeanDefaultUnused, true, true, true, false, 0.7, kTd7AutoFtest, kEasterIncludeTested},
};

bool valid_month(int month) {
    return month >= 1 && month <= 12;
}

tramo_parameter default_parameter() {
    return tramo_parameter{};  // Undefined at 0.0: the pinned Variable.builder() default.
}

const std::vector<tramo_parameter>& request_coefficients(
    const std::vector<tramo_parameter>& coefficients, std::vector<tramo_parameter>& scratch) {
    if (coefficients.empty()) {
        scratch.assign(1, default_parameter());
        return scratch;
    }
    return coefficients;
}

// ModelDescription.sortVariables category index (the six pinned filter
// passes; passes 5/6 sort by the IOutlier position).
int sort_category(const tramo_model_variable& v) {
    switch (v.kind) {
        case tramo_variable_kind::trend_constant:
            return 0;
        case tramo_variable_kind::user_variable:
        case tramo_variable_kind::intervention:
        case tramo_variable_kind::ramp:
            return 1;
        case tramo_variable_kind::trading_days:
        case tramo_variable_kind::stock_trading_days:
        case tramo_variable_kind::length_of_period:
            return 2;
        case tramo_variable_kind::easter_variable:
        case tramo_variable_kind::julian_easter_variable:
            return 3;
        case tramo_variable_kind::additive_outlier:
        case tramo_variable_kind::level_shift:
        case tramo_variable_kind::transitory_change:
        case tramo_variable_kind::periodic_outlier:
            return v.has_ami ? 5 : 4;
    }
    return 6;
}

int position_key(const tramo_model_variable& v) {
    return v.position_year * 100 + v.position_month;
}

}  // namespace

bool tramo_sarima_spec::is_airline() const {
    return phi.size() == 0 && d == 1 && theta.size() == 1 && bphi.size() == 0 && bd == 1 &&
           btheta.size() == 1;
}

bool tramo_sarima_spec::is_undefined() const {
    for (const auto* block : {&phi, &theta, &bphi, &btheta}) {
        for (const auto& p : *block) {
            if (p.kind != tramo_parameter_kind::undefined) return false;
        }
    }
    return true;
}

const tramo_template_fields* tramo_template(const std::string& name) {
    static const std::pair<const char*, int> kNames[14] = {
        {"tr0", 0}, {"tr1", 1}, {"tr2", 2}, {"tr3", 3}, {"tr4", 4}, {"tr5", 5},
        {"trfull", 6}, {"TR0", 0}, {"TR1", 1}, {"TR2", 2}, {"TR3", 3}, {"TR4", 4},
        {"TR5", 5}, {"TRfull", 6},
    };
    for (const auto& entry : kNames) {
        if (name == entry.first) return &kTemplates[entry.second];
    }
    return nullptr;
}

// The pinned initializeVariables: mean, calendar, outliers, users, interventions,
// ramps, in that order. Runs only when RegressionSpec.isUsed() (see
// tramo_regression_used); when skipped, the mean flag keeps its default off and no
// variable registers - the later autoModel overwrite in initializeArima still
// applies.
void initialize_variables(const tramo_build_request& request, tramo_build_result& result) {
    const bool td_registers =
        request.trading_days.used && !request.trading_days.test && !request.trading_days.automatic;
    const bool easter_registers = request.easter.used && !request.easter.test;
    std::vector<tramo_parameter> scratch;

    if (request.mean.used) {
        if (request.mean.coefficient_present &&
            request.mean.coefficient.kind == tramo_parameter_kind::fixed) {
            // The pinned "const" quirk: a fixed mean coefficient becomes a
            // TrendConstant(d, bd, estimation start) named "const" with the
            // Undefined regeffect, and the estimated mean is switched off.
            tramo_model_variable constant;
            constant.name = "const";
            constant.kind = tramo_variable_kind::trend_constant;
            constant.dim = 1;
            constant.has_regeffect = true;
            constant.regeffect = "Undefined";
            constant.coefficients = {request.mean.coefficient};
            constant.constant_d = request.arima.d;
            constant.constant_bd = request.arima.bd;
            constant.constant_ref_year = result.estimation_year;
            constant.constant_ref_month = result.estimation_position + 1;
            result.variables.push_back(std::move(constant));
            result.mean = false;
        } else {
            result.mean = true;
        }
    } else {
        result.mean = false;
    }

    if (td_registers) {
        tramo_model_variable td;
        td.name = "td";
        td.kind = tramo_variable_kind::trading_days;
        td.dim = 1;  // GenericTradingDays.contrasts(TD2): one contrast column.
        td.has_regeffect = true;
        td.regeffect = "CalendarEffect";
        td.coefficients = {default_parameter()};  // getTdCoefficients() default.
        td.td_type = request.trading_days.type;
        td.td_groups = 2;  // DayClustering.TD2.getGroupsCount().
        result.variables.push_back(std::move(td));

        tramo_model_variable lp;
        lp.name = "lp";
        lp.kind = tramo_variable_kind::length_of_period;
        lp.dim = 1;
        lp.has_regeffect = true;
        lp.regeffect = "CalendarEffect";
        lp.coefficients = {default_parameter()};  // getLpCoefficient() default.
        lp.lop = request.trading_days.lop;
        result.variables.push_back(std::move(lp));
    }
    if (easter_registers) {
        tramo_model_variable easter;
        easter.name = "easter";
        easter.kind = tramo_variable_kind::easter_variable;
        easter.dim = 1;
        easter.has_regeffect = true;
        easter.regeffect = "CalendarEffect";
        easter.coefficients = {default_parameter()};  // EasterSpec coefficient default.
        easter.easter_duration = request.easter.duration;
        easter.easter_end_position = -1;  // Standard: endPosition -1, Simple correction.
        result.variables.push_back(std::move(easter));
    }

    for (const auto& outlier : request.outliers) {
        tramo_model_variable variable;
        std::string effect;
        if (outlier.code == "AO") {
            variable.kind = tramo_variable_kind::additive_outlier;
            effect = "Irregular";  // AdditiveOutlierFactory
        } else if (outlier.code == "LS") {
            variable.kind = tramo_variable_kind::level_shift;
            effect = "Trend";  // LevelShiftFactory.FACTORY_ZEROENDED
        } else if (outlier.code == "SO") {
            variable.kind = tramo_variable_kind::periodic_outlier;
            effect = "Seasonal";  // PeriodicOutlierFactory(freq, true)
            variable.period = 12;
        } else if (outlier.code == "TC") {
            variable.kind = tramo_variable_kind::transitory_change;
            effect = "Irregular";  // TransitoryChangeFactory(deltaTC)
            variable.delta = request.delta_tc;
        } else {
            continue;  // pinned exact-case switch: unknown codes are skipped.
        }
        variable.name = outlier.name;
        variable.dim = 1;
        variable.has_regeffect = true;
        variable.regeffect =
            outlier.has_regeffect_override ? outlier.regeffect_override : std::move(effect);
        variable.coefficients = request_coefficients(outlier.coefficients, scratch);
        variable.outlier_code = outlier.code;
        variable.position_year = outlier.position_year;
        variable.position_month = outlier.position_month;
        result.variables.push_back(std::move(variable));
    }

    for (const auto& user : request.user_variables) {
        // TsContextVariable.instantiateFrom(new ModellingContext(), name) is
        // null in this envelope: the pinned warning, and no registration.
        result.warnings.push_back(user.name + ": not found (removed from the model)");
    }

    for (const auto& ivar : request.interventions) {
        tramo_model_variable variable;
        variable.name = ivar.name;
        variable.kind = tramo_variable_kind::intervention;
        variable.dim = 1;
        variable.coefficients = request_coefficients(ivar.coefficients, scratch);
        variable.delta = ivar.delta;
        variable.delta_seasonal = ivar.delta_seasonal;
        variable.sequences = ivar.sequences;
        result.variables.push_back(std::move(variable));
    }
    for (const auto& ramp : request.ramps) {
        tramo_model_variable variable;
        variable.name = ramp.name;
        variable.kind = tramo_variable_kind::ramp;
        variable.dim = 1;
        variable.coefficients = request_coefficients(ramp.coefficients, scratch);
        variable.start_year = ramp.start_year;
        variable.start_month = ramp.start_month;
        variable.end_year = ramp.end_year;
        variable.end_month = ramp.end_month;
        result.variables.push_back(std::move(variable));
    }
}

bool tramo_regression_used(const tramo_build_request& request) {
    // Pinned RegressionSpec.isUsed(): mean, calendar (its used flags carry the
    // spec's equals-vs-DEFAULT state), or any non-empty variable list.
    return request.mean.used || request.trading_days.used || request.easter.used ||
           !request.outliers.empty() || !request.user_variables.empty() ||
           !request.ramps.empty() || !request.interventions.empty();
}

tramo_build_result build_tramo_model(const tramo_build_request& request) {
    // -- envelope validation (fail-closed; every rejection marks a path the
    // double-captured corpus does not qualify) --------------------------------
    if (request.annual_frequency != 12)
        throw std::invalid_argument(
            "the tramo model builder envelope is monthly frequency only (D11)");
    if (request.values.empty())
        throw std::invalid_argument("the tramo model builder requires a non-empty series");
    for (const double value : request.values) {
        if (std::isinf(value))
            throw std::invalid_argument(
                "infinite series values are refused; NaN is the only missing marker");
    }
    bool any_finite = false;
    for (const double value : request.values) any_finite = any_finite || std::isfinite(value);
    if (!any_finite)
        throw std::invalid_argument(
            "the pinned AverageInterpolator leaves an all-missing series undefined");
    if (!valid_month(request.start_month))
        throw std::invalid_argument("series start month must be in 1..12");

    if (!request.using_auto_model && request.arima.period != 0 && request.arima.period != 12)
        throw std::invalid_argument(
            "a declared SARIMA period must be 0 (unset) or 12 in this envelope");

    const bool td_registers =
        request.trading_days.used && !request.trading_days.test && !request.trading_days.automatic;
    if (td_registers) {
        if (request.trading_days.stock)
            throw std::invalid_argument(
                "stock trading days are outside the qualified envelope");
        if (request.trading_days.type != tramo_trading_days_type::td2)
            throw std::invalid_argument(
                "only the TD2 trading-days type is qualified by the corpus");
        if (request.trading_days.lop == tramo_lop::none)
            throw std::invalid_argument(
                "a registered trading-days block without a length-of-period rule would build a "
                "null lp core in the pinned builder; refused");
        if (request.trading_days.lop != tramo_lop::leap_year)
            throw std::invalid_argument(
                "only the LeapYear length-of-period rule is qualified by the corpus");
    }
    const bool easter_registers = request.easter.used && !request.easter.test;
    if (easter_registers) {
        if (request.easter.julian)
            throw std::invalid_argument("julian Easter variables are outside the qualified envelope");
        if (request.easter.type != tramo_easter_type::standard)
            throw std::invalid_argument(
                "only the standard (non-include) Easter type is qualified by the corpus");
        if (request.easter.duration < 0)
            throw std::invalid_argument("Easter duration must be non-negative");
    }
    for (const auto& user : request.user_variables) {
        if (user.resolved)
            throw std::invalid_argument(
                "resolved user variables belong to the 21-I joint surface; refused");
    }
    for (const auto& outlier : request.outliers) {
        if (!valid_month(outlier.position_month))
            throw std::invalid_argument("outlier position month must be in 1..12");
    }
    for (const auto& ramp : request.ramps) {
        if (!valid_month(ramp.start_month) || !valid_month(ramp.end_month))
            throw std::invalid_argument("ramp month must be in 1..12");
    }
    for (const auto& ivar : request.interventions) {
        for (const auto& range : ivar.sequences) {
            if (!valid_month(range.first.second) || !valid_month(range.second.second))
                throw std::invalid_argument("intervention sequence month must be in 1..12");
        }
    }

    tramo_build_result result;

    // -- ModelDescription construction: domain + estimation domain ------------
    const int length = static_cast<int>(request.values.size());
    result.annual_frequency = 12;
    result.domain_year = request.start_year;
    result.domain_position = request.start_month - 1;
    result.domain_length = length;
    result.estimation_year = request.start_year;
    result.estimation_position = request.start_month - 1;
    result.estimation_length = length;

    // -- initializeMissing: the AverageInterpolator.cleanMissings scan,
    // transcribed operation for operation (right-to-left pass, (left + right)
    // / 2.0 for interior runs, left value for trailing runs, first finite
    // value for the leading run; the missing list is then sorted ascending
    // exactly as ModelDescription.interpolate calls Arrays.sort).
    result.interpolated = request.values;
    int next = -1;
    int last = -1;
    for (int i = length - 1; i >= 0; --i) {
        if (!std::isfinite(result.interpolated[static_cast<std::size_t>(i)])) {
            if (next == -1) next = i + 1;
            last = i - 1;
            result.missing.push_back(i);
        } else if (last != -1) {
            const double fill = next < length
                                    ? (result.interpolated[static_cast<std::size_t>(next)] +
                                       result.interpolated[static_cast<std::size_t>(last)]) / 2.0
                                    : result.interpolated[static_cast<std::size_t>(last)];
            for (int j = last + 1; j < next; ++j)
                result.interpolated[static_cast<std::size_t>(j)] = fill;
            next = -1;
            last = -1;
        }
    }
    if (next > 0) {
        for (int j = 0; j < next; ++j)
            result.interpolated[static_cast<std::size_t>(j)] =
                result.interpolated[static_cast<std::size_t>(next)];
    }
    std::sort(result.missing.begin(), result.missing.end());

    // -- initializeTransformation ---------------------------------------------
    result.log_transformation = request.transform_function == tramo_transform_function::log;
    result.preadjustment = request.transform_adjust;

    // -- initializeVariables (pinned RegressionSpec.isUsed gate) ------------------
    if (tramo_regression_used(request)) initialize_variables(request, result);

    // -- initializeArima -------------------------------------------------------
    if (request.using_auto_model) {
        // setAirline(annualFrequency != 1): monthly lands the undefined
        // airline; the annual (m011) branch is outside this envelope.
        result.arima = tramo_sarima_spec{
            12, 1, 1, {}, {default_parameter()}, {}, {default_parameter()}};
        result.mean = true;  // the pinned overwrite AFTER initializeMean
    } else {
        // setArimaSpec(spec.getArima().withPeriod(12)): period 12 keeps the
        // spec (identity), period 0 re-periodizes; the withPeriod(1)
        // seasonal-drop branch is unreachable monthly.
        result.arima = request.arima;
        result.arima.period = 12;
    }

    // -- ModelDescription.regarima: sortVariables, transcribed as the six
    // stable filter passes; the outlier passes additionally sort by position
    // (pinned OUTLIER_COMPARATOR on IOutlier.getPosition). Stream order is
    // the registration order.
    std::stable_sort(result.variables.begin(), result.variables.end(),
                     [](const tramo_model_variable& a, const tramo_model_variable& b) {
                         const int ca = sort_category(a);
                         const int cb = sort_category(b);
                         if (ca != cb) return ca < cb;
                         if (ca == 4 || ca == 5)
                             return position_key(a) < position_key(b);
                         return false;
                     });
    return result;
}

}  // namespace kalos
