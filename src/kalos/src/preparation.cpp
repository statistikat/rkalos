#include "kalos/detail/hostlineage.hpp"
#include "kalos/preparation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace kalos {
namespace {
constexpr auto row_limit = std::numeric_limits<int>::max();
std::int64_t subtract(std::int64_t a, std::int64_t b) {
    if ((b > 0 && a < std::numeric_limits<std::int64_t>::min() + b) ||
        (b < 0 && a > std::numeric_limits<std::int64_t>::max() + b))
        throw std::overflow_error("time-domain index overflow");
    return a - b;
}
TimeDomain domain_at(std::int64_t ordinal, int frequency, std::int64_t count) {
    auto year = ordinal / frequency;
    auto period = ordinal % frequency;
    if (period < 0) { --year; period += frequency; }
    return {year, static_cast<int>(period + 1), frequency, count};
}
void poll(const Interrupt& interrupt, std::size_t i) {
    if (interrupt && (i & 4095U) == 0) interrupt();
}
void names_valid(const std::vector<std::string>& names) {
    std::unordered_set<std::string> seen;
    for (const auto& name : names)
        if (name.empty() || !seen.insert(name).second)
            throw std::invalid_argument("Regressors must have distinct nonempty names.");
}
}

EstimationWindow estimation_window(const TimeDomain& source,
    std::optional<std::int64_t> first, std::optional<std::int64_t> last) {
    if (source.length() < 1 || source.length() > row_limit)
        throw std::invalid_argument("Input domain exceeds the supported matrix row range.");
    const auto start = source.start_ordinal();
    const auto end = start + source.length() - 1;
    const auto a = first.value_or(start), b = last.value_or(end);
    if (a < start || b > end || a > end || b < start)
        throw std::invalid_argument("The estimation window must lie entirely inside the input domain.");
    if (a > b) throw std::invalid_argument("estimation_start must not follow estimation_end.");
    return {domain_at(a, source.frequency(), b - a + 1), static_cast<std::size_t>(a - start)};
}

PreparedRegressors align_regressors(const TimeDomain& observed,
    const std::vector<NamedRegressor>& regressors, std::int64_t backcasts,
    std::int64_t forecasts, const Interrupt& interrupt) {
    if (backcasts < 0 || forecasts < 0 || backcasts > row_limit || forecasts > row_limit)
        throw std::invalid_argument("Horizons must be finite nonnegative integers within the supported range.");
    if (observed.length() < 1 || observed.length() > row_limit - backcasts - forecasts)
        throw std::length_error("Required domain exceeds the supported matrix row range.");
    const auto count = observed.length() + backcasts + forecasts;
    if (regressors.size() > static_cast<std::size_t>(row_limit / count))
        throw std::length_error("Regressor matrix exceeds the supported element range.");
    const auto start = subtract(observed.start_ordinal(), backcasts);
    PreparedRegressors result{observed, domain_at(start, observed.frequency(), count),
        static_cast<std::size_t>(backcasts), {}, {}};
    for (const auto& x : regressors) result.names.push_back(x.name);
    names_valid(result.names);
    result.values.resize(static_cast<std::size_t>(count) * regressors.size());
    for (std::size_t j = 0; j < regressors.size(); ++j) {
        const auto& x = regressors[j];
        const auto label = "Regressor '" + x.name + "'";
        if (x.domain.frequency() != observed.frequency())
            throw std::invalid_argument(label + " frequency does not match the response.");
        if (x.domain.length() != static_cast<std::int64_t>(x.values.size()))
            throw std::invalid_argument(label + " domain does not match its values.");
        if (start < x.domain.start_ordinal() ||
            start + count > x.domain.start_ordinal() + x.domain.length())
            throw std::invalid_argument(label + " does not cover the complete observed/backcast/forecast domain.");
        const auto offset = static_cast<std::size_t>(subtract(start, x.domain.start_ordinal()));
        for (std::size_t i = 0; i < static_cast<std::size_t>(count); ++i) {
            poll(interrupt, i);
            const double value = x.values[offset + i];
            if (!std::isfinite(value))
                throw std::invalid_argument(label + " contains missing or nonfinite values in the required domain.");
            result.values[i + static_cast<std::size_t>(count) * j] = value;
        }
    }
    return result;
}

PreparedRegression prepare_regression(PreparedRegressors regressors,
    const std::vector<NamedCoefficient>& coefficients, int d, int seasonal_d,
    double mean, bool include_mean, bool include_drift, std::size_t source_first,
    const Interrupt& interrupt) {
    if (d < 0 || d > 2 || seasonal_d < 0 || seasonal_d > 2)
        throw std::invalid_argument("Differencing orders must be whole numbers in 0..2.");
    if (include_mean && d + seasonal_d > 0)
        throw std::invalid_argument("include_mean requires a stationary model; use include_drift for a differenced mean.");
    if (include_drift && d + seasonal_d == 0)
        throw std::invalid_argument("include_drift requires differencing; use include_mean for a stationary model.");
    if (include_mean && mean != 0)
        throw std::invalid_argument("Use mean = 0 with include_mean to avoid two level-mean specifications.");
    const auto rows = regressors.required.length(), n = regressors.observed.length();
    if (rows < 1 || rows > row_limit || n < 1 || n > rows ||
        regressors.observed_offset > static_cast<std::size_t>(rows - n) ||
        regressors.required.frequency() != regressors.observed.frequency() ||
        regressors.required.start_ordinal() + static_cast<std::int64_t>(regressors.observed_offset) != regressors.observed.start_ordinal())
        throw std::invalid_argument("Inconsistent prepared regression domain.");
    names_valid(regressors.names);
    std::vector<std::string> generated;
    if (include_mean) generated.push_back("mean");
    if (include_drift) generated.push_back("drift");
    if (regressors.names.size() + generated.size() > 32)
        throw std::invalid_argument("Regression likelihood supports at most 32 design columns.");
    if (regressors.names.size() + generated.size() > static_cast<std::size_t>(row_limit / rows))
        throw std::length_error("Regressor matrix exceeds the supported element range.");
    if (regressors.values.size() != static_cast<std::size_t>(rows) * regressors.names.size())
        throw std::invalid_argument("Prepared regressor values do not match their dimensions.");
    for (std::size_t i = 0; i < regressors.values.size(); ++i) {
        poll(interrupt, i);
        if (!std::isfinite(regressors.values[i])) throw std::invalid_argument("Regressors must be finite.");
    }
    if (source_first > static_cast<std::size_t>(row_limit))
        throw std::length_error("Source row exceeds the supported range.");
    for (const auto& label : generated) {
        if (std::find(regressors.names.begin(), regressors.names.end(), label) != regressors.names.end())
            throw std::invalid_argument("Generated mean/drift names conflict with supplied regressors.");
        regressors.names.push_back(label);
        const auto old_size = regressors.values.size();
        regressors.values.resize(old_size + static_cast<std::size_t>(rows), 1);
        if (label == "drift") {
            const int order = d + seasonal_d;
            double divisor = 1;
            for (int i = 2; i <= order; ++i) divisor *= i;
            divisor *= kalos::detail::glibc_pow(regressors.observed.frequency(), seasonal_d);
            for (std::size_t i = 0; i < static_cast<std::size_t>(rows); ++i) {
                poll(interrupt, i);
                const double time = static_cast<double>(i) + 1 + static_cast<double>(source_first) -
                    static_cast<double>(regressors.observed_offset);
                const double value = kalos::detail::glibc_pow(time, order) / divisor;
                if (!std::isfinite(value)) throw std::overflow_error("Generated drift exceeds the finite numeric range.");
                regressors.values[old_size + i] = value;
            }
        }
    }
    const auto columns = regressors.names.size();
    RegressionDesign design{columns, std::vector<double>(static_cast<std::size_t>(n) * columns),
        std::vector<double>(columns), std::vector<std::uint8_t>(columns)};
    std::unordered_set<std::string> seen;
    for (const auto& coefficient : coefficients) {
        const auto at = std::find(regressors.names.begin(), regressors.names.end(), coefficient.name);
        if (coefficient.name.empty() || !seen.insert(coefficient.name).second || at == regressors.names.end())
            throw std::invalid_argument("coefficients must have distinct names matching the design columns.");
        if (!std::isfinite(coefficient.value))
            throw std::invalid_argument("coefficients must have finite values or explicit free states.");
        const auto index = static_cast<std::size_t>(at - regressors.names.begin());
        design.fixed[index] = coefficient.fixed;
        design.coefficients[index] = coefficient.fixed ? coefficient.value : 0;
    }
    for (std::size_t j = 0; j < columns; ++j)
        for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
            poll(interrupt, i);
            design.values[i + static_cast<std::size_t>(n) * j] =
                regressors.values[regressors.observed_offset + i + static_cast<std::size_t>(rows) * j];
        }
    return {std::move(regressors), std::move(design), std::move(generated)};
}

LikelihoodAdjustment likelihood_adjustment(const MaskedSeries& jacobian,
    Transform transformation, std::size_t degree, std::size_t source_first,
    const Interrupt& interrupt) {
    jacobian.validate();
    if (transformation != Transform::none && transformation != Transform::logarithm)
        throw std::invalid_argument("Unknown transformation.");
    if (source_first > std::numeric_limits<std::size_t>::max() - jacobian.values.size())
        throw std::overflow_error("Source row overflow.");
    LikelihoodAdjustment result{{}, {}, std::vector<double>(jacobian.values.size()), 0};
    std::optional<std::size_t> first;
    long double sum = 0;
    for (std::size_t i = 0; i < jacobian.values.size(); ++i) {
        poll(interrupt, i);
        if (jacobian.missing[i]) continue;
        if (!first) first = i;
        if (transformation == Transform::logarithm && i - *first < degree) continue;
        result.local_rows.push_back(i); result.source_rows.push_back(source_first + i);
        result.terms[i] = jacobian.values[i]; sum += jacobian.values[i];
    }
    result.log_jacobian = static_cast<double>(sum);
    if (!std::isfinite(result.log_jacobian)) throw std::overflow_error("The log-Jacobian correction is nonfinite.");
    return result;
}
}
