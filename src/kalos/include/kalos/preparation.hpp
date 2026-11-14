#ifndef KALOS_PREPARATION_HPP
#define KALOS_PREPARATION_HPP
#include "numerics.hpp"
#include "regression_likelihood.hpp"
#include "transforms.hpp"
#include <optional>
#include <string>

namespace kalos {
struct EstimationWindow {
    TimeDomain domain;
    std::size_t first; // Zero-based position in the original response.
};
// Inclusive endpoints, integer period ordinals. Selection carries no prior state.
EstimationWindow estimation_window(const TimeDomain& source,
    std::optional<std::int64_t> first = {}, std::optional<std::int64_t> last = {});

struct NamedRegressor {
    std::string name;
    TimeDomain domain;
    std::vector<double> values; // Only the requested slice must be finite.
};
struct PreparedRegressors {
    TimeDomain observed, required;
    std::size_t observed_offset;
    std::vector<std::string> names;
    std::vector<double> values; // Column-major over the required domain.
};
PreparedRegressors align_regressors(const TimeDomain& observed,
    const std::vector<NamedRegressor>& regressors,
    std::int64_t backcasts = 0, std::int64_t forecasts = 0,
    const Interrupt& interrupt = {});

struct NamedCoefficient { std::string name; double value; bool fixed; };
struct PreparedRegression {
    PreparedRegressors regressors;
    RegressionDesign design; // Observed rows only, same column ordering.
    std::vector<std::string> generated;
};
// source_first is zero based in the unsliced response. Mean/drift use the same
// origin over observed, backcast and forecast rows and across selected windows.
PreparedRegression prepare_regression(PreparedRegressors regressors,
    const std::vector<NamedCoefficient>& coefficients, int d, int seasonal_d,
    double mean = 0, bool include_mean = false, bool include_drift = false,
    std::size_t source_first = 0, const Interrupt& interrupt = {});

struct LikelihoodAdjustment {
    std::vector<std::size_t> local_rows, source_rows;
    std::vector<double> terms;
    double log_jacobian;
};
// Pinned JD+ calendar-window convention, not a universal improper likelihood.
LikelihoodAdjustment likelihood_adjustment(const MaskedSeries& log_jacobian,
    Transform transformation, std::size_t degree, std::size_t source_first = 0,
    const Interrupt& interrupt = {});
}
#endif
