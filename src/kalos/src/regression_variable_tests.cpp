#include "kalos/regression_variable_tests.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace kalos {

void RegressionTestStatistics::validate() const {
    const std::size_t nx = coefficients.size();
    if (!std::all_of(coefficients.begin(), coefficients.end(),
            [](double v) { return std::isfinite(v); })) {
        throw std::invalid_argument(
            "RegressionTestStatistics: the coefficients must be finite");
    }
    if (unscaled_covariance.size() != nx * nx) {
        throw std::invalid_argument(
            "RegressionTestStatistics: unscaled_covariance must hold "
            "nx * nx values in row-major order");
    }
    if (!std::all_of(unscaled_covariance.begin(), unscaled_covariance.end(),
            [](double v) { return std::isfinite(v); })) {
        throw std::invalid_argument(
            "RegressionTestStatistics: the unscaled covariance must be finite");
    }
    if (!std::isfinite(ssq) || ssq < 0.0) {
        throw std::invalid_argument(
            "RegressionTestStatistics: ssq must be finite and non-negative");
    }
    if (dim < 1) {
        throw std::invalid_argument(
            "RegressionTestStatistics: dim (effective observations) must be "
            "positive");
    }
}

std::vector<double> RegressionTestStatistics::tstats() const {
    // ConcentratedLikelihood.ser(int, boolean) then tstats(int, boolean):
    // ser[j] = sqrt((uvarDiag[j] * ssq) / dof), t[j] = coef[j] / ser[j],
    // with dof = dim for ml = false and no zero guard on the vector form.
    const std::size_t nx = coefficients.size();
    const double dof = static_cast<double>(dim);
    std::vector<double> t(nx);
    for (std::size_t j = 0; j < nx; ++j) {
        const double u = unscaled_covariance[j * nx + j];
        const double ser = std::sqrt((u * ssq) / dof);
        t[j] = coefficients[j] / ser;
    }
    return t;
}

double RegressionTestStatistics::mean_tstat() const {
    if (unscaled_covariance.empty()) {
        throw std::invalid_argument(
            "RegressionTestStatistics: tstat(0, 0, false) needs a mean "
            "column");
    }
    // ConcentratedLikelihood.tstat(i, nobs, ml) / ser(i, nobs, ml) scalar
    // forms with i = 0, nobs = 0, ml = false (dof = dim), including the
    // pinned guards: a zero covariance entry makes ser() NaN, a zero
    // coefficient makes ser() 0 so that tstat() divides 0 by 0.
    const double u = unscaled_covariance[0];
    if (u == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (coefficients[0] == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double ser = std::sqrt((u * ssq) / static_cast<double>(dim));
    return coefficients[0] / ser;
}

double derived_test_tstat(const RegressionTestStatistics& statistics,
                          int first, int dim) {
    // DerivedRegressionTest.accept: vcur = -extract(i0, dim).sum() (a plain
    // left fold from 0.0, DoublesMath.sum = reduce(0.0, Double::sum)); the
    // variance block is extract(i0, dim, i0, dim) whose dense copy is
    // folded row-major over the storage by FastMatrix.sum (isFull branch).
    const std::size_t nx = statistics.coefficients.size();
    double sum = 0.0;
    for (int j = 0; j < dim; ++j) {
        sum += statistics.coefficients[static_cast<std::size_t>(first) + j];
    }
    const double vcur = -sum;
    double block = 0.0;
    for (int r = 0; r < dim; ++r) {
        const std::size_t row =
            (static_cast<std::size_t>(first) + static_cast<std::size_t>(r)) * nx;
        for (int c = 0; c < dim; ++c) {
            block += statistics.unscaled_covariance[row
                                                    + static_cast<std::size_t>(first)
                                                    + static_cast<std::size_t>(c)];
        }
    }
    // ml = true at every live construction -> dof = cl.dim().
    const double var = (block * statistics.ssq)
                       / static_cast<double>(statistics.dim);
    return vcur / std::sqrt(var);
}

bool t_test_accepts(const std::vector<double>& tstats, double tlow,
                    double thigh, int nsig, int first, int dim) {
    // TRegressionTest.accept: high counts |t| >= thigh, the else-if counts
    // the remaining |t| >= tlow; NaN counts in neither. accept =
    // high > 0 || low >= min(nsig, dim).
    int low = 0;
    int high = 0;
    for (int j = 0; j < dim; ++j) {
        const double a = std::fabs(tstats[static_cast<std::size_t>(first) + j]);
        if (a >= thigh) {
            ++high;
        } else if (a >= tlow) {
            ++low;
        }
    }
    return high > 0 || low >= std::min(nsig, dim);
}

void RegressionVariableTestRequest::validate() const {
    statistics.validate();
    if (mean && statistics.coefficients.empty()) {
        // Mean on implies design column 0 (the bytecode would read
        // tvals[-1]; a live description never gets here unestimated).
        throw std::invalid_argument(
            "RegressionVariableTestRequest: a mean-on description has at "
            "least the mean coefficient");
    }
    if (td_joint != 0.0) {
        throw std::invalid_argument(
            "RegressionVariableTestRequest: the joint TD test "
            "(RegressionVariablesTest.Builder.tdJointTest, FRegressionTest) "
            "is dead for every pinned x13 call site and is not qualified");
    }
    if (!(mean_test >= 0.0) || !(td_test >= 0.0) || !(mh_test >= 0.0)
        || !std::isfinite(mean_test) || !std::isfinite(td_test)
        || !std::isfinite(mh_test)) {
        throw std::invalid_argument(
            "RegressionVariableTestRequest: the test cutoffs must be finite "
            "and non-negative (0 disables the test)");
    }
    const std::size_t nx = statistics.coefficients.size();
    for (const auto& variable : variables) {
        if (variable.dim < 1) {
            throw std::invalid_argument(
                "RegressionVariableTestRequest: every variable dim must be "
                "positive");
        }
        if (variable.position < 0
            || static_cast<std::size_t>(variable.position)
                       + static_cast<std::size_t>(variable.dim)
                   > nx) {
            throw std::invalid_argument(
                "RegressionVariableTestRequest: the variable layout must "
                "fit the design columns [0, nx)");
        }
    }
    if (preadjustment < 0 || preadjustment > 2) {
        throw std::invalid_argument(
            "RegressionVariableTestRequest: preadjustment must be 0 (None), "
            "1 (LeapYear) or 2 (GoodFriday)");
    }
}

RegressionVariableTestResult run_regression_variable_tests(
    const RegressionVariableTestRequest& request) {
    request.validate();
    const std::size_t nvars = request.variables.size();
    RegressionVariableTestResult result;
    result.tstats = request.statistics.tstats();
    result.removed.assign(nvars, false);
    result.tested.assign(nvars, false);
    result.accepted.assign(nvars, false);
    result.derived_tested.assign(nvars, false);
    result.derived_tstats.assign(nvars,
                                 std::numeric_limits<double>::quiet_NaN());

    bool changed = false;
    std::vector<std::size_t> removals; // TD + LoP candidates (deferred)
    bool b = false;                    // a TD variable survived
    bool c = false;                    // the LoP variable passed
    bool mean = request.mean;
    int preadjustment = request.preadjustment;

    if (request.td_test > 0.0) {
        // Pass 1, ModellingUtility::isTradingDays: a failed variable gets a
        // second life through the derived joint t-statistic when dim > 1;
        // a survivor (either way) sets b.
        for (std::size_t i = 0; i < nvars; ++i) {
            const auto& variable = request.variables[i];
            if (variable.kind != RegressionVariableKind::trading_days) {
                continue;
            }
            result.tested[i] = true;
            if (t_test_accepts(result.tstats, request.td_test, request.td_test,
                    1, variable.position, variable.dim)) {
                result.accepted[i] = true;
                b = true;
            } else if (variable.dim > 1 && request.derived_test) {
                const double tcur = derived_test_tstat(request.statistics,
                    variable.position, variable.dim);
                result.derived_tested[i] = true;
                result.derived_tstats[i] = tcur;
                if (std::fabs(tcur) > request.td_test) {
                    b = true;
                } else {
                    removals.push_back(i);
                }
            } else {
                removals.push_back(i);
            }
        }
        // Pass 2, isLengthOfPeriod: tested with width 1; a pass sets c.
        for (std::size_t i = 0; i < nvars; ++i) {
            const auto& variable = request.variables[i];
            if (variable.kind != RegressionVariableKind::length_of_period) {
                continue;
            }
            result.tested[i] = true;
            if (t_test_accepts(result.tstats, request.td_test, request.td_test,
                    1, variable.position, 1)) {
                result.accepted[i] = true;
                c = true;
            } else {
                removals.push_back(i);
            }
        }
    }
    if (request.mh_test > 0.0) {
        // Pass 3, isMovingHoliday: failures mark changed and are removed
        // immediately, independently of the TD/LoP gate.
        std::vector<std::size_t> failures;
        for (std::size_t i = 0; i < nvars; ++i) {
            const auto& variable = request.variables[i];
            if (variable.kind != RegressionVariableKind::moving_holiday) {
                continue;
            }
            result.tested[i] = true;
            if (t_test_accepts(result.tstats, request.mh_test, request.mh_test,
                    1, variable.position, 1)) {
                result.accepted[i] = true;
            } else {
                failures.push_back(i);
                changed = true;
            }
        }
        for (const std::size_t i : failures) {
            result.removed[i] = true;
        }
    }
    if (!removals.empty() && !c && !b) {
        // All-fail-or-none gate: the TD/LoP removals fire as a batch only
        // when no tested TD variable survived and the LoP variable passed
        // neither; the preadjustment is reset only here.
        changed = true;
        for (const std::size_t i : removals) {
            result.removed[i] = true;
        }
        if (request.adjusted) {
            preadjustment = 0; // setPreadjustment(LengthOfPeriodType.None)
        }
    }
    if (request.mean_test > 0.0 && mean) {
        if (!t_test_accepts(result.tstats, request.mean_test,
                request.mean_test, 1, 0, 1)) {
            mean = false;
            changed = true;
        }
    }
    result.changed = changed;
    result.mean = mean;
    result.preadjustment = preadjustment;
    return result;
}

void MeanControllerRequest::validate() const {
    statistics.validate();
    if (statistics.coefficients.empty()) {
        throw std::invalid_argument(
            "MeanControllerRequest: the consumed estimation carries the "
            "mean column (coefficient 0 is read by tstat(0, 0, false))");
    }
    if (!(cval > 0.0) || !std::isfinite(cval)) {
        throw std::invalid_argument(
            "MeanControllerRequest: cval must be finite and positive (the "
            "pinned call sites use 1.96 and 1.6)");
    }
}

MeanControllerResult run_mean_controller(const MeanControllerRequest& request) {
    request.validate();
    MeanControllerResult result;
    result.tstat = request.statistics.mean_tstat();
    // sig = |t| > cval (dcmpl ifle: strict, NaN never significant); a
    // decision equal to the current flag is Unchanged without touching the
    // description, otherwise setMean(sig) + clearEstimation.
    const bool sig = std::fabs(result.tstat) > request.cval;
    result.changed = sig != request.mean;
    result.mean = result.changed ? sig : request.mean;
    return result;
}

bool final_mean_rule_adds_mean(bool auto_modelling, bool mean, bool check_mu,
                               double rtval) {
    // Strict rtval0 > 2.5 with the pinned dcmpl ordering (NaN fails).
    return auto_modelling && !mean && check_mu && rtval > kMeanFinalCutoff;
}

} // namespace kalos
