// Implementation of the TRAMO automatic-fit assembly (item 23-F F.4): the
// TramoKernel.process composition of the landed units - the 23-A builder
// mirror and the 23-B transform mirror run LIVE and are cross-checked
// against the captured pre-stage facts, the F.3 round driver runs verbatim,
// the differencing/arma seams re-run their landed stage units live where
// the stage envelope qualifies (header "witness tiers"), and the final
// surface recomputes the D9 regression-block covariance
// V = unscaledCovariance * (ssq / (dim - nx - k)) - the
// ConcentratedLikelihood.covariance(k, true) chain, javap-pinned
// (RegSarimaModel.of passes getArimaSpec().freeParametersCount() and
// iconst_1 for the scaling flag; the chain computes s = ssq / (double)
// (dim - nx - k) once, then FastMatrix.mul(s) elementwise; nx == 0 answers
// FastMatrix.EMPTY) - against the captured coefficientsCovariance upper
// bits. Decoded tables: tmp_opencode/kalos-m5-pinning.md SECTION 32 with
// the header-dated amendments; corpus:
// validation/evidence/native/tramo-fit (Scratch23f, fit-capture-1/2/witness
// byte-identical, 34 identity all-1 rows + the pinned fit_const crash
// witness, which enters no test case).

#include <kalos/tramo_fit.hpp>

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/outlier_scoring.hpp"
#include "mlj_engine.hpp"

#include <kalos/tramo_arma_selection.hpp>
#include <kalos/tramo_differencing_decision.hpp>
#include <kalos/tramo_final_estimation.hpp>
#include <kalos/tramo_outlier_rounds.hpp>
#include <kalos/tramo_regression_tests.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace kalos {
namespace {

void require(bool condition, const std::string& what) {
    if (!condition)
        throw std::invalid_argument("tramo_fit: " + what);
}

bool bit_equal(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

std::vector<double> differenced_response(std::vector<double> y, bool log, int period, int d, int bd) {
    if (log) {
        for (double& v : y) {
            require(v > 0.0, "non-positive value under the log channel");
            v = detail::fdlibm_log(v);
        }
    }
    auto filter = mlj::differencing_filter(period, d, bd);
    if (filter.size() <= 1) return y;
    std::reverse(filter.begin(), filter.end());
    require(y.size() >= filter.size(), "series too short for the differencing filter");
    std::vector<double> result(y.size() - filter.size() + 1);
    mlj::apply_backfilter(filter, y, 0, result);
    return result;
}

std::vector<double> arma_residuals(const tramo_build_result& built, bool log,
                                    const TramoFitSeamInputs& in, double delta_tc) {
    const int period = built.annual_frequency;
    const int d = in.entry_orders.d, bd = in.entry_orders.bd;
    auto y = differenced_response(built.interpolated, log, period, d, bd);
    std::vector<std::vector<double>> columns;
    if (in.entry_mean) columns.emplace_back(y.size(), 1.0);
    require(in.entry_variables.size() == in.entry_outliers.size(), "incomplete ARMA entry design");
    for (std::size_t j = 0; j < in.entry_variables.size(); ++j) {
        const auto& var = in.entry_variables[j];
        require(var.name == in.entry_outliers[j] && var.position >= 0 &&
                    static_cast<std::size_t>(var.position) < built.interpolated.size(),
                "invalid ARMA entry outlier name or position");
        std::vector<double> column(built.interpolated.size());
        const auto type = static_cast<OutlierType>(var.type);
        fill_outlier_column(type, var.position, period, delta_tc, column);
        require(var.column.empty() || var.column == column, "ARMA outlier column disagrees with factory");
        columns.push_back(differenced_response(std::move(column), false, period, d, bd));
    }
    if (columns.empty()) return y;
    const auto beta = detail::pinned_ols_coefficients(y, columns);
    require(beta && beta->size() == columns.size(), "ARMA input regression failed");
    for (std::size_t j = 0; j < columns.size(); ++j)
        for (std::size_t i = 0; i < y.size(); ++i)
            y[i] += -(*beta)[j] * columns[j][i];
    return y;
}

std::string stage_word(TramoKernelStageResult r) {
    switch (r) {
        case TramoKernelStageResult::changed: return "Changed";
        case TramoKernelStageResult::unchanged: return "Unchanged";
        case TramoKernelStageResult::unprocessed: return "Unprocessed";
        case TramoKernelStageResult::failed: return "Failed";
    }
    return "?";
}

std::string orders_word(const TramoOrders& o) {
    return std::to_string(o.p) + "," + std::to_string(o.d) + "," +
           std::to_string(o.q) + "," + std::to_string(o.bp) + "," +
           std::to_string(o.bd) + "," + std::to_string(o.bq);
}

bool orders_equal(const TramoOrders& a, const TramoOrders& b) {
    return orders_word(a) == orders_word(b);
}

struct DetectedOutlier {
    OutlierType type;
    int position;
    std::string name;
};

struct OutlierScanSelection {
    TramoKernelStageResult outcome;
    std::vector<std::string> names;
    std::vector<DetectedOutlier> variables;
};

std::optional<OutlierScanSelection> run_scan(const TramoFitSpec& spec, const tramo_build_result& built,
                                             bool log, const TramoOrders& o, bool mean, double critical_value) {
    TramoOutlierScanRequest request{};
    request.y = built.interpolated;
    request.frequency = built.annual_frequency;
    request.mean = mean;
    request.log = log;
    request.p = o.p;
    request.d = o.d;
    request.q = o.q;
    request.bp = o.bp;
    request.bd = o.bd;
    request.bq = o.bq;
    request.ao = spec.ao;
    request.ls = spec.ls;
    request.tc = spec.tc;
    request.so = spec.so;
    request.tc_rate = spec.delta_tc;
    request.precision = spec.intermediate_precision;
    request.maximum_likelihood = spec.outliers_ml;
    // The pinned builder stores its two default caps in crossed fields.
    request.max_round = 30;
    request.max_outliers = 50;
    request.cv = critical_value;
    request.upper_bound = built.estimation_length;
    if (o.d == 0 && o.bd == 0 && !tramo_stationary_scan_supported(request)) return std::nullopt;
    const auto scan = tramo_outlier_rounds(request);
    const auto outcome = scan.status == TramoOutlierRoundsResult::Status::Failed ? TramoKernelStageResult::failed :
        scan.outliers.empty() ? TramoKernelStageResult::unchanged : TramoKernelStageResult::changed;
    std::vector<OutlierType> factories;
    if (spec.ao) factories.push_back(OutlierType::Additive);
    if (spec.ls) factories.push_back(OutlierType::LevelShift);
    if (spec.tc) factories.push_back(OutlierType::TransitoryChange);
    if (spec.so) factories.push_back(OutlierType::SeasonalOutlier);
    std::vector<std::string> names;
    std::vector<DetectedOutlier> variables;
    for (const auto& selected : scan.outliers) {
        const auto type = factories.at(selected[1]);
        names.push_back(monthly_outlier_name(type, selected[0], built.domain_year, built.domain_position + 1));
        variables.push_back({type, selected[0], names.back()});
    }
    // ModelDescription.sortVariables stably orders automatic outliers by date
    // before estimating the regression. The scan witness uses detection order.
    std::stable_sort(variables.begin(), variables.end(),
                     [](const auto& a, const auto& b) { return a.position < b.position; });
    return OutlierScanSelection{outcome, std::move(names), std::move(variables)};
}

std::optional<std::vector<DetectedOutlier>> check_scan(const TramoFitSpec& spec, const tramo_build_result& built, bool log,
                const TramoFitSeamInputs& in, const TramoKernelSeam& expected) {
    auto scan = run_scan(spec, built, log, in.entry_orders, in.entry_mean, expected.critical_value);
    if (!scan) return std::nullopt;
    require(scan->outcome == expected.result, "live scan outcome disagrees with the module witness");
    require(scan->names == expected.outliers, "live scan selection disagrees with the ordered outlier witness");
    require(orders_equal(in.entry_orders, expected.orders) && expected.mean == in.entry_mean,
            "live scan description disagrees with the module witness");
    return std::move(scan->variables);
}

struct EstimationEntry {
    TramoOrders orders;
    bool mean;
    std::vector<DetectedOutlier> variables;
};

struct FittedOutlierEffects {
    TramoOrders orders;
    bool mean;
    std::vector<double> parameters;
    std::vector<std::string> names;
    std::vector<double> corrected_response;
    TramoFinalEstimatorRunResult fit;
    std::vector<TramoFitVariable> variables;
};

TramoTestStatistics fitted_test_statistics(const FittedOutlierEffects& effects) {
    const auto& fit = effects.fit;
    TramoTestStatistics test;
    test.coefficients = fit.cl_coefs;
    const auto nx = test.coefficients.size();
    require(fit.cl_unscaled_covariance.size() == nx * nx, "invalid fitted regression covariance");
    for (std::size_t i = 0; i < nx; ++i)
        test.cov_diagonal.push_back(fit.cl_unscaled_covariance[i * nx + i]);
    test.ssq = fit.cl_ssq;
    test.dim = fit.cl_dim;
    test.free_parameters = static_cast<int>(effects.parameters.size());
    return test;
}

void check_round_statistics(const FittedOutlierEffects& effects, const TramoKernelStatsRow& row,
                            int frequency, int n, std::size_t seam) {
    auto test = fitted_test_statistics(effects);
    TramoModelStatisticsRequest request;
    request.residuals = effects.fit.cl_residuals;
    request.coefficients = std::move(test.coefficients);
    request.cov_diagonal = std::move(test.cov_diagonal);
    request.ssq = test.ssq;
    request.dim = test.dim;
    request.sigma2 = test.ssq / test.dim;
    request.degrees_of_freedom = test.dim - static_cast<int>(request.coefficients.size());
    request.annual_frequency = frequency;
    request.estimation_domain = n;
    request.ami_outliers = static_cast<int>(effects.names.size());
    request.free_parameters = test.free_parameters;
    const auto stats = tramo_model_statistics(request);
    const std::string stage = " at kernel seam " + std::to_string(seam);
    require(std::equal(row.residuals.begin(), row.residuals.end(), request.residuals.begin(),
                       request.residuals.end(), bit_equal), "fitted likelihood residuals disagree" + stage);
    require(stats.built && row.dim == request.dim && row.free_parameters == request.free_parameters &&
            row.observations == stats.observations_count && row.outliers_count == stats.outliers_count &&
            row.effective_observations == stats.effective_observations_count, "fitted statistics counts disagree" + stage);
    require(bit_equal(row.ljung_box, stats.ljung_box) && bit_equal(row.ljung_box_pvalue, stats.ljung_box_pvalue),
            "fitted Ljung-Box statistics disagree" + stage);
    require(bit_equal(row.bic, stats.bic), "fitted BIC disagrees" + stage);
    require(bit_equal(row.se, stats.se), "fitted standard error disagrees" + stage);
}

bool check_regression_test(const FittedOutlierEffects& effects, bool test_mean,
                           const TramoKernelSeam& seam) {
    TramoRegressionTestRequest request;
    request.mean = effects.mean;
    request.test_mean = test_mean;
    request.mean_threshold = seam.tsig;
    request.statistics = fitted_test_statistics(effects);
    const auto decision = run_tramo_regression_tests(request);
    const auto stage = decision.estimation_cleared ? TramoKernelStageResult::changed
                                                   : TramoKernelStageResult::unchanged;
    require(seam.result == stage && seam.mean == decision.mean, "fitted regression test decision disagrees");
    // The probe asks for tstats after the test: clearing the model also clears its likelihood.
    const auto t = decision.estimation_cleared ? std::vector<double>{}
                                               : tramo_test_tstats(request.statistics, true);
    require(std::equal(t.begin(), t.end(), seam.tstats.begin(), seam.tstats.end(), bit_equal),
            "fitted regression t-statistics disagree");
    return decision.estimation_cleared;
}

FittedOutlierEffects fit_estimation(const TramoFitSpec& spec, const std::vector<double>& y,
                                   EstimationEntry& entry, int pass) {
    TramoFinalEstimatorRunRequest options;
    const auto& o = entry.orders;
    options.orders = {o.p, o.d, o.q, o.bp, o.bd, o.bq};
    options.mean = entry.mean;
    options.eps = spec.precision;
    options.ur = spec.final_unit_root_threshold;
    options.pass = pass;
    options.ami = spec.using_auto_model;
    options.outliers = spec.outliers_used;
    std::vector<std::vector<double>> columns;
    std::vector<std::string> names;
    for (const auto& var : entry.variables) {
        names.push_back(var.name);
        columns.emplace_back(y.size());
        fill_outlier_column(var.type, var.position, spec.annual_frequency, spec.delta_tc, columns.back());
    }
    auto estimate = tramo_final_estimator_run(options, y, columns);
    const auto& t = estimate.orders_after;
    entry.orders = {t.p, t.d, t.q, t.bp, t.bd, t.bq};
    entry.mean = estimate.mean_out;
    FittedOutlierEffects effects{entry.orders, entry.mean, estimate.parameters_after, std::move(names), {}, {}, {}};
    if (!columns.empty() && estimate.cl_coefs.size() == columns.size() + static_cast<std::size_t>(entry.mean)) {
        // RegArimaUtility.regressionEffect accumulates the effects before
        // the differencing module subtracts their sum from the model response.
        std::vector<double> total(y.size(), 0.0);
        for (std::size_t j = 0; j < columns.size(); ++j)
            for (std::size_t i = 0; i < y.size(); ++i)
                total[i] += estimate.cl_coefs[j + entry.mean] * columns[j][i];
        effects.corrected_response = y;
        for (std::size_t i = 0; i < y.size(); ++i) effects.corrected_response[i] -= total[i];
    }
    if (estimate.cl_state != 0) {
        for (std::size_t i = 0; i < columns.size(); ++i) {
            const auto& var = entry.variables[i];
            effects.variables.push_back({var.name, static_cast<TramoOutlierKind>(var.type), var.position,
                estimate.cl_coefs.at(i + entry.mean), std::move(columns[i])});
        }
    }
    effects.fit = std::move(estimate);
    return effects;
}

FittedOutlierEffects check_final_estimation(const TramoFitSpec& spec, const std::vector<double>& y,
                                            EstimationEntry& entry, const TramoKernelSeam& expected) {
    auto effects = fit_estimation(spec, y, entry, expected.pass_at_call);
    const auto& estimate = effects.fit;
    require(estimate.returned == (expected.returned != 0) && estimate.nnsig == expected.nnsig,
            "live final estimation decision disagrees with the module witness");
    require(orders_equal(entry.orders, expected.orders) &&
            effects.names == expected.outliers, "live final estimation specification disagrees with the witness");
    require(std::equal(estimate.parameters_after.begin(), estimate.parameters_after.end(),
                       expected.parameters.begin(), expected.parameters.end(), bit_equal),
            "live final estimation parameters disagree with the module witness");
    return effects;
}

TramoRegressionTestDecision initial_mean_test(const TramoFitSpec& spec, const std::vector<double>& y,
                                              const TramoOrders& orders, bool mean) {
    TramoRegressionTestRequest test;
    test.kind = TramoRegressionTestKind::def;
    test.test_mean = tramo_ami_check_mu(spec.using_auto_model, spec.mean_is_test);
    test.mean = mean;
    if (test.test_mean && mean) {
        // DefaultRegressionTest estimates its own copy with the intermediate
        // precision. The description is still in its post-seasonality orders.
        auto options = spec;
        options.precision = spec.intermediate_precision;
        options.using_auto_model = false;
        EstimationEntry trial{orders, mean, {}};
        auto fit = fit_estimation(options, y, trial, 0);
        require(fit.fit.returned && fit.fit.cl_state != 0, "initial mean-test likelihood unavailable");
        test.statistics = fitted_test_statistics(fit);
    }
    return run_tramo_regression_tests(test);
}

// Compute the first reference model from the dated description. Only
// estimateModel(pass=0) repeats the fast mean test and the model fit.
std::optional<FittedOutlierEffects> fit_reference(const TramoFitSpec& spec, const tramo_build_result& built,
                                                   bool log, const std::vector<double>& y, EstimationEntry entry) {
    if (spec.outliers_used) {
        const double cv = spec.critical_value == 0.0 ? tramo_critical_value(built.estimation_length)
                                                    : spec.critical_value;
        auto selection = run_scan(spec, built, log, entry.orders, entry.mean, cv);
        if (!selection) return std::nullopt;
        entry.variables = std::move(selection->variables);
    }
    for (int retry = 0; ; ++retry) {
        auto fitted = fit_estimation(spec, y, entry, 0);
        if (!fitted.fit.returned || fitted.fit.cl_state == 0) return std::nullopt;
        if (retry == 5) return fitted;
        TramoRegressionTestRequest test;
        test.mean = fitted.mean;
        test.test_mean = tramo_ami_check_mu(spec.using_auto_model, spec.mean_is_test);
        test.mean_threshold = 1.96;
        test.statistics = fitted_test_statistics(fitted);
        const auto decision = run_tramo_regression_tests(test);
        if (!decision.estimation_cleared) return fitted;
        entry.mean = decision.mean;
    }
}

TramoDifferencingDecision difference_orders(const TramoFitSpec& spec, const std::vector<double>& response,
                                            bool log, const TramoOrders& orders, bool mean, bool initial_round,
                                            bool seasonal, int nobs) {
    TramoDifferencingRequest dreq;
    dreq.series = response;
    dreq.log_input = log;
    dreq.period = spec.annual_frequency;
    dreq.d0 = orders.d;
    dreq.bd0 = orders.bd;
    dreq.p0 = orders.p;
    dreq.q0 = orders.q;
    dreq.bp0 = orders.bp;
    dreq.bq0 = orders.bq;
    dreq.mean0 = mean;
    dreq.initial_round = initial_round;
    dreq.eps = spec.intermediate_precision;
    dreq.seasonal = seasonal;
    dreq.ub1 = spec.ub1;
    dreq.ub2 = spec.ub2;
    dreq.cancel = spec.cancel;
    dreq.nobs = nobs;
    return decide_tramo_differencing(dreq);
}

TramoArmaSelectionResult select_orders(TramoArmaSelectionRequest request, const TramoFitSpec& spec,
                                       const TramoOrders& orders, bool seasonal) {
    request.period = spec.annual_frequency;
    request.p0 = orders.p;
    request.d0 = orders.d;
    request.q0 = orders.q;
    request.bp0 = orders.bp;
    request.bd0 = orders.bd;
    request.bq0 = orders.bq;
    request.seasonal = seasonal;
    request.accept_white_noise = spec.accept_airline;
    return select_tramo_arma(request);
}

// The first automatic candidate (pass=1, round=1), before the refStats
// decision. Both the reference and the new design come from the data.
std::optional<FittedOutlierEffects> fit_first_candidate(const TramoFitSpec& spec, const tramo_build_result& built,
        bool log, bool seasonal, const FittedOutlierEffects& reference) {
    const auto& original = reference.orders;
    TramoFitSeamInputs design;
    design.entry_orders = original;
    design.entry_mean = reference.mean;
    design.entry_variables = reference.variables;
    design.entry_outliers = reference.names;
    if (original.d < 2 || original.bd >= 1) {
        const bool raw = reference.names.empty();
        if (!raw && reference.corrected_response.empty()) return std::nullopt;
        const auto d = difference_orders(spec, raw ? built.interpolated : reference.corrected_response,
                                         raw && log, original, reference.mean, true, seasonal, built.estimation_length);
        design.entry_orders = {d.spec_p, d.spec_d, d.spec_q, d.spec_bp, d.spec_bd, d.spec_bq};
        design.entry_mean = d.mean;
        if (d.spec_d + built.annual_frequency * d.spec_bd != original.d + built.annual_frequency * original.bd) {
            design.entry_variables.clear();
            design.entry_outliers.clear();
        }
    }
    TramoArmaSelectionRequest arma;
    arma.residuals = arma_residuals(built, log, design, spec.delta_tc);
    arma.nobs = built.estimation_length;
    const auto selected = select_orders(std::move(arma), spec, design.entry_orders, seasonal);
    EstimationEntry entry{{selected.spec_p, selected.spec_d, selected.spec_q,
                           selected.spec_bp, selected.spec_bd, selected.spec_bq}, design.entry_mean, {}};
    if (spec.outliers_used) {
        const double cv = spec.critical_value == 0.0 ? tramo_critical_value(built.estimation_length)
                                                    : spec.critical_value;
        auto scan = run_scan(spec, built, log, entry.orders, entry.mean, cv);
        if (!scan) return std::nullopt;
        entry.variables = std::move(scan->variables);
    }
    const auto y = differenced_response(built.interpolated, log, built.annual_frequency, 0, 0);
    auto fitted = fit_estimation(spec, y, entry, 1);
    if (!fitted.fit.returned || fitted.fit.cl_state == 0) return std::nullopt;
    return fitted;
}

TramoFitFinalRecord project_fitted_record(const FittedOutlierEffects& effects,
                                         const std::vector<double>& observations, bool log, int frequency) {
    const auto& fit = effects.fit;
    const auto& o = effects.orders;
    TramoFitFinalRecord r;
    r.log = log;
    r.mean = effects.mean;
    r.orders = o;
    r.parameters = effects.parameters;
    r.outlier_names = effects.names;
    r.variables = effects.variables;
    r.coefficients = fit.cl_coefs;
    r.nx = static_cast<int>(r.coefficients.size());
    r.dim = fit.cl_dim;
    r.ssq = r.ls_ssq = fit.cl_ssq;
    r.log_determinant = fit.cl_log_determinant;
    r.unscaled_covariance = fit.cl_unscaled_covariance;
    r.free_arma_parameters = static_cast<int>(r.parameters.size());
    r.scaling_factor = true;
    const auto level = differenced_response(observations, log, frequency, 0, 0);
    r.differenced_response = differenced_response(level, false, frequency, o.d, o.bd);
    require(static_cast<int>(r.differenced_response.size()) == r.dim &&
            r.nx == static_cast<int>(r.variables.size()) + r.mean, "invalid fitted model dimensions");

    // RegArimaUtility.fullResiduals subtracts each fitted regressor on the
    // level scale, differences once, and removes the drift on that domain.
    if (r.nx == 0) {
        r.full_residuals = fit.cl_residuals;
    } else {
        auto corrected = level;
        for (const auto& variable : r.variables)
            for (std::size_t i = 0; i < level.size(); ++i)
                corrected[i] += -variable.coefficient * variable.column[i];
        auto dy = differenced_response(std::move(corrected), false, frequency, o.d, o.bd);
        if (r.mean)
            for (double& v : dy) v -= r.coefficients.front();
        auto residuals = detail::pinned_fast_kalman_residuals(dy, o.p, o.bp, o.q, o.bq, frequency, r.parameters);
        require(residuals && residuals->size() == dy.size(), "final residual filtering failed");
        r.full_residuals = std::move(*residuals);
    }
    r.ls_obs = static_cast<int>(observations.size());
    r.ls_eff = r.dim;
    r.ls_npar = 1 + r.nx + r.free_arma_parameters;
    r.ls_ll = fit.cl_log_likelihood;
    if (log) {
        double sum = 0.0;
        for (std::size_t i = level.size() - r.dim; i < level.size(); ++i) sum += level[i];
        r.ls_tadj = -sum;
    }
    r.ls_all = r.ls_ll + r.ls_tadj;
    const double n = r.dim, k = r.ls_npar, ln = detail::fdlibm_log(n);
    r.ls_aic = -2.0 * (r.ls_all - k);
    r.ls_aicc = -2.0 * (r.ls_all - n * k / (n - k - 1.0));
    r.ls_bic = -2.0 * r.ls_all + k * ln;
    r.ls_bic2 = (-2.0 * r.ls_ll + k * ln) / n;
    r.ls_bicc = detail::fdlibm_log(r.ssq / n) + (k - 1) * ln / n;
    r.ls_hq = -2.0 * (r.ls_all - k * detail::fdlibm_log(ln));
    return r;
}

void check_final_record(const TramoFitFinalRecord& got, const TramoFitFinalRecord& ref) {
    require(orders_equal(got.orders, ref.orders) && got.mean == ref.mean && got.log == ref.log &&
            got.outlier_names == ref.outlier_names && got.dim == ref.dim && got.ls_npar == ref.ls_npar,
            "projected final description disagrees");
    const auto check = [](const auto& a, const auto& b, const char* field) {
        require(std::equal(a.begin(), a.end(), b.begin(), b.end(), bit_equal),
                std::string("projected final ") + field + " disagrees");
    };
    require(got.ls_obs == ref.ls_obs && got.ls_eff == ref.ls_eff && got.nx == ref.nx &&
            got.nmissing == ref.nmissing && got.scaling_factor == ref.scaling_factor &&
            got.free_arma_parameters == ref.free_arma_parameters, "projected final likelihood dimensions disagree");
    check(got.differenced_response, ref.differenced_response, "differenced response");
    require(got.variables.size() == ref.variables.size(), "projected final variable count disagrees");
    for (std::size_t i = 0; i < got.variables.size(); ++i) {
        const auto& actual = got.variables[i];
        const auto& expected = ref.variables[i];
        require(actual.name == expected.name && actual.type == expected.type &&
                actual.position == expected.position && bit_equal(actual.coefficient, expected.coefficient),
                "projected final variable disagrees");
        if (!expected.column.empty()) check(actual.column, expected.column, "variable column");
    }
    check(got.parameters, ref.parameters, "parameters");
    check(got.coefficients, ref.coefficients, "coefficients");
    check(got.unscaled_covariance, ref.unscaled_covariance, "regression covariance");
    check(got.full_residuals, ref.full_residuals, "residuals");
    for (auto member : {&TramoFitFinalRecord::ssq, &TramoFitFinalRecord::log_determinant,
         &TramoFitFinalRecord::ls_ll, &TramoFitFinalRecord::ls_tadj, &TramoFitFinalRecord::ls_all,
         &TramoFitFinalRecord::ls_aic, &TramoFitFinalRecord::ls_aicc, &TramoFitFinalRecord::ls_bic,
         &TramoFitFinalRecord::ls_bicc, &TramoFitFinalRecord::ls_bic2, &TramoFitFinalRecord::ls_hq})
        require(bit_equal(got.*member, ref.*member), "projected final likelihood statistics disagree");
}

void assemble_fitted_result(TramoFitResult& result, const TramoFitFinalRecord& fin,
                            const std::vector<double>* covariance_witness = nullptr) {
    result.ssq = fin.ssq;
    result.cl_dim = fin.dim;
    result.cl_nx = fin.nx;
    result.cl_nmissing = fin.nmissing;
    result.log_determinant = fin.log_determinant;
    auto& covariance = result.covariance;
    covariance.n = static_cast<std::size_t>(fin.nx);
    if (fin.nx > 0) {
        const int df = fin.dim - fin.nx - fin.free_arma_parameters;
        require(df > 0, "non-positive residual degrees of freedom");
        result.covariance_scale = fin.ssq / df;
        covariance.covariance = fin.unscaled_covariance;
        for (double& value : covariance.covariance) value *= result.covariance_scale;
        for (std::size_t i = 0; i < covariance.n; ++i)
            for (std::size_t j = i; j < covariance.n; ++j)
                covariance.covariance_upper.push_back(covariance.covariance[i * covariance.n + j]);
    }
    if (covariance_witness) {
        covariance.captured_upper = *covariance_witness;
        covariance.matches_capture = std::equal(covariance.covariance_upper.begin(), covariance.covariance_upper.end(),
            covariance_witness->begin(), covariance_witness->end(), bit_equal);
        require(covariance.matches_capture, "the recomputed D9 covariance differs from the captured coefficientsCovariance bits");
    }
    result.log = fin.log;
    result.mean = fin.mean;
    result.orders = fin.orders;
    result.parameters = fin.parameters;
    result.variables = fin.variables;
    result.regression_coefficients = fin.coefficients;
    result.ls_obs = fin.ls_obs;
    result.ls_eff = fin.ls_eff;
    result.ls_npar = fin.ls_npar;
    result.ls_ll = fin.ls_ll;
    result.ls_tadj = fin.ls_tadj;
    result.ls_all = fin.ls_all;
    result.ls_ssq = fin.ls_ssq;
    result.ls_aic = fin.ls_aic;
    result.ls_aicc = fin.ls_aicc;
    result.ls_bic = fin.ls_bic;
    result.ls_bicc = fin.ls_bicc;
    result.ls_bic2 = fin.ls_bic2;
    result.ls_hq = fin.ls_hq;
    result.full_residuals = fin.full_residuals;
    result.status = "available";
}

bool outliers_equal(const std::vector<std::string>& a,
                    const std::vector<std::string>& b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

// Seam kinds that mutate the description OUTSIDE the module seams (the
// round-0 restore, the pass-2 gate restore/lastSolution, the pass-1
// difference drop): a module ENTRY state may legitimately differ from the
// previous witnessed post-state across one of them.
bool kernel_transition_between(const std::vector<TramoKernelSeam>& seams,
                               std::size_t from, std::size_t to) {
    for (std::size_t k = from + 1; k <= to; ++k) {
        switch (seams[k].kind) {
            case TramoKernelSeamKind::round0:
            case TramoKernelSeamKind::pass2_gate:
            case TramoKernelSeamKind::differencing_drop:
                return true;
            default:
                break;
        }
    }
    return false;
}

// The module seams carry the post-run description state.
bool carries_state(TramoKernelSeamKind kind) {
    return kind == TramoKernelSeamKind::differencing ||
           kind == TramoKernelSeamKind::arma ||
           kind == TramoKernelSeamKind::scan;
}

const char* module_name(TramoKernelSeamKind kind) {
    switch (kind) {
        case TramoKernelSeamKind::differencing: return "tramo_differencing";
        case TramoKernelSeamKind::arma: return "tramo_arma_selection";
        default: return "tramo_outlier_rounds";
    }
}

}  // namespace

TramoFitResult estimate_tramo_tr1(const std::vector<double>& series, int start_year, int start_month) {
    const auto& preset = *tramo_template("TR1");
    require(series.size() >= 24 && series.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()) &&
            std::all_of(series.begin(), series.end(), [](double v) { return std::isfinite(v); }),
            "TR1 requires a complete monthly series spanning at least two years");
    tramo_build_request build;
    build.values = series;
    build.start_year = start_year;
    build.start_month = start_month;
    build.arima = preset.arima;
    build.mean = preset.mean;
    build.transform_function = preset.transform_function;
    build.using_auto_model = preset.using_auto_model;
    TramoFitResult result;
    result.builder = build_tramo_model(build);

    TramoFitSpec spec;
    spec.using_auto_model = preset.using_auto_model;
    spec.ao = preset.outliers_ao;
    spec.ls = preset.outliers_ls;
    spec.tc = preset.outliers_tc;
    spec.so = preset.outliers_so;
    spec.delta_tc = preset.delta_tc;
    TramoTransformDecisionRequest transform{};
    transform.values = series;
    transform.period = 12;
    transform.mean = result.builder.mean;
    transform.detect_seasonality = preset.using_auto_model;
    transform.log_preference = detail::fdlibm_log(0.95);
    transform.precision = spec.intermediate_precision;
    result.transformation = decide_tramo_transform(transform);
    const bool log = result.transformation.log_transform;
    EstimationEntry entry{{0, 1, 1, 0, 1, 1}, result.builder.mean, {}};
    const auto cv = tramo_critical_value(series.size());
    auto scan = run_scan(spec, result.builder, log, entry.orders, entry.mean, cv);
    require(scan && scan->outcome != TramoKernelStageResult::failed, "TR1 outlier scan failed");
    entry.variables = std::move(scan->variables);
    // Both TR1 passes scan the same fixed airline with the same threshold
    // after removing automatic outliers. The mean and orders are fixed and
    // the non-AMI final estimator returns directly on either pass.
    const auto y = differenced_response(series, log, 12, 0, 0);
    const auto effects = fit_estimation(spec, y, entry, 1);
    require(effects.fit.returned && effects.fit.cl_state != 0, "TR1 final estimation failed");
    const auto record = project_fitted_record(effects, series, log, 12);
    assemble_fitted_result(result, record);
    result.final_model_computed = true;
    return result;
}

TramoFitResult tramo_fit(const TramoFitRequest& request) {
    TramoFitResult result;
    const TramoFitSpec& spec = request.spec;

    // ---- structural validation (throws) ---------------------------------
    require(!request.series.empty(), "empty series");
    require(request.start_month >= 1 && request.start_month <= 12,
            "start_month outside 1..12");

    // ---- fail-closed refusal matrix (before any state is built) ---------
    auto refuse = [&result](const char* reason, std::string note) {
        result.status = std::string("refused:") + reason;
        result.refusal_note = std::move(note);
        return result;
    };
    if (spec.annual_frequency != 12)
        return refuse("d11-frequency",
                      "annual frequency " +
                          std::to_string(spec.annual_frequency) +
                          " is outside the monthly-only envelope (D11)");
    for (const double v : request.series) {
        if (!std::isfinite(v))
            return refuse("missing-values",
                          "the fit envelope carries no interpolation plan "
                          "(qualified at unit level in 23-A/23-B)");
    }
    if (spec.td_automatic || spec.td_test || spec.easter_test ||
        spec.user_variables)
        return refuse("d10-calendar",
                      "the calendar family (trading days automatic/tested, "
                      "Easter tests, user variables) is the joint 21-I "
                      "surface (D10)");
    if (spec.auto_method != tramo_auto_method::unused)
        return refuse("i10-automethod",
                      "the automatic regression methods (FTEST/WALD/BIC/AIC) "
                      "are refused (calendar-roadmap I-10)");
    if (spec.outliers_correction)
        return refuse("oc-precorrection",
                      "the robust outliers pre-correction has never been "
                      "witnessed through the fit seam");
    if (spec.seats_decomposition_mode != -1)
        return refuse("d6-decomposition",
                      "the fit never decomposes (register D6); decomposition "
                      "policy belongs to the 24-D/24-E consumers");

    const TramoKernelRoundsPrestage& prestage = request.kernel.prestage;

    // ---- stage 1: the 23-A builder mirror, cross-checked ----------------
    tramo_build_request breq;
    breq.values = request.series;
    breq.annual_frequency = spec.annual_frequency;
    breq.start_year = request.start_year;
    breq.start_month = request.start_month;
    breq.using_auto_model = spec.using_auto_model;
    breq.transform_function = spec.transform_function;
    breq.transform_adjust = tramo_lop::none;
    breq.arima = spec.fixed_arima;  // read only when !using_auto_model
    breq.mean.used = spec.mean_used;
    breq.mean.trend_constant = spec.mean_used;
    breq.mean.test = spec.mean_is_test;
    breq.delta_tc = spec.delta_tc;
    breq.outliers = request.builder_outliers;
    const tramo_build_result built = build_tramo_model(breq);

    require(static_cast<std::size_t>(built.domain_length) ==
                request.series.size(),
            "builder domain length disagrees with the series");
    require(built.missing.empty(),
            "builder registered missings inside the fit envelope");
    require(built.log_transformation ==
                (spec.transform_function == tramo_transform_function::log),
            "builder log flag disagrees with the transform function");
    {
        // Jar-decoded pre-stage (javap TramoKernel.testSeasonality): the
        // description is rewritten to SarimaOrders.m011(freq) = (0,1,1)(0,0,0)
        // with ctx.seasonal := false only when autoModel is on and the
        // seasonality verdict toInt() < 2 (the fit_gs2 / fit_shortfail shape:
        // prestate orders (0,1,1)(0,0,0) with seasonal off). Every other row
        // keeps the registered orders, so the builder mirror is the authority
        // there; a mismatch is admitted only in the decoded downgrade shape.
        TramoOrders built_orders;
        built_orders.d = built.arima.d;
        built_orders.bd = built.arima.bd;
        built_orders.p = static_cast<int>(built.arima.phi.size());
        built_orders.q = static_cast<int>(built.arima.theta.size());
        built_orders.bp = static_cast<int>(built.arima.bphi.size());
        built_orders.bq = static_cast<int>(built.arima.btheta.size());
        if (!orders_equal(built_orders, prestage.orders)) {
            TramoOrders m011;
            m011.d = 1;
            m011.q = 1;
            require(spec.using_auto_model && !prestage.seasonal &&
                        orders_equal(prestage.orders, m011),
                    "captured prestage orders " + orders_word(prestage.orders) +
                        " match neither the builder registration " +
                        orders_word(built_orders) + " nor the seasonality-"
                        "downgraded (0,1,1)(0,0,0) shape");
        }
    }
    if (prestage.regression_test == "Unprocessed")
        require(built.mean == prestage.mean,
                "builder mean flag disagrees with the captured prestage "
                "mean");
    result.builder = built;

    // ---- stage 2: the 23-B transform mirror, cross-checked --------------
    // testTransformation runs its machine only under TransformationType.Auto
    // with outliersCorrection off (javap, TramoKernel); the +fs row leaves
    // the context defaults in force (seasonal true, log off).
    if (spec.transform_function == tramo_transform_function::automatic) {
        TramoTransformDecisionRequest treq;
        treq.values = built.interpolated;
        treq.period = spec.annual_frequency;
        treq.mean = built.mean;
        treq.detect_seasonality = spec.using_auto_model;
        treq.spec_seasonality = true;  // the ModellingContext default
        // LogSpec.getLogPreference(): the resolved fct through the pinned
        // log (the TramoSpec default is Math.log(0.95)); the resolved line
        // is the authority and the request carries the raw double.
        treq.log_preference = spec.log_preference;
        treq.precision = spec.intermediate_precision;
        result.transformation = decide_tramo_transform(treq);
    } else {
        result.transformation.seasonal = true;
        result.transformation.spec_pre = "(0,1,1)(0,1,1)";
    }
    result.transform_matches_capture =
        result.transformation.seasonal == prestage.seasonal &&
        result.transformation.log_transform == prestage.log;
    require(result.transformation.seasonal == prestage.seasonal,
            "seasonality flag disagrees with the captured prestage");
    require(result.transformation.log_transform == prestage.log,
            "log flag disagrees with the captured prestage");

    // ---- stage 3: the ami() round loop ------------------------------------
    // Jar-decoded ami() branch (javap TramoKernel.ami / isFullySpecified):
    // the round driver is skipped only when the transform is already decided
    // and nothing is automatic - function != Auto && !autoModel &&
    // !outliersDetection && !meanTest && !td.test && !td.automatic &&
    // !easter.test (the last three clauses are vacuous here: the D10 refusal
    // above gates the calendar family out). fit_nomean (auto transform,
    // templates otherwise fixed) therefore runs at least one iterate - the
    // seam stream is cross-checked against the decoded verdict.
    const bool fully_specified =
        spec.transform_function != tramo_transform_function::automatic &&
        !spec.using_auto_model && !spec.outliers_used && !spec.mean_is_test;
    if (fully_specified) {
        require(request.kernel.seams.empty(),
                "a fully-specified request carries round-loop seams");
    } else {
        require(!request.kernel.seams.empty(),
                "a round-driver request carries no round-loop seams");
        result.kernel = tramo_run_kernel_rounds(request.kernel);
    }

    // ---- stage 4: the per-seam live cross-checks --------------------------
    const std::vector<TramoKernelSeam>& seams = request.kernel.seams;
    // The evolving seasonality context (dated 2026-11-23 javap re-decode):
    // the TramoKernel module factories (differencingModule(), armaModule())
    // hand the builder seasonal(ctx.seasonal) - the context at module
    // creation, not the pre-stage value. The seasonality controller seams
    // carry the context mutations (TramoKernelSeam.seasonal_after, from the
    // per-step rsnap witness); the gs2 / shortfail pass-3 modules are created
    // after a "scontroller res Changed" flip back to seasonal.
    const bool has_fixed_regressors = std::any_of(built.variables.begin(), built.variables.end(),
        [](const auto& var) { return var.kind != tramo_variable_kind::trend_constant; });
    std::optional<FittedOutlierEffects> independent_reference;
    if (!fully_specified && !has_fixed_regressors && prestage.outliers.empty()) {
        auto arima = built.arima;
        if (spec.using_auto_model && !result.transformation.seasonal) {
            arima.bd = 0;
            arima.btheta.clear();
        }
        const TramoOrders orders{static_cast<int>(arima.phi.size()), arima.d, static_cast<int>(arima.theta.size()),
            static_cast<int>(arima.bphi.size()), arima.bd, static_cast<int>(arima.btheta.size())};
        const auto y = differenced_response(built.interpolated, result.transformation.log_transform,
                                            built.annual_frequency, 0, 0);
        const auto initial = initial_mean_test(spec, y, orders, built.mean);
        const char* name = initial.result == TramoRegressionTestResult::unprocessed ? "Unprocessed" :
                           initial.result == TramoRegressionTestResult::changed ? "Changed" : "Unchanged";
        require(orders_equal(orders, prestage.orders) && initial.mean == prestage.mean &&
                prestage.regression_test == name, "initial fitted mean test disagrees with the pre-stage");
        result.initial_regression_checked = true;
        if (spec.using_auto_model || spec.outliers_used)
            independent_reference = fit_reference(spec, built, result.transformation.log_transform, y,
                                                  {orders, initial.mean, {}});
    }
    bool ctx_seasonal = prestage.seasonal;
    std::size_t next_input = 0;
    std::optional<EstimationEntry> final_entry;
    std::optional<FittedOutlierEffects> effects, reference_effects;
    bool first_gate = true;
    for (std::size_t s = 0; s < seams.size(); ++s) {
        const TramoKernelSeam& seam = seams[s];
        if (effects && effects->fit.cl_state != 0 &&
            (seam.kind == TramoKernelSeamKind::round0 || seam.kind == TramoKernelSeamKind::pass2_gate)) {
            check_round_statistics(*effects, seam.stats, built.annual_frequency,
                                   static_cast<int>(built.interpolated.size()), s);
            ++result.statistics_checked;
        }
        if (seam.kind == TramoKernelSeamKind::round0) {
            if (independent_reference) {
                const auto& physical = *independent_reference;
                require(effects && orders_equal(effects->orders, physical.orders) &&
                        physical.mean == effects->mean && physical.names == effects->names &&
                        std::equal(physical.parameters.begin(), physical.parameters.end(),
                                   effects->parameters.begin(), effects->parameters.end(), bit_equal),
                        "independent round-zero model disagrees");
                check_round_statistics(physical, seam.stats, built.annual_frequency, built.estimation_length, s);
                effects = std::move(independent_reference);
                result.reference_model_computed = true;
            }
            reference_effects = effects;
        }
        if (seam.kind == TramoKernelSeamKind::pass2_gate && first_gate) {
            first_gate = false;
            if (result.reference_model_computed && reference_effects) {
                auto candidate = fit_first_candidate(spec, built, result.transformation.log_transform,
                                                     result.transformation.seasonal, *reference_effects);
                if (candidate) {
                    require(orders_equal(candidate->orders, seam.arima_orders) && effects &&
                            effects->mean == candidate->mean && effects->names == candidate->names &&
                            std::equal(candidate->parameters.begin(), candidate->parameters.end(),
                                       seam.parameters.begin(), seam.parameters.end(), bit_equal),
                            "independent first automatic orders or parameters disagree");
                    check_round_statistics(*candidate, seam.stats, built.annual_frequency, built.estimation_length, s);
                    effects = std::move(candidate);
                    result.first_candidate_computed = true;
                }
            }
        }
        if (seam.kind == TramoKernelSeamKind::pass2_gate && seam.reverted)
            effects = reference_effects;
        if (seam.kind == TramoKernelSeamKind::estimation && seam.estimation_taken)
            effects.reset();
        if (effects && effects->fit.cl_state != 0 && seam.kind == TramoKernelSeamKind::regression_test) {
            if (check_regression_test(*effects, result.kernel.ami_check_mu, seam))
                effects.reset();
            ++result.regression_tests_checked;
        }
        if (effects && (seam.kind == TramoKernelSeamKind::seasonality_controller ||
                        seam.kind == TramoKernelSeamKind::model_controller) &&
            (seam.parameters != effects->parameters || seam.outliers != effects->names ||
             !orders_equal(seam.orders, effects->orders) || seam.mean != effects->mean))
            effects.reset();
        if (s == 0 && seam.kind == TramoKernelSeamKind::estimation && prestage.outliers.empty()) {
            final_entry = EstimationEntry{prestage.orders, prestage.mean, {}};
            continue;
        }
        if (final_entry && seam.kind == TramoKernelSeamKind::regression_test && seam.retry_site &&
            s > 0 && seams[s - 1].kind == TramoKernelSeamKind::estimation_model) {
            final_entry->mean = seam.mean;
            continue;
        }
        if (seam.kind == TramoKernelSeamKind::seasonality_controller) {
            ctx_seasonal = seam.seasonal_after;
            final_entry.reset();
            continue;
        }
        if (seam.kind == TramoKernelSeamKind::estimation_model) {
            TramoFitSeamCheck check{s, "final_estimator_run", "fact", "entry regression design unavailable"};
            if (!has_fixed_regressors && final_entry) {
                const auto y = differenced_response(built.interpolated, prestage.log, built.annual_frequency, 0, 0);
                effects = check_final_estimation(spec, y, *final_entry, seam);
                check.tier = "live";
                check.note = "match final estimation decision, orders and parameters";
                ++result.seams_live;
            } else {
                effects.reset();
                ++result.seams_fact;
            }
            result.seam_checks.push_back(std::move(check));
            continue;
        }
        // An intervening stage can replace the description. Carry an entry
        // only across its immediately verified predecessor or mean-test retry.
        final_entry.reset();
        if (seam.kind == TramoKernelSeamKind::scan && seam.outliers.empty())
            final_entry = EstimationEntry{seam.orders, seam.mean, {}};
        if (seam.kind != TramoKernelSeamKind::differencing &&
            seam.kind != TramoKernelSeamKind::arma &&
            seam.kind != TramoKernelSeamKind::scan)
            continue;
        while (next_input < request.seam_inputs.size() &&
               request.seam_inputs[next_input].seam < s)
            ++next_input;
        const bool has_input = next_input < request.seam_inputs.size() &&
                               request.seam_inputs[next_input].seam == s;

        TramoFitSeamCheck check;
        check.seam = s;
        check.unit = module_name(seam.kind);
        if (!has_input) {
            check.tier = "fact";
            check.note = "no fitseam entry record";
            result.seam_checks.push_back(check);
            ++result.seams_fact;
            continue;
        }
        const TramoFitSeamInputs& in = request.seam_inputs[next_input];

        // Zero-drift entry chain: the entry state must equal the previous
        // witnessed module post-state unless a kernel-owned transition sits
        // between them.
        for (std::size_t t = s; t-- > 0;) {
            if (!carries_state(seams[t].kind))
                continue;
            if (!kernel_transition_between(seams, t, s)) {
                require(orders_equal(seams[t].orders, in.entry_orders) &&
                            seams[t].mean == in.entry_mean &&
                            outliers_equal(seams[t].outliers,
                                           in.entry_outliers),
                        "seam " + std::to_string(s) +
                            " entry state disagrees with the previous "
                            "witnessed post-state");
            }
            break;
        }

        const bool can_derive_arma = !has_fixed_regressors &&
            (in.entry_outliers.empty() || !in.entry_variables.empty());
        const bool fitted_effects = effects && !effects->corrected_response.empty() &&
            effects->names == in.entry_outliers && orders_equal(effects->orders, in.entry_orders) &&
            effects->mean == in.entry_mean;
        if (seam.kind == TramoKernelSeamKind::differencing && !has_fixed_regressors &&
            (in.entry_outliers.empty() || fitted_effects)) {
            const bool raw = in.entry_outliers.empty();
            const auto d = difference_orders(spec, raw ? built.interpolated : effects->corrected_response,
                                             raw && prestage.log, in.entry_orders, in.entry_mean,
                                             in.round == 1, ctx_seasonal, in.entry_nobs > 0 ? in.entry_nobs : -1);
            const TramoKernelStageResult mapped =
                d.result == TramoDifferencingDecision::Result::Changed
                    ? TramoKernelStageResult::changed
                    : d.result ==
                              TramoDifferencingDecision::Result::Unprocessed
                          ? TramoKernelStageResult::unprocessed
                          : TramoKernelStageResult::unchanged;
            TramoOrders post;
            post.p = d.spec_p;
            post.d = d.spec_d;
            post.q = d.spec_q;
            post.bp = d.spec_bp;
            post.bd = d.spec_bd;
            post.bq = d.spec_bq;
            require(mapped == seam.result,
                    "seam " + std::to_string(s) +
                        ": live differencing result " + stage_word(mapped) +
                        " disagrees with the seam fact " +
                        stage_word(seam.result));
            require(orders_equal(post, seam.orders),
                    "seam " + std::to_string(s) +
                        ": live differencing orders " + orders_word(post) +
                        " disagree with the seam fact " +
                        orders_word(seam.orders));
            require(d.mean == seam.mean,
                    "seam " + std::to_string(s) +
                        ": live differencing mean flag disagrees");
            check.tier = "live";
            check.note = "match " + stage_word(mapped) + " " +
                         orders_word(post);
        } else if (seam.kind == TramoKernelSeamKind::scan && !has_fixed_regressors && in.entry_outliers.empty()) {
            auto selected = check_scan(spec, built, prestage.log, in, seam);
            if (selected) {
                final_entry = EstimationEntry{in.entry_orders, in.entry_mean, std::move(*selected)};
                check.tier = "live";
                check.note = "match " + stage_word(seam.result) + ", ordered outlier selection";
            } else {
                check.tier = "fact";
                check.note = "outside the qualified stationary scan design";
            }
        } else if (seam.kind == TramoKernelSeamKind::arma &&
                   (can_derive_arma || !in.ols_residuals.empty())) {
            TramoArmaSelectionRequest areq;
            if (can_derive_arma) {
                areq.residuals = arma_residuals(built, prestage.log, in, spec.delta_tc);
                areq.nobs = built.estimation_length;
                if (!in.ols_residuals.empty()) {
                    require(areq.nobs == in.nobs && areq.residuals.size() == in.ols_residuals.size() &&
                            std::equal(areq.residuals.begin(), areq.residuals.end(),
                                       in.ols_residuals.begin(), bit_equal),
                            "seam " + std::to_string(s) + ": native ARMA residuals differ from the witness");
                }
            } else {
                areq.residuals = in.ols_residuals;
                areq.nobs = in.nobs;
            }
            const auto a = select_orders(std::move(areq), spec, in.entry_orders, ctx_seasonal);
            const TramoKernelStageResult mapped =
                a.result == TramoArmaSelectionResult::Result::Changed
                    ? TramoKernelStageResult::changed
                    : a.result == TramoArmaSelectionResult::Result::Failed
                          ? TramoKernelStageResult::failed
                          : a.result ==
                                    TramoArmaSelectionResult::Result::
                                        Unprocessed
                                ? TramoKernelStageResult::unprocessed
                                : TramoKernelStageResult::unchanged;
            TramoOrders post;
            post.p = a.spec_p;
            post.d = a.spec_d;
            post.q = a.spec_q;
            post.bp = a.spec_bp;
            post.bd = a.spec_bd;
            post.bq = a.spec_bq;
            require(mapped == seam.result,
                    "seam " + std::to_string(s) +
                        ": live arma result disagrees with the seam fact");
            require(orders_equal(post, seam.orders),
                    "seam " + std::to_string(s) +
                        ": live arma orders disagree with the seam fact");
            check.tier = "live";
            if (in.entry_outliers.empty())
                final_entry = EstimationEntry{post, in.entry_mean, {}};
            check.note = "match " + stage_word(mapped) + " " +
                         orders_word(post);
        } else {
            check.tier = "fact";
            if (seam.kind == TramoKernelSeamKind::scan)
                check.note =
                    "outside the live scan design or differencing envelope";
            else if (seam.kind == TramoKernelSeamKind::arma)
                check.note = "no olsResiduals hand-over captured";
            else
                check.note =
                    "no current fitted outlier effects for the differencing input";
        }
        result.seam_checks.push_back(check);
        if (check.tier == "live")
            ++result.seams_live;
        else
            ++result.seams_fact;
    }

    // ---- stage 5: the final surface and the D9 covariance -----------------
    if (fully_specified && !has_fixed_regressors && built.arima.is_airline() &&
        built.arima.theta.front().kind == tramo_parameter_kind::undefined &&
        built.arima.btheta.front().kind == tramo_parameter_kind::undefined) {
        EstimationEntry entry{{0, 1, 1, 0, 1, 1}, built.mean, {}};
        const auto y = differenced_response(built.interpolated, prestage.log, built.annual_frequency, 0, 0);
        effects = fit_estimation(spec, y, entry, 0);
    }
    std::optional<TramoFitFinalRecord> projected;
    if (effects && effects->fit.cl_state != 0) {
        projected = project_fitted_record(*effects, built.interpolated, prestage.log, built.annual_frequency);
        if (request.verify_final_record) check_final_record(*projected, request.final);
        result.final_model_computed = true;
    }
    require(projected || request.verify_final_record, "final estimator unavailable for native fitted record");
    const auto& fin = projected ? *projected : request.final;
    const auto* covariance_witness = request.verify_final_record ? &request.final.captured_covariance_upper : nullptr;
    require(fin.dim > 0, "final record without a dim");
    require(fin.nx >= 0, "negative nx");
    require(fin.free_arma_parameters == fin.orders.p + fin.orders.q +
                                         fin.orders.bp + fin.orders.bq,
            "free_arma_parameters disagrees with the final orders");
    require(fin.parameters.size() ==
                static_cast<std::size_t>(fin.free_arma_parameters),
            "the final parameter vector length disagrees with k");
    require(fin.coefficients.size() == static_cast<std::size_t>(fin.nx),
            "the coefficient vector length disagrees with nx");
    require(fin.unscaled_covariance.size() ==
                static_cast<std::size_t>(fin.nx) * fin.nx,
            "unscaledCovariance is not nx x nx");
    require(!covariance_witness || covariance_witness->size() ==
                static_cast<std::size_t>(fin.nx) * (fin.nx + 1) / 2,
            "the captured covariance upper triangle has the wrong size");
    require(fin.scaling_factor, "the cl view is not the scaling-factor "
                                "branch covariance(k, true) consumes");
    require(fin.differenced_response.size() ==
                static_cast<std::size_t>(fin.dim),
            "the differenced response length disagrees with cl dim");
    require(fin.variables.size() == fin.outlier_names.size(),
            "the variable records disagree with the outlier names");
    {
        const std::size_t base = fin.mean ? 1 : 0;
        require(base + fin.variables.size() == fin.coefficients.size(),
                "the regression design (mean + outliers) disagrees with nx");
        for (std::size_t i = 0; i < fin.variables.size(); ++i) {
            require(fin.variables[i].name == fin.outlier_names[i],
                    "variable record " + std::to_string(i) +
                        " disagrees with the outlier name list");
            require(bit_equal(fin.variables[i].coefficient,
                              fin.coefficients[base + i]),
                    "the coefficient of variable " +
                        fin.variables[i].name +
                        " disagrees with the cl design slot");
        }
    }

    // The differenced response the landing machine consumed, recomputed
    // natively (pinned log, then the RegArimaUtility.differencingFilter /
    // BackFilter orientation of the 23-C pin: the JVM hands the BackFilter
    // the REVERSED polynomial, applied directly - the (-1)^(d + bd)
    // equivalence of the 23-C note holds entry-wise, but a cancelling sum
    // lands at JVM +0.0 where negating the ascending result would give
    // -0.0 (the fit_fct_air witness at index 99), so the reversed filter is
    // applied as the JVM builds it).
    const auto w = differenced_response(built.interpolated, fin.log, spec.annual_frequency,
                                        fin.orders.d, fin.orders.bd);
    require(w.size() == fin.differenced_response.size(),
            "the recomputed differenced response length disagrees");
    for (std::size_t i = 0; i < w.size(); ++i) {
        require(bit_equal(w[i], fin.differenced_response[i]),
                "the recomputed differenced response differs from the "
                "capture at index " + std::to_string(i));
    }

    assemble_fitted_result(result, fin, covariance_witness);
    return result;
}

}  // namespace kalos
