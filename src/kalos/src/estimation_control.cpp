#include "kalos/estimation_control.hpp"

#include "kalos/arma_optimizer.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/detail/ssq_levmar.hpp"
#include "kalos/hannan_rissanen.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

// Pinned RegSarimaComputer.UR_LIMIT.
constexpr double ur_limit = 0.1;

double unavailable() { return std::numeric_limits<double>::quiet_NaN(); }

const char* stage_name(EstimationStartingPoint point) {
    switch (point) {
    case EstimationStartingPoint::Zero: return "zero";
    case EstimationStartingPoint::Default: return "default";
    case EstimationStartingPoint::HannanRissanen: return "hannan_rissanen";
    default: return "multiple";
    }
}

// The pinned SarimaModel stores the raw coefficients of the polynomial
// 1 + phi(B); the shared Kalos convention is 1 - phi(B) with MA 1 + theta(B)
// on both sides. The unit-root passes are transcribed literally on the
// pinned convention, so every entry crosses this sign flip.
std::array<std::vector<double>, 4> flip_ar(std::array<std::vector<double>, 4> coefficients) {
    for (const std::size_t b : {std::size_t(0), std::size_t(2)}) {
        for (auto& value : coefficients[b]) value = -value;
    }
    return coefficients; // The flip is an involution.
}

// Pinned tryUrpCancelling deflate: par[i] = -acc; acc += old, top-down, for
// the regular AR and regular MA blocks on the pinned arrays.
void deflate_regular(std::array<std::vector<double>, 4>& pinned) {
    for (const std::size_t block : {std::size_t(0), std::size_t(1)}) {
        auto& coefficients = pinned[block];
        double acc = 0.0;
        for (std::size_t i = coefficients.size(); i-- > 0;) {
            const double old = coefficients[i];
            coefficients[i] = -acc;
            acc += old;
        }
    }
}

// Pinned tryUrmCancelling deflate: par[i] = acc; acc = old - acc, top-down,
// for the same (regular) blocks; the seasonal blocks are never touched.
void deflate_seasonal(std::array<std::vector<double>, 4>& pinned) {
    for (const std::size_t block : {std::size_t(0), std::size_t(1)}) {
        auto& coefficients = pinned[block];
        double acc = 0.0;
        for (std::size_t i = coefficients.size(); i-- > 0;) {
            const double old = coefficients[i];
            coefficients[i] = acc;
            acc = old - acc;
        }
    }
}

// The pinned gate accumulators on one regular block: the plain sum for the
// unit-root pass and the toggle sum (-c0 + c1 - c2 + ...) for the -1 pass.
double pinned_gate_sum(const std::vector<double>& coefficients, bool alternating) {
    double acc = 0.0;
    bool negate = true; // The pinned toggle negates the first entry.
    for (const double value : coefficients) {
        acc += (alternating && negate ? -value : value);
        negate = !negate;
    }
    return acc;
}

bool within_limit(double pinned_sum) { return std::abs(1.0 + pinned_sum) <= ur_limit; }

// Solver status names that mark the producing pass as converged: the BFGS
// gradient criterion and the pinned machine's hasConverged stop rule.
bool pass_reports_convergence(const std::string& status) {
    return status == "converged_gradient" || status == "converged_function";
}

// Maps physical coefficients through the qualified parameter map to the
// solver coordinates. Reflection-mode blocks use the map's own inverse
// construction (atanh of the reflection over 1 - margin); coefficient-mode
// blocks pass free entries directly. A non-admissible proposal yields false:
// the pinned raw-coordinate optimizer would throw there and its passes keep
// the incumbent, so the control layer does the same.
bool to_coordinates(const ArmaParameterMap& map,
    const std::array<std::vector<double>, 4>& physical, std::vector<double>& out) {
    out.clear();
    for (std::size_t b = 0; b < 4; ++b) {
        const auto mode = map.modes()[b];
        const auto& block = map.blocks()[b];
        if (mode == "empty" || mode == "fixed") continue;
        if (mode == "coefficient") {
            for (std::size_t j = 0; j < physical[b].size(); ++j) {
                if (!block.fixed[j]) out.push_back(physical[b][j]);
            }
            continue;
        }
        if (mode != "reflection") return false;
        const auto check = arma_admissibility(physical[b], b == 1 || b == 3, map.margin());
        if (check.status != "admissible") return false;
        for (const double reflection : check.reflections) {
            const double coordinate = detail::fdlibm_atanh(reflection / (1 - map.margin()));
            if (!std::isfinite(coordinate)) return false;
            out.push_back(coordinate);
        }
    }
    return true;
}

struct PassOutcome {
    bool available = false;
    std::string optimizer_status;
    std::size_t evaluations = 0;
    double log_likelihood = unavailable();
    std::array<std::vector<double>, 4> coefficients;
    // Item 23-F F.2 landing seam (set when the pinned sum-of-squares path
    // ran and request_landing_point was asked; the FinalEstimator round
    // driver reads it as the machine-failure seam).
    bool landing_derivatives = false;
};

bool pinned_solver(const SarimaEstimationRequest& request) {
    return request.control.solver == EstimationSolver::PinnedLevenbergMarquardt;
}

// Pinned RegArmaSsqFunction with maximumLikelihood(true) over the request
// surface (item 21-D.2). Qualified requests evaluate on the transcribed JVM
// deviances machine (detail::pinned_deviances_evaluate): the residual vector
// is the Evaluation's full-length deviances, the target is ssqErr * factor,
// and the landing value is the concentrated likelihood's logLikelihood, all
// bit-identical to the pinned machine at the differencing-decision seams.
// Requests the machine is not qualified for (missing observations, designs
// that are not the meanCorrection constant, sub-period seasonality) keep the
// 21-A GLS provider, which inverts
// ll = -0.5 * (n * LOGTWOPI + n * (1 + log(ssq/n)) + ldet),
// the pinned ConcentratedLikelihoodWithMissing constructor formula, at the
// surface-identity level the 21-A seams are qualified at. Coordinates are
// the JVM-convention [phi, bphi, theta, btheta] blocks.
class PinnedSsqProvider final : public detail::SsqVectorFunction {
public:
    PinnedSsqProvider(const SarimaEstimationRequest& request, std::size_t p, std::size_t bp,
        std::size_t q, std::size_t bq, const Interrupt& interrupt)
        : request_(request), interrupt_(interrupt), domain_(static_cast<int>(p),
              static_cast<int>(bp), static_cast<int>(q), static_cast<int>(bq), request.gls_sarima_mapping) {
        // Qualify the pinned deviances machine (21-D.2, multi-column
        // generalization item 23-F F.2) for this request: the machine carries
        // no missing observations, needs a seasonal period, and knows the
        // design-free (sfilter), constant-regressor (mfilter, the 1x1
        // machine) and full pivoting-QR multi-column (mfilter, hwp_n)
        // machines. Missing observations and sub-period seasonality keep the
        // 21-A GLS provider.
        if (request_.period < 2 || request_.d < 0 || request_.seasonal_d < 0) {
            return;
        }
        if (!request_.series.missing.empty()) {
            for (const std::uint8_t flag : request_.series.missing) {
                if (flag != 0) return;
            }
        }
        for (const double value : request_.series.values) {
            if (!std::isfinite(value)) return;
        }
        const std::vector<double> kernel = detail::pinned_differencing_polynomial(
            request_.period, request_.d, request_.seasonal_d);
        machine_y_ = detail::pinned_backfilter_apply(kernel, request_.series.values);
        if (machine_y_.empty()) return;
        if (request_.design.columns == 0) {
            machine_usable_ = true;
            return;
        }
        if (!request_.design.fixed.empty()) {
            for (const std::uint8_t flag : request_.design.fixed) {
                if (flag != 0) return; // fixed design entries are not qualified
            }
        }
        const std::size_t rows = request_.series.values.size();
        if (request_.design.values.size() != rows * request_.design.columns) {
            return;
        }
        // The raw-domain columns map to the machine view the pinned
        // ConcentratedLikelihoodComputer sees (mirror of the scan machine's
        // design()): an all-ones column is the meanCorrection constant, 1.0
        // over the differenced domain; every other column carries the
        // unit-root filter.
        std::vector<std::vector<double>> cols;
        cols.reserve(request_.design.columns);
        for (std::size_t c = 0; c < request_.design.columns; ++c) {
            std::vector<double> raw(request_.design.values.begin()
                    + static_cast<std::ptrdiff_t>(c * rows),
                request_.design.values.begin()
                    + static_cast<std::ptrdiff_t>((c + 1) * rows));
            bool constant = true;
            for (const double value : raw) {
                if (value != 1.0) {
                    constant = false;
                    break;
                }
            }
            if (constant) {
                if (request_.design.columns == 1) {
                    machine_design_.assign(machine_y_.size(), 1.0);
                    machine_usable_ = true; // the verified 1x1 constant machine
                    return;
                }
                cols.push_back(std::vector<double>(machine_y_.size(), 1.0));
            } else {
                cols.push_back(detail::pinned_backfilter_apply(kernel, raw));
                if (cols.back().empty()) return;
            }
        }
        if (cols.size() != request_.design.columns) {
            return; // defensive: the column walk above must consume all
        }
        machine_cols_ = std::move(cols);
        machine_multi_ = true;
        machine_usable_ = true;
    }

    const detail::SarimaSsqDomain& domain() const override { return domain_; }

    // True when the qualified deviances machine backs this provider.
    bool machine_active() const { return machine_usable_; }

    std::optional<detail::SsqPoint> ssq_evaluate(const std::vector<double>& x) const override {
        const std::size_t p = domain_.sizes()[0], bp = domain_.sizes()[1];
        const std::size_t q = domain_.sizes()[2], bq = domain_.sizes()[3];
        if (x.size() != p + bp + q + bq) return std::nullopt;
        if (machine_usable_) {
            // x is already the JVM-convention [phi, bphi, theta, btheta]
            // vector the machine runs on; a rejected evaluation mirrors the
            // pinned Evaluation constructor throwing, which the pinned
            // iterateQR catch feeds to the damping tail.
            const auto point = machine_multi_
                ? detail::pinned_deviances_evaluate_multi(
                      machine_y_, p, bp, q, bq, request_.period, machine_cols_, x)
                : detail::pinned_deviances_evaluate(
                      machine_y_, p, bp, q, bq, request_.period, machine_design_, x);
            if (!point || point->e.size() != machine_y_.size() ||
                !std::isfinite(point->ssq_e)) {
                return std::nullopt;
            }
            detail::SsqPoint ssq;
            ssq.e = point->e;
            ssq.ssq_e = point->ssq_e;
            return ssq;
        }
        std::array<std::vector<double>, 4> c;
        c[0].assign(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(p));
        std::for_each(c[0].begin(), c[0].end(), [](double& v) { v = -v; });
        c[2].assign(x.begin() + static_cast<std::ptrdiff_t>(p),
            x.begin() + static_cast<std::ptrdiff_t>(p + bp));
        std::for_each(c[2].begin(), c[2].end(), [](double& v) { v = -v; });
        c[1].assign(x.begin() + static_cast<std::ptrdiff_t>(p + bp),
            x.begin() + static_cast<std::ptrdiff_t>(p + bp + q));
        c[3].assign(x.begin() + static_cast<std::ptrdiff_t>(p + bp + q), x.end());

        std::exception_ptr interrupted;
        const Interrupt checked = [&] {
            if (!interrupt_) return;
            try { interrupt_(); } catch (...) { interrupted = std::current_exception(); throw; }
        };
        RegressionLikelihood fit;
        try {
            fit = sarima_regression_likelihood(request_.series, request_.design, c[0], c[1],
                c[2], c[3], request_.period, request_.d, request_.seasonal_d, 0.0, 1.0, 1e-10,
                checked);
        } catch (const std::invalid_argument&) {
            if (interrupted) std::rethrow_exception(interrupted);
            return std::nullopt; // pinned Evaluation constructor throws
        } catch (const std::runtime_error&) {
            if (interrupted) std::rethrow_exception(interrupted);
            return std::nullopt;
        }
        if (fit.likelihood_status != "available" || fit.profile_status != "available") {
            return std::nullopt;
        }
        const std::size_t n = fit.finite_count;
        const double ssq = fit.sum_squared_standardized;
        const double ll = fit.profiled_log_likelihood;
        if (n == 0 || fit.whitened_residuals.size() != n || !std::isfinite(ssq) ||
            !std::isfinite(ll) || ssq <= 0.0) {
            return std::nullopt;
        }
        // The pinned scaling factor and the two deviances-derived values.
        const double ldet = -2.0 * ll
            - (static_cast<double>(n) * kLogTwoPi
                  + static_cast<double>(n)
                        * (1.0 + detail::fdlibm_log(ssq / static_cast<double>(n))));
        const double factor = detail::fdlibm_exp(ldet / static_cast<double>(n));
        const double root = std::sqrt(factor);
        if (!std::isfinite(factor) || !std::isfinite(root)) return std::nullopt;
        detail::SsqPoint point;
        point.e.reserve(n);
        for (const double value : fit.whitened_residuals) {
            const double scaled = value * root;
            if (!std::isfinite(scaled)) return std::nullopt;
            point.e.push_back(scaled);
        }
        point.ssq_e = ssq * factor; // pinned DefaultLikelihoodEvaluation.deviance
        if (!std::isfinite(point.ssq_e)) return std::nullopt;
        return point;
    }

    // Direct request-surface likelihood at a landed point, bypassing the
    // Kalos map admissibility margin (pinned-domain stabilizations can land
    // inside the pinned limits but outside the Kalos margin).
    std::optional<double> landing_likelihood(const std::array<std::vector<double>, 4>& c,
        const Interrupt& interrupt) const {
        if (machine_usable_) {
            // The machine's ConcentratedLikelihoodWithMissing logLikelihood,
            // the pinned landing value; c is in the Kalos blocks and signs.
            const auto p = domain_.sizes()[0];
            const auto bp = domain_.sizes()[1];
            const auto q = domain_.sizes()[2];
            const auto bq = domain_.sizes()[3];
            if (c[0].size() == p && c[2].size() == bp && c[1].size() == q && c[3].size() == bq) {
                std::vector<double> jvm;
                for (const double v : c[0]) jvm.push_back(-v);
                for (const double v : c[2]) jvm.push_back(-v);
                for (const double v : c[1]) jvm.push_back(v);
                for (const double v : c[3]) jvm.push_back(v);
                const auto point = machine_multi_
                    ? detail::pinned_deviances_evaluate_multi(
                          machine_y_, p, bp, q, bq, request_.period, machine_cols_, jvm)
                    : detail::pinned_deviances_evaluate(
                          machine_y_, p, bp, q, bq, request_.period, machine_design_, jvm);
                if (point && std::isfinite(point->log_likelihood)) {
                    return point->log_likelihood;
                }
                return std::nullopt;
            }
        }
        try {
            const auto fit = sarima_regression_likelihood(request_.series, request_.design,
                c[0], c[1], c[2], c[3], request_.period, request_.d, request_.seasonal_d, 0.0,
                1.0, 1e-10, interrupt);
            if (fit.likelihood_status != "available" || fit.profile_status != "available") {
                return std::nullopt;
            }
            if (!std::isfinite(fit.profiled_log_likelihood)) return std::nullopt;
            return fit.profiled_log_likelihood;
        } catch (const std::invalid_argument&) {
            return std::nullopt;
        } catch (const std::runtime_error&) {
            return std::nullopt;
        }
    }

private:
    // jdplus.api.math.Constants.LOGTWOPI.
    static constexpr double kLogTwoPi = 1.8378770664093453;

    const SarimaEstimationRequest& request_;
    const Interrupt& interrupt_;
    detail::SarimaSsqDomain domain_;
    // Pinned deviances machine setup (21-D.2): the unit-root-filtered
    // response, the machine design (empty or the constant column) and the
    // usability flag; the GLS path stays for everything the machine does not
    // qualify for (missing observations, multi-column or non-constant
    // designs, sub-period seasonality).
    bool machine_usable_ = false;
    bool machine_multi_ = false;
    std::vector<double> machine_y_;
    std::vector<double> machine_design_;
    std::vector<std::vector<double>> machine_cols_;
};

const char* pinned_lm_status(int stop) {
    switch (stop) {
    case 0: return "attempts_exhausted";
    case 1: return "converged_gradient";
    case 2: return "converged_function";
    case 4: return "step_overflow";
    case 5: return "mu_overflow";
    case 6: return "converged_scale";
    default: return "non_finite";
    }
}

// JVM-convention [phi, bphi, theta, btheta] coordinates back to Kalos block
// positions in the JVM sign convention; flip_ar completes the return to the
// Kalos convention.
std::array<std::vector<double>, 4> from_jvm_coordinates(const std::vector<double>& x,
    std::size_t p, std::size_t bp, std::size_t q) {
    std::array<std::vector<double>, 4> c;
    c[0].assign(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(p));
    c[2].assign(x.begin() + static_cast<std::ptrdiff_t>(p),
        x.begin() + static_cast<std::ptrdiff_t>(p + bp));
    c[1].assign(x.begin() + static_cast<std::ptrdiff_t>(p + bp),
        x.begin() + static_cast<std::ptrdiff_t>(p + bp + q));
    c[3].assign(x.begin() + static_cast<std::ptrdiff_t>(p + bp + q), x.end());
    return c;
}

PassOutcome run_pass_pinned_lm(const SarimaEstimationRequest& request,
    const std::array<std::vector<double>, 4>& physical, double precision,
    const Interrupt& interrupt) {
    PassOutcome outcome;
    // Eligibility gate (21-D.2, widened per the M9-H qualification): fully
    // free blocks run the pinned machine directly; partially fixed blocks are
    // handed to the SAME pinned LM machine over the reduced free-coordinate
    // seam (detail::SsqFixedSeam, the pinned SarimaFixedMapping semantics) —
    // the mapping the pinned JVM itself builds for fixed coefficients.
    bool has_fixed = false;
    for (const auto& block : request.blocks) {
        if (std::any_of(block.fixed.begin(), block.fixed.end(),
                [](const std::uint8_t flag) { return flag != 0; })) {
            has_fixed = true;
        }
    }
    const std::size_t p = request.blocks[0].initial.size();
    const std::size_t q = request.blocks[1].initial.size();
    const std::size_t bp = request.blocks[2].initial.size();
    const std::size_t bq = request.blocks[3].initial.size();
    const PinnedSsqProvider provider(request, p, bp, q, bq, interrupt);
    auto jvm = flip_ar(physical); // blockwise JVM sign convention
    if (has_fixed) {
        // Fixed slots carry the request's pinned values, not the starting
        // fill (pinned SarimaFixedMapping.parameters; the Zero/Default fill
        // above overwrites every entry, fixed included).
        for (std::size_t b = 0; b < 4; ++b) {
            for (std::size_t j = 0; j < request.blocks[b].fixed.size(); ++j) {
                if (request.blocks[b].fixed[j] != 0) {
                    jvm[b][j] = (b == 0 || b == 2) ? -request.blocks[b].initial[j]
                                                   : request.blocks[b].initial[j];
                }
            }
        }
    }
    std::vector<double> start;          // JVM block order [phi, bphi, theta, btheta]
    start.insert(start.end(), jvm[0].begin(), jvm[0].end());
    start.insert(start.end(), jvm[2].begin(), jvm[2].end());
    start.insert(start.end(), jvm[1].begin(), jvm[1].end());
    start.insert(start.end(), jvm[3].begin(), jvm[3].end());

    std::exception_ptr interrupted;
    const Interrupt checked = [&] {
        if (!interrupt) return;
        try { interrupt(); } catch (...) { interrupted = std::current_exception(); throw; }
    };
    detail::SsqLmResult lm;
    try {
        if (has_fixed) {
            // Reduced coordinates in the same JVM order as start: mask slots
            // [phi, bphi, theta, btheta] from blocks [ar, sar, ma, sma].
            const std::size_t dim = p + bp + q + bq;
            std::vector<std::uint8_t> mask(dim, 0);
            std::vector<double> pinned(dim, 0.0);
            std::size_t off = 0;
            for (const int b : {0, 2, 1, 3}) {
                const std::size_t size = jvm[b].size();
                for (std::size_t j = 0; j < size; ++j) {
                    mask[off + j] = request.blocks[b].fixed[j];
                    pinned[off + j] = jvm[b][j];
                }
                off += size;
            }
            std::vector<double> reduced;
            for (std::size_t j = 0; j < dim; ++j) {
                if (mask[j] == 0) reduced.push_back(start[j]);
            }
            const detail::SsqFixedSeam seam(
                provider.domain(), std::move(mask), std::move(pinned));
            lm = detail::minimize_ssq_levmar(
                provider, std::move(reduced), precision, checked, &seam);
        } else {
            lm = detail::minimize_ssq_levmar(provider, std::move(start), precision, checked);
        }
    } catch (const std::runtime_error&) {
        if (interrupted) std::rethrow_exception(interrupted);
        outcome.optimizer_status = "start_not_evaluable"; // pinned Evaluation throws
        return outcome;
    }
    const auto landing = flip_ar(from_jvm_coordinates(lm.parameters, p, bp, q));
    outcome.optimizer_status = pinned_lm_status(lm.stop);
    if (request.request_landing_point) {
        outcome.landing_derivatives = true;
    }
    if (provider.machine_active()) {
        // The pinned landing likelihood is the machine's concentrated
        // likelihood at the landing (21-D.2); it is bit-identical to the
        // seam, while the Kalos surface profile differs by rounding.
        const auto fit = provider.landing_likelihood(landing, checked);
        if (!fit) {
            outcome.log_likelihood = unavailable();
            return outcome;
        }
        outcome.available = true;
        outcome.log_likelihood = *fit;
        outcome.coefficients = std::move(landing);
        outcome.evaluations = static_cast<std::size_t>(lm.evaluations);
        return outcome;
    }
    const auto objective = evaluate_arma_physical_objective(request.series, request.design,
        ArmaParameterMap(request.blocks), landing, request.period, request.d,
        request.seasonal_d, 0.0, 1e-10, 0.0, interrupt);
    if (objective.status != "available") {
        // A landing the Kalos map rejects on its admissibility margin (the
        // pinned domain stabilizes with its own limits) keeps the pinned
        // behavior: the likelihood on the request surface decides.
        const auto fit = provider.landing_likelihood(landing, checked);
        if (!fit) {
            outcome.log_likelihood = unavailable();
            return outcome;
        }
        outcome.available = true;
        outcome.log_likelihood = *fit;
        outcome.coefficients = std::move(landing);
        outcome.evaluations = static_cast<std::size_t>(lm.evaluations);
        return outcome;
    }
    outcome.available = true;
    outcome.log_likelihood = objective.adjusted_log_likelihood;
    outcome.coefficients = std::move(landing);
    outcome.evaluations = static_cast<std::size_t>(lm.evaluations);
    return outcome;
}

PassOutcome run_pass(const SarimaEstimationRequest& request, const ArmaParameterMap& map,
    const std::array<std::vector<double>, 4>& physical, const std::vector<double>& coordinates,
    double precision, const Interrupt& interrupt) {
    if (pinned_solver(request)) {
        return run_pass_pinned_lm(request, physical, precision, interrupt);
    }
    PassOutcome outcome;
    OptimizerControl control = request.optimizer;
    control.function_tolerance = precision; // Pinned minimizer function precision.
    control.validate();
    OptimizationResult result = optimize_arma_objective(request.series, request.design, map,
        coordinates, request.period, request.d, request.seasonal_d, 0.0, 1e-10, 0.0, control,
        interrupt);
    outcome.optimizer_status = std::move(result.status);
    outcome.evaluations = static_cast<std::size_t>(result.evaluations);
    const auto point = std::static_pointer_cast<const ArmaOptimizationPoint>(result.final.detail);
    if (!point || result.final.status != "available" ||
        point->objective.status != "available") return outcome;
    outcome.available = true;
    outcome.log_likelihood = point->objective.adjusted_log_likelihood;
    outcome.coefficients = point->objective.parameters.coefficients;
    return outcome;
}

// One pinned unit-root-cancellation pass on the incumbent. The gates and the
// deflation reproduce the pinned bytecode branch by branch; a candidate is
// kept only on a strictly better likelihood and every failure keeps the
// incumbent (the pinned passes swallow all exceptions).
void run_cancelling_pass(std::vector<EstimationPassRecord>& passes,
    EstimationPassRecord record, const SarimaEstimationRequest& request,
    const ArmaParameterMap& map, bool seasonal, const Interrupt& interrupt,
    std::array<std::vector<double>, 4>& current, double& current_ll) {
    // Pinned gates: tryUrpCancelling passes on orders.getDifferenceOrder()
    // != 0 (d + period * bd when period > 1 and bd > 0), tryUrmCancelling on
    // orders.bd != 0; both additionally need a non-empty AR and MA block.
    // The orders are the spec's (spec_d/spec_bd when the response arrives
    // pre-differenced).
    const int spec_d = request.spec_d < 0 ? request.d : request.spec_d;
    const int spec_bd = request.spec_bd < 0 ? request.seasonal_d : request.spec_bd;
    record.gates_passed = !current[0].empty() && !current[1].empty() &&
        (seasonal ? spec_bd != 0
                  : spec_d + (request.period > 1 && spec_bd > 0
                        ? request.period * spec_bd
                        : 0) !=
                        0);
    if (record.gates_passed) {
        const auto pinned = flip_ar(current);
        record.gates_passed = within_limit(pinned_gate_sum(pinned[0], seasonal)) &&
            within_limit(pinned_gate_sum(pinned[1], seasonal));
        if (record.gates_passed) {
            auto deflated = flip_ar(current);
            if (seasonal) deflate_seasonal(deflated); else deflate_regular(deflated);
            record.attempted = true;
            auto candidate_physical = flip_ar(deflated);
            std::vector<double> coordinates;
            if (pinned_solver(request) ||
                to_coordinates(map, candidate_physical, coordinates)) {
                const auto candidate = run_pass(request, map, candidate_physical, coordinates,
                    request.control.eps, interrupt);
                record.optimizer_status = candidate.optimizer_status;
                record.evaluations = candidate.evaluations;
                record.log_likelihood = candidate.log_likelihood;
                if (candidate.available && candidate.log_likelihood > current_ll) {
                    record.accepted = true;
                    current = candidate.coefficients;
                    current_ll = candidate.log_likelihood;
                }
            }
            // An inadmissible deflated start mirrors a pinned exception: the
            // pass is attempted, the candidate is kept out and nothing throws.
        }
    }
    passes.push_back(std::move(record));
}

// One pinned estimate() call: the fixed-precision optimization from the
// start plus the two unit-root-cancellation restart passes on its incumbent.
// Intermediate Multiple-start calls (prefix non-empty) run the identical
// chain at feps; the prefix only names their records.
struct ChainOutcome {
    PassOutcome best; // available holds for a retained candidate.
    bool ok = false;
    std::string producing_status;
    bool producing_converged = false;
    std::vector<EstimationPassRecord> records;
};

ChainOutcome run_chain(const SarimaEstimationRequest& request, const ArmaParameterMap& map,
    const std::array<std::vector<double>, 4>& physical, double precision,
    const std::string& prefix, const Interrupt& interrupt) {
    ChainOutcome chain;
    std::vector<double> coordinates;
    if (!pinned_solver(request) && !to_coordinates(map, physical, coordinates)) {
        throw std::invalid_argument("the pinned starting point is not admissible for this map");
    }
    EstimationPassRecord record{};
    record.stage = prefix + "estimate";
    record.attempted = true;
    record.gates_passed = true;
    const auto initial = run_pass(request, map, physical, coordinates, precision, interrupt);
    record.optimizer_status = initial.optimizer_status;
    record.evaluations = initial.evaluations;
    record.log_likelihood = initial.log_likelihood;
    record.accepted = initial.available; // The retained initial candidate.
    chain.records.push_back(std::move(record));
    chain.best = initial;
    chain.ok = initial.available;
    chain.producing_status = initial.optimizer_status;
    chain.producing_converged = pass_reports_convergence(initial.optimizer_status);
    if (!initial.available) return chain; // No pinned restarts without a base point.
    auto current = initial.coefficients;
    auto current_ll = initial.log_likelihood;
    for (const bool seasonal : {false, true}) {
        EstimationPassRecord pending{};
        pending.stage = prefix + (seasonal ? "ur_seasonal" : "ur_regular");
        run_cancelling_pass(chain.records, std::move(pending),
            request, map, seasonal, interrupt, current, current_ll);
        if (chain.records.back().accepted) {
            chain.producing_converged =
                pass_reports_convergence(chain.records.back().optimizer_status);
            chain.producing_status = chain.records.back().optimizer_status;
        }
    }
    chain.best.coefficients = std::move(current);
    chain.best.log_likelihood = current_ll;
    return chain;
}

// The pinned unit-root polynomial of SarimaUtility.differencingFilter applied
// to one series: (1 - B)^d (1 - B^period)^D through the qualified forward
// differencing, dropping the leading warmup. The pinned filter is that
// polynomial up to a global sign, which cancels in every ARMA quantity the
// Hannan-Rissanen pipeline computes.
std::vector<double> difference_series(std::vector<double> values, int d, int seasonal_d,
    int period) {
    for (int pass = 0; pass < d; ++pass) {
        for (std::size_t t = 0; t + 1 < values.size(); ++t)
            values[t] = values[t + 1] - values[t];
        values.pop_back();
    }
    const std::size_t stride = static_cast<std::size_t>(period);
    for (int pass = 0; pass < seasonal_d; ++pass) {
        for (std::size_t t = 0; t + stride < values.size(); ++t)
            values[t] = values[t + stride] - values[t];
        values.resize(values.size() - stride);
    }
    return values;
}

const char* hr_source_name(HannanRissanenStart::Source source) {
    switch (source) {
    case HannanRissanenStart::Source::Fitted: return "fitted";
    case HannanRissanenStart::Source::ZeroVariance: return "zerovariance";
    default: return "defaulted";
    }
}

} // namespace

void EstimationControl::validate() const {
    if (!std::isfinite(eps) || eps <= 0 || !std::isfinite(feps) || feps <= 0) {
        throw std::invalid_argument("estimation precisions must be finite and positive");
    }
    if (!maximum_likelihood) {
        throw std::invalid_argument(
            "only the maximum-likelihood objective form is qualified for estimation");
    }
}

SarimaEstimationResult estimate_sarima_fixed(
    const SarimaEstimationRequest& request, const Interrupt& interrupt) {
    request.control.validate();
    for (const auto& block : request.blocks) {
        if (block.fixed.size() != block.initial.size()) {
            throw std::invalid_argument(
                "Sarima coefficient states must match their initial values");
        }
    }
    const ArmaParameterMap map(request.blocks);
    std::size_t free_count = 0;
    for (const auto& block : request.blocks) {
        for (const auto flag : block.fixed) free_count += flag ? 0 : 1;
    }
    validate_arma_objective_inputs(request.series, request.design, map, request.period,
        request.d, request.seasonal_d, 0.0, 1e-10, 0.0, interrupt);

    SarimaEstimationResult result;
    if (free_count == 0) { // Pinned parametersCount/dimension guard.
        std::array<std::vector<double>, 4> fixed_coefficients;
        for (std::size_t b = 0; b < 4; ++b) fixed_coefficients[b] = request.blocks[b].initial;
        const auto objective = evaluate_arma_physical_objective(request.series, request.design,
            map, fixed_coefficients, request.period, request.d, request.seasonal_d, 0.0, 1e-10,
            0.0, interrupt);
        result.status = "likelihood_only";
        result.converged = false;
        result.coefficients = objective.parameters.coefficients;
        result.log_likelihood = objective.status == "available" ? objective.adjusted_log_likelihood
                                                                : unavailable();
        return result;
    }

    const bool multiple =
        request.control.starting_point == EstimationStartingPoint::Multiple;
    const bool hr_family = multiple ||
        request.control.starting_point == EstimationStartingPoint::HannanRissanen;

    // The pinned initializer runs on the differenced linear model with all
    // columns estimated and no missing positions; masked series and fixed
    // regression coefficients fail explicitly instead.
    if (hr_family) {
        if (std::any_of(request.series.missing.begin(), request.series.missing.end(),
                [](std::uint8_t flag) { return flag != 0; })) {
            throw std::invalid_argument(std::string("the ") + stage_name(request.control.starting_point) +
                " starting point is only qualified for complete series");
        }
        if (std::any_of(request.design.fixed.begin(), request.design.fixed.end(),
                [](std::uint8_t flag) { return flag != 0; })) {
            throw std::invalid_argument(std::string("the ") + stage_name(request.control.starting_point) +
                " starting point is only qualified for fully free regression designs");
        }
        const std::size_t warmup = static_cast<std::size_t>(request.d) +
            static_cast<std::size_t>(request.seasonal_d) * static_cast<std::size_t>(request.period);
        if (request.series.values.size() <= warmup) {
            throw std::invalid_argument(
                "the response is too short to difference for the Hannan-Rissanen starts");
        }
    }

    // The starting coefficients of the final pass. For the Zero/Default
    // policy the pinned fill overwrites every free entry; fixed entries are
    // preserved exactly as the SarimaFixedMapping seam of the passes.
    std::array<std::vector<double>, 4> physical;
    if (!hr_family) {
        const bool zero_start = request.control.starting_point == EstimationStartingPoint::Zero;
        for (std::size_t b = 0; b < 4; ++b) {
            const double fill = zero_start ? 0.0 : (b == 1 || b == 3 ? -0.2 : 0.1);
            physical[b].assign(request.blocks[b].initial.size(), fill);
        }
    } else {
        // The meanCorrection seam (the pinned airline-model constant column
        // lives on the differenced scale) enters the initializer as that
        // constant, on the single-polynomial differenced response the pinned
        // deviances machine runs on (the JVM reads it from the same
        // differenced model). Every other design keeps the sequential
        // per-pass differencing of the 22-series seam.
        const bool mean_correction = request.design.columns == 1 &&
            request.design.values.size() == request.series.values.size() &&
            std::all_of(request.design.values.begin(), request.design.values.end(),
                [](double value) { return value == 1.0; });
        const auto differenced_response = mean_correction
            ? detail::pinned_backfilter_apply(detail::pinned_differencing_polynomial(
                  request.period, request.d, request.seasonal_d),
                request.series.values)
            : difference_series(
                request.series.values, request.d, request.seasonal_d, request.period);
        RegressionDesign differenced;
        differenced.columns = request.design.columns;
        differenced.coefficients = request.design.coefficients;
        differenced.fixed = request.design.fixed;
        if (mean_correction) {
            differenced.values.assign(differenced_response.size(), 1.0);
        } else {
            for (std::size_t c = 0; c < request.design.columns; ++c) {
                std::vector<double> column(request.design.values.begin() +
                        static_cast<std::ptrdiff_t>(c) * request.series.values.size(),
                    request.design.values.begin() +
                        static_cast<std::ptrdiff_t>(c + 1) * request.series.values.size());
                auto filtered = difference_series(std::move(column), request.d,
                    request.seasonal_d, request.period);
                differenced.values.insert(
                    differenced.values.end(), filtered.begin(), filtered.end());
            }
        }
        const HannanRissanenOrders orders{static_cast<int>(request.blocks[0].initial.size()),
            static_cast<int>(request.blocks[2].initial.size()),
            static_cast<int>(request.blocks[1].initial.size()),
            static_cast<int>(request.blocks[3].initial.size()), request.period};
        const auto start = hannan_rissanen_start(
            differenced_response, differenced, orders, interrupt);
        // The HR unit hands back [ar, sar, ma, sma] (the m_pi extraction
        // order); the request blocks are [ar, ma, sar, sma].
        const std::size_t permutation[4] = {0, 2, 1, 3};
        for (std::size_t b = 0; b < 4; ++b) {
            physical[b] = start.coefficients[permutation[b]];
            for (std::size_t j = 0; j < physical[b].size(); ++j) {
                if (request.blocks[b].fixed[j]) physical[b][j] = request.blocks[b].initial[j];
            }
        }
        EstimationPassRecord initialize_record{};
        initialize_record.stage = "initialize";
        initialize_record.attempted = true;
        initialize_record.gates_passed = true;
        // Provisional: the Multiple selection below rewrites the bit; the
        // HannanRissanen policy always retains the initialized start.
        initialize_record.accepted = !multiple;
        initialize_record.detail = std::string(hr_source_name(start.source)) +
            (start.stabilized ? ",stabilized" : "");
        result.passes.push_back(std::move(initialize_record));
    }
    for (std::size_t b = 0; b < 4; ++b) {
        for (std::size_t j = 0; j < physical[b].size(); ++j) {
            if (request.blocks[b].fixed[j]) physical[b][j] = request.blocks[b].initial[j];
        }
    }

    if (multiple) {
        // Pinned Multiple policy: two internal-precision estimates around
        // the HR start, est2 winning only on a strictly better likelihood
        // (the dcmpg comparison keeps est1 on ties and NaN), and the winner
        // seeding the final eps pass.
        auto first = run_chain(request, map, physical, request.control.feps, "multiple_first.",
            interrupt);
        std::array<std::vector<double>, 4> alternate = first.best.coefficients;
        if (first.ok) {
            const auto& orders = request.blocks;
            if (!orders[0].initial.empty() && !orders[1].initial.empty()) {
                // Pinned restart from a sign-flipped first MA coefficient
                // and zeroed interior AR/MA coefficients (theta(1) <- -/+
                // 0.9 by its sign, phi(1..p) <- 0, theta(2..q) <- 0);
                // seasonal blocks carry over.
                alternate[1][0] = first.best.coefficients[1][0] < 0.0 ? 0.9 : -0.9;
                std::fill(alternate[0].begin(), alternate[0].end(), 0.0);
                for (std::size_t i = 1; i < alternate[1].size(); ++i) alternate[1][i] = 0.0;
            }
            for (std::size_t b = 0; b < 4; ++b) {
                for (std::size_t j = 0; j < alternate[b].size(); ++j) {
                    if (request.blocks[b].fixed[j])
                        alternate[b][j] = request.blocks[b].initial[j];
                }
            }
        }
        auto second = first.ok
            ? run_chain(request, map, alternate, request.control.feps, "multiple_second.",
                  interrupt)
            : ChainOutcome{};
        // Selection. The pinned dcmpg gives est1 on ties and on any NaN;
        // where the pinned would rebuild from an invalid winner, Kalos keeps
        // the best available candidate and falls back to the HR start (a
        // documented fail-closed deviation, see the contract).
        enum { WinnerHr, WinnerFirst, WinnerSecond } winner = WinnerHr;
        if (first.ok) {
            winner = (!second.ok || !(first.best.log_likelihood < second.best.log_likelihood))
                ? WinnerFirst : WinnerSecond;
        } else if (second.ok) {
            winner = WinnerSecond;
        }
        first.records.front().accepted = winner == WinnerFirst;
        if (first.ok) second.records.front().accepted = winner == WinnerSecond;
        result.passes.front().accepted = winner == WinnerHr;
        physical = winner == WinnerFirst ? first.best.coefficients
                 : winner == WinnerSecond ? second.best.coefficients : physical;
        for (auto& record : first.records) result.passes.insert(result.passes.end(), record);
        if (first.ok) {
            for (auto& record : second.records)
                result.passes.insert(result.passes.end(), record);
        }
    }

    auto chain = run_chain(request, map, physical, request.control.eps, "", interrupt);
    for (auto& record : chain.records) result.passes.push_back(std::move(record));
    if (!chain.ok) {
        result.status = chain.producing_status;
        result.converged = false;
        result.coefficients = std::move(physical);
        // No candidate was retained, so there is no likelihood to report. The
        // field has no default initializer: it must be written explicitly, or
        // the result carries the last value this storage held (the M13 stale-
        // likelihood defect: an initial_unavailable chain returned the
        // process's last-computed log-likelihood from the reused stack slot).
        result.log_likelihood = unavailable();
        return result;
    }
    result.coefficients = std::move(chain.best.coefficients);
    result.log_likelihood = chain.best.log_likelihood;
    result.converged = chain.producing_converged;
    result.status = chain.producing_converged ? "converged" : std::move(chain.producing_status);
    return result;
}

SarimaEstimationResult optimize_sarima_from(const SarimaEstimationRequest& request,
    const std::array<std::vector<double>, 4>& start, const Interrupt& interrupt) {
    request.control.validate();
    for (const auto& block : request.blocks) {
        if (block.fixed.size() != block.initial.size()) {
            throw std::invalid_argument(
                "Sarima coefficient states must match their initial values");
        }
    }
    const ArmaParameterMap map(request.blocks);
    std::size_t free_count = 0;
    for (const auto& block : request.blocks) {
        for (const auto flag : block.fixed) free_count += flag ? 0 : 1;
    }
    validate_arma_objective_inputs(request.series, request.design, map, request.period,
        request.d, request.seasonal_d, 0.0, 1e-10, 0.0, interrupt);
    for (std::size_t b = 0; b < 4; ++b) {
        if (start[b].size() != request.blocks[b].initial.size()) {
            throw std::invalid_argument(
                "the warm start must match the request coefficient blocks");
        }
    }

    // The pinned warm start is the incumbent SarimaModel, so fixed entries
    // carry their request value and the free entries the given ones; the map
    // seam keeps them exactly as estimate_sarima_fixed does.
    std::array<std::vector<double>, 4> physical = start;
    for (std::size_t b = 0; b < 4; ++b) {
        for (std::size_t j = 0; j < physical[b].size(); ++j) {
            if (request.blocks[b].fixed[j]) physical[b][j] = request.blocks[b].initial[j];
        }
    }

    SarimaEstimationResult result;
    if (free_count == 0) { // Pinned parametersCount/dimension guard.
        std::array<std::vector<double>, 4> fixed_coefficients;
        for (std::size_t b = 0; b < 4; ++b) fixed_coefficients[b] = request.blocks[b].initial;
        const auto objective = evaluate_arma_physical_objective(request.series, request.design,
            map, fixed_coefficients, request.period, request.d, request.seasonal_d, 0.0, 1e-10,
            0.0, interrupt);
        result.status = "likelihood_only";
        result.converged = false;
        result.coefficients = objective.parameters.coefficients;
        result.log_likelihood = objective.status == "available" ? objective.adjusted_log_likelihood
                                                                : unavailable();
        return result;
    }

    std::vector<double> coordinates;
    if (!pinned_solver(request) && !to_coordinates(map, physical, coordinates)) {
        throw std::invalid_argument("the warm start is not admissible for this map");
    }
    EstimationPassRecord record{};
    record.stage = "estimate";
    record.attempted = true;
    record.gates_passed = true;
    const auto pass = run_pass(request, map, physical, coordinates, request.control.eps,
        interrupt);
    record.optimizer_status = pass.optimizer_status;
    record.evaluations = pass.evaluations;
    record.log_likelihood = pass.log_likelihood;
    record.accepted = pass.available;
    result.passes.push_back(std::move(record));
    if (!pass.available) {
        result.status = pass.optimizer_status;
        result.converged = false;
        result.coefficients = std::move(physical);
        // Same explicit missing likelihood as estimate_sarima_fixed: the
        // field's default is indeterminate and this path retains no candidate.
        result.log_likelihood = unavailable();
        return result;
    }
    result.coefficients = std::move(pass.coefficients);
    result.log_likelihood = pass.log_likelihood;
    result.converged = pass.optimizer_status == "converged_gradient";
    result.status = result.converged ? "converged" : std::move(pass.optimizer_status);
    if (request.request_landing_point && pass.landing_derivatives) {
        result.has_landing_point = true;
    }
    return result;
}

} // namespace kalos
