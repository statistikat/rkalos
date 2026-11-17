#include "kalos/estimation_control.hpp"

#include "kalos/arma_optimizer.hpp"
#include "kalos/hannan_rissanen.hpp"

#include <algorithm>
#include <cmath>
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
            const double coordinate = std::atanh(reflection / (1 - map.margin()));
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
};

PassOutcome run_pass(const SarimaEstimationRequest& request, const ArmaParameterMap& map,
    const std::vector<double>& coordinates, double precision, const Interrupt& interrupt) {
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
    record.gates_passed = !current[0].empty() && !current[1].empty() &&
        (seasonal ? request.seasonal_d : request.d) > 0;
    if (record.gates_passed) {
        const auto pinned = flip_ar(current);
        record.gates_passed = within_limit(pinned_gate_sum(pinned[0], seasonal)) &&
            within_limit(pinned_gate_sum(pinned[1], seasonal));
        if (record.gates_passed) {
            auto deflated = flip_ar(current);
            if (seasonal) deflate_seasonal(deflated); else deflate_regular(deflated);
            record.attempted = true;
            std::vector<double> coordinates;
            if (to_coordinates(map, flip_ar(deflated), coordinates)) {
                const auto candidate = run_pass(request, map, coordinates, request.control.eps,
                    interrupt);
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
    if (!to_coordinates(map, physical, coordinates)) {
        throw std::invalid_argument("the pinned starting point is not admissible for this map");
    }
    EstimationPassRecord record{};
    record.stage = prefix + "estimate";
    record.attempted = true;
    record.gates_passed = true;
    const auto initial = run_pass(request, map, coordinates, precision, interrupt);
    record.optimizer_status = initial.optimizer_status;
    record.evaluations = initial.evaluations;
    record.log_likelihood = initial.log_likelihood;
    record.accepted = initial.available; // The retained initial candidate.
    chain.records.push_back(std::move(record));
    chain.best = initial;
    chain.ok = initial.available;
    chain.producing_status = initial.optimizer_status;
    chain.producing_converged = initial.optimizer_status == "converged_gradient";
    if (!initial.available) return chain; // No pinned restarts without a base point.
    auto current = initial.coefficients;
    auto current_ll = initial.log_likelihood;
    for (const bool seasonal : {false, true}) {
        EstimationPassRecord pending{};
        pending.stage = prefix + (seasonal ? "ur_seasonal" : "ur_regular");
        run_cancelling_pass(chain.records, std::move(pending),
            request, map, seasonal, interrupt, current, current_ll);
        if (chain.records.back().accepted) {
            chain.producing_converged = chain.records.back().optimizer_status == "converged_gradient";
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
        RegressionDesign differenced;
        differenced.columns = request.design.columns;
        differenced.coefficients = request.design.coefficients;
        differenced.fixed = request.design.fixed;
        for (std::size_t c = 0; c < request.design.columns; ++c) {
            std::vector<double> column(request.design.values.begin() +
                    static_cast<std::ptrdiff_t>(c) * request.series.values.size(),
                request.design.values.begin() +
                    static_cast<std::ptrdiff_t>(c + 1) * request.series.values.size());
            auto filtered = difference_series(std::move(column), request.d, request.seasonal_d,
                request.period);
            differenced.values.insert(differenced.values.end(), filtered.begin(), filtered.end());
        }
        const HannanRissanenOrders orders{static_cast<int>(request.blocks[0].initial.size()),
            static_cast<int>(request.blocks[2].initial.size()),
            static_cast<int>(request.blocks[1].initial.size()),
            static_cast<int>(request.blocks[3].initial.size()), request.period};
        const auto start = hannan_rissanen_start(
            difference_series(request.series.values, request.d, request.seasonal_d, request.period),
            differenced, orders, interrupt);
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
        return result; // log_likelihood stays NaN: no candidate was retained.
    }
    result.coefficients = std::move(chain.best.coefficients);
    result.log_likelihood = chain.best.log_likelihood;
    result.converged = chain.producing_converged;
    result.status = chain.producing_converged ? "converged" : std::move(chain.producing_status);
    return result;
}

} // namespace kalos
