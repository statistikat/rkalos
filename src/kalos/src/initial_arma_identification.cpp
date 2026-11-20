#include "kalos/initial_arma_identification.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace kalos {

namespace {

// The ArmaModule grid/panel/refinement scalars, shared by the two public
// request shapes (initial identification and the 21-F select seam). The
// prefix names the request the caller validated.
void validate_module_scalars(const char* prefix, int period, int models_count,
    int max_p, int max_q, int max_bp, int max_bq, double eps) {
    if (period != 4 && period != 12) {
        throw std::invalid_argument(
            std::string(prefix) + ": the period must be 4 or 12 (the "
            "qualified estimation-seam domain)");
    }
    if (models_count < 1) {
        throw std::invalid_argument(std::string(prefix)
            + ": models_count must be >= 1 (the pinned estimations array "
              "holds nmod slots)");
    }
    if (max_p < 0 || max_q < 0 || max_bp < 0 || max_bq < 0) {
        throw std::invalid_argument(std::string(prefix)
            + ": the grid bounds must not be negative");
    }
    if (!std::isfinite(eps) || !(eps > 0.0)) {
        throw std::invalid_argument(std::string(prefix)
            + ": the estimation precision must be finite and positive");
    }
}

} // namespace

void InitialArmaIdentificationRequest::validate() const {
    if (series.empty()) {
        throw std::invalid_argument(
            "InitialArmaIdentificationRequest: the series must not be empty");
    }
    if (!std::all_of(series.begin(), series.end(),
            [](double v) { return std::isfinite(v); })) {
        throw std::invalid_argument(
            "InitialArmaIdentificationRequest: the series must be finite");
    }
    for (const auto& column : regressors) {
        if (column.size() != series.size()) {
            throw std::invalid_argument(
                "InitialArmaIdentificationRequest: every regressor must have "
                "series.size() rows");
        }
        if (!std::all_of(column.begin(), column.end(),
                [](double v) { return std::isfinite(v); })) {
            throw std::invalid_argument(
                "InitialArmaIdentificationRequest: every regressor must be finite");
        }
    }
    if (d < 0 || seasonal_d < 0) {
        throw std::invalid_argument(
            "InitialArmaIdentificationRequest: the differencing orders must not "
            "be negative");
    }
    validate_module_scalars("InitialArmaIdentificationRequest", period,
        models_count, max_p, max_q, max_bp, max_bq, eps);
}

void ArmaSelectRequest::validate() const {
    if (cursor.empty()) {
        throw std::invalid_argument(
            "ArmaSelectRequest: the cursor must not be empty (the pinned "
            "module never selects over an empty sequence)");
    }
    if (!std::all_of(cursor.begin(), cursor.end(),
            [](double v) { return std::isfinite(v); })) {
        throw std::invalid_argument(
            "ArmaSelectRequest: the cursor must be finite");
    }
    if (d < 0 || seasonal_d < 0) {
        throw std::invalid_argument(
            "ArmaSelectRequest: the differencing orders must not be negative");
    }
    validate_module_scalars(
        "ArmaSelectRequest", period, models_count, max_p, max_q, max_bp, max_bq, eps);
}

namespace {

// Pinned literal: the RegArmaBic NO_BIC sentinel (the est == null branch that
// fills it is dead code in the pinned 3.9.0 bytecode; the sort's drop test and
// the refinement's tail trim still key on the literal).
constexpr double kNoBic = 99999.0;

void require(bool condition, const char* where, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(std::string(where) + ": " + message);
    }
}

// Double.compare(a, b): the total order where -0.0 < +0.0 and NaN is the
// greatest element (RegArmaBic.compareTo drives the stable TimSort).
int double_compare(double a, double b) {
    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    if (a == b) {
        if (std::signbit(a) == std::signbit(b)) {
            return 0;
        }
        return std::signbit(a) ? -1 : 1;
    }
    const bool an = std::isnan(a);
    const bool bn = std::isnan(b);
    if (an && bn) {
        return 0;
    }
    return an ? 1 : -1;
}

// SarmaOrders: value equality covers all five fields (the pinned equals).
struct SarmaOrders {
    int period = 0;
    int p = 0;
    int q = 0;
    int bp = 0;
    int bq = 0;
    bool equals(const SarmaOrders& other) const {
        return period == other.period && p == other.p && q == other.q
            && bp == other.bp && bq == other.bq;
    }
    int parameters_count() const { return p + bp + q + bq; }
    std::string text() const {
        return "(" + std::to_string(p) + "," + std::to_string(q) + ")("
            + std::to_string(bp) + "," + std::to_string(bq) + ")";
    }
};

// One RegArmaBic evaluation as the sort/merge/refinement machinery uses it.
struct RegArmaBic {
    SarmaOrders spec;
    bool active = false; // the pinned array slot was filled (null check)
    int stage = 0;
    int slot = 0;
    std::string initialSource;
    std::string fitStatus;
    double logLikelihood = std::numeric_limits<double>::quiet_NaN();
    double bic = kNoBic;
    std::array<std::vector<double>, 4> coefficients; // Kalos blocks [ar, ma, sar, sma]
};

bool all_finite(const std::vector<double>& v) {
    return std::all_of(v.begin(), v.end(),
        [](double x) { return std::isfinite(x); });
}

// HannanRissanenOrders(p, bp, q, bq, period) fit through the qualified
// estimate_sarima_fixed with the pinned sum-of-squares solver: the HR
// initialization plus one eps-precision pinned-Levenberg-Marquardt chain over
// the qualified deviances machine, on the cursor with NO regression design
// (the pinned candidate model is RegArimaModel.builder().y(cursor).arima(...),
// no meanCorrection) and d = bd = 0 - exactly GlsSarimaComputer.process,
// which carries no starting-point policy and no unit-root restarts. The
// pinned machine has no null return path: any numerical failure throws,
// mirrored here as std::runtime_error (the ProcessingResult.Failed
// equivalence, never a silent fallback).
RegArmaBic evaluate_candidate(const SarmaOrders& orders, const std::vector<double>& cursor,
    double eps, int stage, int slot, const Interrupt& interrupt, const char* where) {
    const std::size_t n = cursor.size();
    SarimaEstimationRequest request;
    request.series.values = cursor;
    request.series.missing.assign(n, 0);
    // No design at all: the mean lives in the cursor upstream (the OLS
    // residual cursor already carries the mean correction), and the pinned
    // candidate fits a design-free ARMA. The qualified seam runs the
    // design-free (sfilter) deviances machine here.
    request.design.columns = 0;
    request.period = orders.period;
    request.d = 0;
    request.seasonal_d = 0;
    request.blocks[0].initial.assign(static_cast<std::size_t>(orders.p), 0.0);
    request.blocks[1].initial.assign(static_cast<std::size_t>(orders.q), 0.0);
    request.blocks[2].initial.assign(static_cast<std::size_t>(orders.bp), 0.0);
    request.blocks[3].initial.assign(static_cast<std::size_t>(orders.bq), 0.0);
    for (auto& block : request.blocks) {
        block.fixed.assign(block.initial.size(), 0);
    }
    request.control = EstimationControl{eps, 1.0e-4, true, false, true, true,
        EstimationStartingPoint::HannanRissanen, EstimationSolver::PinnedLevenbergMarquardt};

    const SarimaEstimationResult result = estimate_sarima_fixed(request, interrupt);
    require(result.status != "start_not_evaluable", where,
        "candidate " + orders.text()
            + ": the pinned machine rejects the start point (the pinned "
              "RegArmaSsqFunction evaluation throws before minimization, the "
              "ProcessingResult.Failed path)");
    require(std::isfinite(result.log_likelihood), where,
        "candidate " + orders.text()
            + ": the pinned machine rejects the landing point (the pinned "
              "DEFAULT_COMPUTER.compute at the landing throws, the "
              "ProcessingResult.Failed path)");
    for (const auto& block : result.coefficients) {
        require(all_finite(block), where,
            "candidate " + orders.text() + ": the landing carries non-finite coefficients");
    }

    // BIC2 = ((-2 logLikelihood) + k log n) / n with the FdLibm log, n the
    // cursor length (no differencing, no missing) and k = (p+bp+q+bq) + 1:
    // the ARMA parameters plus the innovation variance and NO mean - the
    // candidate model has no regressors (RegArimaEstimation.statistics plus
    // LikelihoodStatistics.calc; the +1 is the variance degree of freedom).
    // The ll is the machine's concentrated likelihood at the landing, which
    // the seam evaluates bit-identically (21-D.2).
    RegArmaBic evaluation;
    evaluation.spec = orders;
    evaluation.active = true;
    evaluation.stage = stage;
    evaluation.slot = slot;
    for (const auto& pass : result.passes) {
        if (pass.stage == "initialize") {
            evaluation.initialSource = pass.detail;
        }
    }
    evaluation.fitStatus = result.status;
    evaluation.logLikelihood = result.log_likelihood;
    evaluation.coefficients = result.coefficients;
    const double nd = static_cast<double>(n);
    const double k = static_cast<double>(orders.parameters_count() + 1);
    evaluation.bic = ((-2.0 * result.log_likelihood) + k * detail::fdlibm_log(nd)) / nd;
    return evaluation;
}

// RegArmaBic.sort: evaluate all grid candidates in grid order (every
// evaluation recorded, including those the sort drops), drop the NO_BIC
// sentinel (a NaN bic is kept - the pinned dcmpl test), then sort ascending
// in the Double.compare total order with insertion order preserved on bit-
// equal BIC (stable TimSort).
std::vector<RegArmaBic> sort_candidates(const std::vector<SarmaOrders>& grid,
    const std::vector<double>& cursor, double eps, int stage,
    std::vector<RegArmaBic>& evaluations, const Interrupt& interrupt, const char* where) {
    std::vector<RegArmaBic> sorted;
    int slot = 0;
    for (const auto& orders : grid) {
        RegArmaBic evaluation =
            evaluate_candidate(orders, cursor, eps, stage, slot++, interrupt, where);
        evaluations.push_back(evaluation);
        if (evaluation.bic != kNoBic) { // NaN kept
            sorted.push_back(std::move(evaluation));
        }
    }
    std::stable_sort(sorted.begin(), sorted.end(),
        [](const RegArmaBic& x, const RegArmaBic& y) {
            return double_compare(x.bic, y.bic) < 0;
        });
    return sorted;
}

// RegArmaBic.mergeInto: merge the ascending list into the fixed nmod-slot
// estimations array: an empty slot takes the entry, a duplicate specification
// is dropped, a better entry inserts with an upward shift that evicts the
// worst tail. The merge stops when the sorted list or the array is done.
void merge_into(
    const std::vector<RegArmaBic>& sorted, std::vector<RegArmaBic>& estimations) {
    const std::size_t n = estimations.size();
    std::size_t j = 0;
    const std::size_t m = std::min(sorted.size(), n);
    for (std::size_t i = 0; i < m && j < n; ++i) {
        const double bic = sorted[i].bic;
        for (std::size_t k = j; k < n;) {
            if (!estimations[k].active) { // the pinned null slot
                estimations[k] = sorted[i];
                j = k + 1;
                break;
            }
            if (estimations[k].spec.equals(sorted[i].spec)) {
                j = k + 1; // duplicate specification: dropped
                break;
            }
            if (estimations[k].bic > bic) { // dcmpl: NaN entries keep scanning
                for (std::size_t l = n - 1; l > k; --l) {
                    estimations[l] = estimations[l - 1];
                }
                estimations[k] = sorted[i];
                j = k + 1;
                break;
            }
            ++k;
        }
    }
}

// The pinned ArmaModule.select(cursor, period, d, bd): the three-stage grid
// (bp, bq) / (p, q) / (bp, bq) with the mixed gate and the bd clamp, the
// sort/merge panel and the eleven-rule BIC-margin refinement. The request is
// validated by the public entry point; this engine assumes validity.
ArmaSelectResult select_engine(
    const ArmaSelectRequest& request, const Interrupt& interrupt, const char* where) {
    const std::vector<double>& cursor = request.cursor;

    // ArmaModule.select(cursor, period, d, bd). The request domain is
    // period in {4, 12}; the pinned period == 1 shortcuts stay unreachable.
    const int period = request.period;
    const int nmod = request.models_count;
    std::vector<RegArmaBic> estimations(static_cast<std::size_t>(nmod));
    std::vector<RegArmaBic> evaluations;

    SarmaOrders tmpl;
    tmpl.period = period;
    tmpl.p = period == 2 ? 1 : 3; // the pinned template (3, 0) for period 4/12
    tmpl.q = 0;

    const auto mixed_ok = [request](int a, int b) {
        return request.mixed || a == 0 || b == 0;
    };

    // Stage 1: the (bp, bq) grid on the (3, 0) regular part; its winner's
    // (bp, bq) is all stage 2 carries - the stage-1 winner itself (p = 3
    // exceeds the stage-2 grid bound) is never merged.
    SarmaOrders best = tmpl;
    {
        std::vector<SarmaOrders> grid;
        for (int bp = 0; bp <= request.max_bp; ++bp) {
            for (int bq = 0; bq <= request.max_bq; ++bq) {
                if (!mixed_ok(bp, bq)) {
                    continue;
                }
                SarmaOrders orders = tmpl;
                orders.bp = bp;
                orders.bq = bq;
                grid.push_back(orders);
            }
        }
        const std::vector<RegArmaBic> sorted =
            sort_candidates(grid, cursor, request.eps, 1, evaluations, interrupt, where);
        require(!sorted.empty(), where, "the seasonal grid produced no candidate");
        best = sorted[0].spec;
    }

    // Stage 2: the (p, q) grid on the stage-1 winner's (bp, bq); its winner
    // is adopted and the sorted list is merged into the fresh array (the
    // first merge degenerates to copying the leading min(m, nmod) entries).
    {
        std::vector<SarmaOrders> grid;
        for (int p = 0; p <= request.max_p; ++p) {
            for (int q = 0; q <= request.max_q; ++q) {
                if (!mixed_ok(p, q)) {
                    continue;
                }
                SarmaOrders orders = best;
                orders.p = p;
                orders.q = q;
                grid.push_back(orders);
            }
        }
        const std::vector<RegArmaBic> sorted =
            sort_candidates(grid, cursor, request.eps, 2, evaluations, interrupt, where);
        require(!sorted.empty(), where, "the regular grid produced no candidate");
        best = sorted[0].spec;
        merge_into(sorted, estimations);
    }

    // Stage 3: the (bp, bq) refinement grid on the stage-2 winner; the
    // winner is adopted only through the merged list and the refinement
    // scan. The stage-3 bp clamp keys on the seasonal differencing order
    // (bd == 1), never on the stored, never-read acceptWhiteNoise field.
    {
        int m1 = request.max_bp;
        if (request.seasonal_d == 1) {
            m1 = 0;
        }
        std::vector<SarmaOrders> grid;
        for (int bp = 0; bp <= m1; ++bp) {
            for (int bq = 0; bq <= request.max_bq; ++bq) {
                if (!mixed_ok(bp, bq)) {
                    continue;
                }
                SarmaOrders orders = best;
                orders.bp = bp;
                orders.bq = bq;
                grid.push_back(orders);
            }
        }
        const std::vector<RegArmaBic> sorted =
            sort_candidates(grid, cursor, request.eps, 3, evaluations, interrupt, where);
        require(!sorted.empty(), where,
            "the seasonal refinement grid produced no candidate");
        merge_into(sorted, estimations);
    }

    // All nmod slots must be filled: the pinned dead tail-trim then NPEs on
    // an empty slot, which the module turns into ProcessingResult.Failed.
    for (const auto& entry : estimations) {
        require(entry.active, where,
            "the merged list cannot fill all models_count slots (the pinned "
            "refinement would NPE on the empty slot, the "
            "ProcessingResult.Failed path)");
    }
    const int n = nmod;

    // select(cursor, d, bd): the refinement scan over the merged list. The
    // polarity measures count the unit roots the spec differencing implies
    // against MA cancellation; the threshold band keys on the merged list's
    // BIC spread and shrinks by each accepted move's gain.
    SarmaOrders incumbent = estimations[0].spec;
    int npq = incumbent.p + incumbent.q;
    int nbpq = incumbent.bp + incumbent.bq;
    int dpq = std::abs(incumbent.p + request.d - incumbent.q);
    int dbpq = std::abs(incumbent.bp + request.seasonal_d - incumbent.bq);
    double dif = estimations[static_cast<std::size_t>(n - 1)].bic - estimations[0].bic;
    if (dif < 0.003) { // dcmpg: a NaN spread takes 1.0
        dif = 0.0625;
    } else if (dif < 0.03) {
        dif = 0.25;
    } else {
        dif = 1.0;
    }
    double lim1 = 0.01 * dif;
    double lim2 = 0.0025 * dif;
    double lim3 = 0.0075 * dif;
    int id = 0;
    int cnt = 0;
    std::vector<ArmaSelectionStep> steps;

    for (int i = 1; i < n; ++i) {
        if (interrupt) {
            interrupt();
        }
        const SarmaOrders& cur = estimations[static_cast<std::size_t>(i)].spec;
        const int cpq = cur.p + cur.q;
        const int cbpq = cur.bp + cur.bq;
        const int cdpq = std::abs(cur.p + request.d - cur.q);
        const int cdbpq = std::abs(cur.bp + request.seasonal_d - cur.bq);
        // curdif is measured against the INCUMBENT, not the list head (the
        // pinned local 26 accumulates from the incumbent's BIC).
        const double curdif = estimations[static_cast<std::size_t>(i)].bic
            - estimations[static_cast<std::size_t>(id)].bic;
        int rule = 0;
        const bool polar = cdpq < dpq || cdbpq < dbpq;
        if (request.balanced && polar && npq == cpq && nbpq == cbpq && curdif <= lim1) {
            rule = 1;
        } else if (request.balanced && cdpq < dpq && cpq <= npq && cbpq == nbpq
            && curdif <= lim2 && cur.p > 0 && cur.q > 0) {
            rule = 2;
        } else if (request.balanced
            && ((cdpq == 0 && cdpq < dpq && request.d > 0)
                || (cdbpq == 0 && cdbpq < dbpq && request.seasonal_d > 0))
            && npq == cpq && nbpq == cbpq && curdif <= lim1) {
            rule = 3;
        } else if (request.balanced && cdpq == 0 && cdbpq == 0 && curdif < lim2) {
            rule = 4;
        } else if (request.balanced && cpq > npq && cdpq == 0 && cbpq == nbpq
            && curdif < lim2) {
            rule = 5;
        } else if (request.balanced && cbpq > nbpq && cdbpq == 0 && cpq == npq
            && curdif < lim2) {
            rule = 6;
        } else if (request.balanced && cbpq < nbpq && cbpq > 0 && cpq == npq
            && cdbpq == 0 && curdif < lim2) {
            rule = 7;
        } else if (i == 1 && npq == 0 && cpq == 1 && cbpq == nbpq && curdif < lim2) {
            rule = 8;
        } else if (cpq < npq && cpq > 0 && cbpq == nbpq && curdif < lim2) {
            rule = 9;
        } else if (cbpq < nbpq && cbpq > 0 && cpq == npq && curdif < lim2) {
            rule = 10;
        } else if (cur.p < incumbent.p && cur.q == incumbent.q && cpq > 0
            && cbpq == nbpq && curdif < lim3) {
            rule = 11;
        }
        if (rule != 0) {
            ++cnt;
            const double gain =
                estimations[static_cast<std::size_t>(i)].bic - estimations[0].bic;
            lim1 -= gain;
            lim2 -= gain;
            lim3 -= gain;
            npq = cpq;
            nbpq = cbpq;
            dpq = cdpq;
            dbpq = cdbpq;
            id = i;
            incumbent = cur;
        }
        steps.push_back(ArmaSelectionStep{i, rule, id, cnt});
    }

    // A zero-parameter incumbent escapes to the next merged entry whenever
    // one exists behind it.
    SarmaOrders selected = estimations[static_cast<std::size_t>(id)].spec;
    if (selected.parameters_count() == 0 && id <= nmod - 2) {
        selected = estimations[static_cast<std::size_t>(id + 1)].spec;
    }

    ArmaSelectResult out;
    out.p = selected.p;
    out.q = selected.q;
    out.bp = selected.bp;
    out.bq = selected.bq;
    for (const auto& evaluation : evaluations) {
        out.evaluations.push_back(ArmaCandidateEvaluation{evaluation.stage, evaluation.slot,
            evaluation.spec.p, evaluation.spec.q, evaluation.spec.bp, evaluation.spec.bq,
            evaluation.initialSource, evaluation.fitStatus, evaluation.logLikelihood,
            evaluation.bic, evaluation.coefficients});
    }
    for (const auto& entry : estimations) {
        out.preferred_models.push_back(ArmaPreferredModel{entry.spec.p, entry.spec.q,
            entry.spec.bp, entry.spec.bq, entry.bic, entry.logLikelihood,
            entry.coefficients});
    }
    out.refinement = steps;
    return out;
}

} // namespace

ArmaSelectResult select_arma_orders(
    const ArmaSelectRequest& request, const Interrupt& interrupt) {
    request.validate();
    if (interrupt) {
        interrupt();
    }
    return select_engine(request, interrupt, "select_arma_orders");
}

InitialArmaIdentification identify_initial_arma(
    const InitialArmaIdentificationRequest& request, const Interrupt& interrupt) {
    request.validate();
    if (interrupt) {
        interrupt();
    }

    // RegArimaUtility.olsResiduals(RegArimaModel.of(model)): the response is
    // BackFilter(spec.doDifferencing())-filtered and the design is built on
    // the differenced domain - a fresh ones column (not a differenced
    // constant) under mean_correction, then the BackFilter'd user columns in
    // request order. The cursor is the OLS residual of that design.
    const std::vector<double> filter =
        detail::pinned_differencing_polynomial(request.period, request.d,
            request.seasonal_d);
    const auto backfilter = [&filter](const std::vector<double>& x) {
        // The qualified seam returns the filtered vector already shortened by
        // the filter degree (an empty result marks a filter not shorter than
        // the series, where the pinned machine fails).
        return filter.empty() ? x : detail::pinned_backfilter_apply(filter, x);
    };
    const std::vector<double> y = backfilter(request.series);
    if (y.empty()) {
        throw std::runtime_error(
            "identify_initial_arma: not enough observations after differencing "
            "(the pinned EcoException)");
    }
    std::vector<std::vector<double>> columns;
    if (request.mean_correction) {
        columns.push_back(std::vector<double>(y.size(), 1.0));
    }
    for (const auto& regressor : request.regressors) {
        columns.push_back(backfilter(regressor));
    }
    if (interrupt) {
        interrupt();
    }

    std::vector<double> cursor;
    if (columns.empty()) { // LinearModel without regressors: the response
        cursor = y;
    } else {
        double norm2 = 0.0;
        for (const double v : y) {
            norm2 += v * v;
        }
        if (norm2 < 1.1102230246251565e-16) {
            throw std::runtime_error(
                "identify_initial_arma: the OLS cursor response is degenerate "
                "(norm2 below the pinned epsilon, the pinned null-results NPE)");
        }
        std::vector<double> beta;
        if (columns.size() == 1) {
            beta = {detail::pinned_ols_1col(columns[0], y)};
        } else {
            // Pinned QRLeastSquaresSolver.robustLeastSquares: HouseholderWith-
            // Pivoting plus leastSquares at rcond 1.0E-9. Kalos's DGELSY seam
            // (a documented deviation, same qualification class as the 21-A
            // Yule-Walker seam): rank-deficient designs take the minimum-norm
            // solution where the JVM solve truncates the pivoted system.
            DenseMatrix design(y.size(), columns.size());
            for (std::size_t j = 0; j < columns.size(); ++j) {
                for (std::size_t i = 0; i < y.size(); ++i) {
                    design(i, j) = columns[j][i];
                }
            }
            DenseMatrix response(y.size(), 1);
            for (std::size_t i = 0; i < y.size(); ++i) {
                response(i, 0) = y[i];
            }
            const LeastSquaresResult solution = qr_least_squares(design, response, 1.0e-9);
            beta.reserve(columns.size());
            for (std::size_t j = 0; j < columns.size(); ++j) {
                beta.push_back(solution.coefficients(j, 0));
            }
        }
        // LinearModel.calcResiduals: r = y; per column ASCENDING over the
        // whole vector r += (-b_j) * col_j (the DataBlock.addAY per-column
        // accumulation, not per-row dot products).
        cursor = y;
        for (std::size_t j = 0; j < columns.size(); ++j) {
            const double bj = beta[j];
            for (std::size_t i = 0; i < cursor.size(); ++i) {
                cursor[i] += (-bj) * columns[j][i];
            }
        }
    }

    ArmaSelectRequest select;
    select.cursor = cursor;
    select.period = request.period;
    select.d = request.d;
    select.seasonal_d = request.seasonal_d;
    select.models_count = request.models_count;
    select.max_p = request.max_p;
    select.max_q = request.max_q;
    select.max_bp = request.max_bp;
    select.max_bq = request.max_bq;
    select.balanced = request.balanced;
    select.mixed = request.mixed;
    select.accept_white_noise = request.accept_white_noise;
    select.eps = request.eps;
    // select.validate() holds no new information past the validated request
    // (the cursor is built from finite inputs by the differencing/OLS sweep);
    // the engine runs directly, keeping the pinned-abort messages under the
    // identify_initial_arma prefix.
    ArmaSelectResult result =
        select_engine(select, interrupt, "identify_initial_arma");

    InitialArmaIdentification out;
    out.p = result.p;
    out.q = result.q;
    out.bp = result.bp;
    out.bq = result.bq;
    out.residuals = std::move(cursor);
    out.evaluations = std::move(result.evaluations);
    out.preferred_models = std::move(result.preferred_models);
    out.refinement = std::move(result.refinement);
    return out;
}

} // namespace kalos
