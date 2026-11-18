#include "kalos/arma_refinement.hpp"

#include "kalos/detail/ssq_deviances.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace kalos {

void ArmaRefinementRequest::validate() const {
    if (series.empty()) {
        throw std::invalid_argument(
            "ArmaRefinementRequest: the series must not be empty");
    }
    if (!std::all_of(series.begin(), series.end(),
            [](double v) { return std::isfinite(v); })) {
        throw std::invalid_argument(
            "ArmaRefinementRequest: the series must be finite");
    }
    for (const auto& column : regressors) {
        if (column.size() != series.size()) {
            throw std::invalid_argument(
                "ArmaRefinementRequest: every regressor must have series.size() rows");
        }
        if (!std::all_of(column.begin(), column.end(),
                [](double v) { return std::isfinite(v); })) {
            throw std::invalid_argument(
                "ArmaRefinementRequest: every regressor must be finite");
        }
    }
    const std::size_t expected =
        regressors.size() + (mean ? 1u : 0u); // cl.coefficients() shape
    if (regression_coefficients.size() != expected) {
        throw std::invalid_argument(
            "ArmaRefinementRequest: regression_coefficients must hold the mean "
            "coefficient (when the mean flag is on) plus one value per "
            "regressor");
    }
    if (!std::all_of(regression_coefficients.begin(), regression_coefficients.end(),
            [](double v) { return std::isfinite(v); })) {
        throw std::invalid_argument(
            "ArmaRefinementRequest: the regression coefficients must be finite");
    }
    if (p < 0 || q < 0 || bp < 0 || bq < 0) {
        throw std::invalid_argument(
            "ArmaRefinementRequest: the current orders must not be negative");
    }
    if (d < 0 || seasonal_d < 0) {
        throw std::invalid_argument(
            "ArmaRefinementRequest: the current differencing orders must not "
            "be negative");
    }
    // Compose the two inner validators instead of duplicating their rules:
    // the DifferencingDecisionRequest and ArmaSelectRequest scalar checks
    // with placeholder response data (the real linearized response is
    // re-validated by decide_differencing itself inside the guarded region).
    DifferencingDecisionRequest dif;
    dif.series.assign(2, 1.0);
    dif.period = period;
    dif.max_d = max_d;
    dif.max_seasonal_d = max_seasonal_d;
    dif.ub1 = ub1;
    dif.ub2 = ub2;
    dif.cancel = cancel;
    dif.precision = precision;
    dif.validate();
    ArmaSelectRequest arma;
    arma.cursor.assign(2, 1.0);
    arma.period = period;
    arma.d = 0;
    arma.seasonal_d = 0;
    arma.models_count = models_count;
    arma.max_p = max_p;
    arma.max_q = max_q;
    arma.max_bp = max_bp;
    arma.max_bq = max_bq;
    arma.eps = eps;
    arma.validate();
}

namespace {

// The AutoModellingModule re-differencing polynomial: UnitRoots.D(d) times
// UnitRoots.D(period, bd). The one-argument UnitRoots.D(n) is the SINGLE tap
// 1 - B^n (bytecode: fromDegree(n), coeff[0] = 1, coeff[n] = -1), NOT the
// d-fold (1 - B)^d the SarimaUtility.differencingFilter path builds via
// D(1, d); at d = 2 (reachable: the decoder maxd is 2) the filters genuinely
// differ, (1 - B^2) vs (1 - B)^2. D(period, bd) is the bd-fold power
// (1 - B^period)^bd; setRoots there only caches the roots, so the pinned
// coefficients are the exact integer convolutions below.
std::vector<double> automod_unitroot_polynomial(int period, int d, int bd) {
    auto convolve = [](const std::vector<double>& a, const std::vector<double>& b) {
        std::vector<double> c(a.size() + b.size() - 1, 0.0);
        for (std::size_t i = 0; i < a.size(); ++i) {
            for (std::size_t j = 0; j < b.size(); ++j) {
                c[i + j] += a[i] * b[j];
            }
        }
        return c; // integer taps: exact; no trailing zeros can appear here
    };
    std::vector<double> seasonal{1.0};
    for (int i = 0; i < bd; ++i) {
        std::vector<double> tap(static_cast<std::size_t>(period) + 1, 0.0);
        tap[0] = 1.0;
        tap[static_cast<std::size_t>(period)] = -1.0;
        seasonal = convolve(seasonal, tap);
    }
    std::vector<double> regular{1.0};
    if (d > 0) {
        regular.assign(static_cast<std::size_t>(d) + 1, 0.0);
        regular[0] = 1.0;
        regular[static_cast<std::size_t>(d)] = -1.0;
    }
    return convolve(regular, seasonal); // Polynomial.times order: D(d) first
}

} // namespace

ArmaRefinement refine_arma_model(
    const ArmaRefinementRequest& request, const Interrupt& interrupt) {
    request.validate();
    if (interrupt) {
        interrupt();
    }

    ArmaRefinement out;

    // The pinned module catches RuntimeException around the whole guarded
    // region (linearizedData .. SarimaOrders.of), so a runtime_error from
    // any inner qualified unit becomes the airline result there. The
    // Kalos-side interrupt contract ("interrupt exceptions always propagate")
    // is kept by arming on the first callback firing: any exception that
    // arrives after the caller's callback has thrown is rethrown instead of
    // being converted into the pinned fallback.
    bool interrupt_fired = false;
    const Interrupt guarded = [&interrupt, &interrupt_fired]() {
        if (!interrupt) {
            return;
        }
        try {
            interrupt();
        } catch (...) {
            interrupt_fired = true;
            throw;
        }
    };

    const int period = request.period;
    int d = 0;
    int bd = 0;
    bool mean = false;
    std::optional<ArmaSelectResult> selected;
    const char* stage = "linearized";

    try {
        // RegArimaUtility.linearizedData(model, cl): z = y; missing
        // corrections (no missing in the qualified domain - preparation is
        // upstream, like every earlier item-21 seam); then, when the cl
        // carries coefficients, per user column ASCENDING
        // z[i] = z[i] - coefficient * column[i] (the DataBlock.apply
        // element loop with the captured-coefficient lambda). The mean
        // coefficient is SKIPPED, never removed: linearized data keeps the
        // mean in the series.
        std::vector<double> lin = request.series;
        if (!request.regression_coefficients.empty()) {
            std::size_t position = request.mean ? 1u : 0u;
            for (const auto& column : request.regressors) {
                const double coefficient = request.regression_coefficients[position++];
                for (std::size_t i = 0; i < lin.size(); ++i) {
                    lin[i] = lin[i] - (coefficient * column[i]);
                }
            }
        }
        out.linearized = lin;

        // DifferencingModule.process(lin, period) - the 21-D unit, verbatim.
        // The pinned "== false" remark branch is dead code in 3.9.0 (calc()
        // always returns true; every failure throws, landing in the outer
        // catch), so the transcription has no (1, 1, false) path: a throw
        // here is the airline fallback, exactly as the JVM reaches it.
        stage = "differencing";
        DifferencingDecisionRequest dif;
        dif.series = lin;
        dif.period = period;
        dif.max_d = request.max_d;
        dif.max_seasonal_d = request.max_seasonal_d;
        dif.ub1 = request.ub1;
        dif.ub2 = request.ub2;
        dif.cancel = request.cancel;
        dif.precision = request.precision;
        out.differencing = decide_differencing(dif, guarded);
        d = out.differencing->d;
        bd = out.differencing->seasonal_d;
        mean = out.differencing->mean_correction;

        // if (d != 0 || bd != 0) z = BackFilter(D(d) * D(period, bd)).apply(lin)
        // else z = lin. The apply is the qualified pinned BackFilter (exact
        // in 21-E); a response not longer than the filter degree is the
        // pinned NegativeArraySizeException/failed-candidate path, which the
        // outer catch turns into the airline spec.
        stage = "filter";
        if (d == 0 && bd == 0) {
            out.re_differenced = lin;
        } else {
            out.re_differenced = detail::pinned_backfilter_apply(
                automod_unitroot_polynomial(period, d, bd), lin);
            if (out.re_differenced.empty()) {
                throw std::runtime_error(
                    "refine_arma_model: the linearized response is not longer "
                    "than the re-differencing filter (the pinned negative-size "
                    "array, the airline path)");
            }
        }

        // ArmaModule.select(z, period, d, bd) + getPreferedModels(): the 21-E
        // machine on the filtered sequence. Pinned candidate/panel aborts are
        // runtime_errors landing in the same airline catch (the module
        // survives where ArmaModule.process would report Failed).
        stage = "select";
        ArmaSelectRequest arma;
        arma.cursor = out.re_differenced;
        arma.period = period;
        arma.d = d;
        arma.seasonal_d = bd;
        arma.models_count = request.models_count;
        arma.max_p = request.max_p;
        arma.max_q = request.max_q;
        arma.max_bp = request.max_bp;
        arma.max_bq = request.max_bq;
        arma.balanced = request.balanced;
        arma.mixed = request.mixed;
        arma.accept_white_noise = request.accept_white_noise;
        arma.eps = request.eps;
        selected = select_arma_orders(arma, guarded);
    } catch (const std::runtime_error&) {
        if (interrupt_fired) {
            throw;
        }
        // SarimaOrders.airline(period); mean = false; remark
        // "default model selected". period > 1 always holds in the
        // qualified domain {4, 12}; the period == 1 branch (bd = bq = 0) is
        // transcribed for fidelity.
        out.fallback = true;
        out.fallback_stage = std::string(stage);
        out.p = 0;
        out.q = 1;
        out.bp = 0;
        d = 1;
        mean = false;
        if (period > 1) {
            out.bq = 1;
            bd = 1;
        } else {
            out.bq = 0;
            bd = 0;
        }
    }

    if (!out.fallback) {
        out.p = selected->p;
        out.q = selected->q;
        out.bp = selected->bp;
        out.bq = selected->bq;
        out.evaluations = std::move(selected->evaluations);
        out.preferred_models = std::move(selected->preferred_models);
        out.refinement = std::move(selected->refinement);
    }
    out.d = d;
    out.seasonal_d = bd;
    out.mean = mean;

    // The pinned comparison block: SarimaOrders.equals covers period, p, d,
    // q, bp, bd, bq (the result period is desc.getAnnualFrequency(), so the
    // six order fields decide), then the separate mean-flag test.
    out.specification_changed = out.p != request.p || out.q != request.q
        || out.bp != request.bp || out.bq != request.bq || out.d != request.d
        || out.seasonal_d != request.seasonal_d;
    out.mean_changed = out.mean != request.mean;
    out.changed = out.specification_changed || out.mean_changed;
    return out;
}

} // namespace kalos
