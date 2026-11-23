// Implementation of the TramoKernel round-driver seam (item 23-F.3):
// the ami() ten-iterate loop with the guards, the change accumulator,
// the estimateModel() retry ladder, the round-0 reference snapshot, the
// refStats gate call, the seasonality-controller pass3 reset and the
// controller sweep, over the ordered seam facts of the request. The heavy
// stages ride as facts; the kernel machines are the landed 23-E entries.
// Decoded tables: tmp_opencode/kalos-m5-pinning.md SECTIONS 25.1, 26.4,
// 30.8; corpus: validation/evidence/native/tramo-kernel-rounds
// (Scratch23r manual driver, 35/35 identity).

#include <kalos/tramo_kernel_rounds.hpp>
#include <kalos/tramo_transform_decision.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace kalos {
namespace {

void require(bool condition, const std::string& what) {
    if (!condition)
        throw std::invalid_argument("tramo kernel rounds: " + what);
}

void require_finite(double v, const std::string& what) {
    if (!std::isfinite(v))
        throw std::invalid_argument("tramo kernel rounds: a non-finite " +
                                    what + " is outside the envelope");
}

int parameter_count(const TramoOrders& o) {
    return o.p + o.q + o.bp + o.bq;
}

int difference_order(const TramoOrders& o, int period) {
    return o.d + period * o.bd;
}

bool orders_equal(const TramoOrders& a, const TramoOrders& b) {
    return a.p == b.p && a.d == b.d && a.q == b.q && a.bp == b.bp &&
           a.bd == b.bd && a.bq == b.bq;
}

// desc.arima().phi(1) / bphi(1) over the carried parameter vector
// [phi(p)..., bphi(bp)..., theta(q)..., btheta(bq)...] (the p2dbg pars
// line). The shape guards read the coefficient only when the matching
// order is 1; no kr corpus case reaches a guard (all 34 gate echoes
// keep), so this indexing is decode-qualified, not corpus-witnessed.
double arima_phi(int before, const std::vector<double>& pars) {
    if (before < 0 || static_cast<std::size_t>(before) >= pars.size())
        return 0.0;
    return pars[static_cast<std::size_t>(before)];
}

}  // namespace

TramoKernelRoundsRun tramo_run_kernel_rounds(
    const TramoKernelRoundsRequest& request) {
    const TramoKernelRoundsSpec& spec = request.spec;
    if (spec.annual_frequency != 12)
        throw std::invalid_argument(
            "the tramo kernel rounds envelope is monthly frequency only (D11)");
    if (spec.observations <= 0)
        throw std::invalid_argument(
            "the tramo kernel rounds request needs a positive observation count");
    if (spec.controller_count < 0)
        throw std::invalid_argument(
            "the tramo kernel rounds controller count must be non-negative");
    if (request.seams.empty())
        throw std::invalid_argument(
            "the tramo kernel rounds request carries no seam facts");
    require_finite(spec.pcr, "pcr");
    require_finite(spec.reduce_va, "reduceVa");
    require_finite(spec.tsig, "tsig");
    require_finite(spec.critical_value, "critical value");

    TramoKernelRoundsRun run;
    run.ami_check_mu =
        tramo_ami_check_mu(spec.using_auto_model, spec.mean_is_test);

    // --- initProcessing(regarima.getActualObservationsCount()) ---
    double curva = spec.critical_value;
    // TramoUtility.calcCv runs ONLY when outlier detection is on and
    // dcmpl(curva, 0.0) == 0, i.e. the user left the critical value at
    // 0.0 (SECTION 25 correction; a NaN curva skips, dcmpl(NaN, 0) != 0).
    if (spec.outliers_detection && curva == 0.0)
        curva = tramo_critical_value(static_cast<std::size_t>(spec.observations));
    double pcr = spec.pcr;
    int pass = 0, round = 0;
    bool pass3 = false;
    bool need_out = spec.outliers_detection;
    bool need_auto = false;  // initProcessing sets the field false
    bool ref_auto = false, ref_stats = false;
    // refAuto's description state + refStats' row (the gate operands).
    TramoOrders ref_orders = request.prestage.orders;
    bool ref_mean = request.prestage.mean;
    std::vector<std::string> ref_outlier_names = request.prestage.outliers;
    double ref_pvalue = 0.0, ref_se = 0.0;

    // The working description state (pre-stage hand-over).
    TramoOrders orders = request.prestage.orders;
    bool mean = request.prestage.mean;
    std::vector<std::string> outliers = request.prestage.outliers;
    bool seasonal = request.prestage.seasonal;
    int original_seasonality = request.prestage.original_seasonality_test;
    // needEstimation(): TramoModelBuilder.build() leaves the model
    // unestimated, so the first iterate() estimates (kr witness).
    bool estimation_absent = true;

    run.initial_curva = curva;
    run.initial_pcr = pcr;

    std::vector<TramoKernelLogEvent>& events = run.events;
    events.push_back({TramoKernelLogKind::push, 0, "tramo"});
    events.push_back(
        {TramoKernelLogKind::push, 1, "automatic model identification"});
    for (const TramoKernelLogEvent& e : request.prestage.events)
        events.push_back(e);

    const std::vector<TramoKernelSeam>& seams = request.seams;
    std::size_t next = 0;

    auto expect = [&](TramoKernelSeamKind kind,
                      const std::string& where) -> const TramoKernelSeam& {
        require(next < seams.size(),
                "the facts end where " + where + " expects another seam");
        require(seams[next].kind == kind,
                where + " disagrees with the seam facts at fact " +
                    std::to_string(next));
        return seams[next++];
    };
    auto guard_facts = [](const TramoKernelSeam& f) {
        require_finite(f.critical_value, "curva fact");
        require_finite(f.intermediate_precision, "precision fact");
        require_finite(f.tsig, "tsig fact");
        require_finite(f.returned_pcr, "pcr fact");
        require_finite(f.stats.ljung_box_pvalue, "stats lbp");
        require_finite(f.stats.se, "stats se");
        require_finite(f.stats.ljung_box, "stats lb");
        require_finite(f.stats.bic, "stats bic");
        for (const double v : f.parameters)
            require_finite(v, "arima parameter");
        for (const double v : f.tstats)
            require_finite(v, "t-statistic");
    };

    // The needXXX desc predicates (bytecode 1251/1311/1300).
    const auto round2_guard = [&]() { return round == 2 && outliers.empty(); };
    const auto need_differencing = [&]() {
        if (!need_auto || round2_guard())
            return false;
        return orders.d < 2 || orders.bd >= 1;
    };
    const auto need_auto_desc = [&]() { return need_auto && !round2_guard(); };
    const auto need_outliers_desc = [&]() {
        return spec.outliers_detection && need_out;
    };

    auto snapshot_now = [&]() {
        run.snapshots.push_back({pass, round, pass3, need_out, need_auto,
                                 curva, pcr, ref_auto, ref_stats, seasonal,
                                 original_seasonality});
    };
    auto push_decision = [&](TramoKernelSeamKind kind) {
        run.decisions.push_back(TramoKernelSeamDecision{});
        run.decisions.back().kind = kind;
    };

    bool converged = false;
    for (int iteration = 0; iteration < 10 && !converged; ++iteration) {
        // ---------------- iterate(m) ----------------
        events.push_back({TramoKernelLogKind::step, 2,
                          "round " + std::to_string(pass)});
        snapshot_now();
        run.iterations = iteration + 1;

        // 1. The estimate site (m.needEstimation()).
        {
            const TramoKernelSeam& fact =
                expect(TramoKernelSeamKind::estimation, "the estimate site");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::estimation);
            run.decisions.back().flag = estimation_absent ? 1 : 0;
            require(fact.estimation_taken == (estimation_absent),
                    "the estimate site disagrees with the estimation state "
                    "at fact " + std::to_string(next - 1));
            if (fact.estimation_taken) {
                require(fact.intermediate_precision ==
                            spec.intermediate_precision,
                        "the estimate site consumed a precision the spec "
                        "does not carry");
                estimation_absent = false;
            }
        }

        bool changed = false;
        const TramoOrders orders0 = orders;
        const bool mean0 = mean;

        // 2. Differencing (changed := res == Changed - an assignment).
        if (need_differencing()) {
            const TramoKernelSeam& fact = expect(
                TramoKernelSeamKind::differencing, "the differencing block");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::differencing);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            changed = fact.result == TramoKernelStageResult::changed;
            orders = fact.orders;
            mean = fact.mean;
            outliers = fact.outliers;
            // pass-1 difference-order change: the ami outliers go.
            if (pass == 1 &&
                difference_order(fact.orders, spec.annual_frequency) != difference_order(orders0, spec.annual_frequency)) {
                const TramoKernelSeam& drop = expect(
                    TramoKernelSeamKind::differencing_drop,
                    "the pass-1 differencing drop");
                for (const TramoKernelLogEvent& e : drop.events)
                    events.push_back(e);
                push_decision(TramoKernelSeamKind::differencing_drop);
                outliers.clear();
                estimation_absent = true;  // m.clearEstimation()
            }
        }

        // 3. Arma (changed := !(spec.equals(orders0) && mean == mean0)).
        if (need_auto_desc()) {
            const TramoKernelSeam& fact =
                expect(TramoKernelSeamKind::arma, "the arma block");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::arma);
            // The pass >= 3 zero-parameter rewrite runs BEFORE the probe
            // dumped the arma fact, so the fact state is post-rewrite.
            if (parameter_count(fact.orders) == 0 && pass >= 3) {
                const TramoKernelSeam& q1 = expect(
                    TramoKernelSeamKind::arma_q1, "the pass-3 q rewrite");
                guard_facts(q1);
                push_decision(TramoKernelSeamKind::arma_q1);
                TramoOrders rewritten = fact.orders;
                rewritten.q = 1;
                require(orders_equal(rewritten, q1.orders),
                        "the pass-3 q rewrite disagrees with its fact");
            }
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            changed = !(orders_equal(fact.orders, orders0) && fact.mean == mean0);
            orders = fact.orders;
            mean = fact.mean;
        }

        // 4. Outliers (changed |= res == Changed).
        if (need_outliers_desc()) {
            if (!outliers.empty()) {
                // desc.removeVariable(isOutlier(v, true)) cleared the
                // estimation: every kernel round feeds the scan a clean
                // design (26.5).
                estimation_absent = true;
                outliers.clear();
            }
            const TramoKernelSeam& fact =
                expect(TramoKernelSeamKind::scan, "the outliers block");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::scan);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            require(fact.critical_value == curva,
                    "the outliers block ran with a critical value the "
                    "kernel does not carry");
            changed = changed || fact.result == TramoKernelStageResult::changed;
            orders = fact.orders;
            mean = fact.mean;
            outliers = fact.outliers;
        }

        // 5. estimateModel(m): the ladder over the estmodel/retry facts,
        // cross-checked against the landed 23-E retry machine.
        bool ladder_ok = false;
        bool reduce_fired = false;
        double reduce_curva = 0.0;
        const double curva_at_ladder = curva;
        std::vector<bool> estimation_results;
        std::vector<bool> test_results;
        int retry = 0;
        for (;;) {
            const TramoKernelSeam& fact = expect(
                TramoKernelSeamKind::estimation_model, "estimateModel");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::estimation_model);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            require(fact.pass_at_call == pass,
                    "the FinalEstimator was built at a pass the kernel "
                    "does not carry");
            require(fact.returned == 0 || fact.returned == 1,
                    "the estimateModel fact carries an impossible return");
            estimation_results.push_back(fact.returned == 1);
            orders = fact.orders;
            outliers = fact.outliers;
            // FinalEstimator.estimate ASSIGNS its estimation before any
            // gate (30.1); the corpus never reaches the failed-attempt
            // path, where the common-roots rewrite leaves no estimation
            // (documented, replay-qualified by the desync throw).
            estimation_absent = fact.returned != 1;
            if (fact.returned == 0) {
                if (pass == 1 && need_outliers_desc() && need_auto_desc()) {
                    const TramoKernelSeam& rv = expect(
                        TramoKernelSeamKind::estimation_reduce,
                        "the estimateModel failure-arm reduceVa");
                    guard_facts(rv);
                    push_decision(TramoKernelSeamKind::estimation_reduce);
                    const TramoReduceVa reduce =
                        tramo_reduce_va(curva, spec.reduce_va);
                    run.decisions.back().value = reduce.curva;
                    require(rv.critical_value == reduce.curva,
                            "the estimateModel reduceVa disagrees with the "
                            "landed machine");
                    for (const TramoKernelLogEvent& e : rv.events)
                        events.push_back(e);
                    curva = reduce.curva;
                    reduce_fired = true;
                    reduce_curva = reduce.curva;
                }
                break;  // estimateModel returned false
            }
            if (pass3 || pass != 0) {
                ladder_ok = true;
                break;
            }
            if (retry >= 5) {
                ladder_ok = true;  // the bytecode returns true at i >= 5
                break;
            }
            ++retry;
            const TramoKernelSeam& t = expect(
                TramoKernelSeamKind::regression_test,
                "the estimateModel mean test");
            guard_facts(t);
            push_decision(TramoKernelSeamKind::regression_test);
            for (const TramoKernelLogEvent& e : t.events)
                events.push_back(e);
            require(t.retry_site && t.retry_index == retry,
                    "the mean-test site disagrees with the ladder");
            require(t.tsig == 1.96,
                    "the ladder mean test ran with a threshold other than "
                    "the pinned 1.96");
            mean = t.mean;
            const bool ok = t.result == TramoKernelStageResult::unchanged;
            test_results.push_back(ok);
            if (ok) {
                ladder_ok = true;
                break;
            }
            // failure: back to L0 (the for-counter increments at the
            // next loop head, so five failures still complete).
        }
        {
            TramoEstimateModelRetryRequest retry_request;
            retry_request.pass = pass;
            retry_request.pass3 = pass3;
            retry_request.need_outliers_desc = need_outliers_desc();
            retry_request.need_auto_desc = need_auto_desc();
            retry_request.curva = curva_at_ladder;
            retry_request.reduce_va = spec.reduce_va;
            retry_request.estimation_results = estimation_results;
            retry_request.test_results = test_results;
            const TramoEstimateModelRetry machine =
                tramo_estimate_model_retry(retry_request);
            require(machine.completed == ladder_ok,
                    "the estimateModel ladder disagrees with the landed "
                    "retry machine");
            require(machine.estimations_consumed == estimation_results.size(),
                    "the estimateModel ladder consumed a different number "
                    "of estimations than the landed retry machine");
            require(machine.tests_consumed == test_results.size(),
                    "the estimateModel ladder consumed a different number "
                    "of mean tests than the landed retry machine");
            require(machine.reduce.changed == reduce_fired,
                    "the estimateModel reduceVa arm disagrees with the "
                    "landed retry machine");
            if (reduce_fired)
                require(machine.reduce.curva == reduce_curva,
                        "the estimateModel reduceVa value disagrees with "
                        "the landed machine");
        }
        if (!ladder_ok) {
            need_out = spec.outliers_detection;
            need_auto = spec.using_auto_model;
            ++round;
            ++pass;
            const TramoKernelSeam& end = expect(
                TramoKernelSeamKind::iterate_end, "the estimateModel end");
            push_decision(TramoKernelSeamKind::iterate_end);
            run.decisions.back().flag =
                static_cast<int>(TramoKernelIterateEnd::estimation_fail);
            require(end.end == TramoKernelIterateEnd::estimation_fail,
                    "the iterate end should be estmodelfail");
            continue;  // iterate returned false
        }

        // 6. The round-0 reference block (the acceptAirline rule).
        if (round == 0 &&
            (spec.using_auto_model || spec.outliers_detection)) {
            need_out = spec.outliers_detection;
            need_auto = spec.using_auto_model;
            ++round;
            ++pass;
            const TramoKernelSeam& fact =
                expect(TramoKernelSeamKind::round0, "the round-0 block");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::round0);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            ref_auto = true;
            ref_stats = true;
            ref_orders = orders;
            ref_mean = mean;
            ref_outlier_names = outliers;
            ref_pvalue = fact.stats.ljung_box_pvalue;
            ref_se = fact.stats.se;
            run.decisions.back().flag =
                tramo_round0_accepts(fact.stats.ljung_box_pvalue, spec.pcr,
                                     spec.accept_airline)
                    ? 1
                    : 0;
            require(fact.accept == (run.decisions.back().flag == 1),
                    "the round-0 accept disagrees with the landed rule");
            if (run.decisions.back().flag == 1)
                converged = true;  // iterate returned true (no end fact)
            continue;              // iterate returned (accept)
        }

        // 7. The refStats gate (pass2(!changed, m) - the negated
        // accumulator is the operand; the continue arm reads it back).
        if (pass <= 3 && !pass3 && spec.using_auto_model) {
            const TramoKernelSeam& fact =
                expect(TramoKernelSeamKind::pass2_gate, "the refStats gate");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::pass2_gate);
            TramoKernelSeamDecision& d = run.decisions.back();
            TramoRefStatsGateRequest gate;
            gate.pass = pass;
            gate.round = round;
            gate.pcr = pcr;
            gate.changed = !changed;
            gate.seasonal = seasonal;
            gate.curva = curva;
            gate.reduce_va = spec.reduce_va;
            gate.outliers_detection = spec.outliers_detection;
            // The rules read desc.arima().orders() and the arima
            // coefficients (p2dbg), not the specification.
            gate.orders = fact.arima_orders;
            gate.mean = mean;
            gate.phi1 = fact.arima_orders.p >= 1
                            ? arima_phi(0, fact.parameters)
                            : 0.0;
            gate.bphi1 = fact.arima_orders.bp >= 1
                             ? arima_phi(fact.arima_orders.p,
                                         fact.parameters)
                             : 0.0;
            gate.ljung_box_complement = 1.0 - fact.stats.ljung_box_pvalue;
            gate.standard_error = fact.stats.se;
            gate.current_outliers = static_cast<int>(outliers.size());
            require(gate.current_outliers == fact.outlier_regressors,
                    "the current outlier tally disagrees with the "
                    "countRegressors echo (dim-1 outlier variables only)");
            const int ref_count = static_cast<int>(ref_outlier_names.size());
            if (ref_auto) {
                require(ref_count == fact.reference_regressors,
                        "the reference outlier tally disagrees with the "
                        "countRegressors echo");
                gate.reference_outliers = ref_count;
                gate.reference_orders = ref_orders;
                gate.reference_mean = ref_mean;
                gate.reference_ljung_box_complement = 1.0 - ref_pvalue;
                gate.reference_standard_error = ref_se;
            } else {
                gate.reference_outliers = -1;
            }
            const TramoRefStatsGateResult g = tramo_ref_stats_gate(gate);

            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            d.effective_ljung_box_complement = g.effective_ljung_box_complement;
            d.current_standard_error = gate.standard_error;
            d.reference_ljung_box_complement =
                gate.reference_ljung_box_complement;
            d.reference_standard_error = gate.reference_standard_error;
            d.current_regressors = gate.current_outliers;
            d.reference_regressors = gate.reference_outliers;
            d.reverted = g.reverted;
            d.skipped_rules = gate.reference_outliers >= 0 &&
                              gate.reference_outliers > gate.current_outliers;
            d.accepted = g.accepted;
            d.remark = g.remark;
            d.last_solution = g.last_solution;
            d.pass_after = g.pass;
            d.round_after = g.round;
            d.pcr_after = g.pcr;
            d.curva_after = g.curva;
            d.need_auto_modelling_after = g.need_auto_modelling;
            d.need_outliers_after = g.need_outliers;
            d.orders_after = g.orders_after;
            d.mean_after = g.mean_after;
            d.outliers_after = g.outliers_after;

            // JVM-echoed verdict cross-checks (the rule id rides as
            // data; the landed 23-E corpus owns the A..F attribution).
            require(fact.reverted == g.reverted,
                    "the gate revert disagrees with the landed machine");
            require(d.skipped_rules == (fact.rule == "skip-nrnb"),
                    "the skip-nrnb echo disagrees with nr > nb");
            if (fact.rule != "first") {
                // The JVM gate-line operands, bit-for-bit: the lbp/se
                // pair the return test saw (post-restore on a revert).
                require(fact.gate_ljung_box_complement ==
                            g.effective_ljung_box_complement,
                        "the gate lbp echo disagrees with the landed "
                        "machine");
                require(fact.gate_standard_error ==
                            (g.reverted ? gate.reference_standard_error
                                        : gate.standard_error),
                        "the gate se echo disagrees with the landed "
                        "machine");
            }

            // The gate's own log events (in bytecode order: the remark
            // after reduceVa, then the last-chance warning).
            if (g.remark)
                events.push_back({TramoKernelLogKind::remark, 2,
                                  "reduction of the critical value for "
                                  "outliers detection"});

            // State transitions.
            pass = g.pass;
            round = g.round;
            pcr = g.pcr;
            curva = g.curva;
            if (!g.accepted) {
                need_out = g.need_outliers;
                need_auto = g.need_auto_modelling;
            }
            if (g.reverted) {
                // restore(m): m.set(copyOf(refAuto.desc), refAuto.est).
                // The corpus never reaches this arm (all 34 echoes keep);
                // the restore is qualified by the landed gate corpus.
                orders = ref_orders;
                mean = ref_mean;
                outliers = ref_outlier_names;
            } else {
                if (g.last_solution) {
                    events.push_back({TramoKernelLogKind::warning, 2,
                                      "last chance model"});
                    // lastSolution rewrites the SPECIFICATION (not the
                    // arima operand), clears the outliers and pins round
                    // to 1; the kr_q_air pass-3 re-estimation witness
                    // invalidates the estimation.
                    orders = tramo_last_solution_orders(orders, seasonal);
                    outliers.clear();
                    estimation_absent = true;
                }
                // refAuto := copyOf(m); refStats := stats.
                ref_orders = orders;
                ref_mean = mean;
                ref_outlier_names = outliers;
                ref_pvalue = fact.stats.ljung_box_pvalue;
                ref_se = fact.stats.se;
                ref_auto = true;
                ref_stats = true;
            }

            const TramoKernelSeam& ret = expect(
                TramoKernelSeamKind::pass2_return, "the gate return");
            guard_facts(ret);
            push_decision(TramoKernelSeamKind::pass2_return);
            run.decisions.back().flag = g.accepted ? 1 : 0;
            run.decisions.back().value = g.pcr;
            require(ret.returned == (g.accepted ? 1 : 0),
                    "the gate return disagrees with the landed machine");
            require(ret.returned_pcr == g.pcr,
                    "the evolved pcr disagrees with the landed ladder");
            if (!g.accepted)
                continue;  // iterate returned false (no end fact)
        }

        // 8. The raw-tsig mean test.
        {
            const TramoKernelSeam& fact = expect(
                TramoKernelSeamKind::regression_test, "the iterate mean test");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::regression_test);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            require(!fact.retry_site,
                    "an iterate mean test must not be a retry site");
            require(fact.tsig == spec.tsig,
                    "the iterate mean test ran with a threshold the spec "
                    "does not carry");
            mean = fact.mean;
            if (fact.result != TramoKernelStageResult::unchanged) {
                pass = 4;
                need_auto = false;
                need_out = false;
                const TramoKernelSeam& end = expect(
                    TramoKernelSeamKind::iterate_end, "the mean-test end");
                push_decision(TramoKernelSeamKind::iterate_end);
                run.decisions.back().flag =
                    static_cast<int>(TramoKernelIterateEnd::regression_fail);
                require(end.end == TramoKernelIterateEnd::regression_fail,
                        "the iterate end should be regfail");
                continue;  // iterate returned false
            }
        }

        // 9. The seasonality controller (the pass3 wiring).
        if (spec.using_auto_model) {
            const TramoKernelSeam& fact = expect(
                TramoKernelSeamKind::seasonality_controller,
                "the seasonality controller");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::seasonality_controller);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            orders = fact.orders;
            mean = fact.mean;
            outliers = fact.outliers;
            seasonal = fact.seasonal_after;
            original_seasonality = fact.original_seasonality_after;
            if (fact.result == TramoKernelStageResult::changed && !pass3) {
                pass3 = true;
                pass = 1;
                need_auto = true;
                need_out = spec.outliers_detection;
                const TramoKernelSeam& end = expect(
                    TramoKernelSeamKind::iterate_end, "the pass3 reset");
                push_decision(TramoKernelSeamKind::iterate_end);
                run.decisions.back().flag =
                    static_cast<int>(TramoKernelIterateEnd::pass3_reset);
                require(end.end == TramoKernelIterateEnd::pass3_reset,
                        "the iterate end should be pass3reset");
                continue;  // iterate returned false
            }
        }

        // 10. The controller sweep (the kernel ignores the results).
        for (int c = 0; c < spec.controller_count; ++c) {
            const TramoKernelSeam& fact = expect(
                TramoKernelSeamKind::model_controller, "the controller sweep");
            guard_facts(fact);
            push_decision(TramoKernelSeamKind::model_controller);
            for (const TramoKernelLogEvent& e : fact.events)
                events.push_back(e);
            orders = fact.orders;
            mean = fact.mean;
            outliers = fact.outliers;
        }

        {
            const TramoKernelSeam& end =
                expect(TramoKernelSeamKind::iterate_end, "the iterate tail");
            push_decision(TramoKernelSeamKind::iterate_end);
            run.decisions.back().flag =
                static_cast<int>(TramoKernelIterateEnd::done);
            require(end.end == TramoKernelIterateEnd::done,
                    "the iterate end should be done");
        }
        converged = true;
    }

    events.push_back({TramoKernelLogKind::pop, 1, "-"});
    events.push_back({TramoKernelLogKind::pop, 0, "-"});
    run.converged = converged;
    require(next == seams.size(),
            "seam facts are left over after the round loop (fact " +
                std::to_string(next) + " of " + std::to_string(seams.size()) +
                ")");

    run.final_state = {pass, round, pass3, need_out, need_auto, curva, pcr,
                       ref_auto, ref_stats, seasonal, original_seasonality};
    return run;
}

}  // namespace kalos
