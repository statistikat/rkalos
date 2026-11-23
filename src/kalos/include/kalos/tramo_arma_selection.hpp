#ifndef KALOS_TRAMO_ARMA_SELECTION_HPP
#define KALOS_TRAMO_ARMA_SELECTION_HPP

// Transcription of the arma selection of the live tramo.internal.ArmaModule
// of the pinned jdplus 3.9.0 stack (item 23-D of
// docs/tramo-seats-roadmap.md, docs/contracts/tramo-arma-selection-contract.md
// carries the coordinate map). The unit owns the whole module stage: the
// comespa size gate with its default-model/airline branches, the calcmaxspec
// candidate envelope, the three-stage ArmaModelSelector (stage-1 seasonal
// AR block grid, stage-2 ARMA grid on the preferred model, stage-3 seasonal
// grid) whose per-candidate fits run the qualified HannanRissanen pipeline
// seam and the pinned FastKalmanFilter.fastProcessing fast BIC
// log(ssq/n) + k log(n)/n - no maximum-likelihood machine is involved (the
// force-load witness pins HannanRissanen/FastKalmanFilter loads against zero
// RegSarimaComputer/GlsSarimaComputer/HannanRissanenInitializer/
// LevenbergMarquardtMinimizer loads), the mergeInto panel merge into the
// five-row hrModels panel, and the eleven-gate select() tie-break over that
// panel.
//
// Scope envelope (register D11): monthly data only. Requests with period !=
// 12, an empty residual block or non-finite residuals are refused with
// std::invalid_argument; the unit never falls back to the X-13 selection.
//
// The residuals are RegArimaUtility.olsResiduals(regarima()) as the module
// reads them - request data owned by the upstream regression stage (the
// JVM-captured golden vectors pin the hand-over; the quarterly-shape and the
// degenerate-constant captures pin the refusal and failure semantics). The
// request also carries regarima().getObservationsCount() separately: the
// comespa gate sizes on the observations, not on the residual block.

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace kalos {

struct TramoArmaSelectionRequest {
    // RegArimaUtility.olsResiduals(desc.regarima()): the mean/log-corrected,
    // regression- and outlier-adjusted series the selector consumes. Must be
    // complete and finite.
    std::vector<double> residuals;
    // desc.regarima().getObservationsCount() behind the comespa gate.
    int nobs = 0;
    // SarimaOrders.getPeriod() of the entry specification (12; D11).
    int period = 12;

    // Entry specification of the description.
    int p0 = 0, d0 = 0, q0 = 0, bp0 = 0, bd0 = 0, bq0 = 0;

    // ArmaModule builder pins (TR defaults): seasonal = true, wn = false.
    // seasonal routes the stage-1/stage-3 grids and the select() return;
    // accept_white_noise is the selector acceptwn flag.
    bool seasonal = true;
    bool accept_white_noise = false;
};

struct TramoSarmaOrders {
    int p = 0, q = 0, bp = 0, bq = 0;
};

// One row of the CInfo panel (ArmaModelSelector.gePreferredModels, five
// slots): trailing nulls stay present = false with a NaN BIC, exactly as
// CInfo.models()/bic() expose them.
struct TramoArmaPanelRow {
    bool present = false;
    TramoSarmaOrders orders;
    double bic = std::numeric_limits<double>::quiet_NaN();
};

struct TramoArmaSelectionResult {
    enum class Result {
        Changed,      // specification replaced by SarimaOrders.of(sel, d, bd)
        Unchanged,    // selection equals the stationary part of the entry
        Unprocessed,  // comespa gate with an already-airline entry
        Failed,       // RuntimeException path: "arma selection failed"
    };
    Result result = Result::Unchanged;

    // Which remark fired, if any: default_model is the pcnt == 0 remark
    // "default model selected (not enough obs.)", selection_failed the
    // RuntimeException catch remark "arma selection failed".
    bool default_model = false;
    bool selection_failed = false;

    // Description specification after the module (entry orders unless the
    // stage changed them).
    int spec_p = 0, spec_d = 0, spec_q = 0;
    int spec_bp = 0, spec_bd = 0, spec_bq = 0;

    // CInfo of the "selected model" info line; absent on the gate and the
    // failure paths (the JVM never reaches the log line there).
    bool has_selection = false;
    TramoSarmaOrders selection;
    std::vector<TramoArmaPanelRow> panel; // five rows, JVM order

    // ProcessingLog sequence of the module section: "push arma selection",
    // optional remark, optional "info selected model (p,q)(bp,bq)", "pop".
    std::vector<std::string> events;
};

// Runs the live ArmaModule.process(RegSarimaModelling) stage. Throws
// std::invalid_argument outside the ported domain (see the scope envelope).
// Pinned numerical failures of the candidate fits (HannanRissanen pipeline
// exceptions, fast-filter aborts, the empty-panel null selection) follow the
// JVM RuntimeException path into the Failed result and never throw.
TramoArmaSelectionResult select_tramo_arma(const TramoArmaSelectionRequest& request);

} // namespace kalos

#endif // KALOS_TRAMO_ARMA_SELECTION_HPP
