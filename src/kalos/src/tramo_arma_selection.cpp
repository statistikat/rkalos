// Item 23-D transcription of the arma selection of the live
// tramo.internal.ArmaModule over a RegSarimaModelling (docs/contracts/
// tramo-arma-selection-contract.md carries the coordinate map). The stage
// constants are the pinned literals: the panel capacity NMOD = 5 (the
// ArmaModelSelector builder default modelsCount), the comespa descent from
// maxInic(period) = 3, the fast-BIC tie-break tiers {0.003 -> 0.0625,
// 0.03 -> 0.25, else 1.0} scaling the accept budgets {0.01, 0.0025, 0.0075}
// and the per-accept budget relief a -= bic - bic0.
#include "kalos/tramo_arma_selection.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/filter_utility.hpp"
#include "kalos/detail/hannan_rissanen_pipeline.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/hannan_rissanen.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

// JVM RuntimeException inside the module try block; the catch turns any of
// them into the "arma selection failed" remark and the Failed result.
struct SelectionAbort {
    std::string message;
};

// SarmaOrders over the fixed monthly period (the selection-space orders:
// p, q, bp, bq; SarimaOrders.equals reduces to these four once the period,
// d and bd agree, which they do inside the selector).
struct SarmaOrders {
    int p = 0, q = 0, bp = 0, bq = 0;

    int parameters_count() const { return p + bp + q + bq; }
    bool equals(const SarmaOrders& o) const {
        return p == o.p && q == o.q && bp == o.bp && bq == o.bq;
    }
};

std::string to_string(const SarmaOrders& o) {
    return "(" + std::to_string(o.p) + "," + std::to_string(o.q) + ")("
        + std::to_string(o.bp) + "," + std::to_string(o.bq) + ")";
}

// Orders of the ladder gates (SarimaOrders with its own period, as
// TramoUtility.autlar reads sarimaOrders.getPeriod() internally).
struct Orders {
    int p = 0, d = 0, q = 0, bp = 0, bq = 0, bd = 0;
    int period = 12;
};

// TramoUtility.autlar(n, orders) (same transcription as the item 23-B scan
// gate and the item 23-C ladder; note the JVM integer log argument n2 * n2
// and the integer divisions).
int tramo_autlar(int n, const Orders& o) {
    const int d_total = o.d + o.period * o.bd;
    const int q_total = o.q + o.period * o.bq;
    const int p_total = o.p + o.period * o.bp;
    const int n2 = n - d_total;
    int k = static_cast<int>(std::log(static_cast<double>(n2 * n2)));
    const int kmax = std::max(2 * q_total, k);
    if (kmax > k) {
        k = kmax;
    }
    if (!(k < n2)) {
        k = n2 - n2 / 4;
    }
    if (k > 50) {
        k = 50;
    }
    const int tmp = o.p + o.bp * (1 + o.p) + o.q + o.bq * (1 + o.q);
    return (n2 - k) - std::max(p_total, q_total) - tmp;
}

// ArmaModule.maxInic(freq).
int max_inic(int period) {
    if (period == 2) {
        return 1;
    }
    if (period == 3) {
        return 2;
    }
    return 3;
}

// ArmaModule.checkmaxspec: the autlar envelope candidate. m = 1..3 take
// p = q = m with bq = 1 and bp = 1 only on an undifferenced seasonal;
// m = 4 takes (3,3)(2,2) under seasonality.
Orders checkmaxspec(int period, int m, int d, int bd, bool seasonal) {
    Orders o;
    o.period = period;
    o.d = d;
    if (seasonal) {
        o.bd = bd;
    }
    switch (m) {
    case 1:
    case 2:
    case 3:
        o.p = m;
        o.q = m;
        if (seasonal) {
            if (bd == 0) {
                o.bp = 1;
            }
            o.bq = 1;
        }
        break;
    case 4:
        o.p = 3;
        o.q = 3;
        if (seasonal) {
            o.bp = 2;
            o.bq = 2;
        }
        break;
    default:
        break;
    }
    return o;
}

// ArmaModule.calcmaxspec: the candidate envelope the module builder reads.
// m = 1..3 take p = q = m with bp = bq = 1 under seasonality (unconditional
// on bd); m = 4 takes (3,3)(2,2) unconditionally.
Orders calcmaxspec(int period, int m, int d, int bd, bool seasonal) {
    Orders o;
    o.period = period;
    o.d = d;
    if (seasonal) {
        o.bd = bd;
    }
    switch (m) {
    case 1:
    case 2:
    case 3:
        o.p = m;
        o.q = m;
        if (seasonal) {
            o.bp = 1;
            o.bq = 1;
        }
        break;
    case 4:
        o.p = 3;
        o.q = 3;
        o.bp = 2;
        o.bq = 2;
        break;
    default:
        break;
    }
    return o;
}

// ArmaModule.checkespa(freq, n, m, d, bd, seasonal). The second gate is the
// bare autlar sign on the p = n - p - period * bp length (the pinned bytecode
// keeps no subtraction; the n = 56 flat capture discriminates).
bool checkespa(int period, int n, int m, int d, int bd, bool seasonal) {
    Orders spec = checkmaxspec(period, m, d, bd, seasonal);
    if (tramo_autlar(n, spec) < 0) {
        return false;
    }
    const int p = n - spec.p - period * spec.bp;
    spec.p = 0;
    spec.bp = 0;
    return tramo_autlar(p, spec) >= 0;
}

// ArmaModule.comespa: descent from maxInic; monthly reaches m = 3 or 2
// (the loop never checks m = 1), else 0.
int comespa(int period, int n, int i0, int d, int bd, bool seasonal) {
    for (int i = i0; i > 1; --i) {
        if (checkespa(period, n, i, d, bd, seasonal)) {
            return i;
        }
    }
    return 0;
}

// SarimaOrders.setDefault(seasonal): airline with no regular AR; the
// seasonal pair only on a seasonal period.
void set_default(Orders& o, bool seasonal) {
    o.p = 0;
    o.d = 1;
    o.q = 1;
    o.bp = 0;
    if (o.period > 1 && seasonal) {
        o.bd = 1;
        o.bq = 1;
    }
}

// SarimaOrders.isAirline(seasonal).
bool is_airline(const Orders& o, bool seasonal) {
    const bool base = o.p == 0 && o.q == 1 && o.d == 1;
    if (!base || !seasonal) {
        return base;
    }
    return o.bp == 0 && o.bq == 1 && o.bd == 1;
}

// ---------------------------------------------------------------------------
// HannanRissanen candidate fits and the fast-BIC panel machine
// ---------------------------------------------------------------------------

// SarimaModel parameter blocks in JVM order [phi | bphi | theta | btheta].
struct SarModel {
    std::vector<double> phi, bphi, th, bth;

    std::vector<double> parameters() const {
        std::vector<double> p;
        for (const std::vector<double>* b : {&phi, &bphi, &th, &bth}) {
            p.insert(p.end(), b->begin(), b->end());
        }
        return p;
    }
};

// FilterUtility.checkStability block step (same transcription as the item
// 23-C ladder): trailing zeros trimmed, then the roots of the coefficient
// polynomial strictly inside the unit circle.
bool block_stable(const std::vector<double>& block) {
    std::vector<double> trimmed = block;
    while (!trimmed.empty() && trimmed.back() == 0.0) {
        trimmed.pop_back();
    }
    return trimmed.empty() || detail::check_roots(trimmed.data(), trimmed.size(), 1.0);
}

// SarimaModel.isStable(true): all four blocks inside the unit circle.
bool model_stable(const SarModel& m) {
    return block_stable(m.phi) && block_stable(m.bphi) && block_stable(m.th)
        && block_stable(m.bth);
}

// ArmaModelSelector$FastBIC over a fitted candidate: orders (getSpecification
// = the model orders, the HR pipeline never rewrites the requested spec),
// the parameters count behind getParametersCount() and the fast BIC
// FastKalmanFilter.fastProcessing(csr, k) = log(ssq/n) + k log(n)/n with the
// DataBlock.ssq() forward fold over the finite innovations.
struct FastBic {
    SarmaOrders orders;
    int count = 0;
    double bic = 0.0;
};

FastBic make_fast_bic(const std::vector<double>& csr, const SarModel& model,
    const SarmaOrders& orders) {
    FastBic fb;
    fb.orders = orders;
    fb.count = orders.parameters_count();
    const auto innovations = detail::pinned_fast_filter_residuals(csr,
        static_cast<std::size_t>(orders.p), static_cast<std::size_t>(orders.bp),
        static_cast<std::size_t>(orders.q), static_cast<std::size_t>(orders.bq),
        12, model.parameters());
    if (!innovations) {
        // The JVM aborts the machine here (short series index failure, NaN h
        // ArimaException, non-stationary acgf guard); the RuntimeException
        // escapes sort() into the module catch.
        throw SelectionAbort{"fastFilter aborted on " + to_string(orders)};
    }
    double ssq = 0.0; // DataBlock.ssq(): finite-guarded forward fold
    for (const double v : *innovations) {
        if (std::isfinite(v)) {
            ssq += v * v;
        }
    }
    const double n = static_cast<double>(innovations->size());
    fb.bic = detail::fdlibm_log(ssq / n)
        + static_cast<double>(fb.count) * detail::fdlibm_log(n) / n;
    return fb;
}

// Double.compare: the natural order behind List.sort(naturalOrder) and
// FastBIC.compareTo (NaN greatest, -0.0 below +0.0).
int double_compare(double x, double y) {
    if (std::isnan(x)) {
        return std::isnan(y) ? 0 : 1;
    }
    if (std::isnan(y)) {
        return -1;
    }
    if (x < y) {
        return -1;
    }
    if (x > y) {
        return 1;
    }
    if (std::signbit(x) == std::signbit(y)) {
        return 0;
    }
    return std::signbit(x) ? -1 : 1;
}

// ArmaModelSelector.sort(csr, orders): HannanRissanen.process over every
// candidate, keep the stable fits, sort ascending by BIC (the TimSort is
// stable). HannanRissanen.process returns true on every path, so a pipeline
// exception is the only per-candidate failure and it escapes to the module
// catch; the JVM builder defaults match the qualified pipeline seam.
std::vector<FastBic> sort_candidates(const std::vector<double>& csr,
    const std::vector<SarmaOrders>& grid) {
    std::vector<FastBic> kept;
    for (const SarmaOrders& cand : grid) {
        const auto blocks = detail::hannan_rissanen_pipeline(csr,
            HannanRissanenOrders{cand.p, cand.bp, cand.q, cand.bq, 12});
        SarModel model;
        model.phi = blocks[0];
        model.bphi = blocks[1];
        model.th = blocks[2];
        model.bth = blocks[3];
        model.phi.resize(static_cast<std::size_t>(cand.p), 0.0);
        model.bphi.resize(static_cast<std::size_t>(cand.bp), 0.0);
        model.th.resize(static_cast<std::size_t>(cand.q), 0.0);
        model.bth.resize(static_cast<std::size_t>(cand.bq), 0.0);
        if (!model_stable(model)) {
            continue;
        }
        kept.push_back(make_fast_bic(csr, model, cand));
    }
    std::stable_sort(kept.begin(), kept.end(),
        [](const FastBic& a, const FastBic& b) { return double_compare(a.bic, b.bic) < 0; });
    return kept;
}

// ArmaModelSelector.getPreferredSpecification: with a single row or an
// acceptwn selector the best row; else the first row with parameters. An
// all-white-noise multi-row panel trips the JVM array-index failure.
SarmaOrders get_preferred_specification(const std::vector<FastBic>& tmp, bool acceptwn) {
    if (tmp.size() == 1 || acceptwn) {
        return tmp[0].orders;
    }
    std::size_t i = 0;
    while (i < tmp.size() && tmp[i].count == 0) {
        ++i;
    }
    if (i == tmp.size()) {
        throw SelectionAbort{"white-noise-only panel (JVM array index failure)"};
    }
    return tmp[i].orders;
}

// The five-slot hrModels panel (builder default modelsCount = NMOD = 5).
constexpr int kNmod = 5;
using Panel = std::array<std::optional<FastBic>, kNmod>;

// ArmaModelSelector.mergeInto(src, dst): insert each src row at the first
// dst slot at or after the watermark that is empty (take it), duplicates the
// existing row (skip it), or is beaten on BIC (shift the tail down and take
// the slot). A scanned slot with a better BIC just advances k; exhausting
// dst leaves the watermark untouched for the next src row.
void merge_into(const std::vector<FastBic>& src, Panel& dst) {
    const int m = static_cast<int>(dst.size());
    const int n = std::min<int>(static_cast<int>(src.size()), m);
    int j = 0;
    for (int i = 0; i < n && j < m; ++i) {
        for (int k = j; k < m; ++k) {
            if (!dst[k]) {
                dst[k] = src[i];
                j = k + 1;
                break;
            }
            if (dst[k]->orders.equals(src[i].orders)) {
                j = k + 1;
                break;
            }
            if (dst[k]->bic - src[i].bic > 0.0) { // dcmpl ifle continues on <=
                for (int l = m - 1; l > k; --l) {
                    dst[l] = dst[l - 1];
                }
                dst[k] = src[i];
                j = k + 1;
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ArmaModelSelector
// ---------------------------------------------------------------------------

class ArmaModelSelector {
public:
    ArmaModelSelector(int models_count, bool acceptwn, int max_p, int max_q, int max_bp,
        int max_bq)
        : nmodels_(models_count), acceptwn_(acceptwn), maxP_(max_p), maxQ_(max_q),
          maxBp_(max_bp), maxBq_(max_bq) {}

    const Panel& hr_models() const { return hr_; }

    // process(csr, freq, d, bd, seasonal); nullopt mirrors the JVM null
    // return (the module dereferences it, so it fails closed as Failed).
    std::optional<SarmaOrders> process(const std::vector<double>& csr, bool seasonal, int d,
        int bd) {
        hr_ = Panel{};
        // The preferred base object the JVM mutates and clones (p = 3, q = 0
        // after the stage-1 preamble; its bp/bq tail is never read back).
        SarmaOrders base;
        base.p = 3;
        base.q = 0;
        std::optional<SarmaOrders> preferred;

        if (seasonal) {
            auto grid = seasonal_grid(3, 0);
            auto tmp = sort_candidates(csr, grid);
            if (tmp.empty()) {
                for (SarmaOrders& cand : grid) {
                    cand.p = 1; // the setP(1) re-entry on the same candidates
                }
                tmp = sort_candidates(csr, grid);
                if (tmp.empty()) {
                    return std::nullopt;
                }
            }
            preferred = get_preferred_specification(tmp, acceptwn_);
            if (base.p <= maxP_) { // the mutable base still carries p = 3
                merge_into(tmp, hr_);
            }
        } else {
            preferred = base;
        }

        {
            std::vector<SarmaOrders> grid;
            for (int p = 0; p <= maxP_; ++p) {
                for (int q = 0; q <= maxQ_; ++q) {
                    SarmaOrders cand = *preferred;
                    cand.p = p;
                    cand.q = q;
                    grid.push_back(cand);
                }
            }
            auto tmp = sort_candidates(csr, grid);
            if (tmp.empty()) {
                return std::nullopt;
            }
            preferred = get_preferred_specification(tmp, acceptwn_);
            merge_into(tmp, hr_);
        }

        if (seasonal) {
            auto grid = seasonal_grid(preferred->p, preferred->q);
            auto tmp = sort_candidates(csr, grid);
            if (tmp.empty()) {
                return std::nullopt;
            }
            merge_into(tmp, hr_);
        }

        if (!seasonal) {
            if (!hr_[0]) {
                throw SelectionAbort{"empty hrModels panel (JVM null dereference)"};
            }
            if (hr_[1] && hr_[0]->count == 0 && !acceptwn_) {
                return hr_[1]->orders; // white noise on top: take the runner-up
            }
            return hr_[0]->orders;
        }
        return select(d, bd);
    }

private:
    std::vector<SarmaOrders> seasonal_grid(int fixed_p, int fixed_q) const {
        std::vector<SarmaOrders> grid;
        for (int bp = 0; bp <= maxBp_; ++bp) {
            for (int bq = 0; bq <= maxBq_; ++bq) {
                SarmaOrders cand;
                cand.p = fixed_p;
                cand.q = fixed_q;
                cand.bp = bp;
                cand.bq = bq;
                grid.push_back(cand);
            }
        }
        return grid;
    }

    // select(d, bd): the eleven accept gates over the merged panel, all
    // measured against the current best row with the tier-scaled budgets
    // relieved by the BIC improvement on every accept.
    std::optional<SarmaOrders> select(int p, int q) {
        int n = nmodels_;
        for (;;) {
            if (n < 1) {
                throw SelectionAbort{"all-null hrModels (JVM array index failure)"};
            }
            if (!hr_[n - 1]) {
                --n;
                continue;
            }
            break; // FastBIC.arma is never null
        }
        if (n == 0) {
            return std::nullopt;
        }
        if (n == 1) {
            return hr_[0]->orders;
        }
        SarmaOrders cur = hr_[0]->orders;
        int p0 = cur.p + cur.q;
        int q0 = cur.bp + cur.bq;
        int dist_p = std::abs(cur.p + p - cur.q);
        int dist_q = std::abs(cur.bp + q - cur.bq);
        const double diff = hr_[n - 1]->bic - hr_[0]->bic;
        const double alpha = diff < 0.003 ? 0.0625 : (diff < 0.03 ? 0.25 : 1.0);
        double a = 0.01 * alpha;
        double b = 0.0025 * alpha;
        double c = 0.0075 * alpha;
        int best = 0;
        for (int i = 1; i < n; ++i) {
            const SarmaOrders& o = hr_[i]->orders;
            const int p1 = o.p + o.q;
            const int q1 = o.bp + o.bq;
            const int dp = std::abs(o.p + p - o.q);
            const int dq = std::abs(o.bp + q - o.bq);
            const double bic = hr_[i]->bic - hr_[best]->bic;
            const bool accept
                = ((dp < dist_p || dq < dist_q) && p0 == p1 && q0 == q1 && bic <= a)
                || (dp < dist_p && p1 <= p0 && q1 == q0 && bic <= b && o.p > 0 && o.q > 0)
                || (((dp == 0 && dp < dist_p && p > 0) || (dq == 0 && dq < dist_q && q > 0))
                       && p0 == p1 && q0 == q1 && bic <= a)
                || (dp == 0 && dq == 0 && bic < b)
                || (p1 > p0 && dp == 0 && q1 == q0 && bic < b)
                || (i == 1 && p0 == 0 && p1 == 1 && q1 == q0 && bic < b)
                || (q1 > q0 && dq == 0 && p1 == p0 && bic < b)
                || (p1 < p0 && p1 > 0 && q1 == q0 && bic < b)
                || (o.p < cur.p && o.q == cur.q && p1 > 0 && q1 == q0 && bic < c)
                || (q1 < q0 && q1 > 0 && p1 == p0 && dq == 0 && bic < b)
                || (q1 < q0 && q1 > 0 && p1 == p0 && bic < b);
            if (accept) {
                const double relief = hr_[i]->bic - hr_[0]->bic;
                a -= relief;
                b -= relief;
                c -= relief;
                p0 = p1;
                q0 = q1;
                dist_p = dp;
                dist_q = dq;
                best = i;
                cur = o;
            }
        }
        if (cur.parameters_count() == 0 && best < nmodels_ - 1) {
            return hr_[best + 1]->orders;
        }
        return hr_[best]->orders;
    }

    int nmodels_;
    bool acceptwn_;
    int maxP_, maxQ_, maxBp_, maxBq_;
    Panel hr_;
};

} // namespace

TramoArmaSelectionResult select_tramo_arma(const TramoArmaSelectionRequest& request) {
    if (request.period != 12) {
        throw std::invalid_argument(
            "select_tramo_arma: only monthly series are in the ported domain");
    }
    if (request.residuals.empty()) {
        throw std::invalid_argument("select_tramo_arma: empty residual block");
    }
    for (const double v : request.residuals) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("select_tramo_arma: non-finite residual");
        }
    }

    TramoArmaSelectionResult result;
    Orders spec;
    spec.period = request.period;
    spec.p = request.p0;
    spec.d = request.d0;
    spec.q = request.q0;
    spec.bp = request.bp0;
    spec.bd = request.bd0;
    spec.bq = request.bq0;

    result.events.push_back("push arma selection");
    try {
        const int pcnt = comespa(spec.period, request.nobs, max_inic(spec.period), spec.d,
            spec.bd, request.seasonal);
        if (pcnt == 0) {
            result.default_model = true;
            result.events.push_back("remark default model selected (not enough obs.)");
            if (!is_airline(spec, request.seasonal)) {
                set_default(spec, request.seasonal);
                result.result = TramoArmaSelectionResult::Result::Changed;
            } else {
                result.result = TramoArmaSelectionResult::Result::Unprocessed;
            }
        } else {
            const Orders maxspec = calcmaxspec(
                spec.period, pcnt, spec.d, spec.bd, request.seasonal);
            ArmaModelSelector selector(
                kNmod, request.accept_white_noise, maxspec.p, maxspec.q, maxspec.bp,
                maxspec.bq);
            const auto selection = selector.process(
                request.residuals, request.seasonal, maxspec.d, maxspec.bd);
            if (!selection) {
                // log.info("selected model " + sel.toString()) dereferences
                // the null selection (NullPointerException into the catch).
                throw SelectionAbort{"null selection (JVM null dereference)"};
            }
            result.has_selection = true;
            result.selection = TramoSarmaOrders{selection->p, selection->q, selection->bp,
                selection->bq};
            for (const auto& row : selector.hr_models()) {
                TramoArmaPanelRow out;
                if (row) {
                    out.present = true;
                    out.orders = TramoSarmaOrders{
                        row->orders.p, row->orders.q, row->orders.bp, row->orders.bq};
                    out.bic = row->bic;
                }
                result.panel.push_back(out);
            }
            result.events.push_back("info selected model: " + to_string(*selection));
            const SarmaOrders stationary{spec.p, spec.q, spec.bp, spec.bq};
            if (selection->equals(stationary)) {
                result.result = TramoArmaSelectionResult::Result::Unchanged;
            } else {
                // SarimaOrders.of(sel, spec.getD(), spec.getBd()).
                spec.p = selection->p;
                spec.q = selection->q;
                spec.bp = selection->bp;
                spec.bq = selection->bq;
                result.result = TramoArmaSelectionResult::Result::Changed;
            }
        }
    } catch (const SelectionAbort&) {
        result.selection_failed = true;
        result.result = TramoArmaSelectionResult::Result::Failed;
        result.has_selection = false;
        result.panel.clear();
        result.events.push_back("remark arma selection failed");
    } catch (const std::exception&) {
        // Dated 2026-11-23 (javap evidence): ArmaModelSelector carries no
        // Exception table at all (no per-candidate try in sort()/Arma()),
        // and ArmaModule.process catches the escape into
        // ProcessingResult.Failed. A HannanRissanen pipeline RuntimeException
        // (the pinned short-window design class) therefore kills the module
        // as Failed exactly like a SelectionAbort - it is never a throw to
        // the caller (fit_shortfail witness: the JVM survives short residual
        // blocks through the comespa/airline entry gates, not a selector-
        // internal candidate catch).
        result.selection_failed = true;
        result.result = TramoArmaSelectionResult::Result::Failed;
        result.has_selection = false;
        result.panel.clear();
        result.events.push_back("remark arma selection failed");
    }
    result.events.push_back("pop");

    result.spec_p = spec.p;
    result.spec_d = spec.d;
    result.spec_q = spec.q;
    result.spec_bp = spec.bp;
    result.spec_bd = spec.bd;
    result.spec_bq = spec.bq;
    return result;
}

} // namespace kalos
