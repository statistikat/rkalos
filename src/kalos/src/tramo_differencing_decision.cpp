// Item 23-C transcription of the live round-1 differencing decision,
// tramo.internal.DifferencingModule (docs/contracts/
// tramo-differencing-decision-contract.md carries the coordinate map). The
// ladder constants are the pinned literals: MAXD = 2, MAXBD = 1, the
// near-unit-root tolerance 1.02, the out-of-range modulus 1.12, the root
// windows 0.05 / 0.02, the din source 1.005 - ub2, the per-iteration
// cancellation decay c -= 0.002, the default fills -0.1 (AR) / -0.2 (MA) and
// the isMeanCorrection length table {<=80: 1.96, <=155: 1.98, <=230: 2.1,
// <=320: 2.3, else 2.5}.
#include "kalos/tramo_differencing_decision.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/filter_utility.hpp"
#include "kalos/detail/hannan_rissanen_pipeline.hpp"
#include "kalos/detail/sarima_stabilization.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/linear_algebra.hpp"

#include "mlj_engine.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

struct SarModel {
    std::vector<double> phi, bphi, th, bth;

    double phi1() const { return phi.empty() ? 0.0 : phi[0]; }
    double bphi1() const { return bphi.empty() ? 0.0 : bphi[0]; }
    double th1() const { return th.empty() ? 0.0 : th[0]; }
    double bth1() const { return bth.empty() ? 0.0 : bth[0]; }

    // JVM SarimaModel.parameters(): [phi | bphi | theta | btheta].
    std::vector<double> parameters() const {
        std::vector<double> p;
        for (const std::vector<double>* b : {&phi, &bphi, &th, &bth}) {
            p.insert(p.end(), b->begin(), b->end());
        }
        return p;
    }
};

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

// SarimaMapping.stabilize block step at a declared length (the pinned model
// keeps its order even when the stabilizer drops trailing zeros).
void stabilize_block(std::vector<double>& block) {
    if (block.empty()) {
        return;
    }
    std::vector<double> work = block;
    detail::stabilize_coefficient_block(work);
    work.resize(block.size(), 0.0);
    block = std::move(work);
}

void stabilize_model(SarModel& m) {
    stabilize_block(m.phi);
    stabilize_block(m.bphi);
    stabilize_block(m.th);
    stabilize_block(m.bth);
}

// Pinned 1 + c(1) ... value of a coefficient block summed at B = 1 / B = -1
// (the tryUrpCancelling / tryUrmCancelling gate polynomials).
double block_at_one(const std::vector<double>& block) {
    double s = 1.0;
    for (const double v : block) {
        s += v;
    }
    return s;
}

double block_at_minus_one(const std::vector<double>& block) {
    double s = 1.0;
    bool minus = false;
    for (const double v : block) {
        s += minus ? -v : v;
        minus = !minus;
    }
    return s;
}

double data_block_average(const std::vector<double>& x) {
    double s = 0.0;
    for (const double v : x) {
        s += v;
    }
    return s / static_cast<double>(x.size());
}

struct Orders {
    int p = 0, d = 0, q = 0, bp = 0, bq = 0, bd = 0;
    int period = 12;
};

// Pinned internal failure of the ladder (TramoException in the JVM); the
// process() catch turns any of them into the airline fallback.
struct TramoLadderError {
    std::string message;
};

// Warm maximum-likelihood re-estimation of initstep:
// TramoUtility.processor(true, eps).optimize(regarima, null), i.e. the
// GlsSarimaComputer ML pass from the stabilized landing over the raw
// response with the meanCorrection constant column. The qualified seam is
// the 23-B warm Levenberg-Marquardt machine; a landing that reaches the
// private optimize tail's unit-root-cancellation restarts fails closed.
void optimize_landing(const Orders& o, const SarModel& in, const std::vector<double>& y,
    double eps, SarModel& out) {
    SarimaEstimationRequest warm;
    warm.series.values = y;
    warm.series.missing.assign(y.size(), 0);
    warm.design.columns = 1;
    warm.design.values.assign(y.size(), 1.0);
    warm.design.coefficients = {0.0};
    warm.design.fixed = {0};
    warm.blocks[0].initial.assign(static_cast<std::size_t>(o.p), 0.0);
    warm.blocks[0].fixed.assign(static_cast<std::size_t>(o.p), 0);
    warm.blocks[1].initial.assign(static_cast<std::size_t>(o.q), 0.0);
    warm.blocks[1].fixed.assign(static_cast<std::size_t>(o.q), 0);
    warm.blocks[2].initial.assign(static_cast<std::size_t>(o.bp), 0.0);
    warm.blocks[2].fixed.assign(static_cast<std::size_t>(o.bp), 0);
    warm.blocks[3].initial.assign(static_cast<std::size_t>(o.bq), 0.0);
    warm.blocks[3].fixed.assign(static_cast<std::size_t>(o.bq), 0);
    warm.period = o.period;
    warm.d = o.d;
    warm.seasonal_d = o.bd;
    warm.control = EstimationControl{eps, 1e-4, true, false, true, true,
        EstimationStartingPoint::HannanRissanen,
        EstimationSolver::PinnedLevenbergMarquardt};

    std::array<std::vector<double>, 4> start;
    start[0].reserve(in.phi.size());
    for (const double v : in.phi) {
        start[0].push_back(-v);
    }
    start[1] = in.th;
    start[2].reserve(in.bphi.size());
    for (const double v : in.bphi) {
        start[2].push_back(-v);
    }
    start[3] = in.bth;

    SarimaEstimationResult est;
    try {
        est = optimize_sarima_from(warm, start);
    } catch (const std::exception&) {
        throw TramoLadderError{"Non convergence in ESPDIF"};
    }
    SarModel result;
    for (const double v : est.coefficients[0]) {
        result.phi.push_back(-v);
    }
    for (const double v : est.coefficients[2]) {
        result.bphi.push_back(-v);
    }
    result.th = est.coefficients[1];
    result.bth = est.coefficients[3];
    // The pinned optimize returns null only on a thrown machine failure
    // (mirrored by the catch above); function-tolerance convergence counts,
    // matching the item 23-B qualification of the same seam. Non-finite
    // landings follow the pinned Non-convergence path.
    for (const auto& block : est.coefficients) {
        for (const double v : block) {
            if (!std::isfinite(v)) {
                throw TramoLadderError{"Non convergence in ESPDIF"};
            }
        }
    }

    // Fail-closed seam guard: the private optimize tail runs finalProcessing
    // (tryUrpCancelling, then tryUrmCancelling) on this landing. Their gates
    // (orders and the |1 + phi(B)|, |1 + theta(B)| threshold 0.1) are not
    // part of the qualified warm machine; a reachable gate throws instead of
    // silently skipping the pinned restart attempts. The JVM reads the
    // SarimaModel parameters, so the AR blocks enter sign-flipped.
    const bool urp = o.p > 0 && o.q > 0 && o.d > 0
        && std::fabs(block_at_one(result.phi)) <= 0.1
        && std::fabs(block_at_one(result.th)) <= 0.1;
    const bool urm = o.bp > 0 && o.bq > 0 && o.bd > 0
        && std::fabs(block_at_minus_one(result.bphi)) <= 0.1
        && std::fabs(block_at_minus_one(result.bth)) <= 0.1;
    if (urp || urm) {
        throw std::runtime_error("decide_tramo_differencing: the maximum-likelihood "
                                 "landing reaches the unit-root cancellation restarts, "
                                 "outside the qualified warm seam");
    }

    // est.getModel().arima().stationaryTransformation().getStationaryModel():
    // the stable representative of the landing (root inversion per block,
    // the same routine SarimaMapping.stabilize drives).
    stabilize_model(result);
    out = std::move(result);
}

struct Ladder {
    const TramoDifferencingRequest& rq;
    const int maxd = 2;
    int maxbd = 0;

    std::vector<double> y; // module field y (the request series)
    Orders spec;
    SarModel last;
    bool has_last = false;

    double rmax = 0.0, rsmax = 0.0, c = 0.0, din = 0.0, tmean = 0.0;
    int iter = 0;
    bool ml = false, useml = false, mlused = false;

    explicit Ladder(const TramoDifferencingRequest& request)
        : rq(request), maxbd(request.seasonal ? 1 : 0), ml(request.ml) {}

    // TramoUtility.processor(...)-driven warm ML and the raw HannanRissanen
    // step of the ladder.
    void initstep(bool b) {
        if (spec.d == 0 && spec.bd == 0 && b) {
            spec.p = spec.period == 2 ? 1 : 2;
            spec.q = 0;
            spec.bq = 0;
            if (rq.seasonal) {
                spec.bp = 1;
            }
        } else {
            spec.p = 1;
            spec.q = 1;
            if (maxbd > 0) {
                spec.bp = 1;
                spec.bq = 1;
            }
        }

        const auto filter = mlj::differencing_filter(spec.period, spec.d, spec.bd);
        std::vector<double> w;
        if (filter.size() > 1) {
            w.assign(y.size() - (filter.size() - 1), 0.0);
            mlj::apply_backfilter(filter, y, 0, w);
            // Orientation pin: RegArimaUtility.differencingFilter hands the
            // BackFilter the reversed polynomial while mlj::apply_backfilter
            // reads its coefficients ascending, so the JVM correlation runs
            // on the mlj correlation times (-1)^(d + bd) (each reversed
            // (1-B) / (1-B^s) factor contributes one sign; verified element
            // wise against the live BackFilter.apply on both a bd = 1 and a
            // d = bd = 1 case).
            if ((spec.d + spec.bd) % 2 != 0) {
                for (double& v : w) {
                    v = -v;
                }
            }
        } else {
            w = y;
        }
        const double mu = data_block_average(w);
        for (double& v : w) {
            v -= mu;
        }

        // A raw HannanRissanen.process exception exits the entire module.
        // In particular, a short wide design fails at the pinned rawR solve.
        std::array<std::vector<double>, 4> blk;
        try {
            blk = detail::hannan_rissanen_pipeline(
                w, HannanRissanenOrders{spec.p, spec.bp, spec.q, spec.bq, spec.period});
        } catch (const std::runtime_error& err) {
            throw TramoLadderError{err.what()};
        }
        last.phi = blk[0];
        last.bphi = blk[1];
        last.th = blk[2];
        last.bth = blk[3];
        last.phi.resize(static_cast<std::size_t>(spec.p), 0.0);
        last.bphi.resize(static_cast<std::size_t>(spec.bp), 0.0);
        last.th.resize(static_cast<std::size_t>(spec.q), 0.0);
        last.bth.resize(static_cast<std::size_t>(spec.bq), 0.0);
        has_last = true;
        bool err = false;
        if (b && !model_stable(last)) {
            if (spec.p > 1) {
                err = true;
            } else if (spec.p == 1 && std::fabs(last.phi1()) > 1.02) {
                err = true;
            } else if (spec.bp != 1) {
                stabilize_model(last);
            } else if (std::fabs(last.bphi1()) > 1.02) {
                err = true;
            } else {
                stabilize_model(last);
            }
        }
        if (err) {
            // SarimaModel.builder(doStationary).setDefault(): -0.1 AR,
            // -0.2 MA fills at the stationary orders.
            last = SarModel{};
            last.phi.assign(static_cast<std::size_t>(spec.p), -0.1);
            last.bphi.assign(static_cast<std::size_t>(spec.bp), -0.1);
            last.th.assign(static_cast<std::size_t>(spec.q), -0.2);
            last.bth.assign(static_cast<std::size_t>(spec.bq), -0.2);
            has_last = true;
        }

        if (err || ml || useml) {
            stabilize_model(last);
            SarModel ml_model;
            optimize_landing(spec, last, y, rq.eps, ml_model);
            last = std::move(ml_model);
            has_last = true;
            mlused = true;
        } else {
            mlused = false;
        }
        useml = false;
    }

    // searchur over the characteristic roots of the mirrored AR block: the
    // count of near-positive unit roots and the largest remaining near-real
    // modulus (stored into rmax / rsmax).
    int searchur(const std::vector<double>& block, double ub, bool regular) {
        // getRegularAR().mirror().roots(): the full polynomial
        // [1, phi1, ..., phip] reversed, i.e. the block coefficients in
        // descending order with the implicit constant 1 as leading term.
        std::vector<double> mirrored(block.rbegin(), block.rend());
        mirrored.push_back(1.0);
        int count = 0;
        double max = 0.0;
        for (const auto& root : detail::jvm_closed_roots(mirrored)) {
            const double im = std::fabs(root.imag);
            const double modulus = detail::jvm_complex_abs(root.real, root.imag);
            if (modulus >= ub && im <= 0.05 && root.real > 0.0) {
                ++count;
            } else if (im <= 0.02 && root.real > 0.0 && modulus > max) {
                max = modulus;
            }
        }
        if (regular) {
            rmax = max;
        } else {
            rsmax = max;
        }
        return count;
    }

    void step0() {
        initstep(true);
        if (spec.d != 0 || spec.bd != 0) {
            rmax = last.phi1();
            if (maxbd > 0) {
                rsmax = last.bphi1();
            }
        }
        spec.d += searchur(last.phi, rq.ub1, true);
        if (maxbd > 0) {
            spec.bd += searchur(last.bphi, rq.ub1, false);
        }
    }

    int maincondition() {
        const double phi1 = last.phi1();
        const double theta1 = last.th1();
        const double bphi1 = maxbd > 0 ? last.bphi1() : 0.0;
        const double btheta1 = maxbd > 0 ? last.bth1() : 0.0;

        c -= 0.002;
        din = 1.005 - rq.ub2;
        int m = 0;
        if (std::fabs(phi1 + 1.0) <= din) {
            if (-phi1 > 1.02) {
                m = 1;
                useml = true;
            } else if (std::fabs(phi1 - theta1) > c) {
                ++m;
                spec.d += 1;
            }
        } else if (std::fabs(phi1) > 1.12) {
            m = 1;
            useml = true;
        }
        if (maxbd > 0) {
            if (std::fabs(bphi1 + 1.0) <= din) {
                if (-bphi1 > 1.02) {
                    useml = true;
                    m = 1;
                } else if (spec.bd == 0 && std::fabs(bphi1 - btheta1) > c) {
                    ++m;
                    spec.bd += 1;
                    if (useml) {
                        --m;
                        useml = false;
                    }
                }
            } else if (std::fabs(bphi1) > 1.12) {
                m = 1;
                useml = true;
            }
        }
        return m;
    }

    int cond1(int m) {
        if (spec.d + spec.bd == 0) {
            return m;
        }
        const double phi1 = last.phi1();
        const double theta1 = last.th1();
        const double bphi1 = maxbd > 0 ? last.bphi1() : 0.0;
        const double btheta1 = maxbd > 0 ? last.bth1() : 0.0;

        if (std::fabs(phi1 - theta1) >= c && (maxbd <= 0 || std::fabs(bphi1 - btheta1) < c)) {
            if (rmax < 0.9 && rsmax < 0.9) {
                return m;
            }
            if (useml && m == 1) {
                useml = false;
            } else {
                ++m;
            }
            if (rmax > rsmax) {
                spec.d += 1;
            } else {
                spec.bd += 1;
            }
        }
        return m;
    }

    int finalcond(int m) {
        if (m == 2) {
            spec.d -= 1;
            spec.bd -= 1;
            if (mlused) {
                if (last.phi1() < last.bphi1()) {
                    spec.d += 1;
                } else {
                    spec.bd += 1;
                }
            } else if (rmax > rsmax) {
                if (rmax > 0.0) {
                    spec.d += 1;
                }
            } else if (rsmax > 0.0) {
                spec.bd += 1;
            }
        }
        if (spec.d > maxd) {
            spec.d = maxd;
            m = 0;
        }
        if (spec.bd > maxbd) {
            spec.bd = maxbd;
            m = 0;
        }
        return m;
    }

    bool nextstep() {
        initstep(false);
        int m = maincondition();
        if (iter == 0) {
            m = cond1(m);
        }
        return finalcond(m) != 0;
    }

    void computeTMean() {
        std::vector<double> w;
        if (spec.d == 0 && spec.bd == 0) {
            w = y;
        } else {
            if (!has_last) {
                throw TramoLadderError{
                    "Failure in the identification of the differencing orders"};
            }
            stabilize_model(last);
            const auto filter = mlj::differencing_filter(spec.period, spec.d, spec.bd);
            w.assign(y.size() - (filter.size() - 1), 0.0);
            mlj::apply_backfilter(filter, y, 0, w);
            if ((spec.d + spec.bd) % 2 != 0) {
                for (double& v : w) {
                    v = -v; // same parity orientation pin as initstep
                }
            }
            const auto innovations = detail::pinned_fast_filter_residuals(w,
                last.phi.size(), last.bphi.size(), last.th.size(), last.bth.size(),
                spec.period, last.parameters());
            if (!innovations) {
                throw TramoLadderError{"Non convergence in ESPDIF"};
            }
            w = *innovations;
        }
        double s = 0.0, ss = 0.0;
        for (const double v : w) {
            s += v;
            ss += v * v;
        }
        const double n = static_cast<double>(w.size());
        tmean = s / std::sqrt((ss * n - s * s) / n);
    }

    bool calc() {
        c = rq.cancel;
        useml = false;
        mlused = false;
        rmax = 0.0;
        rsmax = 0.0;
        din = 0.0;
        step0();
        iter = 0;
        while (nextstep()) {
            if (iter >= 5) {
                break;
            }
            ++iter;
        }
        computeTMean();
        return true;
    }

    // isMeanCorrection: the pinned length table on the response length.
    bool mean_correction() const {
        const int n = static_cast<int>(y.size());
        double threshold = 2.5;
        if (n <= 80) {
            threshold = 1.96;
        } else if (n <= 155) {
            threshold = 1.98;
        } else if (n <= 230) {
            threshold = 2.1;
        } else if (n <= 320) {
            threshold = 2.3;
        }
        return std::fabs(tmean) > threshold;
    }
};

// TramoUtility.autlar(n, orders) (same transcription as the item 23-B scan
// gate; note the JVM integer log argument n2 * n2 and the integer
// divisions). The orders enter with their own period, as the JVM reads
// sarimaOrders.getPeriod() internally.
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

// comespd(period, n, seasonal): the ladder runs only when the autlar
// heuristics still leave usable data on the (0,2,0)(0, seasonal ? 1 : 0)
// probe model.
bool comespd(int period, int n, bool seasonal) {
    Orders o;
    o.period = period;
    o.d = 2;
    o.bd = seasonal ? 1 : 0;
    return tramo_autlar(n, o) >= 0;
}

} // namespace

TramoDifferencingDecision decide_tramo_differencing(
    const TramoDifferencingRequest& request) {
    std::vector<double> work_series;
    if (request.series.empty()) {
        throw std::invalid_argument("decide_tramo_differencing: empty series");
    }
    for (const double v : request.series) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("decide_tramo_differencing: non-finite value");
        }
    }
    if (request.period < 2) {
        throw std::invalid_argument("decide_tramo_differencing: period below 2");
    }
    if (request.log_input) {
        // The description's log channel (regarima().y): the pinned fdlibm
        // natural log applied per observation, the same application as
        // transform_series(Transform::logarithm).
        for (std::size_t i = 0; i < request.series.size(); ++i) {
            const double v = request.series[i];
            if (v <= 0.0) {
                throw std::invalid_argument(
                    "decide_tramo_differencing: non-positive value under the log channel");
            }
            work_series.push_back(detail::fdlibm_log(v));
        }
    } else {
        work_series = request.series;
    }

    TramoDifferencingDecision decision;
    decision.mean = request.mean0;
    decision.ml = request.ml;
    auto& events = decision.events;
    events.push_back("push differencing selection");

    auto apply_airline = [&] {
        const bool regular = request.p0 == 0 && request.q0 == 1 && request.d0 == 1;
        const bool seasonal = request.bp0 == 0 && request.bq0 == 1 && request.bd0 == 1;
        if (regular && (!request.seasonal || seasonal)) {
            decision.result = TramoDifferencingDecision::Result::Unprocessed;
            decision.spec_p = request.p0;
            decision.spec_d = request.d0;
            decision.spec_q = request.q0;
            decision.spec_bp = request.bp0;
            decision.spec_bd = request.bd0;
            decision.spec_bq = request.bq0;
        } else {
            decision.airline_applied = true;
            decision.result = TramoDifferencingDecision::Result::Changed;
            decision.spec_p = 0;
            decision.spec_d = 1;
            decision.spec_q = 1;
            decision.spec_bp = 0;
            decision.spec_bd = request.seasonal ? 1 : 0;
            decision.spec_bq = request.seasonal ? 1 : 0;
            decision.mean = false;
        }
        events.push_back("pop");
    };
    const int nobs = request.nobs > 0 ? request.nobs
                                      : static_cast<int>(request.series.size());
    if (!comespd(request.period, nobs, request.seasonal)) {
        events.push_back("remark default model selected (not enough obs.)");
        decision.default_model = true;
        apply_airline();
        return decision;
    }

    Ladder ladder(request);
    bool failed = false;
    try {
        ladder.y = work_series;
        ladder.spec = Orders{};
        ladder.spec.period = request.period;
        ladder.spec.d = request.initial_round ? 0 : request.d0;
        ladder.spec.bd = request.initial_round ? 0 : request.bd0;
        ladder.calc();
    } catch (const std::runtime_error&) {
        throw; // qualified-seam errors keep propagating
    } catch (const TramoLadderError&) {
        failed = true;
    }

    if (failed) {
        events.push_back("remark differencing selection failed");
        decision.selection_failed = true;
        apply_airline();
        return decision;
    }

    decision.mean_correction = ladder.mean_correction();
    decision.tmean = ladder.tmean;
    decision.rmax = ladder.rmax;
    decision.rsmax = ladder.rsmax;
    decision.c = ladder.c;
    decision.din = ladder.din;
    decision.iter = ladder.iter;
    decision.useml = ladder.useml;
    decision.mlused = ladder.mlused;
    decision.has_last_model = ladder.has_last;
    decision.last_model = ladder.last.parameters();
    decision.mspec_p = ladder.spec.p;
    decision.mspec_d = ladder.spec.d;
    decision.mspec_q = ladder.spec.q;
    decision.mspec_bp = ladder.spec.bp;
    decision.mspec_bd = ladder.spec.bd;
    decision.mspec_bq = ladder.spec.bq;

    // Description update: the module spec replaces the entry specification
    // when the differencing orders moved, and the mean follows
    // isMeanCorrection(); either change clears the estimation (the driver
    // seam) and marks Changed.
    bool changed = false;
    decision.spec_p = request.p0;
    decision.spec_d = request.d0;
    decision.spec_q = request.q0;
    decision.spec_bp = request.bp0;
    decision.spec_bd = request.bd0;
    decision.spec_bq = request.bq0;
    // Dated 2026-11-23 (javap DifferencingModule.process(RegSarimaModelling)
    // offsets 293-316 + fit-tr3_air witness): the Changed comparison is
    // against the ENTRY description orders, not the (initial-round zeroed)
    // ladder start - a ladder that climbs from 0,0 back onto an entry
    // 1,1/0,1 leaves the description untouched (Unchanged). The 23-C battery
    // always entered at 0,0, so both readings agreed there; the fit corpus
    // supplies the discriminating rows.
    if (ladder.spec.d != request.d0 || ladder.spec.bd != request.bd0) {
        changed = true;
        decision.spec_p = ladder.spec.p;
        decision.spec_d = ladder.spec.d;
        decision.spec_q = ladder.spec.q;
        decision.spec_bp = ladder.spec.bp;
        decision.spec_bd = ladder.spec.bd;
        decision.spec_bq = ladder.spec.bq;
    }
    if (decision.mean_correction != request.mean0) {
        changed = true;
        decision.mean = decision.mean_correction;
    }

    events.push_back("info diffInfo");
    events.push_back("pop");
    decision.result = changed ? TramoDifferencingDecision::Result::Changed
                              : TramoDifferencingDecision::Result::Unchanged;
    return decision;
}

} // namespace kalos
