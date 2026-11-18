#include "kalos/differencing_decision.hpp"

#include "kalos/detail/hannan_rissanen_pipeline.hpp"
#include "kalos/detail/sarima_stabilization.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

// Pinned FilterUtility.checkStability(DoubleSeq): Schur-type deflation of the
// Toeplitz system for the polynomial 1 + b(B). The recursion is sign-
// equivariant under negation of the block, so it may run on either sign
// convention. Pinned NaN behavior: |NaN| >= 1 is false and the recursion
// continues with NaN coefficients.
bool check_stability(const std::vector<double>& block) {
    const std::size_t m = block.size();
    if (m == 0) return true;
    if (m == 1) return std::abs(block[0]) < 1.0;
    std::vector<double> arr = block;
    std::vector<double> c(m, 0.0);
    std::vector<double> d(m, 0.0);
    for (std::size_t i = m; i-- > 0;) {
        c[i] = arr[i];
        // Pinned dcmpl/iflt: a NaN coefficient is not >= 1 and the pinned
        // recursion continues with it.
        if (std::abs(c[i]) >= 1.0) return false;
        for (std::size_t j = 0; j < i; ++j) d[j] = arr[i - 1 - j];
        const double den = 1.0 - c[i] * c[i];
        for (std::size_t k = 0; k < i / 2; ++k) {
            arr[k] = (arr[k] - c[i] * d[k]) / den;
            arr[i - 1 - k] = (arr[i - 1 - k] - c[i] * d[i - 1 - k]) / den;
        }
        if (i % 2 != 0) arr[i / 2] = d[i / 2] / (1.0 + c[i]);
    }
    return true;
}

// Pinned SarimaModel.isStable(invertible): checkStability over the regular
// AR, the seasonal AR and, when invertible is set, both MA blocks.
struct JvmModel { // Pinned SarimaModel: blocks of 1 + phi1(B), 1 + th(B).
    std::vector<double> phi, bphi, th, bth;
    double phi_at(std::size_t i) const { return phi.at(i); }
};

bool model_stable(const JvmModel& m, bool invertible) {
    if (!check_stability(m.phi)) return false;
    if (!check_stability(m.bphi)) return false;
    if (invertible && !check_stability(m.th)) return false;
    if (invertible && !check_stability(m.bth)) return false;
    return true;
}

// Pinned SarimaMapping.stabilize(true, orders, block) over the model
// parameter layout [phi, bphi, th, bth] with whole seasonal blocks (the
// differencing call site always has bp = bq = 1 here). Returns whether any
// block was rewritten.
bool stabilize_jvm_model(JvmModel& m) {
    bool changed = detail::stabilize_coefficient_block(m.phi);
    changed |= detail::stabilize_coefficient_block(m.bphi);
    changed |= detail::stabilize_coefficient_block(m.th);
    changed |= detail::stabilize_coefficient_block(m.bth);
    return changed;
}

// Pinned ComplexType.abs(re, im): scale-aware Euclidean norm.
double complex_abs(double re, double im) {
    if (re == 0.0 && im == 0.0) return 0.0;
    const double a = std::abs(re);
    const double b = std::abs(im);
    const double hi = std::max(a, b);
    const double lo = std::min(a, b);
    if (lo == 0.0) return hi;
    const double ratio = lo / hi;
    return hi * std::sqrt(1.0 + ratio * ratio);
}

// Roots of Polynomial.mirror() of the AR polynomial 1 + a(B). Pinned chain:
// getRegularAR()/getSeasonalAR() build [1, a...] with Polynomial.ofInternal
// (NO zero trimming) and mirror() reverses to [a..., 1] with a trailing one
// that Coefficients.ofInternal can never trim, so the mirrored degree is
// exactly the block size. Pinned FastEigenValuesSolver: degree 1 is monic,
// degree 2 the textbook quadratic with the complex pair stored as
// conjugates. Degrees above 2 would enter the pinned QR machine; the
// differencing call site never exceeds degree 2 (p <= 2, bp <= 1).
std::vector<std::pair<double, double>> mirrored_ar_roots(const std::vector<double>& block) {
    std::vector<double> poly; // Ascending coefficients of the mirrored polynomial.
    poly.reserve(block.size() + 1);
    for (std::size_t i = block.size(); i-- > 0;) poly.push_back(block[i]);
    poly.push_back(1.0);
    if (poly.size() == 1) return {};                 // Pinned: constant, no roots.
    if (poly.size() == 2) return {{-poly[0] / poly[1], 0.0}};
    if (poly.size() > 3) throw std::logic_error("differencing decision: mirrored AR degree above 2 is unreachable");
    const double p2 = poly[2], p1 = poly[1], p0 = poly[0];
    const double d8 = 2.0 * p2;
    const double disc = p1 * p1 - 4.0 * p2 * p0;
    if (disc < 0.0) {
        const double im = std::sqrt(-disc);
        return {{-p1 / d8, im / d8}, {-p1 / d8, -(im / d8)}};
    }
    const double root = std::sqrt(disc);
    return {{(-p1 + root) / d8, 0.0}, {(-p1 - root) / d8, 0.0}};
}

// Pinned UnitRoots.D products: (1 - B)^d * (1 - B^period)^bd as an
// ascending-coefficient polynomial and its degree.
std::vector<double> differencing_polynomial(int period, int d, int bd) {
    return detail::pinned_differencing_polynomial(period, d, bd);
}

// Pinned BackFilter.apply (shared transcription in detail::
// pinned_backfilter_apply): the filter acts through the reversed weight
// array, accumulating over k ascending and skipping zero weights. A filter
// longer than the series reproduces the pinned negative-size array failure.
std::vector<double> filter_apply(const std::vector<double>& f, const std::vector<double>& x) {
    std::vector<double> out = detail::pinned_backfilter_apply(f, x);
    if (out.empty()) {
        throw std::runtime_error(
            "differencing decision: the differencing filter is not shorter than the series");
    }
    return out;
}

// Pinned DataBlock.average() on a finite block: naive left-to-right sum / n.
double block_mean(const std::vector<double>& x) {
    double sum = 0.0;
    for (const double value : x) sum += value;
    return sum / static_cast<double>(x.size());
}

bool all_finite(const std::array<std::vector<double>, 4>& blocks) {
    for (const auto& block : blocks) {
        for (const double value : block) {
            if (!std::isfinite(value)) return false;
        }
    }
    return true;
}

struct Spec { // Pinned SarimaOrders fields the machine reads and writes.
    int p = 0, q = 0, bp = 0, bq = 0, d = 0, bd = 0;
};

class Machine {
public:
    Machine(const DifferencingDecisionRequest& request, const Interrupt& interrupt)
        : x_(request.series), period_(request.period), maxd_(request.max_d),
          maxbd_(request.max_seasonal_d), ub1_(request.ub1), ub2_(request.ub2),
          c_(request.cancel), cancel_(request.cancel), eps_(request.precision),
          interrupt_(interrupt) {}

    DifferencingDecision run() {
        useml_ = false;
        mlused_ = false;
        rmax_ = rsmax_ = tmean_ = 0.0;
        mu_ = false;
        DifferencingDecision decision;
        step0(decision);
        iter_ = 0;
        int nextstep_calls = 0;
        while (true) {
            const bool more = nextstep(decision);
            ++nextstep_calls;
            if (!more) break;
            if (iter_ >= 5) break;
            ++iter_;
        }
        calc_mean();
        decision.d = spec_.d;
        decision.seasonal_d = spec_.bd;
        decision.mean_correction = mu_;
        decision.t_mean = tmean_;
        decision.ml_used = mlused_;
        decision.nextstep_calls = nextstep_calls;
        return decision;
    }

private:
    bool has_seas() const { return maxbd_ > 0 && period_ > 1; }

    double phi1() const { return last_.phi.at(0); }
    double theta1() const { return last_.th.at(0); }
    double bphi1() const { return last_.bphi.at(0); }
    double btheta1() const { return last_.bth.at(0); }

    // The two estimator seams work in the Kalos convention [ar, ma, sar,
    // sma] with AR 1 - phi1(B); the machine lives in the JVM convention.
    std::array<std::vector<double>, 4> to_kalos() const {
        auto negated = [](const std::vector<double>& v) {
            std::vector<double> out(v.size());
            for (std::size_t i = 0; i < v.size(); ++i) out[i] = -v[i];
            return out;
        };
        return {negated(last_.phi), last_.th, negated(last_.bphi), last_.bth};
    }

    void from_kalos(const std::array<std::vector<double>, 4>& c) {
        auto negated = [](const std::vector<double>& v) {
            std::vector<double> out(v.size());
            for (std::size_t i = 0; i < v.size(); ++i) out[i] = -v[i];
            return out;
        };
        last_.phi = negated(c[0]);
        last_.th = c[1];
        last_.bphi = negated(c[2]);
        last_.bth = c[3];
    }

    void fill_model_sizes() {
        last_.phi.assign(static_cast<std::size_t>(spec_.p), 0.0);
        last_.bphi.assign(static_cast<std::size_t>(spec_.bp), 0.0);
        last_.th.assign(static_cast<std::size_t>(spec_.q), 0.0);
        last_.bth.assign(static_cast<std::size_t>(spec_.bq), 0.0);
    }

    void default_fill_model() {
        // Pinned SarimaModel.Builder.setDefault() in the JVM convention:
        // phi = bphi = -0.1, th = bth = -0.2.
        std::fill(last_.phi.begin(), last_.phi.end(), -0.1);
        std::fill(last_.bphi.begin(), last_.bphi.end(), -0.1);
        std::fill(last_.th.begin(), last_.th.end(), -0.2);
        std::fill(last_.bth.begin(), last_.bth.end(), -0.2);
    }

    void initstep(bool arfirst, DifferencingIterationRecord& record) {
        if (spec_.d == 0 && spec_.bd == 0 && arfirst) {
            spec_.p = period_ == 2 ? 1 : 2;
            spec_.q = 0;
            spec_.bq = 0;
            if (has_seas()) spec_.bp = 1;
        } else {
            spec_.p = 1;
            spec_.q = 1;
            if (period_ > 1) {
                spec_.bp = 1;
                spec_.bq = 1;
            }
        }
        record.p = spec_.p;
        record.q = spec_.q;
        record.bp = spec_.bp;
        record.bq = spec_.bq;

        const auto filter = differencing_polynomial(period_, spec_.d, spec_.bd);
        // The raw filtered series is what the pinned RegArmaModel.differencedModel()
        // carries (y = unit-root filter applied, no centering); the HannanRissanen
        // pipeline consumes the centered copy below.
        std::vector<double> z_filtered = filter_apply(filter, x_);
        std::vector<double> z = z_filtered;
        const double mean = block_mean(z);
        for (double& value : z) value -= mean;

        fill_model_sizes();
        bool failure = false;
        const std::array<std::vector<double>, 4> raw =
            detail::hannan_rissanen_pipeline(z, {spec_.p, spec_.bp, spec_.q, spec_.bq, period_},
                interrupt_);
        record.hr_source = "fitted";
        last_.phi = raw[0];
        last_.bphi = raw[1];
        last_.th = raw[2];
        last_.bth = raw[3];
        if (arfirst && !model_stable(last_, true)) {
            // Pinned dcmpl/ifgt: a NaN coefficient takes the failure branch,
            // so the gates transcribe as the NaN-true !(<=) form.
            if (spec_.p > 1 || (spec_.p == 1 && !(std::abs(phi1()) <= 1.02))) {
                failure = true;
            } else if (spec_.bp == 1 && !(std::abs(bphi1()) <= 1.02)) {
                failure = true;
            } else {
                record.hr_stabilized = stabilize_jvm_model(last_);
            }
        }
        record.failure = failure;
        if (failure) {
            fill_model_sizes();
            default_fill_model();
            record.hr_source = "defaultfill";
        }
        if (failure || useml_) { // Pinned module-level ml is false.
            stabilize_jvm_model(last_);
            record.trial_fit = true;
            SarimaEstimationRequest trial;
            // Pinned call site: RegArimaModel(y = x, no regressors, no mean
            // correction) goes through IRegArimaComputer.optimize, which works
            // on RegArmaModel.differencedModel(): the unit-root-filtered
            // response, the stationary part of the model and no regressors
            // (unlike the nstMean seam below, this machine carries no mean
            // column; corpus-qualified through the pinned final-parameter
            // records, c[20]).
            trial.series.values = std::move(z_filtered);
            trial.series.missing.assign(trial.series.values.size(), 0);
            trial.design.columns = 0;
            trial.period = period_;
            trial.d = 0;
            trial.seasonal_d = 0;
            const std::array<std::vector<double>, 4> start = to_kalos();
            trial.blocks[0].initial = start[0];
            trial.blocks[1].initial = start[1];
            trial.blocks[2].initial = start[2];
            trial.blocks[3].initial = start[3];
            for (auto& block : trial.blocks) block.fixed.assign(block.initial.size(), 0);
            trial.control = EstimationControl{eps_, 1e-4, true, false, true, true,
                EstimationStartingPoint::Default, EstimationSolver::PinnedLevenbergMarquardt};
            const SarimaEstimationResult est = optimize_sarima_from(trial, start, interrupt_);
            record.trial_status = est.status;
            if (!all_finite(est.coefficients)) {
                throw std::runtime_error("differencing decision: non convergence in IDDIF");
            }
            from_kalos(est.coefficients);
            mlused_ = true;
        } else {
            mlused_ = false;
        }
        useml_ = false;
    }

    void step0(DifferencingDecision& decision) {
        DifferencingIterationRecord record;
        record.iteration = 0;
        record.cancel_threshold = c_;
        initstep(true, record);
        if (spec_.d != 0 || spec_.bd != 0) {
            rmax_ = phi1();
            if (has_seas()) rsmax_ = bphi1();
        }
        const auto regular_roots = mirrored_ar_roots(last_.phi);
        spec_.d += searchur(regular_roots, ub1_, true);
        if (has_seas()) {
            const auto seasonal_roots = mirrored_ar_roots(last_.bphi);
            spec_.bd += searchur(seasonal_roots, ub1_, false);
        }
        record.rmax = rmax_;
        record.rsmax = rsmax_;
        record.useml = useml_;
        record.mlused = mlused_;
        record.phi = last_.phi.empty() ? not_a_number : phi1();
        record.theta = last_.th.empty() ? not_a_number : theta1();
        record.bphi = last_.bphi.empty() ? not_a_number : bphi1();
        record.btheta = last_.bth.empty() ? not_a_number : btheta1();
        record.d = spec_.d;
        record.seasonal_d = spec_.bd;
        decision.trace.push_back(record);
    }

    int searchur(const std::vector<std::pair<double, double>>& roots, double ub, bool regular) {
        int count = 0;
        double maxmod = 0.0;
        for (const auto& root : roots) {
            const double im = std::abs(root.second);
            const double mod = complex_abs(root.first, root.second);
            if (mod >= ub && im <= 0.05 && root.first > 0.0) {
                ++count;
            } else if (im <= 0.02 && root.first > 0.0 && mod > maxmod) {
                maxmod = mod;
            }
        }
        if (regular) rmax_ = maxmod;
        else rsmax_ = maxmod;
        return count;
    }

    int maincondition() {
        const double phi = last_.phi.empty() ? 0.0 : phi1();
        const double theta = last_.th.empty() ? 0.0 : theta1();
        const double bphi = !has_seas() || last_.bphi.empty() ? 0.0 : bphi1();
        const double btheta = !has_seas() || last_.bth.empty() ? 0.0 : btheta1();
        c_ -= 0.002;
        double uc;
        // Pinned dcmpl/iflt: a NaN ub2 takes the 1.005 branch, so the
        // comparison transcribes as the NaN-false form.
        if (!mlused_ && ub2_ >= 0.869) uc = 0.136;
        else uc = 1.005 - ub2_;
        int res = 0;
        // Pinned dcmpg/ifgt: the inside test is <= (NaN goes to the else branch).
        if (std::abs(phi + 1.0) <= uc) {
            if (-phi > 1.02) {
                res = 1;
                useml_ = true;
            } else if (std::abs(phi - theta) > c_) {
                res = 1;
                ++spec_.d;
            }
        } else if (std::abs(phi) > 1.12) {
            res = 1;
            useml_ = true;
        }
        if (has_seas()) {
            // Pinned dcmpg/ifgt: inside is <= 0.19 (NaN goes to the else branch).
            if (std::abs(bphi + 1.0) <= 0.19) {
                if (-bphi > 1.02) {
                    useml_ = true;
                    res = 1;
                } else if (spec_.bd == 0 && std::abs(bphi - btheta) > c_) {
                    ++res;
                    ++spec_.bd;
                    if (useml_) {
                        --res;
                        useml_ = false;
                    }
                }
            } else if (std::abs(bphi) > 1.12) {
                res = 1;
                useml_ = true;
            }
        }
        return res;
    }

    int cond1(int n) {
        if (spec_.d + spec_.bd != 0) return n;
        const double phi = phi1();
        const double theta = theta1();
        const double bphi = has_seas() ? bphi1() : 0.0;
        const double btheta = has_seas() ? btheta1() : 0.0;
        const auto increment = [&]() {
            if (useml_ && n == 1) useml_ = false;
            else ++n;
            if (rmax_ > rsmax_) ++spec_.d;
            else ++spec_.bd;
            return n;
        };
        // Pinned fallthrough: the cancellation branch continues into the
        // plain tests when neither rmax nor rsmax reaches 0.9. The
        // cancellation entry and the near tests are dcmpg-based (NaN-false
        // <= and <); the far tests route NaN to the far side (!(<=)).
        if (std::abs(phi - theta) < c_ || (has_seas() && std::abs(bphi - btheta) < c_)) {
            if (rmax_ >= 0.9 || rsmax_ >= 0.9) return increment();
        }
        if (std::abs(phi + 1.0) <= 0.15 || (has_seas() && std::abs(bphi + 1.0) <= 0.16)) {
            if (rmax_ >= 0.9 || rsmax_ >= 0.88) return increment();
        }
        if (!(std::abs(phi + 1.0) <= 0.16) && (!has_seas() || !(std::abs(bphi + 1.0) <= 0.17))) {
            return n;
        }
        if (rmax_ >= 0.91 || rsmax_ >= 0.89) return increment();
        return n;
    }

    void allcond() {
        if (spec_.d + spec_.bd != 0) return;
        const double phi = phi1();
        // Pinned dcmpg directions: the plain test is <= (NaN skips), the
        // seasonal guard enters on the NaN-false <= pair.
        if (!has_seas()) {
            if (std::abs(phi + 1.0) <= 0.15) ++spec_.d;
            return;
        }
        const double bphi = bphi1();
        if (!(std::abs(phi + 1.0) <= 0.15) && !(std::abs(bphi + 1.0) <= 0.25)) return;
        rmax_ = -phi;
        rsmax_ = -bphi;
        if (phi < bphi) ++spec_.d;
        else ++spec_.bd;
    }

    int finalcond(int n) {
        if (n == 2) {
            --spec_.d;
            --spec_.bd;
            if (mlused_) {
                if (phi1() < bphi1()) ++spec_.d;
                else ++spec_.bd;
            } else if (rmax_ > rsmax_) {
                if (rmax_ > 0.0) ++spec_.d;
            } else if (rsmax_ > 0.0) {
                ++spec_.bd;
            }
        }
        if (spec_.d > maxd_) {
            spec_.d = maxd_;
            n = 0;
        }
        if (spec_.bd > maxbd_) {
            spec_.bd = maxbd_;
            n = 0;
        }
        return n;
    }

    bool nextstep(DifferencingDecision& decision) {
        DifferencingIterationRecord record;
        record.iteration = iter_ + 1;
        record.cancel_threshold = c_;
        initstep(false, record);
        // maincondition reads the parameters initstep produced; snapshot
        // them before the tests can change the orders.
        record.phi = phi1();
        record.theta = theta1();
        record.bphi = has_seas() ? bphi1() : 0.0;
        record.btheta = has_seas() ? btheta1() : 0.0;
        const int main = maincondition();
        record.maincondition = main;
        record.cancel_threshold = c_; // Post-decrement c_ as the tests used it.
        int res = main;
        if (iter_ == 0) res = cond1(res);
        record.carried = res;
        allcond();
        const int settled = finalcond(res);
        record.rmax = rmax_;
        record.rsmax = rsmax_;
        record.useml = useml_;
        record.mlused = mlused_;
        record.d = spec_.d;
        record.seasonal_d = spec_.bd;
        decision.trace.push_back(record);
        return settled != 0;
    }

    void calc_mean() {
        const int order = spec_.d + (period_ > 1 && spec_.bd > 0 ? period_ * spec_.bd : 0);
        if (order == 0) { // Pinned stMean().
            double s = 0.0, q = 0.0;
            for (const double value : x_) {
                s += value;
                q += value * value;
            }
            tmean_ = s / std::sqrt(q);
            const double th = x_.size() > 200 ? 2.55 : (x_.size() > 80 ? 2.0 : 1.96);
            mu_ = std::abs(tmean_) > th;
            return;
        }
        // Pinned nstMean(): X13Utility.processor(true, eps) == RegSarimaComputer
        // {ml = true, precision = eps, startingPoint = Multiple, feps = 1e-4,
        // mt = false, cdf = true, fast = true} over the defaulted model with a
        // mean correction; tstat of the concentrated likelihood's first
        // coefficient. The pinned meanCorrection (RegArimaModel.isMean() with
        // xCount == 0) is an ARIMA-level Griffies mean, so every pinned
        // computation happens on RegArmaModel.differencedModel(): the
        // unit-root-filtered response, the stationary model and one constant
        // regressor column (differencedModel bytecode: isMean writes the
        // first column with 1.0). Passing the filtered response with zero
        // differencing and a constant design is the identical machine; the
        // mean column is written fresh on the differenced domain (raw ones
        // differenced would be zeros), which is why the spec differencing is
        // carried in spec_d/spec_bd rather than on d/seasonal_d. A drift
        // column instead would change the optimizer's search surface (and is
        // unidentifiable at d = 2).
        const std::vector<double> filter =
            differencing_polynomial(period_, spec_.d, spec_.bd);
        const std::vector<double> filtered = filter_apply(filter, x_);
        SarimaEstimationRequest test;
        test.series.values = filtered;
        test.series.missing.assign(test.series.values.size(), 0);
        test.design.columns = 1;
        test.design.values.assign(test.series.values.size(), 1.0);
        test.design.coefficients.assign(1, 0.0);
        test.design.fixed.assign(1, 0);
        test.blocks[0].initial.assign(static_cast<std::size_t>(spec_.p), 0.1);
        test.blocks[0].fixed.assign(static_cast<std::size_t>(spec_.p), 0);
        test.blocks[1].initial.assign(static_cast<std::size_t>(spec_.q), -0.2);
        test.blocks[1].fixed.assign(static_cast<std::size_t>(spec_.q), 0);
        test.blocks[2].initial.assign(static_cast<std::size_t>(spec_.bp), 0.1);
        test.blocks[2].fixed.assign(static_cast<std::size_t>(spec_.bp), 0);
        test.blocks[3].initial.assign(static_cast<std::size_t>(spec_.bq), -0.2);
        test.blocks[3].fixed.assign(static_cast<std::size_t>(spec_.bq), 0);
        test.period = period_;
        test.d = 0;
        test.seasonal_d = 0;
        test.spec_d = spec_.d;
        test.spec_bd = spec_.bd;
        test.control = EstimationControl{eps_, 1e-4, true, false, true, true,
            EstimationStartingPoint::Multiple, EstimationSolver::PinnedLevenbergMarquardt};
        const SarimaEstimationResult est = estimate_sarima_fixed(test, interrupt_);
        if (!all_finite(est.coefficients)) return; // Pinned null-estimation return.
        // The pinned tstat comes from the same machine that produced the
        // estimation: the ConcentratedLikelihoodWithMissing of the landed
        // coordinates, coefficient 0 over its standard error (21-D.2).
        {
            std::vector<double> jvm;
            for (const double v : est.coefficients[0]) jvm.push_back(-v);
            for (const double v : est.coefficients[2]) jvm.push_back(-v);
            for (const double v : est.coefficients[1]) jvm.push_back(v);
            for (const double v : est.coefficients[3]) jvm.push_back(v);
            const auto point = detail::pinned_deviances_evaluate(filtered,
                static_cast<std::size_t>(spec_.p), static_cast<std::size_t>(spec_.bp),
                static_cast<std::size_t>(spec_.q), static_cast<std::size_t>(spec_.bq), period_,
                std::vector<double>(filtered.size(), 1.0), jvm);
            if (!point) return; // Pinned machine failure: the nstMean catch
                                // swallows it and the mean stays untested.
            if (point->rank == 0) {
                // Pinned IndexOutOfBounds on the coefficient-less concentrated
                // likelihood of a rank-deficient design.
                throw std::runtime_error(
                    "differencing decision: the mean regressor is not identifiable after "
                    "differencing");
            }
            // ConcentratedLikelihood.tstat mirrors the uc == 0 NaN and the
            // zero-coefficient shortcuts, so the division never sees a zero
            // standard error.
            tmean_ = point->tstat;
        }
        const std::size_t n = x_.size();
        const double th = n <= 80 ? 1.96
                          : n <= 155 ? 1.98
                          : n <= 230 ? 2.1
                          : n <= 350 ? 2.3
                                     : 2.5;
        mu_ = std::abs(tmean_) > th;
    }

    std::vector<double> x_;
    int period_;
    int maxd_, maxbd_;
    double ub1_, ub2_;
    double c_, cancel_, eps_;
    const Interrupt& interrupt_;
    Spec spec_;
    JvmModel last_;
    double rmax_ = 0.0, rsmax_ = 0.0;
    double tmean_ = 0.0;
    bool mu_ = false;
    bool useml_ = false, mlused_ = false;
    int iter_ = 0;
};

} // namespace

void DifferencingDecisionRequest::validate() const {
    if (series.empty()) throw std::invalid_argument("differencing decision: empty series");
    for (const double value : series) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("differencing decision: series must be finite throughout");
        }
    }
    // Qualified domain: every decision path (trial fit, stMean, nstMean)
    // enters the SARIMA estimation seams, which are qualified for periods 4
    // and 12 only. The pinned JVM machine itself accepts any period; the
    // restriction is the Kalos estimator domain, documented in the contract.
    if (period != 4 && period != 12) {
        throw std::invalid_argument("differencing decision: period must be 4 or 12");
    }
    if (max_d < 0 || max_seasonal_d < 0) {
        throw std::invalid_argument("differencing decision: differencing clamps must be non-negative");
    }
    for (const double value : {ub1, ub2, cancel}) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("differencing decision: thresholds must be finite");
        }
    }
    if (!std::isfinite(precision) || precision <= 0.0) {
        throw std::invalid_argument("differencing decision: precision must be finite and positive");
    }
}

DifferencingDecision decide_differencing(const DifferencingDecisionRequest& request,
    const Interrupt& interrupt) {
    request.validate();
    Machine machine(request, interrupt);
    return machine.run();
}

} // namespace kalos
