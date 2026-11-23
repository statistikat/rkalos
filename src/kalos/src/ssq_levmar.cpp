#include "kalos/detail/sarima_stabilization.hpp"
#include "kalos/detail/filter_utility.hpp"
#include "kalos/detail/hostlineage.hpp"
#include "kalos/detail/ssq_levmar.hpp"

#include "kalos/linear_algebra.hpp" // polynomial_roots

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>

namespace kalos::detail {

// FilterUtility.checkRoots (declared in kalos/detail/filter_utility.hpp);
// visible to the stabilization helpers below and to the 21-H final-estimator unit.
bool check_roots(const double* c, std::size_t length, double limit) {
    std::size_t k = length;
    while (k > 0 && c[k - 1] == 0.0) --k;
    if (k == 0) return true;
    if (k == 1) return std::abs(c[0]) < limit;
    if (k == 2) {
        const double disc = c[0] * c[0] - 4.0 * c[1];
        if (disc > 0.0) {
            const double sq = std::sqrt(disc);
            const double r1 = (c[0] * -1.0 + sq) / (2.0 * c[1]);
            const double r2 = (c[0] * -1.0 - sq) / (2.0 * c[1]);
            const double i1 = 1.0 / std::abs(r1);
            const double i2 = 1.0 / std::abs(r2);
            return i1 < limit && i2 < limit;
        }
        return std::sqrt(c[1]) < limit;
    }
    // The pinned >= 3 path delegates to checkRoots(Complex[], limit), whose
    // guard returns true for a negative limit (and for a null root array,
    // which the deterministic port cannot produce).
    if (limit < 0.0) return true;
    std::vector<double> polynomial(k + 1, 0.0);
    polynomial[0] = 1.0;
    for (std::size_t i = 0; i < k; ++i) polynomial[i + 1] = c[i];
    for (const auto& root : polynomial_roots(polynomial)) {
        const double modulus = kalos::detail::glibc_hypot(root.real, root.imag);
        if (!(1.0 / modulus < limit)) return false;
    }
    return true;
}

namespace {

// Pinned SarimaMapping2.STEP and the of() default limits.
constexpr double kStep = 1.4901161028591027e-08; // Math.sqrt(2.220446E-16)
constexpr double kArLimit = 0.99999;
constexpr double kMaLimit = 1.0;
// Pinned LmBuilder defaults.
constexpr double kGradientPrecision = 1e-15;
constexpr double kParameterPrecision = 1e-15;
constexpr double kTau = 0.001;
constexpr int kMaxIter = 200;
constexpr double kAcceptance = 1e-4;
constexpr double kMuFloor = 0.3333333333333333;
// Pinned SsqNumericalDerivatives.g_nsteps.
constexpr int kEpsilonDeflations = 2;
// jdplus.toolkit.base.api.math.Constants, used by the pinned Reflector.
constexpr double kEpsilon = 1.1102230246251565e-16;   // Constants.getEpsilon()
constexpr double kSafeMin = 2.2250738585072014e-308;  // Constants.getSafeMin()

// Pinned DataPointer.norm2 (LAPACK dnrm2 scaling, order-sensitive).
double data_norm2(const double* p, std::size_t n) {
    if (n < 1) return 0.0;
    if (n == 1) return std::abs(p[0]);
    double scale = 0.0, ssq = 1.0;
    for (std::size_t t = 0; t < n; ++t) {
        const double xi = p[t];
        if (xi == 0.0) continue;
        const double a = std::abs(xi);
        if (std::isnan(scale) || scale < a) { // pinned dcmpg: NaN falls through
            const double r = scale / a;
            ssq = 1.0 + ssq * r * r;
            scale = a;
        } else {
            const double r = a / scale;
            ssq += r * r;
        }
    }
    return scale * std::sqrt(ssq);
}
} // namespace

// A faithful transcription of the pinned least-squares path of iterateQR:
// Householder2.decompose (Reflector.larfg/lapply on the shared storage) and
// QRLeastSquaresSolver.fastLeastSquares with UpperTriangularMatrix.rank(R, 0)
// (the count of nonzero |R(j,j)|, NaN counting as nonzero) followed by
// solveUx on the rank x rank block. The unpivoted reflector chain decides the
// split between the range and the near-null directions of a rank-deficient
// Jacobian differently from a pivoted LAPACK solver, and near ARMA
// cancellations that split steers the whole LM trajectory, so the pinned
// arithmetic is reproduced operation by operation (mirroring IEEE
// comparisons with the JVM dcmp/dcmpl/dcmpg branch tables where the pinned
// code uses them). singular mirrors the pinned m_err_sing throw, which the
// pinned iterateQR catch (Exception) feeds into the damping tail.
PinnedLsSolution pinned_least_squares(std::vector<double> qr, const std::vector<double>& e,
    std::size_t m, std::size_t dim) {
    PinnedLsSolution out;
    std::vector<double> beta(dim, 0.0);
    std::vector<double> b(m, 0.0);
    std::copy(e.begin(), e.end(), b.begin());
    const double s0 = kSafeMin / kEpsilon;

    // Householder2.householder: one Reflector per column on the shared array.
    for (std::size_t i = 0; i < dim; ++i) {
        const std::size_t xstart = i * (m + 1);
        const std::size_t n = m - i;
        double alpha = qr[xstart];
        double norm = data_norm2(&qr[xstart], n);
        int count = 0;
        if (norm != 0.0 && std::isless(norm, s0)) { // dcmpg ifge: NaN skips
            while (true) {
                const double s = 1.0 / s0;
                for (std::size_t t = 1; t < n; ++t) qr[xstart + t] *= s;
                alpha *= s;
                norm *= s;
                if (!std::isless(norm, s0)) { // dcmpg ifge: NaN recomputes
                    norm = data_norm2(&qr[xstart], n);
                    break;
                }
                ++count;
                if (count >= 4) {
                    norm = data_norm2(&qr[xstart], n);
                    break;
                }
            }
        }
        if (std::isless(alpha, 0.0)) norm = -norm; // dcmpg ifge: NaN keeps
        for (int t = 0; t < count; ++t) norm *= s0;
        double b_ = norm == 0.0 ? 0.0 : norm / (norm + alpha);
        for (std::size_t t = 1; t < n; ++t) qr[xstart + t] /= norm;
        qr[xstart] = 1.0 + alpha / norm;
        alpha = -norm;
        // Reflector.lapply on the trailing columns (window rows i..m-1,
        // columns i+1..dim-1); dcmpl ifeq: a NaN beta still applies.
        if (dim - i > 1 && b_ != 0.0) {
            for (std::size_t c = i + 1; c < dim; ++c) {
                const std::size_t pos = i + c * m;
                double w = qr[pos] / b_;
                for (std::size_t t = 1; t < n; ++t) w += qr[pos + t] * qr[xstart + t];
                if (w != 0.0) {
                    qr[pos] -= w;
                    w = -w * b_;
                    for (std::size_t t = 1; t < n; ++t) qr[pos + t] += w * qr[xstart + t];
                }
            }
        }
        beta[i] = b_;
        qr[xstart] = alpha; // R(j,j), the rawR() diagonal
    }

    // QRDecomposition.applyQt on the right-hand side [e; 0].
    for (std::size_t j = 0; j < dim; ++j) {
        if (beta[j] == 0.0) continue;
        const std::size_t xstart = j * (m + 1);
        double w = b[j] / beta[j];
        for (std::size_t t = j + 1; t < m; ++t) w += b[t] * qr[xstart + 1 + (t - j - 1)];
        if (w != 0.0) {
            b[j] -= w;
            w = -w * beta[j];
            for (std::size_t t = j + 1; t < m; ++t) {
                b[t] += w * qr[xstart + 1 + (t - j - 1)];
            }
        }
    }

    // UpperTriangularMatrix.rank(R, 0): the count of nonzero diagonal
    // entries; the fastRank predicate is dcmpg-based so a NaN counts.
    int rank = 0;
    for (std::size_t j = 0; j < dim; ++j) {
        const double d = qr[j * (m + 1)];
        if (!(std::abs(d) <= 0.0)) ++rank;
    }
    out.rank = rank;

    // solveUx on the rank x rank leading block; a zero pivot mirrors the
    // pinned m_err_sing throw (the |x| >= tol test is vacuous at tol = 0).
    for (std::size_t k = rank; k-- > 0;) {
        const double ui = qr[k * m + k];
        if (ui == 0.0) { // dcmpg ifgt: NaN is invertible, 0 is not
            out.singular = true;
            return out;
        }
        const double q = b[k] / ui;
        b[k] = q;
        for (std::size_t t = k; t-- > 0;) b[t] -= q * qr[k * m + t];
    }
    out.x.assign(dim, 0.0);
    std::copy_n(b.begin(), static_cast<std::size_t>(rank), out.x.begin());
    out.residual_tail.assign(b.begin() + static_cast<std::ptrdiff_t>(rank), b.end());
    return out;
}

namespace {

// Pinned FilterUtility.checkRoots(DoubleSeq, limit) for coefficient blocks of
// the polynomial 1 + c(B): trailing-zero trim, the |c| < limit case, the
// quadratic with its real and complex branch and the general root scan.

std::complex<double> pinned_inverse(const std::complex<double>& z) {
    // Pinned Complex.inv: plain IEEE component arithmetic, no guards; a zero
    // root produces the same infinities or NaNs as the pinned JVM.
    const double norm = z.real() * z.real() + z.imag() * z.imag();
    return {z.real() / norm, -z.imag() / norm};
}

// Pinned FilterUtility.stabilize(DoubleSeq.Mutable, limit): the single
// coefficient clamp or inversion and the root-inverting multi-coefficient
// reconstruction of stabilize(Polynomial, limit) followed by fromComplexRoots
// and the constant-term normalization. IEEE arithmetic only; like the pinned
// routine it never reports numerical dead ends.
bool stabilize_block(double* c, std::size_t length, double limit) {
    std::size_t k = length;
    while (k > 0 && c[k - 1] == 0.0) --k;
    if (k == 0) return false;
    if (k == 1) {
        const double magnitude = std::abs(c[0]);
        if (magnitude < limit) return false;
        if (magnitude > 1.0 / limit) {
            c[0] = 1.0 / c[0];
        } else {
            c[0] = c[0] > 0.0 ? limit - 1e-9 : -limit + 1e-9;
        }
        return true;
    }
    std::vector<double> polynomial(k + 1, 0.0);
    polynomial[0] = 1.0;
    for (std::size_t i = 0; i < k; ++i) polynomial[i + 1] = c[i];
    const auto found = polynomial_roots(polynomial);
    std::vector<std::complex<double>> roots;
    roots.reserve(found.size());
    for (const auto& root : found) roots.emplace_back(root.real, root.imag);
    bool changed = false;
    for (auto& root : roots) {
        const double modulus = std::abs(root);
        if (!(modulus < 1.0 / limit)) continue;
        if (modulus < 1.0) root = pinned_inverse(root);
        if (modulus > limit) root = root / (limit - 1e-9);
        changed = true;
    }
    if (!changed) return false;
    // Pinned Polynomial.fromComplexRoots (the z-roots recurrence); the caller
    // reads the non-constant coefficients through the constant-term view.
    std::vector<std::complex<double>> z(roots.size() + 1, {0.0, 0.0});
    z[0] = {1.0, 0.0};
    for (std::size_t i = 0; i < roots.size(); ++i) {
        z[i + 1] = z[i];
        for (std::size_t j = i; j >= 1; --j) z[j] = z[j - 1] - z[j] * roots[i];
        z[0] = z[0] * (-roots[i]);
    }
    for (std::size_t i = 0; i < k; ++i) c[i] = z[i + 1].real() / z[0].real();
    return true;
}

} // namespace

SarimaSsqDomain::SarimaSsqDomain(int p, int bp, int q, int bq, bool gls_mapping)
    : sizes_{static_cast<std::size_t>(p), static_cast<std::size_t>(bp),
          static_cast<std::size_t>(q), static_cast<std::size_t>(bq)},
      offsets_{0}, gls_mapping_(gls_mapping),
      limits_{gls_mapping ? 1.0 : kArLimit, gls_mapping ? 1.0 : kArLimit, kMaLimit, kMaLimit} {
    for (const std::size_t size : sizes_) offsets_.push_back(offsets_.back() + size);
}

double SarimaSsqDomain::epsilon(const std::vector<double>& x, std::size_t j) const {
    const double value = x[j];
    // Pinned dcmpg: a NaN coordinate takes the negative (second) branch.
    return value < 0.0 ? kStep * std::max(1.0, std::abs(value))
                       : -kStep * std::max(1.0, std::abs(value));
}

bool SarimaSsqDomain::check_boundaries(const std::vector<double>& x) const {
    for (std::size_t b = 0; b < 4; ++b) {
        if (sizes_[b] > 0 &&
            !check_roots(x.data() + offsets_[b], sizes_[b], limits_[b])) {
            return false;
        }
    }
    return true;
}

SsqParamValidation SarimaSsqDomain::validate(std::vector<double>& x) const {
    if (x.size() != dim()) return SsqParamValidation::Invalid;
    bool changed = false;
    for (std::size_t b = 0; b < 4; ++b) {
        if (sizes_[b] == 0) continue;
        if (gls_mapping_) {
            auto first = x.begin() + offsets_[b];
            if (check_roots(x.data() + offsets_[b], sizes_[b], 1.0)) continue;
            std::vector<double> block(first, first + sizes_[b]);
            if (stabilize_coefficient_block(block)) {
                std::copy(block.begin(), block.end(), first);
                changed = true;
            }
        } else if (stabilize_block(x.data() + offsets_[b], sizes_[b], limits_[b])) {
            changed = true;
        }
    }
    return changed ? SsqParamValidation::Changed : SsqParamValidation::Valid;
}

namespace {

struct Derivatives {
    std::vector<std::vector<double>> de; // dim columns of length nE
    std::vector<double> gradient;        // dim entries
    std::vector<bool> degenerate;        // epsp[j] == 0: the pinned NaN column
};

// Pinned SsqNumericalDerivatives.calcgrad with m_sym = m_mt = false: the
// per-coordinate checkepsilon deflation, then err(j, eps) in pinned order.
// The degenerate columns (eps == 0, so err returns the current vector and the
// division 0/0 produces NaN) are flagged: the pinned machine pushes NaN into
// the QR, loses the rank there and takes the rejection path.
Derivatives numerical_derivatives(const SsqVectorFunction& fn, const std::vector<double>& x,
    const std::vector<double>& ecur, const Interrupt& interrupt) {
    const SarimaSsqDomain& domain = fn.domain();
    const std::size_t dim = domain.dim();
    const std::size_t nE = ecur.size();
    Derivatives out;
    out.de.assign(dim, {});
    out.gradient.assign(dim, std::numeric_limits<double>::quiet_NaN());
    out.degenerate.assign(dim, false);
    std::vector<double> epsilon(dim, 0.0);
    std::vector<std::optional<SsqPoint>> plus(dim);
    for (std::size_t j = 0; j < dim; ++j) {
        if (interrupt) interrupt();
        epsilon[j] = domain.epsilon(x, j);
        if (epsilon[j] == 0.0) continue;
        // Pinned checkepsilon: keep the forward step when its point is
        // inside; else deflate by halves (at most g_nsteps), then negate,
        // deflate by halves again and give up with eps = 0.
        std::vector<double> probe = x;
        const double x0 = probe[j];
        const double original = epsilon[j];
        probe[j] += original;
        if (domain.check_boundaries(probe)) continue;
        double eps = original;
        bool inside = false;
        for (int steps = 1; steps <= kEpsilonDeflations; ++steps) {
            eps /= 2.0;
            probe[j] = x0 + eps;
            inside = domain.check_boundaries(probe);
            if (inside) break;
        }
        if (inside) {
            epsilon[j] = eps;
            continue;
        }
        eps = -original;
        probe[j] = x0 + eps;
        if (domain.check_boundaries(probe)) {
            epsilon[j] = eps;
            continue;
        }
        inside = false;
        for (int steps = 1; steps <= kEpsilonDeflations; ++steps) {
            eps /= 2.0;
            probe[j] = x0 + eps;
            inside = domain.check_boundaries(probe);
            if (inside) break;
        }
        epsilon[j] = inside ? eps : 0.0;
    }
    for (std::size_t j = 0; j < dim; ++j) {
        if (epsilon[j] == 0.0) {
            out.degenerate[j] = true; // grad[j] keeps its NaN: 0/0
            continue;
        }
        std::vector<double> probe = x;
        probe[j] += epsilon[j];
        plus[j] = fn.ssq_evaluate(probe); // pinned err: failures give ecur
        const std::vector<double>& ep =
            plus[j] && !plus[j]->e.empty() ? plus[j]->e : ecur;
        out.de[j].resize(nE);
        double acc = 0.0;
        bool finite = true;
        for (std::size_t i = 0; i < nE; ++i) {
            out.de[j][i] = (ep[i] - ecur[i]) / epsilon[j];
            acc += ecur[i] * out.de[j][i];
            finite = finite && std::isfinite(out.de[j][i]);
        }
        out.gradient[j] = 2.0 * acc;
        if (!finite) out.degenerate[j] = true;
    }
    return out;
}

struct LmState {
    const SsqVectorFunction* fn;
    std::vector<double> x;
    std::vector<double> e;
    double objective;
    double scale;
    double scale2;
    int stop = 0;
    int evaluations = 0;
    double mu = 0.0;
    long long nu = 4;
};

// Pinned iterateQR body; returns whether the outer loop continues.
bool iterate_qr(LmState& s, double eps3, const Interrupt& interrupt) {
    const SarimaSsqDomain& domain = s.fn->domain();
    const std::size_t dim = domain.dim();
    if (!std::isfinite(s.objective)) {
        s.stop = 7;
        return false;
    }
    if (!(s.objective > eps3 * s.scale2)) {
        s.stop = 6;
        return false;
    }
    const std::size_t nE = s.e.size();
    Derivatives derivs;
    try {
        derivs = numerical_derivatives(*s.fn, s.x, s.e, interrupt);
    } catch (const std::exception&) {
        if (interrupt) interrupt(); // rethrows an interruption unchanged
        return false;              // pinned catch (Exception) -> return false
    }
    // Pinned DoublesMath.normInf: seed abs(g0), then keep abs(gi) only when
    // strictly greater (dcmpl/iflt), so a NaN element never wins and a NaN
    // seed is never dislodged.
    double gnorm = 0.0;
    for (std::size_t j = 0; j < derivs.gradient.size(); ++j) {
        const double magnitude = std::abs(derivs.gradient[j]);
        if (j == 0 || magnitude > gnorm) gnorm = magnitude;
    }
    // Pinned dcmpg/ifgt: the loop continues when gnorm > eps1 * scale or the
    // norm is NaN (dcmpg pushes greater); a non-NaN value at or below the
    // bound stops. A NaN gradient surfaces later as stop 7 on the
    // non-finite least-squares right-hand side.
    if (gnorm <= kGradientPrecision * s.scale) {
        s.stop = 1;
        return false;
    }
    const bool degenerate = std::any_of(derivs.degenerate.begin(), derivs.degenerate.end(),
        [](const bool flag) { return flag; });
    for (int attempt = 0; attempt < 100; ++attempt) {
        if (interrupt) interrupt();
        // Stacked system [J; sqrt(mu) I] x = [e; 0]; step = -x.
        bool reject = degenerate;
        std::vector<double> step(dim, 0.0);
        if (!reject) {
            // Pinned stack: FastMatrix(nE + dim, dim) column-major, the top
            // nE x dim block carrying the Jacobian and sqrt(mu) on the
            // diagonal of the bottom block; the right-hand side is the e
            // vector padded with dim zeros.
            std::vector<double> a((nE + dim) * dim, 0.0);
            for (std::size_t j = 0; j < dim; ++j) {
                for (std::size_t i = 0; i < nE; ++i) {
                    a[j * (nE + dim) + i] = derivs.de[j][i];
                }
            }
            if (s.mu > 0.0) {
                const double root = std::sqrt(s.mu);
                for (std::size_t j = 0; j < dim; ++j) a[j * (nE + dim) + nE + j] = root;
            }
            std::vector<double> rhs(nE + dim, 0.0);
            std::copy(s.e.begin(), s.e.end(), rhs.begin());
            try {
                const auto solution = pinned_least_squares(std::move(a), rhs, nE + dim, dim);
                // Pinned: a zero pivot throws m_err_sing into the
                // catch (Exception) below and rank != dim branches to the
                // same damping tail.
                if (solution.singular || static_cast<std::size_t>(solution.rank) != dim) {
                    reject = true;
                } else {
                    for (std::size_t j = 0; j < dim; ++j) step[j] = -solution.x[j];
                }
            } catch (const std::exception&) {
                if (interrupt) interrupt();
                reject = true;
            }
        }
        if (!reject) {
            double stepSsq = 0.0;
            for (const double value : step) stepSsq += value * value;
            if (!std::isfinite(stepSsq)) {
                s.stop = 7;
                return false;
            }
            double xSsq = 0.0;
            for (const double value : s.x) xSsq += value * value;
            if (!(stepSsq > kParameterPrecision * (xSsq + kParameterPrecision))) {
                s.stop = 2;
                return false;
            }
            if (!(stepSsq < (xSsq + kParameterPrecision) / 1.0000000000000001E-34)) {
                s.stop = 4;
                return false;
            }
            std::vector<double> candidate = s.x;
            for (std::size_t j = 0; j < dim; ++j) candidate[j] += step[j];
            const auto validation = domain.validate(candidate);
            std::optional<SsqPoint> point;
            if (validation != SsqParamValidation::Invalid) {
                point = s.fn->ssq_evaluate(candidate); // pinned evaluate after validate
            }
            if (point) ++s.evaluations;
            const double change = point ? s.objective - point->ssq_e
                                        : std::numeric_limits<double>::quiet_NaN();
            if (point && change > 0.0) {
                if (validation == SsqParamValidation::Changed) {
                    s.x = std::move(candidate);
                    s.e = std::move(point->e);
                    s.objective = point->ssq_e;
                    return true;
                }
                double rho = std::numeric_limits<double>::quiet_NaN();
                {
                    double pred = 0.0;
                    for (std::size_t j = 0; j < dim; ++j) {
                        pred += (derivs.gradient[j] - 2.0 * s.mu * step[j]) * step[j];
                    }
                    pred = -pred / 2.0;
                    rho = change / pred;
                }
                if (rho > kAcceptance) {
                    const double u0 = 2.0 * rho - 1.0;
                    const double damping = 1.0 - u0 * u0 * u0;
                    if (s.mu != 0.0) {
                        s.mu *= damping < kMuFloor ? kMuFloor : damping;
                    }
                    s.nu = 4;
                    const bool done = !(change > eps3 * s.scale2);
                    s.x = std::move(candidate);
                    s.e = std::move(point->e);
                    s.objective = point->ssq_e;
                    if (done) {
                        s.stop = 2;
                        return false;
                    }
                    return true;
                }
            }
        }
        // Pinned rejection tail: the first failure arms mu = tau * scale,
        // later failures scale by nu and quadruple nu (nu << 2); a
        // non-growing nu stops.
        if (s.mu == 0.0) {
            s.mu = kTau * s.scale;
        } else {
            s.mu *= static_cast<double>(s.nu);
        }
        const long long quadrupled = s.nu << 2;
        if (!(quadrupled > s.nu)) {
            s.stop = 5;
            return false;
        }
        s.nu = quadrupled;
    }
    return false; // pinned 100 damping attempts exhausted
}

} // namespace

SsqLmResult minimize_ssq_levmar(const SsqVectorFunction& fn, std::vector<double> start,
    double fn_precision, const Interrupt& interrupt) {
    const double eps3 = fn_precision;
    auto point = fn.ssq_evaluate(start);
    if (!point) {
        throw std::runtime_error(
            "pinned sum-of-squares minimizer: the starting point is not evaluable");
    }
    LmState s;
    s.fn = &fn;
    s.x = std::move(start);
    s.e = std::move(point->e);
    s.objective = point->ssq_e;
    s.scale2 = s.objective;
    s.scale = std::sqrt(s.objective);
    s.evaluations = 1;
    // Pinned calc(): the loop compares the OLD iter against maxIter
    // (dup_x1/iadd/if_icmpge), so iterate() runs at most maxIter times and an
    // exhausted run leaves iter == maxIter + 1 (JVM: iter 201, 200 calls).
    int iterations = 0;
    for (;;) {
        const int iter_before = iterations;
        iterations = iter_before + 1;
        if (iter_before >= kMaxIter) break;
        if (!iterate_qr(s, eps3, interrupt)) break;
    }
    SsqLmResult result;
    result.parameters = std::move(s.x);
    result.e = std::move(s.e);
    result.objective = s.objective;
    result.iterations = iterations;
    result.evaluations = s.evaluations;
    result.stop = s.stop;
    result.converged = s.stop != 7 && s.stop != 4 && iterations < kMaxIter;
    result.has_converged = s.stop == 2;
    return result;
}

} // namespace kalos::detail
