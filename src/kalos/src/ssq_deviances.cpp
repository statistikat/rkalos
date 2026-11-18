// Pinned jdplus 3.9.0 deviances machine for the differencing-decision
// estimation seams (item 21-D.2). See include/kalos/detail/ssq_deviances.hpp
// for the scope; every stage below mirrors one pinned class:
//
//   Polynomial / Polynomials / Coefficients.getUsedCoefficients  (poly_*)
//   RationalFunction.prepare                                     (rf_coefficients)
//   AutoCovarianceFunction + defaultComputer                     (acf_values)
//   LULinearSystemSolver2                                        (lu_solve)
//   DeterminantalTerm                                            (det_term)
//   internal KalmanFilter (prepare/tx/calcC/sfilter/mfilter/tlast)(kalman)
//   HouseholderWithPivoting + Reflector.larfg (single column)     (hwp_1col)
//   QRDecomposition.applyQt + QRLeastSquaresSolver.leastSquares   (least_squares)
//   ConcentratedLikelihoodComputer.process + CLWM builder/ctor    (evaluate)
//
// Comparison flavours follow the pinned bytecode: a `dcmpl`-guarded test
// treats NaN as "less", a `dcmpg`-guarded test treats NaN as "greater" and
// `dcmp`-guarded tests treat NaN as unordered-equal; the transcriptions below
// reproduce each branch table explicitly.

#include "kalos/detail/ssq_deviances.hpp"

#include "kalos/detail/fdlibm.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

namespace kalos::detail {
namespace {

// jdplus.toolkit.base.api.math.Constants.
constexpr double kEpsilon = 1.1102230246251565e-16;   // getEpsilon()
constexpr double kSafeMin = 2.2250738585072014e-308;  // getSafeMin()
constexpr double kMinValue = 4.9406564584124654e-324; // Double.MIN_VALUE
constexpr double kLogTwoPi = 1.8378770664093453;      // Constants.LOGTWOPI
constexpr double kLog2 = 0.6931471805599453;          // DeterminantalTerm.LOG2
// Coefficients.getUsedCoefficients tolerance behind Polynomial.ofInternal.
constexpr double kTrimTol = 1e-9;
// The steady-state scan bound in KalmanFilter.calcC/sfilter.
constexpr double kSteadyTol = 1e-15;
// LULinearSystemSolver2(1e-13) singularity bound and the
// ConcentratedLikelihoodComputer rcond handed to QRLeastSquaresSolver.
constexpr double kLUEps = 1e-13;
constexpr double kRcond = 1e-12;

// Polynomial.ofInternal: Coefficients.getUsedCoefficients(c, 1e-9) strips
// trailing entries while !(|c| > 1e-9); the dcmpg guard treats a NaN trailing
// coefficient as stripped.
void trim_used(std::vector<double>& c) {
    while (!c.empty() && !(std::abs(c.back()) > kTrimTol)) {
        c.pop_back();
    }
}

double poly_get(const std::vector<double>& c, std::size_t i) {
    return i < c.size() ? c[i] : 0.0;
}

// Polynomial.isIdentity: a single 1.0 coefficient; Polynomial.isZero: a
// single 0.0 coefficient.
bool poly_is_one(const std::vector<double>& c) {
    return c.size() == 1 && c[0] == 1.0; // dcmpl ifne: NaN is not the identity
}
bool poly_is_zero(const std::vector<double>& c) {
    return c.size() == 1 && c[0] == 0.0;
}

// internal.toolkit...Polynomials.times: skip zero operands, out[i + j] +=
// a[i] * b[j] with i ascending outside; the zero-element skips (dcmpl ifeq)
// keep untouched products out of the accumulator, preserving signed zeros.
std::vector<double> poly_times(const std::vector<double>& a, const std::vector<double>& b) {
    if (poly_is_zero(a) || poly_is_zero(b)) {
        return {0.0};
    }
    if (poly_is_one(b)) {
        return a;
    }
    if (poly_is_one(a)) {
        return b;
    }
    std::vector<double> out(a.size() + b.size() - 1, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0.0) {
            continue;
        }
        for (std::size_t j = 0; j < b.size(); ++j) {
            if (b[j] == 0.0) {
                continue;
            }
            out[i + j] += a[i] * b[j];
        }
    }
    return out;
}

// The pinned JVM evaluates the estimation machine with the default floating
// point environment: flush-to-zero and denormals-are-zero are clear. BLAS
// back ends leave those MXCSR bits set (and the C library's exp and log set
// flush-to-zero themselves), which changes results by an ulp whenever a
// denormal appears in an intermediate scaling sequence. The machine therefore
// restores the default environment for the duration of an evaluation and puts
// the caller's state back afterwards. The instruction encoding matches the
// 32-bit MXCSR operand of ldmxcsr and stmxcsr; x86-64 always assembles it.
namespace {

#if defined(__x86_64__) || defined(__i386__)
inline unsigned mxcsr_get() {
    unsigned value = 0;
    asm volatile("stmxcsr %0" : "=m"(value));
    return value;
}
inline void mxcsr_set(unsigned value) { asm volatile("ldmxcsr %0" :: "m"(value)); }

class default_fp_environment {
    unsigned saved_ = 0;

public:
    default_fp_environment() {
        saved_ = mxcsr_get();
        // Bit 6 (0x40, flush-to-zero) and bit 1 (0x2, denormals-are-zero)
        // cleared; every other control bit is left as the caller set it.
        mxcsr_set(saved_ & ~0x42u);
    }
    ~default_fp_environment() { mxcsr_set(saved_); }
    default_fp_environment(const default_fp_environment&) = delete;
    default_fp_environment& operator=(const default_fp_environment&) = delete;
};
#else
// On other architectures the denormal mode is not observable through MXCSR;
// the JVM semantics hold to the extent the platform provides them.
class default_fp_environment {};
#endif

} // namespace

// Polynomial.factor(coef, s): [1, -coef * B^s]; coef == 0 gives ONE. The
// result goes through ofInternal (trailing trim).
std::vector<double> poly_factor(double coef, int s) {
    if (coef == 0.0) { // dcmpl ifne: NaN factors
        return {1.0};
    }
    std::vector<double> c(static_cast<std::size_t>(s) + 1, 0.0);
    c[0] = 1.0;
    c[static_cast<std::size_t>(s)] = -coef;
    trim_used(c);
    return c;
}

// RationalFunction.of(num, den).coefficients(k): the power-series quotient
// coefficients [0..k). The first prepare(n) call with n = k - 1 sizes the
// cache at max(n, maxdeg) or ((n - 1) / 32 + 1) * 32 entries (Java integer
// division) before the literal recurrence.
std::vector<double> rf_coefficients(const std::vector<double>& num, const std::vector<double>& den,
    std::size_t k) {
    if (k == 0) {
        return {};
    }
    const std::size_t dendeg = den.size() - 1;
    const std::size_t maxdeg = std::max(num.size() - 1, dendeg);
    std::size_t n = k - 1;
    if (n < maxdeg) {
        n = maxdeg;
    } else {
        n = static_cast<std::size_t>(
                (static_cast<std::int64_t>(n) - 1) / 32 + 1) * 32; // Java / truncates
    }
    std::vector<double> c(n + 1, 0.0);
    const double d0 = poly_get(den, 0);
    c[0] = poly_get(num, 0) / d0;
    for (std::size_t i = 1; i <= n; ++i) {
        double s = 0.0;
        const std::size_t jmax = std::min(i, dendeg);
        for (std::size_t j = 1; j <= jmax; ++j) {
            s += poly_get(den, j) * c[i - j];
        }
        const double q = poly_get(num, i); // 0 beyond the numerator degree
        c[i] = (q - s) / d0;
    }
    c.resize(std::min(k, c.size()));
    return c;
}

// LULinearSystemSolver2.solve(A, b) with the pinned column bookkeeping: the
// m x (m + 1) workspace is the column-major storage of A followed by b, the
// pivot search at step k scans row k over the unassigned columns (storage
// index k + m * col), and the solution of column col is written from row
// ip[col] of the augmented block. False mirrors the pinned MatrixException.
bool lu_solve(const std::vector<double>& a, std::size_t m, const double* b, double* sol) {
    const std::size_t nn = m + 1;
    std::vector<double> ab(m * nn, 0.0);
    for (std::size_t i = 0; i < m * m; ++i) {
        ab[i] = a[i];
    }
    for (std::size_t i = 0; i < m; ++i) {
        ab[m * m + i] = b[i];
    }
    std::vector<int> ip(nn, -1);
    for (std::size_t k = 0; k < m; ++k) {
        double max = kMinValue;
        std::int64_t imax = -1;
        for (std::size_t r = 0; r < m; ++r) {
            if (ip[r] != -1) {
                continue;
            }
            const double v = std::abs(ab[k + r * m]);
            if (v > max) { // dcmpl ifle: NaN never wins
                imax = static_cast<std::int64_t>(r);
                max = v;
            }
        }
        if (imax < 0) {
            return false;
        }
        ip[static_cast<std::size_t>(imax)] = static_cast<int>(k);
        const double piv = ab[k + static_cast<std::size_t>(imax) * m];
        if (std::abs(piv) < kLUEps) { // dcmpg: a NaN pivot passes
            return false;
        }
        const double inv = 1.0 / piv;
        for (std::size_t i = 0; i < m; ++i) { // Gauss-Jordan: all other rows
            if (i == k) {
                continue;
            }
            const double u = ab[static_cast<std::size_t>(imax) * m + i];
            if (u == 0.0) { // dcmpl ifeq: a NaN factor eliminates
                continue;
            }
            const double f = inv * u;
            for (std::size_t col = 0; col < nn; ++col) {
                if (ip[col] == -1) {
                    ab[col * m + i] -= f * ab[col * m + k];
                }
            }
        }
        for (std::size_t col = 0; col < nn; ++col) {
            if (ip[col] == -1) {
                ab[col * m + k] *= inv;
            }
        }
    }
    for (std::size_t col = 0; col < m; ++col) {
        if (ip[col] != -1) {
            sol[col] = ab[static_cast<std::size_t>(ip[col]) + m * m];
        }
    }
    return true;
}

// AutoCovarianceFunction.ac: the sma (ardeg == 0) and defaultComputer
// branches, then the 36-aligned prepare growth and values(nb) copy that the
// KalmanFilter consumes. nullopt mirrors "acgf of non stationary model".
std::optional<std::vector<double>> acf_values(    const std::vector<double>& ar, const std::vector<double>& ma, std::size_t nb_want) {    const std::size_t ardeg = ar.size() - 1;
    std::vector<double> ac;
    if (ardeg == 0) {
        // SymmetricFilter.convolutionOf(BackFilter(ma), 1): c[k] accumulates
        // w[u] * w[u + k] in the pinned u-major sweep; ivar enters as 1.
        std::vector<double> c(ma.size(), 0.0);
        for (std::size_t u = 0; u < ma.size(); ++u) {
            for (std::size_t v = u; v < ma.size(); ++v) {
                c[v - u] += ma[u] * ma[v];
            }
        }
        ac = std::move(c);
    } else {
        const std::size_t madeg = ma.size() - 1;
        std::size_t nb = nb_want;
        const std::size_t m6 = std::max(ardeg, madeg) + 1;
        if (nb < m6) {
            nb = m6;
        }
        ac.assign(nb + 1, 0.0);
        const auto psi = rf_coefficients(ma, ar, madeg + 1);
        for (std::size_t i = 0; i <= madeg; ++i) {
            double s = 0.0;
            for (std::size_t j = i; j <= madeg; ++j) {
                s += poly_get(ma, j) * psi[j - i];
            }
            ac[i] = s;
        }
        std::vector<double> a(m6 * m6, 0.0); // column-major
        for (std::size_t i = 0; i < m6; ++i) {
            for (std::size_t j = 0; j <= ardeg; ++j) {
                const double c = poly_get(ar, j);
                if (c != 0.0) {
                    const std::size_t col = i >= j ? i - j : j - i;
                    a[col * m6 + i] += c;
                }
            }
        }
        std::vector<double> sol(m6, 0.0);
        if (!lu_solve(a, m6, ac.data(), sol.data())) {
            return std::nullopt;
        }
        for (std::size_t i = 0; i < m6; ++i) {
            ac[i] = sol[i];
        }
        for (std::size_t i = m6; i <= nb; ++i) {
            double s = 0.0;
            for (std::size_t j = 1; j <= ardeg; ++j) {
                s += poly_get(ar, j) * ac[i - j];
            }
            ac[i] = -s;
        }
        return std::optional<std::vector<double>>(std::move(ac));
    }
    // The sma branch: prepare() grows the cache to the 36-aligned length and
    // the outer ac() loop extends it to that length with the AR recursion
    // (ardeg == 0 here, so every extended entry is -0.0).
    std::size_t nb = nb_want;
    if (nb == 0) {
        nb = 36;
    } else if (const std::size_t rem = nb % 36; rem != 0) {
        nb = nb + 36 - rem;
    }
    if (ac.size() <= nb) {
        const std::size_t old = ac.size();
        ac.resize(nb + 1);
        for (std::size_t i = old; i <= nb; ++i) {
            double s = 0.0;
            for (std::size_t j = 1; j <= ardeg; ++j) {
                s += poly_get(ar, j) * ac[i - j];
            }
            ac[i] = -s;
        }
    }
    ac.resize(std::min(nb_want, ac.size()));
    return std::optional<std::vector<double>>(std::move(ac));
}

// DeterminantalTerm: the scaled detman / detcar accumulator. The JVM throws
// EcoException on a non-finite or non-positive variance; callers catch that
// through std::-domain failure via the returned flag.
struct det_term {
    double detman = 1.0;
    double detcar = 0.0;
    bool failed = false;

    void add(double x) {
        if (!std::isfinite(x) || failed) {
            failed = true; // "Invalid variance"
            return;
        }
        if (!(x > 0.0)) { // dcmpg ifgt: NaN throws, already caught
            failed = true; // "Negative variance"
            return;
        }
        detman *= x;
        while (detman >= 1.0) {
            detman *= 0.0625;
            detcar += 4.0;
        }
        if (detman != 0.0) {
            while (detman <= 0.0625) {
                detman *= 16.0;
                detcar -= 4.0;
            }
        }
    }

    double log_determinant() const {
        // DeterminantalTerm.getLogDeterminant: Math.log(detman) + detcar*LOG2.
        return fdlibm_log(detman) + detcar * kLog2;
    }
};

// Pinned DataPointer.norm2 (LAPACK dnrm2 scaling, order-sensitive), used by
// Reflector.larfg; the dcmpg-guarded scale comparison treats NaN specially.
double data_norm2(const double* p, std::size_t n) {
    if (n < 1) {
        return 0.0;
    }
    if (n == 1) {
        return std::abs(p[0]);
    }
    double scale = 0.0, ssq = 1.0;
    for (std::size_t t = 0; t < n; ++t) {
        const double xi = p[t];
        if (xi == 0.0) {
            continue;
        }
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

// internal KalmanFilter over the pinned stationary model. prepare() runs tx
// and, for the multi-use (design present) machine, calcC; sfilter runs the
// design-free machine and re-derives the state pair like the JVM.
struct kalman {
    bool multi_use = false;
    std::vector<double> phi; // getAr() coefficients, raw JVM signs
    std::vector<double> c0;  // the acf window, tx-transformed in place
    std::vector<double> C;   // calcC state rows (multi-use only)
    std::vector<double> s;   // innovation standard deviations
    double var = 1.0;
    double h0 = 0.0;
    double ldet = std::numeric_limits<double>::quiet_NaN();
    std::size_t dim = 0;
    std::size_t n = 0;
    std::size_t steady_pos = 0;
    bool failed = false; // "invalid ARIMA" / det_term failure

    // KalmanFilter.tlast(w) = -sum_{j=1}^{|phi|-1} phi[j] * w[dim - j], the
    // dsub chain from +0.0.
    double tlast(const double* w) const {
        double s = 0.0;
        for (std::size_t j = 1; j < phi.size(); ++j) {
            s -= phi[j] * w[dim - j];
        }
        return s;
    }

    void prepare(std::vector<double> arfull, std::vector<double> acf, std::size_t n_,
        std::size_t dim_) {
        phi = std::move(arfull);
        dim = dim_;
        n = n_;
        c0 = std::move(acf);
        h0 = c0[0];
        // tx(): the last entry becomes the -sum chain over the original
        // window, then the window shifts left.
        double last = 0.0;
        for (std::size_t j = 1; j < phi.size(); ++j) {
            last -= phi[j] * c0[dim - j];
        }
        for (std::size_t j = 1; j < dim; ++j) {
            c0[j - 1] = c0[j];
        }
        c0[dim - 1] = last;
        if (multi_use) {
            calc_c();
        }
    }

    void calc_c() {
        std::vector<double> cc = c0;
        C.assign(dim * n, 0.0);
        for (std::size_t j = 0; j < dim; ++j) {
            C[j] = c0[j];
        }
        s.assign(n, 0.0);
        double h = h0;
        det_term dterm;
        steady_pos = n - 1;
        const std::size_t kmax = dim - 1;
        std::size_t rowbase = 0;
        std::size_t nextbase = dim;
        for (std::size_t i = 0; i < n; ++i) {
            if (std::isnan(h) || h < 0.0) { // dcmpg ifge: NaN threw earlier
                failed = true;
                return;
            }
            dterm.add(h);
            s[i] = std::sqrt(h);
            if (i < steady_pos) {
                const double c_0 = cc[0];
                const double r = c_0 / h;
                const double tl = tlast(cc.data());
                const double last_c = C[rowbase + kmax];
                for (std::size_t j = 0; j < kmax; ++j) {
                    const double nxt = cc[j + 1];
                    if (r != 0.0) { // dcmpl ifeq: a NaN r eliminates
                        cc[j] = nxt - C[rowbase + j] * r;
                        C[nextbase + j] = C[rowbase + j] - r * nxt;
                    } else {
                        cc[j] = nxt;
                        C[nextbase + j] = C[rowbase + j];
                    }
                }
                cc[kmax] = tl - r * last_c;
                C[nextbase + kmax] = C[rowbase + kmax] - r * tl;
                h = h - c_0 * r;
                if (h < var) { // dcmpg: a NaN h keeps NaN
                    h = var;
                }
                std::size_t j = 0;
                for (; j < dim; ++j) {
                    if (std::abs(cc[j]) > kSteadyTol) { // dcmpl ifle: NaN keeps
                        // scanning, matching the pinned loop
                        break;
                    }
                }
                if (j == dim) {
                    steady_pos = i + 1;
                }
            }
            rowbase = nextbase;
            nextbase += dim;
        }
        if (dterm.failed) {
            failed = true;
            return;
        }
        ldet = dterm.log_determinant();
    }

    // sfilter: the design-free filter with the (cc, cc2) state pair; sets
    // ldet from its own DeterminantalTerm run, exactly like the JVM.
    std::optional<std::vector<double>> sfilter(const std::vector<double>& y) {
        det_term dterm;
        std::vector<double> cc = c0;
        std::vector<double> cc2 = cc;
        double h = h0;
        std::vector<double> w(dim, 0.0);
        std::vector<double> out(n, 0.0);
        steady_pos = n - 1;
        const std::size_t kmax = dim - 1;
        for (std::size_t i = 0; i < n; ++i) {
            if (std::isnan(h) || h < 0.0) {
                failed = true;
                return std::nullopt;
            }
            dterm.add(h);
            if (dterm.failed) {
                failed = true;
                return std::nullopt;
            }
            const double root = std::sqrt(h);
            const double innov = y[i] - w[0];
            out[i] = innov / root;
            const double tl = tlast(w.data());
            if (innov != 0.0) { // dcmpl ifeq: a NaN innovation updates
                const double r = innov / h;
                for (std::size_t j = 0; j < kmax; ++j) {
                    w[j] = w[j + 1] + cc[j] * r;
                }
                w[kmax] = tl + cc[kmax] * r;
            } else {
                for (std::size_t j = 0; j < kmax; ++j) {
                    w[j] = w[j + 1];
                }
                w[kmax] = tl;
            }
            if (i < steady_pos) {
                const double c_0 = cc2[0];
                const double r2 = c_0 / h;
                const double tl2 = tlast(cc2.data());
                const double last_c = cc[kmax];
                for (std::size_t j = 0; j < kmax; ++j) {
                    const double nxt = cc2[j + 1];
                    if (r2 != 0.0) {
                        cc2[j] = nxt - cc[j] * r2;
                        cc[j] = cc[j] - r2 * nxt;
                    } else {
                        cc2[j] = nxt;
                    }
                }
                cc2[kmax] = tl2 - r2 * last_c;
                cc[kmax] = cc[kmax] - r2 * tl2;
                h = h - c_0 * r2;
                if (h < var) {
                    h = var;
                }
                std::size_t j = 0;
                for (; j < dim; ++j) {
                    if (std::abs(cc2[j]) > kSteadyTol) {
                        break;
                    }
                }
                if (j == dim) {
                    steady_pos = i + 1;
                }
            }
        }
        if (dterm.failed) {
            failed = true;
            return std::nullopt;
        }
        ldet = dterm.log_determinant();
        return std::optional<std::vector<double>>(std::move(out));
    }

    // mfilter: the multi-use whitening with the pinned sticky "seen a
    // non-zero response" flag.
    std::vector<double> mfilter(const std::vector<double>& y) const {
        std::vector<double> w(dim, 0.0);
        std::vector<double> out(n, 0.0);
        bool nonzero_seen = false;
        const std::size_t kmax = dim - 1;
        std::size_t crow = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const double yi = y[i];
            if (yi != 0.0) { // dcmpl ifeq: NaN counts as non-zero
                nonzero_seen = true;
            }
            const double si = s[i];
            const double r = (yi - w[0]) / si;
            out[i] = r;
            if (nonzero_seen) {
                const double tl = tlast(w.data());
                if (r != 0.0) { // dcmpl ifeq: a NaN r updates
                    const double k = r / si;
                    for (std::size_t j = 0; j < kmax; ++j) {
                        w[j] = w[j + 1] + C[crow + j] * k;
                    }
                    w[kmax] = tl + C[crow + kmax] * k;
                } else {
                    for (std::size_t j = 0; j < kmax; ++j) {
                        w[j] = w[j + 1];
                    }
                    w[kmax] = tl;
                }
            }
            if (i < steady_pos) {
                crow += dim;
            }
        }
        return out;
    }
};

// HouseholderWithPivoting over a single column: the fixed-block branch
// (nfixed = 1) and the dqrdc2 branch (nfixed = 0) both run exactly one
// Reflector.larfg over the whole column, so they produce identical storage;
// only the pivot differs (null vs [0]), which the least-squares tail treats
// identically for a single column.
struct hwp_1col {
    std::vector<double> v; // v[0] is the R diagonal, v[1..] the reflector
    double beta = 0.0;

    explicit hwp_1col(const std::vector<double>& col) : v(col) {
        const std::size_t m = v.size();
        const double s0 = kSafeMin / kEpsilon;
        double alpha = v[0];
        double norm = data_norm2(v.data(), m);
        int count = 0;
        if (norm != 0.0 && std::isless(norm, s0)) { // dcmpg ifge: NaN skips
            while (true) {
                const double sc = 1.0 / s0;
                for (std::size_t t = 1; t < m; ++t) {
                    v[t] *= sc;
                }
                alpha *= sc;
                norm *= sc;
                if (!std::isless(norm, s0)) { // dcmpg ifge: NaN recomputes
                    norm = data_norm2(v.data(), m);
                    break;
                }
                ++count;
                if (count >= 4) {
                    norm = data_norm2(v.data(), m);
                    break;
                }
            }
        }
        if (std::isless(alpha, 0.0)) { // dcmpg ifge: NaN keeps the sign
            norm = -norm;
        }
        for (int t = 0; t < count; ++t) {
            norm *= s0;
        }
        beta = norm == 0.0 ? 0.0 : norm / (norm + alpha);
        for (std::size_t t = 1; t < m; ++t) {
            v[t] /= norm; // the stored reflector tail; the pivot is implicit 1
        }
        v[0] = -norm; // Reflector.alpha, the rawR() diagonal after decompose
    }
};

// QRDecomposition.applyQt on the (single-reflector) whitened right-hand
// side, then UpperTriangularMatrix.rank / solveUx with the pinned rcond:
// fastRank counts |R(j,j)| > rcond and a NaN diagonal does not count (the
// lambda$fastRank$1 dcmpl returns false for NaN).
struct ls_solution {
    int rank = 0;
    double bb = 0.0;
    std::vector<double> e;
    double ssq_err = 0.0;
    bool singular = false;
};

ls_solution least_squares(const hwp_1col& hwp, const std::vector<double>& w,
    double rcond = kRcond) {
    ls_solution out;
    const double r00 = hwp.v[0];
    if (std::isgreaterequal(std::abs(r00), rcond)) { // NaN does not count
        out.rank = 1;
    }
    std::vector<double> b = w;
    const std::size_t m = b.size();
    if (hwp.beta != 0.0) { // dcmpl ifeq: a NaN beta still applies
        double t = b[0] / hwp.beta;
        for (std::size_t i = 1; i < m; ++i) {
            t += b[i] * hwp.v[i];
        }
        if (t != 0.0) {
            b[0] -= t;
            const double c = -t * hwp.beta;
            for (std::size_t i = 1; i < m; ++i) {
                b[i] += c * hwp.v[i];
            }
        }
    }
    if (out.rank == 1) {
        if (r00 == 0.0) { // the pinned m_err_sing mirror
            out.singular = true;
            return out;
        }
        out.bb = b[0] / r00;
    }
    // DataBlock.ssq: ascending, skipping the non-finite coordinates.
    double acc = 0.0;
    for (std::size_t i = static_cast<std::size_t>(out.rank); i < m; ++i) {
        if (std::isfinite(b[i])) {
            acc += b[i] * b[i];
        }
    }
    out.e.assign(b.begin() + out.rank, b.end());
    out.ssq_err = acc;
    return out;
}

// SarimaModel.getRegularAR/getRegularMA: [1, c...] through ofInternal.
std::vector<double> raw_poly(const std::vector<double>& c) {
    std::vector<double> p;
    p.reserve(c.size() + 1);
    p.push_back(1.0);
    p.insert(p.end(), c.begin(), c.end());
    trim_used(p);
    return p;
}

// SarimaModel.seasonalAR/seasonalMA: ONE for no seasonals, factor(-b, s) for
// a single seasonal (the JVM stores the raw sign, factor flips it back), the
// lag-padded ofInternal array otherwise.
std::vector<double> seasonal_poly(const std::vector<double>& b, int s) {
    if (b.empty()) {
        return {1.0};
    }
    if (b.size() == 1) {
        return poly_factor(-b[0], s);
    }
    std::vector<double> p(static_cast<std::size_t>(s) * b.size() + 1, 0.0);
    p[0] = 1.0;
    for (std::size_t j = 0; j < b.size(); ++j) {
        p[(j + 1) * static_cast<std::size_t>(s)] = b[j];
    }
    trim_used(p);
    return p;
}

} // namespace

std::vector<double> pinned_differencing_polynomial(int period, int d, int bd) {
    std::vector<double> f{1.0};
    for (int i = 0; i < d; ++i) {
        std::vector<double> next(f.size() + 1, 0.0);
        for (std::size_t j = 0; j < f.size(); ++j) {
            next[j] += f[j];
            next[j + 1] -= f[j];
        }
        f = std::move(next);
    }
    for (int i = 0; i < bd; ++i) {
        std::vector<double> next(f.size() + static_cast<std::size_t>(period), 0.0);
        for (std::size_t j = 0; j < f.size(); ++j) {
            next[j] += f[j];
            next[j + static_cast<std::size_t>(period)] -= f[j];
        }
        f = std::move(next);
    }
    return f;
}

std::vector<double> pinned_backfilter_apply(
    const std::vector<double>& f, const std::vector<double>& x) {
    const std::size_t deg = f.size() - 1;
    if (deg >= x.size()) {
        return {};
    }
    std::vector<double> out(x.size() - deg, 0.0);
    for (std::size_t t = deg; t < x.size(); ++t) {
        double acc = 0.0;
        for (std::size_t k = 0; k <= deg; ++k) {
            const double w = f[deg - k];
            if (w != 0.0) { // the pinned apply skips zero weights
                acc += w * x[t - deg + k];
            }
        }
        out[t - deg] = acc;
    }
    return out;
}

double pinned_ols_1col(const std::vector<double>& col, const std::vector<double>& y) {
    // Pinned OlsComputer with one regressor column: HouseholderWithPivot
    // (nfixed = 0) plus QRLeastSquaresSolver.leastSquares at the
    // robustLeastSquares rcond 1e-9. A rank-0 solution hands back the
    // zero-padded getB(), i.e. a zero coefficient.
    const hwp_1col hwp(col);
    return least_squares(hwp, y, 1e-9).bb;
}

std::optional<PinnedMachinePoint> pinned_deviances_evaluate(
    const std::vector<double>& y, std::size_t p, std::size_t bp, std::size_t q, std::size_t bq,
    int period, const std::vector<double>& design, const std::vector<double>& jvm) {
    if (jvm.size() != p + bp + q + bq || period < 2 || y.empty()) {
        return std::nullopt;
    }
    if (!design.empty() && design.size() != y.size()) {
        return std::nullopt;
    }
    // JVM semantics: the machine runs with the default floating point
    // environment regardless of what the caller's BLAS back end left behind.
    const default_fp_environment fp_guard;
    for (const double v : jvm) {
        if (!std::isfinite(v)) {
            return std::nullopt; // the pinned domain never emits them; the JVM
                                 // machine would carry NaN into a throw
        }
    }
    auto split = [&](std::size_t off, std::size_t len) {
        return std::vector<double>(
            jvm.begin() + static_cast<std::ptrdiff_t>(off),
            jvm.begin() + static_cast<std::ptrdiff_t>(off + len));
    };
    const auto phi = split(0, p);
    const auto bphi = split(p, bp);
    const auto th = split(p + bp, q);
    const auto bth = split(p + bp + q, bq);

    // Stationary model polynomials (the unit roots are already differenced
    // out of y): getStationaryAr = regAR.times(seasAR), getMa likewise.
    const auto reg_ar = raw_poly(phi);
    const auto seas_ar = seasonal_poly(bphi, period);
    const auto ar = poly_times(reg_ar, seas_ar);
    const auto reg_ma = raw_poly(th);
    const auto seas_ma = seasonal_poly(bth, period);
    const auto ma = poly_times(reg_ma, seas_ma);

    // KalmanFilter dimensions: getArOrder / getMaOrder over the stationary
    // spec counts (d = bd = 0 after the unit-root filter).
    const std::size_t ar_order = p + static_cast<std::size_t>(period) * bp;
    const std::size_t ma_order = q + static_cast<std::size_t>(period) * bq;
    const std::size_t dim = std::max(ar_order, ma_order + 1);

    auto acf = acf_values(ar, ma, dim);
    if (!acf || acf->size() < dim) {
        return std::nullopt; // "acgf of non stationary model"
    }
    acf->resize(dim);

    const bool multi_use = !design.empty();
    kalman kf;
    kf.multi_use = multi_use;
    kf.prepare(ar, std::move(*acf), y.size(), dim);
    if (kf.failed) {
        return std::nullopt;
    }

    PinnedMachinePoint pt;
    std::vector<double> w;
    if (multi_use) {
        w = kf.mfilter(y);
    } else {
        auto filtered = kf.sfilter(y);
        if (!filtered || kf.failed) {
            return std::nullopt;
        }
        w = std::move(*filtered);
    }
    pt.ldet = kf.ldet;

    std::vector<double> xx;
    ls_solution sol;
    if (multi_use) {
        xx = kf.mfilter(design);
        hwp_1col hwp(xx);
        sol = least_squares(hwp, w);
        if (sol.singular) {
            return std::nullopt;
        }
        if (sol.rank == 0) {
            // ConcentratedLikelihoodComputer: the rank-0 branch profiles no
            // regressor and uses w.ssq() with the raw residuals.
            double wssq = 0.0;
            for (const double v : w) {
                if (std::isfinite(v)) {
                    wssq += v * v;
                }
            }
            pt.ssq_err = wssq;
            pt.rank = 0;
        } else {
            pt.ssq_err = sol.ssq_err;
            pt.rank = 1;
            pt.coefficient = sol.bb;
            const double inv = 1.0 / hwp.v[0];
            pt.uc00 = inv * inv; // UpperTriangularMatrix.inverse on a 1x1
                                 // block, then UUt(inv)
        }
    } else {
        double wssq = 0.0;
        for (const double v : w) {
            if (std::isfinite(v)) {
                wssq += v * v;
            }
        }
        pt.ssq_err = wssq;
        pt.rank = 0;
    }

    const double n = static_cast<double>(y.size());
    // Likelihood.factor(): Math.exp(logDeterminant()/dim()).
    pt.factor = fdlibm_exp(pt.ldet / n);
    if (!std::isfinite(pt.factor)) {
        return std::nullopt;
    }
    pt.ssq_e = pt.ssq_err * pt.factor;
    pt.log_likelihood =
        -0.5 * (n * kLogTwoPi + n * (1.0 + fdlibm_log(pt.ssq_err / n)) + pt.ldet);

    // ConcentratedLikelihoodComputer.process residual assembly: the full
    // residuals w - sum b_j xx_j (design present, rank 1) or w, then
    // deviances(): select(isFinite) and map(x -> x * sqrt(factor)) with the
    // pinned factor == 1.0 shortcut.
    std::vector<double> res = w;
    if (multi_use && sol.rank == 1) {
        for (std::size_t i = 0; i < res.size(); ++i) {
            res[i] += -sol.bb * xx[i]; // the pinned dneg + addAY order
        }
    }
    const double root = std::sqrt(pt.factor);
    pt.e.clear();
    pt.e.reserve(res.size());
    for (const double v : res) {
        if (!std::isfinite(v)) {
            continue; // DoubleSeq.select(isFinite)
        }
        pt.e.push_back(pt.factor == 1.0 ? v : v * root); // dcmpl ifne shortcut
    }

    if (multi_use && sol.rank == 1) {
        // ConcentratedLikelihood.tstat(0, 0, false): the ser order is
        // sqrt((uc * ssq) / dim) with the uc == 0 and coefficient == 0
        // shortcuts.
        if (pt.uc00 == 0.0) {
            pt.tstat = pt.coefficient / std::numeric_limits<double>::quiet_NaN();
        } else if (pt.coefficient == 0.0) {
            pt.tstat = 0.0;
        } else {
            const double ser = std::sqrt((pt.uc00 * pt.ssq_err) / n);
            pt.tstat = pt.coefficient / ser;
        }
    }
    return pt;
}

} // namespace kalos::detail
