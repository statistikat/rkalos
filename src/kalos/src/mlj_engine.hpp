#ifndef KALOS_MLJ_ENGINE_HPP
#define KALOS_MLJ_ENGINE_HPP

// Shared ML-residual machine extracted verbatim (item 21-H, pure move) from
// src/outlier_scan_engine.cpp: the pinned AutoCovarianceFunction computer,
// the ModifiedLjungBox / MaLjungBox filter, the rational-function expansion
// and the pinned triangular-solver / Cholesky leaves they run on. The
// 21-C scan engine and the 21-H accept-model machine (src/arma_diagnostics
// .cpp) build on the identical code; semantics are qualified by the
// 21-C corpus sweep (docs/contracts/outlier-scoring-contract.md).

#include "outlier_detail.hpp"

#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace kalos {
namespace mlj {

using outlier_detail::dcmpl;
using outlier_detail::dcmpg;
using outlier_detail::robust_lu_solve;

inline double dot(const double* a, const double* b, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

inline double ssq(const double* a, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * a[i];
    return s;
}

// RationalFunction long division: expansion of num / den, pinned recursion
// order (the cache-block split in the bytecode is bookkeeping; every branch
// computes this recursion). den[0] must be nonzero (monic polynomials).
inline std::vector<double> rational_coefficients(
    const std::vector<double>& num, const std::vector<double>& den, std::size_t count) {
    std::vector<double> a(count, 0.0);
    if (count == 0) return a;
    a[0] = num[0] / den[0];
    for (std::size_t k = 1; k < count; ++k) {
        double s = 0.0;
        for (std::size_t j = 1; j <= k && j < den.size(); ++j) s += den[j] * a[k - j];
        const double n_k = k < num.size() ? num[k] : 0.0;
        a[k] = (n_k - s) / den[0];
    }
    return a;
}

// ---- AutoCovarianceFunction (default linear-system computer) -------------------

// Pinned triangular-solver and Cholesky leaves (defined with the other leaves
// below); declared up front so the ML-residual machines can use them.
inline void symmetric_lcholesky(std::vector<double>& a, std::size_t n);
inline void solve_lx(const std::vector<double>& a, std::size_t n, double* x, std::size_t xlen);
inline void solve_lx(const std::vector<double>& a, std::size_t n, std::vector<double>& x);
inline void solve_xl(const std::vector<double>& a, std::size_t n, std::vector<double>& x);
inline void solve_xu(const std::vector<double>& u, std::size_t n, std::vector<double>& x);
inline void solve_ux(const std::vector<double>& u, std::size_t n, std::vector<double>& x);
inline std::vector<double> mat_mul(const std::vector<double>& a, const std::vector<double>& b,
    std::size_t ar, std::size_t ac, std::size_t bc);

class Acvf {
public:
    Acvf(std::vector<double> ma, std::vector<double> ar) : ma_(std::move(ma)), ar_(std::move(ar)) {}

    // values(n): gamma[0..n) with the pinned prepare() rounding to multiples
    // of 36 and the AR-recursion extension inside.
    std::vector<double> values(std::size_t n) {
        prepare(n);
        std::vector<double> out(n, 0.0);
        for (std::size_t i = 0; i < n && i < cache_.size(); ++i) out[i] = cache_[i];
        return out;
    }

private:
    void prepare(std::size_t n) {
        std::size_t nn = n;
        if (nn == 0) {
            nn = 36;
        } else {
            const std::size_t rem = nn % 36;
            if (rem != 0) nn = nn - rem + 36;
        }
        if (!cache_.empty() && cache_.size() > nn) return;
        cache_ = extend(nn);
    }

    std::vector<double> extend(std::size_t n) {
        const std::size_t p = ar_.size() - 1;
        const std::size_t q = ma_.size() - 1;
        std::size_t m = std::max(p, q) + 1;
        if (n < m) n = m;
        std::vector<double> res(n + 1, 0.0);
        const auto psi = rational_coefficients(ma_, ar_, q + 1);
        std::vector<double> A(m * m, 0.0);
        for (std::size_t i = 0; i <= q; ++i) {
            double s = 0.0;
            for (std::size_t j = i; j <= q; ++j) s += ma_[j] * psi[j - i];
            res[i] = s;
        }
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t jj = 0; jj <= p; ++jj) {
                const double a = ar_[jj];
                if (dcmpl(a, 0.0) == 0) continue;
                const std::size_t col = i > jj ? i - jj : jj - i;
                A[col * m + i] += a; // add(row = i, col = |i - jj|, a)
            }
        }
        std::vector<double> rhs(res.begin(), res.begin() + static_cast<std::ptrdiff_t>(m));
        try {
            robust_lu_solve(A, m, rhs);
        } catch (const std::exception&) {
            throw std::runtime_error("acgf of non stationary model");
        }
        for (std::size_t i = 0; i < m; ++i) res[i] = rhs[i];
        for (std::size_t i = m; i <= n; ++i) {
            double s = 0.0;
            for (std::size_t j = 1; j <= p; ++j) s += ar_[j] * res[i - j];
            res[i] = -s;
        }
        // ivar == 1.0 always here (SarimaModel constant) — pinned skip.
        return res;
    }

    std::vector<double> ma_;
    std::vector<double> ar_;
    std::vector<double> cache_;
};

// ---- ModifiedLjungBox / MaLjungBox (ML residuals) -------------------------------

class MaLjungBox {
public:
    void prepare(const std::vector<double>& ma, std::size_t n) {
        ma_ = ma;
        n_ = n;
        q_ = ma.size() - 1;
        // The JVM bounds semantics of calcg(int): psi has m_n coefficients
        // and the pinned column recursion reads psi[m_n - j] for every
        // j < q, so a ma degree beyond n + 1 raises an
        // ArrayIndexOutOfBoundsException in the JVM (propagating out of
        // prepare). Transcribe the trap instead of performing the read
        // (21-H.1: an out-of-range read here is undefined behaviour).
        if (q_ >= n_ + 2) {
            throw std::runtime_error(
                "MaLjungBox: ma degree exceeds the ml residual window");
        }
        v1_.assign(q_ * q_, 0.0);
        if (q_ > 0) {
            for (std::size_t k = 0; k < q_; ++k) { // diagonal + subDiagonal(k)
                const std::size_t len = q_ - k;
                const double value = ma_[q_ - k];
                for (std::size_t c = 0; c < len; ++c) {
                    v1_[(c + k) * q_ + c] = value; // (row = c + k, col = c)
                }
            }
        }
        calc_g();
        // X = XtSX(G, V1) = V1^T G V1 (tmp = G*V1, lower-triangle dots,
        // fromLower), then diag += 1 and the symmetric Cholesky.
        std::vector<double> tmp = mat_mul(g_, v1_, q_, q_, q_);
        x_.assign(q_ * q_, 0.0);
        for (std::size_t c = 0; c < q_; ++c) {
            for (std::size_t r = c; r < q_; ++r) {
                x_[c * q_ + r] = dot(&v1_[r * q_], &tmp[c * q_], q_);
            }
        }
        for (std::size_t c = 0; c < q_; ++c) {
            for (std::size_t r = 0; r < c; ++r) x_[c * q_ + r] = x_[r * q_ + c];
        }
        for (std::size_t k = 0; k < q_; ++k) x_[k * q_ + k] += 1.0;
        symmetric_lcholesky(x_, q_);
    }

    // filter(a, out): out[0..q) = u, out[q..q + len) = a0 - z.
    void filter(const std::vector<double>& a, double* out) const {
        const std::size_t len = a.size();
        std::vector<double> a0 = a;
        rma(a0);
        std::vector<double> g = calc_g_values(a0);
        u_.assign(q_, 0.0);
        for (std::size_t i = 0; i < q_; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j <= i; ++j) s += ma_[q_ - i + j] * g[j];
            u_[i] = s;
        }
        std::vector<double> u = u_;
        solve_lx(x_, q_, u);    // LowerTriangularMatrix.solveLx, tol 0
        solve_xl(x_, q_, u);    // LowerTriangularMatrix.solvexL, tol 0
        std::vector<double> z(len, 0.0);
        for (std::size_t i = 0; i < q_; ++i) {
            double s = 0.0;
            for (std::size_t j = i; j < q_; ++j) s += ma_[q_ + i - j] * u[j];
            z[i] = s;
        }
        rma(z);
        for (std::size_t i = 0; i < q_; ++i) out[i] = u[i];
        for (std::size_t i = 0; i < len; ++i) out[q_ + i] = a0[i] - z[i];
    }

    std::size_t warmup() const { return q_; }

private:
    // rma: recursive MA filtering, in place, pinned two-phase loops.
    void rma(std::vector<double>& a) const {
        if (q_ == 0) return;
        for (std::size_t i = 1; i < q_ && i < a.size(); ++i) {
            double s = 0.0;
            for (std::size_t j = 1; j <= i; ++j) s += ma_[j] * a[i - j];
            a[i] -= s;
        }
        for (std::size_t i = q_; i < a.size(); ++i) {
            double s = 0.0;
            for (std::size_t j = 1; j <= q_; ++j) s += ma_[j] * a[i - j];
            a[i] -= s;
        }
    }

    // calcg(int): the G matrix.
    void calc_g() {
        g_.assign(q_ * q_, 0.0);
        if (q_ == 0) return;
        const std::vector<double> one{1.0};
        const auto psi = rational_coefficients(one, ma_, n_);
        for (std::size_t i = 0; i < q_; ++i) {
            double s = 0.0;
            for (std::size_t j = i; j < n_; ++j) s += psi[j] * psi[j - i];
            g_[0 * q_ + i] = s; // set(row = i, col = 0, s)
        }
        for (std::size_t i = 1; i < q_; ++i) {
            for (std::size_t j = i; j < q_; ++j) {
                // Pinned column(j).set(i, column(j-1).get(i-1) - psi[n-i]*psi[n-j]):
                // the cell is OVERWRITTEN from the previous (already updated)
                // column, not decremented - there is no dup2 read of the
                // current cell in the pinned bytecode.
                g_[i * q_ + j] = g_[(i - 1) * q_ + (j - 1)] - psi[n_ - j] * psi[n_ - i];
            }
        }
        for (std::size_t c = 0; c < q_; ++c) {
            for (std::size_t r = 0; r < c; ++r) g_[c * q_ + r] = g_[r * q_ + c];
        }
    }

    // calcg(double[]): the g recursion on a data block.
    std::vector<double> calc_g_values(const std::vector<double>& x) const {
        std::vector<double> g = x;
        if (q_ > 0) {
            const std::size_t n = n_;
            if (n >= 2) {
                for (std::size_t jj = n - 2; jj + 1 >= 1; --jj) {
                    double acc = 0.0;
                    for (std::size_t i = 1, j = jj + 1; i <= q_ && j < n; ++i, ++j) {
                        acc += ma_[i] * g[j];
                    }
                    g[jj] -= acc;
                    if (jj == 0) break;
                }
            }
        }
        return g;
    }

    std::vector<double> ma_;
    std::size_t n_ = 0;
    std::size_t q_ = 0;
    std::vector<double> g_;
    std::vector<double> x_;
    std::vector<double> v1_;
    mutable std::vector<double> u_;
};

class ModifiedLjungBox {
public:
    ModifiedLjungBox(std::vector<double> ar, std::vector<double> ma, Acvf& acvf,
        std::vector<double> psi, std::size_t n)
        : ar_(std::move(ar)), ma_(std::move(ma)), acvf_(acvf), psi_(std::move(psi)), n_(n) {
        p_ = ar_.size() - 1;
        q_ = ma_.size() - 1;
        // The JVM operand n - p is passed as a Java int: a negative value
        // raises a NegativeArraySizeException inside MaLjungBoxFilter.
        // Transcribe the trap (size_t arithmetic would otherwise wrap).
        if (n_ < p_) {
            throw std::runtime_error(
                "ModifiedLjungBox: window below the ar degree");
        }
        if (q_ > 0) malb_.prepare(ma_, n_ - p_);
        if (p_ == 0) return;
        const auto gamma = acvf_.values(p_);
        l_.assign(p_ * p_, 0.0);
        for (std::size_t k = 0; k < p_; ++k) { // diagonal + subDiagonal(k)
            const std::size_t len = p_ - k;
            for (std::size_t c = 0; c < len; ++c) l_[(c + k) * p_ + c] = gamma[k];
        }
        if (q_ > 0) {
            const std::size_t rows_g = n_ - p_;
            std::vector<double> G(rows_g * p_, 0.0);
            c_.assign((rows_g + q_) * p_, 0.0);
            for (std::size_t j = 0; j < p_; ++j) {
                std::vector<double> g(rows_g, 0.0);
                for (std::size_t i = 0; i < q_; ++i) {
                    const std::size_t k = p_ + i - j;
                    if (k > q_) break;
                    double s = 0.0;
                    for (std::size_t t = k; t <= q_; ++t) s += ma_[t] * psi_[t - k];
                    g[i] = s;
                }
                malb_.filter(g, &c_[j * (rows_g + q_)]);
            }
            for (std::size_t j = 0; j < p_; ++j) {
                for (std::size_t i = 0; i <= j; ++i) {
                    const double d = dot(&c_[j * (rows_g + q_)], &c_[i * (rows_g + q_)],
                        rows_g + q_);
                    l_[j * p_ + i] -= d; // add(row = i, col = j, -d)
                }
            }
        }
        for (std::size_t c = 0; c < p_; ++c) {
            for (std::size_t r = 0; r < c; ++r) l_[c * p_ + r] = l_[r * p_ + c]; // fromUpper
        }
        symmetric_lcholesky(l_, p_);
    }

    // ML residuals of w: warmup = n + q; the caller keeps the last w entries.
    std::vector<double> residuals(const std::vector<double>& w) const {
        std::vector<double> out(n_ + q_, 0.0);
        apply(w, out.data());
        return std::vector<double>(out.begin() + static_cast<std::ptrdiff_t>(q_), out.end());
    }

public:
    // apply is public here as it is in the JVM (public ArmaFilter method);
    // residuals() remains the convenience wrapper.
    void apply(const std::vector<double>& src, double* dest) const {
        if (p_ > 0) {
            std::vector<double> a(n_ - p_);
            for (std::size_t j = 0; j < a.size(); ++j) a[j] = src[p_ + j];
            for (std::size_t k = 1; k <= p_; ++k) {
                const double c = ar_[k];
                for (std::size_t j = 0; j < a.size(); ++j) a[j] += c * src[p_ + j - k];
            }
            if (q_ > 0) {
                malb_.filter(a, dest + p_);
            } else {
                for (std::size_t j = 0; j < a.size(); ++j) dest[p_ + j] = a[j];
            }
            if (!c_.empty()) {
                const std::size_t tail = n_ + q_ - p_;
                for (std::size_t i = 0; i < p_; ++i) {
                    dest[i] = src[i] - dot(&c_[i * tail], dest + p_, tail);
                }
            } else {
                for (std::size_t i = 0; i < p_; ++i) dest[i] = src[i];
            }
            solve_lx(l_, p_, dest, p_);
        } else {
            if (q_ > 0) {
                malb_.filter(src, dest);
            } else {
                for (std::size_t j = 0; j < n_; ++j) dest[j] = src[j];
            }
        }
    }

private:
    std::vector<double> ar_;
    std::vector<double> ma_;
    Acvf& acvf_;
    std::vector<double> psi_;
    std::size_t n_;
    std::size_t p_ = 0;
    std::size_t q_ = 0;
    MaLjungBox malb_;
    std::vector<double> l_;
    std::vector<double> c_;
};

// ---- Triangular solvers and symmetric Cholesky (pinned leaves) -------------------

// SymmetricMatrix.lcholesky via the provider (lapack.Cholesky, tol 0):
// lower triangle in place; diag <= -tol throws m_err_chol, |diag| < tol sets
// 0 (impossible at tol 0 unless exactly 0, which the first branch catches).
inline void symmetric_lcholesky(std::vector<double>& a, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double s = a[j * n + i];
            for (std::size_t k = 0; k < j; ++k) s -= a[k * n + i] * a[k * n + j];
            if (i == j) {
                if (dcmpg(s, 0.0) <= 0) throw std::runtime_error("m_err_chol");
                a[j * n + i] = std::sqrt(s);
            } else {
                a[j * n + i] = s / a[j * n + j];
            }
        }
    }
}

// LowerTriangularMatrix.solveLx, tol 0, pinned eliminate-after-divide loop:
// q = x[i] / diag; x[i] = q; then x[i+1..] -= q * (column i below the diag).
// Singular (|diag| not > 0, so 0 or NaN): x[i] = 0 only when |x| < tol (never
// at tol 0), else the pinned MatrixException m_err_sing.
inline void solve_lx(const std::vector<double>& a, std::size_t n, double* x, std::size_t xlen) {
    (void)xlen;
    for (std::size_t i = 0, pos = 0; i < n; ++i, pos += n + 1) {
        const double xi = x[i];
        const double d = a[pos];
        if (!(dcmpg(std::fabs(d), 0.0) > 0)) {
            if (dcmpl(std::fabs(xi), 0.0) < 0) {
                x[i] = 0.0;
                continue;
            }
            throw std::runtime_error("m_err_sing");
        }
        const double q = xi / d;
        x[i] = q;
        std::size_t p = pos + 1;
        for (std::size_t j = i + 1; j < n; ++j, ++p) x[j] -= q * a[p];
    }
}

inline void solve_lx(const std::vector<double>& a, std::size_t n, std::vector<double>& x) {
    solve_lx(a, n, x.data(), x.size());
}

// LowerTriangularMatrix.solvexL (solveLtx), tol 0, pinned loop: i from n-1
// down; s = x[i] - sum_{j>i} L[j][i] * x[j] (column i below the diag, pinned
// product order storage[p] * a[j]); x[i] = s / diag with the same singular
// mapping as solveLx on the reduced s.
inline void solve_xl(const std::vector<double>& a, std::size_t n, std::vector<double>& x) {
    for (std::size_t ii = n; ii-- > 0;) {
        const std::size_t pos = ii * n + ii;
        double s = x[ii];
        std::size_t p = pos + 1;
        for (std::size_t j = ii + 1; j < n; ++j, ++p) s -= a[p] * x[j];
        const double d = a[pos];
        if (!(dcmpg(std::fabs(d), 0.0) > 0)) {
            if (dcmpl(std::fabs(s), 0.0) < 0) {
                x[ii] = 0.0;
                continue;
            }
            throw std::runtime_error("m_err_sing");
        }
        x[ii] = s / d;
    }
}

// UpperTriangularMatrix.solvexU (solveUtx), tol 0: x <- x * inv(U).
// Pinned: j = 0..n-1 ascending, x[j] = (x[j] - sum_{i<j} x[i] * U[i][j]) / U[j][j].
inline void solve_xu(const std::vector<double>& u, std::size_t n, std::vector<double>& x) {
    for (std::size_t j = 0; j < n; ++j) {
        const double d = u[j * n + j];
        double s = x[j];
        if (!(dcmpg(std::fabs(d), 0.0) > 0)) {
            if (dcmpl(std::fabs(s), 0.0) < 0) {
                x[j] = 0.0;
                continue;
            }
            throw std::runtime_error("m_err_sing");
        }
        for (std::size_t i = 0; i < j; ++i) s -= x[i] * u[j * n + i];
        x[j] = s / d;
    }
}

// UpperTriangularMatrix.solveUx (tol 0), back substitution; x is the leading
// n entries (pinned works in place on the qy buffer).
inline void solve_ux(const std::vector<double>& u, std::size_t n, std::vector<double>& x) {
    for (std::size_t ii = n; ii-- > 0;) {
        const double d = u[ii * n + ii];
        double v = x[ii];
        if (!(dcmpg(std::fabs(d), 0.0) > 0)) {
            if (dcmpl(std::fabs(v), 0.0) < 0) {
                x[ii] = 0.0;
                continue;
            }
            throw std::runtime_error("m_err_sing");
        }
        v = v / d;
        x[ii] = v;
        for (std::size_t j = ii; j-- > 0;) x[j] -= v * u[ii * n + j];
    }
}

// GeneralMatrix.AB (dense, column-major): transcribed dot accumulation with
// the ikj column order of the pinned BLAS kernel.
inline std::vector<double> mat_mul(const std::vector<double>& a, const std::vector<double>& b,
    std::size_t ar, std::size_t ac, std::size_t bc) {
    std::vector<double> r(ar * bc, 0.0);
    for (std::size_t c = 0; c < bc; ++c) {
        for (std::size_t k = 0; k < ac; ++k) {
            const double bk = b[c * ac + k];
            if (bk == 0.0) continue;
            for (std::size_t i = 0; i < ar; ++i) r[c * ar + i] += a[k * ar + i] * bk;
        }
    }
    return r;
}

// Polynomial.times: r[i + j] += a[i] * b[j], i then j ascending.
inline std::vector<double> poly_times(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> r(a.size() + b.size() - 1, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] += a[i] * b[j];
        }
    }
    return r;
}

// UnitRoots.D(freq, times) = (1 - B^freq)^times, ascending coefficients.
inline std::vector<double> unit_roots_d(int freq, int times) {
    if (times == 0) return std::vector<double>{1.0};
    std::vector<double> root(static_cast<std::size_t>(freq) + 1, 0.0);
    root[0] = 1.0;
    root[static_cast<std::size_t>(freq)] = -1.0;
    std::vector<double> acc = root; // Pinned: acc = D(freq); acc = acc.times(D) x times-1.
    for (int k = 1; k < times; ++k) acc = poly_times(acc, root);
    return acc;
}

// SarimaUtility.differencingFilter: (1 - B)^d convolved with (1 - B^P)^D.
inline std::vector<double> differencing_filter(int period, int d, int seasonal_d) {
    std::vector<double> poly;
    if (d > 0) poly = unit_roots_d(1, d);
    if (seasonal_d > 0) {
        const auto seasonal = unit_roots_d(period, seasonal_d);
        poly = poly.empty() ? seasonal : poly_times(poly, seasonal);
    }
    if (poly.empty()) poly.push_back(1.0);
    return poly;
}

// BackFilter.apply: dest[i] = sum_k coef[k] * src[i + k], k ascending.
inline void apply_backfilter(const std::vector<double>& coef, const std::vector<double>& src,
    std::size_t src_offset, std::vector<double>& dest) {
    for (std::size_t i = 0; i < dest.size(); ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k < coef.size(); ++k) {
            s += coef[k] * src[src_offset + i + k];
        }
        dest[i] = s;
    }
}

// SarimaModel.getStationaryRepresentation() polynomials from the kalos
// [phi, theta, bphi, btheta] blocks: AR sign flip, regular x seasonal
// convolution and the trailing-zero strip that normalizes the degree.
// Moved verbatim from src/outlier_scan_engine.cpp (item 21-H).
inline void expand_sarima_polynomials(const std::vector<double>& phi,
    const std::vector<double>& theta, const std::vector<double>& bphi,
    const std::vector<double>& btheta, int period, std::vector<double>& ar_full,
    std::vector<double>& ma_full) {
    std::vector<double> ar_reg = phi;
    std::vector<double> sar_reg = bphi;
    for (auto& v : ar_reg) v = -v;
    for (auto& v : sar_reg) v = -v;
    ar_full = std::vector<double>{1.0};
    ma_full = std::vector<double>{1.0};
    {
        std::vector<double> reg_ar{1.0};
        reg_ar.insert(reg_ar.end(), ar_reg.begin(), ar_reg.end());
        std::vector<double> seas_ar{1.0};
        for (std::size_t i = 0; i < sar_reg.size(); ++i) {
            seas_ar.resize((i + 1) * static_cast<std::size_t>(period) + 1, 0.0);
            seas_ar[(i + 1) * static_cast<std::size_t>(period)] = sar_reg[i];
        }
        ar_full = poly_times(reg_ar, seas_ar);
        std::vector<double> reg_ma{1.0};
        reg_ma.insert(reg_ma.end(), theta.begin(), theta.end());
        std::vector<double> seas_ma{1.0};
        for (std::size_t i = 0; i < btheta.size(); ++i) {
            seas_ma.resize((i + 1) * static_cast<std::size_t>(period) + 1, 0.0);
            seas_ma[(i + 1) * static_cast<std::size_t>(period)] = btheta[i];
        }
        ma_full = poly_times(reg_ma, seas_ma);
        // Polynomial degree normalization: the pinned Polynomial strips
        // trailing zeros, so caller-passed zero AR/MA coefficients (a block
        // of size b with all-zero coefficients) shrink the expanded degree,
        // exactly as the pinned SarimaModel.getStationaryRepresentation().
        const auto strip = [](std::vector<double>& poly) {
            while (poly.size() > 1 && poly.back() == 0.0) poly.pop_back();
        };
        strip(ar_full);
        strip(ma_full);
    }
}


}  // namespace mlj
}  // namespace kalos

#endif
