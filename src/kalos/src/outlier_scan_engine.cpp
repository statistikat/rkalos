// Self-contained scan engine for the 21-C outlier decision unit: the pinned
// ExactSingleOutlierDetector.calc() pipeline (differenced model, Ansley
// whitening, robust pivoting-QR initialization, Mad(50) scale over the ML
// residuals, per-factory T/coefficient tables and the searchMax winner).
// Every loop, comparison and constant is the literal transcription of the
// jdplus 3.9.0 bytecode qualified in docs/contracts/outlier-scoring-contract.md.
// The engine never allocates across scans and uses no BLAS/LAPACK.

#include "kalos/outlier_scoring.hpp"

#include "outlier_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace kalos {

namespace {

using outlier_detail::dcmpl;
using outlier_detail::dcmpg;
using outlier_detail::java_max;
using outlier_detail::robust_lu_solve;

constexpr double kEpsilon = 1.1102230246251565e-16; // Constants.getEpsilon().

// ---- Polynomial / filter helpers ----------------------------------------------

// Polynomial.times: r[i + j] += a[i] * b[j], i then j ascending.
std::vector<double> poly_times(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> r(a.size() + b.size() - 1, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] += a[i] * b[j];
        }
    }
    return r;
}

// UnitRoots.D(freq, times) = (1 - B^freq)^times, ascending coefficients.
std::vector<double> unit_roots_d(int freq, int times) {
    if (times == 0) return std::vector<double>{1.0};
    std::vector<double> root(static_cast<std::size_t>(freq) + 1, 0.0);
    root[0] = 1.0;
    root[static_cast<std::size_t>(freq)] = -1.0;
    std::vector<double> acc = root; // Pinned: acc = D(freq); acc = acc.times(D) x times-1.
    for (int k = 1; k < times; ++k) acc = poly_times(acc, root);
    return acc;
}

// SarimaUtility.differencingFilter: (1 - B)^d convolved with (1 - B^P)^D.
std::vector<double> differencing_filter(int period, int d, int seasonal_d) {
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
void apply_backfilter(const std::vector<double>& coef, const std::vector<double>& src,
    std::size_t src_offset, std::vector<double>& dest) {
    for (std::size_t i = 0; i < dest.size(); ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k < coef.size(); ++k) {
            s += coef[k] * src[src_offset + i + k];
        }
        dest[i] = s;
    }
}

double dot(const double* a, const double* b, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

double ssq(const double* a, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * a[i];
    return s;
}

// RationalFunction long division: expansion of num / den, pinned recursion
// order (the cache-block split in the bytecode is bookkeeping; every branch
// computes this recursion). den[0] must be nonzero (monic polynomials).
std::vector<double> rational_coefficients(
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
void symmetric_lcholesky(std::vector<double>& a, std::size_t n);
void solve_lx(const std::vector<double>& a, std::size_t n, double* x, std::size_t xlen);
void solve_lx(const std::vector<double>& a, std::size_t n, std::vector<double>& x);
void solve_xl(const std::vector<double>& a, std::size_t n, std::vector<double>& x);
void solve_xu(const std::vector<double>& u, std::size_t n, std::vector<double>& x);
void solve_ux(const std::vector<double>& u, std::size_t n, std::vector<double>& x);
std::vector<double> mat_mul(const std::vector<double>& a, const std::vector<double>& b,
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

private:
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
void symmetric_lcholesky(std::vector<double>& a, std::size_t n) {
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
void solve_lx(const std::vector<double>& a, std::size_t n, double* x, std::size_t xlen) {
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

void solve_lx(const std::vector<double>& a, std::size_t n, std::vector<double>& x) {
    solve_lx(a, n, x.data(), x.size());
}

// LowerTriangularMatrix.solvexL (solveLtx), tol 0, pinned loop: i from n-1
// down; s = x[i] - sum_{j>i} L[j][i] * x[j] (column i below the diag, pinned
// product order storage[p] * a[j]); x[i] = s / diag with the same singular
// mapping as solveLx on the reduced s.
void solve_xl(const std::vector<double>& a, std::size_t n, std::vector<double>& x) {
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
void solve_xu(const std::vector<double>& u, std::size_t n, std::vector<double>& x) {
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
void solve_ux(const std::vector<double>& u, std::size_t n, std::vector<double>& x) {
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
std::vector<double> mat_mul(const std::vector<double>& a, const std::vector<double>& b,
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

// ---- AnsleyFilter ---------------------------------------------------------------

class Ansley {
public:
    // prepare on the stationary arma (raw polynomials, leading 1) and length n.
    void prepare(std::vector<double> ar, std::vector<double> ma, std::size_t n, Acvf& acvf) {
        ar_ = std::move(ar);
        ma_ = std::move(ma);
        n_ = n;
        p_ = ar_.size() - 1;
        q_ = ma_.size() - 1;
        var_ = 1.0; // SarimaModel.getInnovationVariance()
        trivial_ = (p_ == 0 && q_ == 0);
        if (trivial_) return;
        mm_ = std::max(p_, q_ + 1);
        std::vector<double> mabuf(mm_, 0.0);
        std::vector<double> acvf_vals;
        if (p_ > 0) {
            acvf_vals = acvf.values(mm_);
            const auto psi = rational_coefficients(ma_, ar_, q_);
            for (std::size_t i = 1; i <= q_; ++i) {
                double t = ma_[i];
                for (std::size_t j = i + 1; j <= q_; ++j) t += ma_[j] * psi[j - i];
                mabuf[i] = t * var_;
            }
        }
        l_.assign(mm_ * n_, 0.0);
        for (std::size_t c = 0; c < p_; ++c) {
            double* col = &l_[c * mm_];
            for (std::size_t i = 0; i < p_ - c; ++i) col[i] = acvf_vals[i];
            for (std::size_t i = p_ - c; i < mm_; ++i) col[i] = mabuf[i];
        }
        // W = extract(0, q+1, p, n-p): rows 0..q, cols p..n-1 get the
        // MA-convolution weight of the row (overwrite; the block was zero).
        for (std::size_t r = 0; r <= q_ && r < mm_; ++r) {
            double w = 0.0;
            for (std::size_t i = 0; i + r <= q_; ++i) w += ma_[i] * ma_[i + r];
            if (var_ != 1.0) w *= var_; // dcmpl(var, 1) != 0; var == 1 skips.
            if (dcmpl(w, 0.0) == 0) continue;
            for (std::size_t c = p_; c < n_; ++c) l_[c * mm_ + r] = w;
        }
        lcholesky();
    }

    std::vector<double> apply(const std::vector<double>& y) const {
        std::vector<double> a = y;
        if (trivial_) {
            if (var_ != 1.0) {
                const double s = std::sqrt(var_);
                for (auto& v : a) v /= s;
            }
            return a;
        }
        for (std::size_t i = n_ - 1; i >= p_; --i) {
            double s = 0.0;
            for (std::size_t k = 1; k <= p_; ++k) s += ar_[k] * a[i - k];
            a[i] += s;
        }
        rsolve(a);
        return a;
    }

private:
    // Pinned lcholesky() on the banded storage (column c = rows c..c+mm-1 of
    // the implicit lower-triangular factor): band recursion of section
    // F_AnsleyFilter.lcholesky.
    void lcholesky() {
        if (mm_ == 1) {
            for (auto& v : l_) {
                if (dcmpg(v, 0.0) <= 0) throw std::runtime_error("m_err_chol");
                v = std::sqrt(v);
            }
            return;
        }
        const std::size_t rr = mm_ - 1;
        const std::size_t rr2 = rr * rr;
        long long ptr = 0;
        for (std::size_t i = 0; i < n_; ++i) {
            double diag = l_[static_cast<std::size_t>(ptr)];
            long long lim = ptr - static_cast<long long>(rr2);
            if (lim < 0) lim = 0;
            long long k = ptr - static_cast<long long>(rr);
            while (k >= lim) {
                const double v = l_[static_cast<std::size_t>(k)];
                if (dcmpl(v, 0.0) != 0) diag -= v * v;
                k -= static_cast<long long>(rr);
            }
            if (dcmpg(diag, 0.0) <= 0) throw std::runtime_error("m_err_chol");
            diag = std::sqrt(diag);
            l_[static_cast<std::size_t>(ptr)] = diag;
            k = ptr - static_cast<long long>(rr);
            long long jj = static_cast<long long>(i) + static_cast<long long>(mm_) - 1;
            while (k >= lim) {
                const double v = l_[static_cast<std::size_t>(k)];
                if (dcmpl(v, 0.0) != 0) {
                    // Pinned bound: while (s < k + cnt), i.e. cnt - 1 updates
                    // (the pinned loop starts at s = k + 1).
                    const std::size_t cnt =
                        std::min(static_cast<std::size_t>(jj), n_) - i;
                    long long d = ptr + 1;
                    long long s = k + 1;
                    for (std::size_t t = 0; t + 1 < cnt; ++t, ++d, ++s) {
                        l_[static_cast<std::size_t>(d)] -=
                            v * l_[static_cast<std::size_t>(s)];
                    }
                }
                k -= static_cast<long long>(rr);
                --jj;
            }
            const long long end = static_cast<long long>(mm_) * (static_cast<long long>(i) + 1);
            for (long long d = ptr + 1; d < end; ++d) {
                l_[static_cast<std::size_t>(d)] /= diag;
            }
            ptr += static_cast<long long>(mm_);
        }
    }

    // rsolve: forward substitution with the banded factor (first nonzero).
    void rsolve(std::vector<double>& x) const {
        std::size_t f = 0;
        while (f < n_ && x[f] == 0.0) ++f;
        std::size_t k = f * mm_;
        for (std::size_t i = f; i < n_; ++i, k += mm_) {
            const double v = x[i] / l_[k];
            const std::size_t cnt = std::min(mm_, n_ - i);
            for (std::size_t t = 1; t < cnt; ++t) x[i + t] -= v * l_[k + t];
            x[i] = v;
        }
    }

    std::vector<double> ar_;
    std::vector<double> ma_;
    std::vector<double> l_;
    std::size_t n_ = 0;
    std::size_t p_ = 0;
    std::size_t q_ = 0;
    std::size_t mm_ = 0;
    double var_ = 1.0;
    bool trivial_ = false;
};

// ---- Mad(50) ---------------------------------------------------------------------

double mad50(const std::vector<double>& x) {
    std::vector<double> a(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) a[i] = std::fabs(x[i]);
    // Arrays.sort places NaN last; std::sort with NaN is UB — map NaN to +inf,
    // the same tail position. (Qualified fixtures are NaN-free; tiered.)
    for (auto& v : a) {
        if (std::isnan(v)) v = std::numeric_limits<double>::infinity();
    }
    std::sort(a.begin(), a.end());
    const double q = (static_cast<double>(a.size()) + 1.0) * 50.0 / 100.0;
    const int i = static_cast<int>(q);
    const double f = q - static_cast<double>(i);
    double m;
    if (dcmpg(f, 1.0e-9) < 0) {
        m = a[static_cast<std::size_t>(i) - 1];
    } else {
        m = a[static_cast<std::size_t>(i) - 1] * (1.0 - f) + a[static_cast<std::size_t>(i)] * f;
    }
    return m / standard_normal_lower_quantile(0.5 + 0.005 * 50.0); // Phi^-1(0.75)
}

// ---- HouseholderWithPivoting (full pivoting, nfixed = 0) --------------------------

struct Qr {
    std::size_t m = 0;
    std::size_t n = 0;
    std::size_t kmax = 0;                 // min(m, n)
    std::vector<double> qr;               // m x n storage, column-major
    std::vector<double> beta;             // kmax entries
    std::vector<int> pivot;               // kmax entries
};

double fast_norm2(const double* p, std::size_t n) {
    return std::sqrt(ssq(p, n));
}

// Safe dnrm2-style norm (scaling by max |x|); the pinned DataPointer.norm2.
double safe_norm2(const double* p, std::size_t n) {
    double scale = 0.0;
    for (std::size_t i = 0; i < n; ++i) scale = std::max(scale, std::fabs(p[i]));
    if (std::isnan(scale)) return std::nan("");
    if (scale == 0.0) return 0.0;
    if (scale == std::numeric_limits<double>::infinity()) return scale;
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double t = p[i] / scale;
        s += t * t;
    }
    return scale * std::sqrt(s);
}

struct Reflector {
    double beta = 0.0;
    double alpha = 0.0;
};

// Reflector.set(ptr, len).larfg() on the storage window [ptr, ptr + len).
// Pinned: norm = DataPointer(x).norm2(n) over the FULL window (x0 included),
// the underflow scaling multiplies the tail (n-1 entries) plus x0 and always
// recomputes the full norm on either loop exit, then beta = norm/(norm + x0)
// (0 only when norm == 0 exactly), v /= norm, x0 <- 1 + x0/norm, alpha = -norm.
Reflector larfg(std::vector<double>& storage, std::size_t ptr, std::size_t len) {
    Reflector refl;
    if (len == 0) return refl;
    double x0 = storage[ptr];
    double* v = storage.data() + ptr + 1;
    const std::size_t vn = len - 1;
    double norm = safe_norm2(storage.data() + ptr, len);
    const double safmin = std::numeric_limits<double>::min();
    const double underlim = safmin / kEpsilon;
    int scales = 0;
    if (dcmpl(norm, 0.0) != 0 && dcmpg(norm, underlim) < 0) {
        const double scale = 1.0 / underlim;
        while (true) {
            for (std::size_t i = 0; i < vn; ++i) v[i] *= scale;
            x0 *= scale;
            norm *= scale;
            if (dcmpg(norm, underlim) >= 0) break;
            if (++scales >= 4) break;
        }
        norm = safe_norm2(storage.data() + ptr, len); // Pinned recompute at both exits.
    }
    if (!(dcmpg(x0, 0.0) >= 0)) norm = -norm;
    for (int i = 0; i < scales; ++i) norm *= underlim;
    refl.beta = dcmpl(norm, 0.0) == 0 ? 0.0 : norm / (norm + x0);
    for (std::size_t i = 0; i < vn; ++i) v[i] /= norm;
    storage[ptr] = 1.0 + x0 / norm;
    refl.alpha = -norm;
    return refl;
}

// Reflector.lapply over the trailing window: dest column = rows [ptr_start..).
// dest points at the reflector's diag row for each column; v = storage[ptr+1..].
void lapply(const std::vector<double>& storage, std::size_t ptr, std::size_t len,
    const Reflector& refl, std::vector<double>& matrix, std::size_t matrix_rows,
    std::size_t first_col, std::size_t last_col, std::size_t row0) {
    const double* v = storage.data() + ptr + 1;
    for (std::size_t col = first_col; col <= last_col; ++col) {
        double* dest = matrix.data() + col * matrix_rows + row0;
        double tau = dest[0] / refl.beta;
        for (std::size_t i = 1; i < len; ++i) tau += dest[i] * v[i - 1];
        if (dcmpl(tau, 0.0) != 0) {
            dest[0] -= tau;
            const double t = -tau * refl.beta;
            for (std::size_t i = 1; i < len; ++i) dest[i] += t * v[i - 1];
        }
    }
}

// QRLeastSquaresSolver.robustLeastSquares(y, X): full-pivoting decompose,
// tol 1e-9 least squares. Returns rank, b (original order), pivot, R
// (kmax x kmax upper) and qy (Q^T y).
struct QrSolution {
    std::size_t rank;
    std::vector<double> b;
    std::vector<double> u; // rank-relevant rawR, kmax x kmax column-major
    std::vector<double> qy;
    std::vector<int> pivot;
    double ssq_err;
};

QrSolution robust_least_squares(
    const std::vector<double>& y, const std::vector<double>& x, std::size_t m, std::size_t n) {
    Qr q;
    q.m = m;
    q.n = n;
    q.qr = x; // column-major m x n (deep copy of the design)
    q.beta.assign(n, 0.0);
    q.kmax = std::min(m, n);
    q.pivot.assign(q.kmax, 0);
    std::vector<double> cnorm(q.kmax, 0.0);
    std::vector<double> wnorm(q.kmax, 0.0);
    for (std::size_t i = 0; i < q.kmax; ++i) q.pivot[i] = static_cast<int>(i);
    for (std::size_t j = 0; j < q.kmax; ++j) {
        const double nv = fast_norm2(&q.qr[j * m], m);
        cnorm[j] = nv;
        wnorm[j] = nv;
    }
    const double lim = std::sqrt(kEpsilon);
    std::size_t ptr = 0;
    for (std::size_t k = 0; k < q.kmax; ++k) {
        std::size_t sel = k;
        for (std::size_t j = k + 1; j < q.kmax; ++j) {
            if (dcmpl(wnorm[j], wnorm[sel]) > 0) sel = j;
        }
        if (sel != k) {
            for (std::size_t i = 0; i < m; ++i) std::swap(q.qr[k * m + i], q.qr[sel * m + i]);
            std::swap(q.pivot[k], q.pivot[sel]);
            wnorm[sel] = wnorm[k];
            cnorm[sel] = cnorm[k];
        }
        const Reflector refl = larfg(q.qr, ptr, m - k);
        const std::size_t rem_cols = n > k + 1 ? n - k - 1 : 0;
        if (rem_cols > 0 && dcmpl(refl.beta, 0.0) != 0) {
            // Pinned Reflector.lapply: each trailing window column starts AT
            // the reflector's first row (bhshrink only drops the left
            // column), so dest[0] is row k and the vector runs to row m-1.
            lapply(q.qr, ptr, m - k, refl, q.qr, m, k + 1, n - 1, k);
        }
        q.beta[k] = refl.beta;
        q.qr[ptr] = refl.alpha;
        std::size_t col_start = ptr + m; // storage pos of R(k, k+1)
        for (std::size_t jp = k + 1; jp < q.kmax; ++jp) {
            if (dcmpl(wnorm[jp], 0.0) != 0) {
                double t = q.qr[col_start] / wnorm[jp];
                t = java_max(1.0 - t * t, 0.0); // Pinned max(1 - t*t, 0.0); NaN propagates.
                const double ratio = wnorm[jp] / cnorm[jp];
                const double upd = t * ratio * ratio;
                if (dcmpg(upd, lim) > 0) { // Pinned ifgt after dcmpg: NaN takes
                    wnorm[jp] *= std::sqrt(t); // the cheap multiply.
                } else {
                    const double nv = fast_norm2(q.qr.data() + col_start + 1, m - k - 1);
                    wnorm[jp] = nv;
                    cnorm[jp] = nv;
                }
            }
            col_start += m;
        }
        ptr += m + 1;
    }

    // leastSquares(qr, y, tol = 1e-9).
    QrSolution sol;
    sol.qy = y;
    // applyQt: forward reflector order, skip beta == 0 cells.
    {
        std::size_t pos = 0;
        for (std::size_t k = 0; k < n; ++k) {
            if (dcmpl(q.beta[k], 0.0) != 0) {
                double tau = sol.qy[k] / q.beta[k];
                for (std::size_t i = k + 1; i < m; ++i) {
                    tau += sol.qy[i] * q.qr[pos + (i - k)];
                }
                if (dcmpl(tau, 0.0) != 0) {
                    sol.qy[k] -= tau;
                    const double t = -tau * q.beta[k];
                    for (std::size_t i = k + 1; i < m; ++i) {
                        sol.qy[i] += t * q.qr[pos + (i - k)];
                    }
                }
            }
            pos += m + 1;
        }
    }
    // rawR: extract(0, n, 0, n) upper part.
    std::vector<double> u(n * n, 0.0);
    for (std::size_t c = 0; c < n; ++c) {
        for (std::size_t r = 0; r <= c && r < m; ++r) u[c * n + r] = q.qr[c * m + r];
        for (std::size_t r = c + 1; r < n && r < m; ++r) u[c * n + r] = 0.0;
    }
    sol.u = u;
    // fastRank: count |diag| > 1e-9 (dcmpl: a NaN diagonal entry does not count).
    sol.rank = 0;
    const std::size_t r = std::min(q.kmax, n);
    for (std::size_t k = 0; k < r; ++k) {
        if (dcmpl(std::fabs(u[k * n + k]), 1.0e-9) > 0) ++sol.rank;
    }
    sol.ssq_err = ssq(sol.qy.data() + r, m - r);
    solve_ux(u, r, sol.qy); // tol 0, leading r entries
    sol.b.assign(n, 0.0);
    for (std::size_t k = 0; k < r; ++k) sol.b[static_cast<std::size_t>(q.pivot[k])] = sol.qy[k];
    sol.pivot = q.pivot;
    return sol;
}

} // namespace

// ---- scan_outlier_candidates -------------------------------------------------------

OutlierScanResult scan_outlier_candidates(
    const OutlierScanRequest& request, const Interrupt& interrupt) {
    OutlierScanResult result;
    result.available = false;
    const std::size_t rows = request.series.values.size();
    result.rows = rows;
    const std::size_t nf = request.types.size();
    if (request.types.empty()) {
        throw std::invalid_argument("the scan requires at least one candidate type");
    }
    for (std::size_t i = 0; i < nf; ++i) {
        for (std::size_t j = i + 1; j < nf; ++j) {
            if (request.types[i] == request.types[j]) {
                throw std::invalid_argument("duplicate outlier candidate type");
            }
        }
    }
    if (request.period <= 1 && std::find(request.types.begin(), request.types.end(),
                                    OutlierType::SeasonalOutlier) != request.types.end()) {
        throw std::invalid_argument("the seasonal factory requires period > 1");
    }
    if (std::any_of(request.series.missing.begin(), request.series.missing.end(),
            [](std::uint8_t flag) { return flag != 0; })) {
        throw std::invalid_argument("the scan requires a complete response");
    }
    if (request.d < 0 || request.seasonal_d < 0 || request.period < 1) {
        throw std::invalid_argument("invalid differencing orders");
    }
    const std::size_t nx = request.design.columns;
    if (request.design.values.size() != nx * rows
        || request.design.coefficients.size() != nx || request.design.fixed.size() != nx) {
        throw std::invalid_argument("inconsistent regression design");
    }
    if (std::any_of(request.design.fixed.begin(), request.design.fixed.end(),
            [](std::uint8_t flag) { return flag != 0; })) {
        throw std::invalid_argument("the scan requires a fully free design");
    }
    if (std::any_of(request.design.values.begin(), request.design.values.end(),
            [](double v) { return !std::isfinite(v); })) {
        throw std::invalid_argument("the scan design must carry finite values");
    }
    if (request.excluded.size() != rows * nf) {
        throw std::invalid_argument("the exclusion table must cover rows x types");
    }
    if (request.upper_bound < request.lower_bound || request.upper_bound > rows) {
        throw std::invalid_argument("scan bounds must satisfy l <= r <= length");
    }
    if (std::isnan(request.transitory_rate)) {
        throw std::invalid_argument("the transitory rate must not be NaN");
    }
    result.t_table.assign(rows * nf, 0.0);
    result.coefficient_table.assign(rows * nf, 0.0);
    result.excluded_after = request.excluded;

    if (request.upper_bound <= request.lower_bound) {
        result.status = "empty_bounds"; // Pinned calc(): ubound <= lbound -> false.
        return result;
    }
    if (interrupt) interrupt();

    // Kalos blocks -> pinned raw polynomials (AR sign flip, MA identity).
    std::vector<double> ar_reg = request.arma[0];
    std::vector<double> sar_reg = request.arma[2];
    for (auto& v : ar_reg) v = -v;
    for (auto& v : sar_reg) v = -v;
    std::vector<double> ar_full{1.0};
    std::vector<double> ma_full{1.0};
    {
        std::vector<double> reg_ar{1.0};
        reg_ar.insert(reg_ar.end(), ar_reg.begin(), ar_reg.end());
        std::vector<double> seas_ar{1.0};
        for (std::size_t i = 0; i < sar_reg.size(); ++i) {
            seas_ar.resize((i + 1) * static_cast<std::size_t>(request.period) + 1, 0.0);
            seas_ar[(i + 1) * static_cast<std::size_t>(request.period)] = sar_reg[i];
        }
        ar_full = poly_times(reg_ar, seas_ar);
        std::vector<double> reg_ma{1.0};
        reg_ma.insert(reg_ma.end(), request.arma[1].begin(), request.arma[1].end());
        std::vector<double> seas_ma{1.0};
        for (std::size_t i = 0; i < request.arma[3].size(); ++i) {
            seas_ma.resize((i + 1) * static_cast<std::size_t>(request.period) + 1, 0.0);
            seas_ma[(i + 1) * static_cast<std::size_t>(request.period)] = request.arma[3][i];
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
    const std::vector<double> ns =
        differencing_filter(request.period, request.d, request.seasonal_d);
    const std::size_t w = ns.size() - 1;
    if (rows <= w) {
        result.status = "too_short"; // Pinned "Not enough observations" -> calc false.
        return result;
    }
    const std::size_t m = rows - w;

    try {
        Acvf acvf(ma_full, ar_full);
        const auto psi = rational_coefficients(ma_full, ar_full, ma_full.size() - 1);
        Ansley ansley;
        ansley.prepare(ar_full, ma_full, m, acvf);
        // ML-residual machine over the differenced length.
        const ModifiedLjungBox mlj(ar_full, ma_full, acvf, psi, m);

        std::vector<double> ydiff(m);
        apply_backfilter(ns, request.series.values, 0, ydiff);

        std::vector<double> yl = ansley.apply(ydiff);

        const std::size_t allowed_types = nf;
        auto zone_allows = [&](OutlierType type, std::size_t p) {
            const OutlierZone z = outlier_excluding_zone(type, request.period);
            return p >= z.start && p + z.end < rows;
        };

        std::vector<double> xdiff;   // m x nx column-major
        std::vector<double> xl;      // whitened design
        std::vector<double> u;       // rawR (nx x nx)
        std::vector<int> pivot;
        std::vector<double> z;
        double mad = 0.0;

        if (nx > 0) {
            xdiff.resize(m * nx);
            xl.resize(m * nx);
            for (std::size_t c = 0; c < nx; ++c) {
                std::vector<double> col(m);
                apply_backfilter(ns, request.design.values, c * rows, col);
                for (std::size_t i = 0; i < m; ++i) xdiff[c * m + i] = col[i];
                const auto wl = ansley.apply(col);
                for (std::size_t i = 0; i < m; ++i) xl[c * m + i] = wl[i];
            }
            const auto sol = robust_least_squares(yl, xl, m, nx);
            if (sol.rank != nx) {
                result.status = "rank_deficient";
                return result;
            }
            u = sol.u;
            pivot = sol.pivot;
            z.assign(nx, 0.0);
            for (std::size_t k = 0; k < nx; ++k) {
                z[k] = dot(yl.data(), xl.data() + static_cast<std::size_t>(pivot[k]) * m, m);
            }
            solve_xu(u, nx, z);
            std::vector<double> wvec = ydiff;
            for (std::size_t k = 0; k < nx; ++k) {
                const double b = sol.b[k];
                if (dcmpl(b, 0.0) == 0) continue;
                for (std::size_t i = 0; i < m; ++i) wvec[i] -= b * xdiff[k * m + i];
            }
            mad = mad50(mlj.residuals(wvec));
        } else {
            mad = mad50(mlj.residuals(ydiff));
        }
        if (interrupt) interrupt();

        // processOutlier for each factory.
        std::vector<double> buf(2 * rows, 0.0);
        std::vector<double> xlbuf(2 * rows - w, 0.0);
        std::vector<double> v(nx, 0.0);
        for (std::size_t j = 0; j < allowed_types; ++j) {
            buf.assign(2 * rows, 0.0);
            fill_outlier_column(request.types[j], rows, request.period,
                request.transitory_rate, buf);
            apply_backfilter(ns, buf, 0, xlbuf);
            for (std::size_t p = request.lower_bound; p < request.upper_bound; ++p) {
                const std::size_t cell = p + j * rows;
                if (!zone_allows(request.types[j], p) || request.excluded[cell] != 0) {
                    result.t_table[cell] = std::nan("");
                    result.coefficient_table[cell] = std::nan("");
                    continue;
                }
                // win[i] = xlbuf[rows - p + i], i < m.
                std::vector<double> win(m);
                const double* base = xlbuf.data() + (rows - p);
                std::copy(base, base + m, win.begin());
                const std::vector<double> rvec = ansley.apply(win);
                const double ssqr = ssq(rvec.data(), m);
                const double num0 = dot(rvec.data(), yl.data(), m);
                if (nx > 0) {
                    for (std::size_t k = 0; k < nx; ++k) {
                        v[k] = dot(rvec.data(),
                            xl.data() + static_cast<std::size_t>(pivot[k]) * m, m);
                    }
                    solve_xu(u, nx, v);
                    const double den = ssqr - ssq(v.data(), nx);
                    if (dcmpg(den, 0.0) > 0) {
                        const double num = num0 - dot(z.data(), v.data(), nx);
                        result.t_table[cell] = num / (std::sqrt(den) * mad);
                        result.coefficient_table[cell] = num / den;
                    } else {
                        result.excluded_after[cell] = 1; // Pinned exclude; T stays 0.
                    }
                } else {
                    if (dcmpg(ssqr, 0.0) > 0) {
                        result.t_table[cell] = num0 / (std::sqrt(ssqr) * mad);
                        result.coefficient_table[cell] = num0 / ssqr;
                    } else {
                        result.excluded_after[cell] = 1;
                        result.t_table[cell] = std::nan("");
                        result.coefficient_table[cell] = std::nan("");
                    }
                }
            }
            if (interrupt) interrupt();
        }

        // searchMax: column-major scan, weights all 1.0.
        double best = 0.0;
        long long idx = -1;
        long long k = 0;
        for (std::size_t j = 0; j < nf; ++j) {
            for (std::size_t i = 0; i < rows; ++i, ++k) {
                const double v2 = std::fabs(result.t_table[static_cast<std::size_t>(k)]);
                if (dcmpl(v2, best) > 0) {
                    idx = k;
                    best = v2;
                }
            }
            if (idx == -1) {
                // Pinned oMax escape: getMaxTStat would throw; the driver
                // treats the sentinel as "stop" (see the contract).
                result.available = true;
                result.status = "ok";
                result.max_position = OutlierScanResult::max_type_sentinel;
                result.max_type = OutlierScanResult::max_type_sentinel;
                result.max_t_stat = 0.0;
                return result;
            }
        }
        result.max_position = static_cast<std::size_t>(idx) % rows;
        result.max_type = static_cast<std::size_t>(idx) / rows;
        result.max_t_stat = result.t_table[static_cast<std::size_t>(idx)];
        result.available = true;
        result.status = "ok";
        return result;
    } catch (const std::invalid_argument&) {
        throw; // Request-level problems stay errors (fail-closed).
    } catch (const std::exception&) {
        // Pinned calc() catches every exception and returns false.
        result.status = "whitening_failed";
        result.available = false;
        return result;
    }
}

} // namespace kalos
