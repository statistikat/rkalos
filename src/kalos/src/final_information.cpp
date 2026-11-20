// Item 21-H.0: transcription of the final-information producer tail
// (RegArmaSsqFunction.ssqEvaluate -> LinearModel.calcResiduals ->
// RegArmaFunction.evaluate -> NumericalDerivatives(point, false, false)
// -> RegArmaProcessor.compute scalings -> SymmetricMatrix.inverse),
// docs/contracts/final-information-contract.md. Comparison flavours
// (dcmp/dcmpl/dcmpg) and operation order follow the pinned bytecode
// literally; the inverse kernels keep the JVM column-major storage layout
// (row stride 1, column stride n) so every transcribed index expression
// matches the original bytecode offsets. The public matrices are row-major
// flats, converted at the boundary.
#include "kalos/final_information.hpp"

#include "kalos/detail/filter_utility.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace kalos {

double final_information_step() {
    // SarimaMapping2: STEP = Math.sqrt(2.220446E-16) (static field read by
    // the builder default and by epsilon()).
    return std::sqrt(2.220446E-16);
}

namespace {

// Java double comparison flavours.
int dcmpg(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return 1;
    }
    return a > b ? 1 : (a < b ? -1 : 0);
}

int dcmpl(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return -1;
    }
    return a > b ? 1 : (a < b ? -1 : 0);
}

// java.lang.Math.max.
double java_max(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (a > b) {
        return a;
    }
    if (b > a) {
        return b;
    }
    // equal magnitudes (or both zeros): Math.max(-0.0, +0.0) == +0.0.
    if (std::signbit(a) && !std::signbit(b)) {
        return b;
    }
    return a;
}

// MatrixException of the pinned stack, raised by the Cholesky and
// triangular-solve kernels and caught by symmetric_inverse.
struct matrix_exception {
};

// lapack.Cholesky.lcholesky(m, tol) on the column-major n x n storage m.
// ip walks the diagonal, ir walks row i from column 0 (k += column stride),
// col_end ends the current column.
void lcholesky(std::vector<double>& m, int n, double tol) {
    if (n == 0) {
        return;
    }
    if (n == 1) {
        const double v = m[0];
        if (dcmpg(v, -tol) <= 0) { // !(v > -tol); NaN keeps going
            throw matrix_exception{};
        }
        m[0] = dcmpg(v, tol) < 0 ? 0.0 : std::sqrt(v);
        return;
    }
    for (int i = 0; i < n; ++i) {
        const int ir = i;
        const int ip = i * (n + 1);
        const int col_end = (i + 1) * n;
        double s = m[ip];
        for (int k = ir; k != ip; k += n) { // row i, columns 0..i-1
            const double v = m[k];
            s -= v * v;
        }
        if (dcmpg(s, -tol) <= 0) { // !(s > -tol); NaN keeps going
            throw matrix_exception{};
        }
        if (dcmpg(s, tol) < 0) { // s < tol (unreachable at tol = 0)
            m[ip] = 0.0;
            for (int k = ir; k != ip; k += n) {
                const double v = m[k];
                if (dcmpl(v, 0.0) != 0) { // v != 0, NaN included
                    int q = k + 1;
                    int r = ip + 1;
                    while (r < col_end) {
                        m[r] -= v * m[q];
                        ++q;
                        ++r;
                    }
                }
            }
            for (int r = ip + 1; r < col_end; ++r) {
                if (dcmpl(std::fabs(m[r]), tol) > 0) {
                    throw matrix_exception{};
                }
                m[r] = 0.0;
            }
            continue;
        }
        s = std::sqrt(s);
        m[ip] = s;
        for (int k = ir; k != ip; k += n) { // A[i][j] != 0 updates column i
            const double v = m[k];
            if (dcmpl(v, 0.0) != 0) {
                int q = k + 1;
                int r = ip + 1;
                while (r < col_end) {
                    m[r] -= v * m[q];
                    ++q;
                    ++r;
                }
            }
        }
        const double inv = 1.0 / s; // multiply by 1/L[i][i], never divide
        for (int r = ip + 1; r < col_end; ++r) {
            m[r] *= inv;
        }
    }
}

// LowerTriangularMatrix.inverse(l): identity(n) then solveLX(l, ., 0).
std::vector<double> lower_inverse(const std::vector<double>& l, int n) {
    for (int i = 0; i < n; ++i) {
        if (dcmpl(l[i * (n + 1)], 0.0) == 0) { // diagonal anyMatch(d == 0)
            throw matrix_exception{};
        }
    }
    std::vector<double> x(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        x[i * (n + 1)] = 1.0; // FastMatrix.identity
    }
    for (int off = 0; off < n * n; off += n) { // solveLX per column
        int dcur = 0;
        for (int k = off; k != off + n; ++k, dcur += n + 1) {
            const double xv = x[k];
            const double d = l[dcur];
            if (dcmpg(std::fabs(d), 0.0) <= 0) { // |d| <= 0 (NaN divides)
                if (dcmpl(std::fabs(xv), 0.0) >= 0) { // |xv| >= 0, NaN excluded
                    throw matrix_exception{};
                }
                x[k] = 0.0;
            } else {
                const double r = xv / d;
                x[k] = r;
                int jj = k + 1;
                int dd = dcur + 1;
                while (jj != off + n) {
                    x[jj] -= r * l[dd]; // L below the diagonal in column k-off
                    ++jj;
                    ++dd;
                }
            }
        }
    }
    return x;
}

// SymmetricMatrix.LtL(l) into fresh storage, then fromLower.
std::vector<double> lt_l(const std::vector<double>& l, int n) {
    std::vector<double> out(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {   // written column
        for (int j = i; j < n; ++j) { // written row
            double dot = 0.0;
            int a = i * n + j; // l[j..][i]
            int b = j * n + j; // l[j..][j]
            const int end = (i + 1) * n;
            while (a < end) {
                dot += l[a] * l[b];
                ++a;
                ++b;
            }
            out[i * n + j] = dot;
        }
    }
    if (n != 1) { // fromLower: mirror the computed half across the diagonal
        for (int col = 0; col < n; ++col) {
            for (int row = col + 1; row < n; ++row) {
                out[row * n + col] = out[col * n + row];
            }
        }
    }
    return out;
}

// CroutDoolittle.decompose(a, Constants.getEpsilon()) plus
// LUDecomposition.solve(identity), both on column-major storage; x receives
// the inverse.
void crout_inverse(const std::vector<double>& a_in, std::vector<double>& x_out, int n) {
    const double tol = 1.1102230246251565E-16; // Constants.getEpsilon()
    std::vector<double> a = a_in;              // decompose mutates its argument
    std::vector<int> pivot;                    // null until the first swap
    std::vector<double> v(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            v[j] = a[k * n + j]; // current column k
        }
        for (int i = 0; i < n; ++i) {
            const int mlim = i < k ? i : k;
            double s = 0.0;
            for (int m = 0; m < mlim; ++m) {
                s += a[i + m * n] * v[m];
            }
            v[i] -= s; // v keeps the final column entries for later rows
            a[k * n + i] = v[i];
        }
        int imax = k;
        for (int i = k + 1; i < n; ++i) {
            if (dcmpl(std::fabs(v[i]), std::fabs(v[imax])) > 0) {
                imax = i;
            }
        }
        if (imax != k) {
            if (pivot.empty()) {
                pivot.resize(static_cast<std::size_t>(n));
                for (int i = 0; i < n; ++i) {
                    pivot[i] = i;
                }
            }
            for (int t = 0; t < n; ++t) { // swap rows imax and k
                const double tmp = a[t * n + imax];
                a[t * n + imax] = a[t * n + k];
                a[t * n + k] = tmp;
            }
            const int ptmp = pivot[imax];
            pivot[imax] = pivot[k];
            pivot[k] = ptmp;
        }
        const double p = a[k * n + k];
        if (dcmpg(std::fabs(p), tol) < 0) { // NaN pivots keep going
            throw matrix_exception{};
        }
        for (int i = k + 1; i < n; ++i) {
            a[k * n + i] = a[k * n + i] / p;
        }
    }
    x_out.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        x_out[i * (n + 1)] = 1.0; // FastMatrix.identity(n)
    }
    if (!pivot.empty()) { // swap every rhs column through the pivot
        std::vector<double> tmp(static_cast<std::size_t>(n));
        for (int col = 0; col < n; ++col) {
            for (int t = 0; t < n; ++t) {
                tmp[t] = x_out[col * n + t];
            }
            for (int t = 0; t < n; ++t) {
                x_out[col * n + t] = tmp[pivot[t]];
            }
        }
    }
    for (int off = 0; off < n * n; off += n) { // solveL1X, no singularity guard
        for (int k = off; k != off + n; ++k) {
            const double xv = x_out[k];
            int r = k + 1;
            int dcur = (k - off) * (n + 1) + 1; // dd advances with i*(colInc+1)
            while (r != off + n) {
                x_out[r] -= xv * a[dcur];
                ++r;
                ++dcur;
            }
        }
    }
    for (int off = 0; off < n * n; off += n) { // solveUX (tol 0)
        int dcur = (n - 1) * (n + 1);
        for (int k = off + n - 1; k >= off; --k, dcur -= n + 1) {
            const double xv = x_out[k];
            const double d = a[dcur];
            if (dcmpg(std::fabs(d), 0.0) <= 0) {
                if (dcmpl(std::fabs(xv), 0.0) >= 0) {
                    throw matrix_exception{};
                }
                x_out[k] = 0.0;
            } else {
                const double r = xv / d;
                x_out[k] = r;
                int jj = k - 1;
                int dd = dcur - 1;
                while (jj >= off) {
                    x_out[jj] -= r * a[dd];
                    --jj;
                    --dd;
                }
            }
        }
    }
}

// SymmetricMatrix.inverse on the column-major storage of a symmetric
// matrix. a_row_major is the caller's information; the pinned code mutates
// neither (Cholesky runs on a deepClone, decompose on a toArray copy), and
// the unit works on an internal copy so result.information is untouched.
std::vector<double> symmetric_inverse(const std::vector<double>& a_row_major, int n) {
    std::vector<double> a(static_cast<std::size_t>(n) * n);
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            a[c * n + r] = a_row_major[r * n + c];
        }
    }
    std::vector<double> out; // column-major result
    const std::vector<double> original = a;
    try {
        // The pinned inverse deep-clones its argument for the Cholesky
        // branch and the MatrixException catch decomposes the ORIGINAL
        // matrix (bytecode aload_0), so the Crout input is never the
        // partially factorized clone.
        lcholesky(a, n, 0.0);
        const std::vector<double> li = lower_inverse(a, n);
        out = lt_l(li, n);
    } catch (const matrix_exception&) {
        crout_inverse(original, out, n); // rethrows on the pinned singular path
    }
    std::vector<double> row_major(static_cast<std::size_t>(n) * n);
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            row_major[r * n + c] = out[c * n + r];
        }
    }
    return row_major;
}

} // namespace

final_information_result final_information(const final_information_request& req) {
    final_information_result res;
    const int p = req.orders.p;
    const int bp = req.orders.bp;
    const int q = req.orders.q;
    const int bq = req.orders.bq;
    if (p < 0 || bp < 0 || q < 0 || bq < 0) {
        throw std::invalid_argument("negative sarima orders");
    }
    const int dim = p + bp + q + bq;
    if (static_cast<int>(req.point.size()) != dim) {
        throw std::invalid_argument("point must cover the declared orders");
    }
    const std::vector<double>& point = req.point;

    // SarimaMapping2.checkBoundaries: FilterUtility.checkRoots per block
    // (AR blocks against arLimit; MA blocks against maLimit only when
    // checkAll), first failing block returns false.
    auto check_boundaries = [&](const std::vector<double>& x) {
        int c = 0;
        if (p > 0 && !detail::check_roots(&x[c], static_cast<std::size_t>(p), req.ar_limit)) {
            return false;
        }
        c += p;
        if (bp > 0 && !detail::check_roots(&x[c], static_cast<std::size_t>(bp), req.ar_limit)) {
            return false;
        }
        c += bp;
        if (req.check_all) {
            if (q > 0 && !detail::check_roots(&x[c], static_cast<std::size_t>(q), req.ma_limit)) {
                return false;
            }
            c += q;
            if (bq > 0 && !detail::check_roots(&x[c], static_cast<std::size_t>(bq), req.ma_limit)) {
                return false;
            }
        }
        return true;
    };

    // NumericalDerivatives.newval: copy the point, add the offsets,
    // evaluate; any JVM failure surfaces as fx.
    auto newval = [&](int i, double d) {
        std::vector<double> x = point;
        x[i] += d;
        const std::optional<double> v = req.deviance(x);
        return v.has_value() ? *v : res.fx;
    };
    auto newval2 = [&](int i, int j, double di, double dj) {
        std::vector<double> x = point;
        x[i] += di;
        x[j] += dj;
        const std::optional<double> v = req.deviance(x);
        return v.has_value() ? *v : res.fx;
    };

    // The pinned tail evaluates the deviance at the final point first
    // (RegArmaSsqFunction.ssqEvaluate / Evaluation.getSsqE); a failure there
    // aborts the round before any stencil work.
    const std::optional<double> base = req.deviance(point);
    if (!base.has_value()) {
        res.status = final_information_status::final_point_rejected;
        return res;
    }
    res.fx = *base;

    // NumericalDerivatives(point, false, false): per coordinate the domain
    // epsilon, checkepsilon, then the eager forward node fp[i].
    std::vector<double> fp(static_cast<std::size_t>(dim));
    res.eps.resize(static_cast<std::size_t>(dim));
    for (int i = 0; i < dim; ++i) {
        // SarimaMapping2.epsilon(x, i), dcmpg against 0 (NaN takes the
        // negative branch, Math.max(1, NaN) = NaN).
        const double v = point[i];
        double e = dcmpg(v, 0.0) < 0 ? req.epsilon * java_max(1.0, -v)
                                     : (-req.epsilon) * java_max(1.0, v);
        // NumericalDerivatives.checkepsilon over a mutable copy of x.
        if (e != 0.0) { // dcmpl(e, 0) != 0; NaN enters
            std::vector<double> x = point;
            const double base_i = x[i];
            const double e_orig = e;
            x[i] = base_i + e;
            if (!check_boundaries(x)) {
                int c = 0;
                bool ok = false;
                for (;;) {
                    e = e / 2.0;
                    x[i] = base_i + e;
                    ++c;
                    if (c > 2) {
                        break;
                    }
                    if (check_boundaries(x)) {
                        ok = true;
                        break;
                    }
                }
                if (!ok) {
                    // flip; eps[i] still carries the original epsilon
                    e = -e_orig;
                    x[i] = base_i + e;
                    if (check_boundaries(x)) {
                        ok = true;
                    } else {
                        c = 0;
                        for (;;) {
                            e = e / 2.0;
                            x[i] = base_i + e;
                            ++c;
                            if (c > 2) {
                                break;
                            }
                            if (check_boundaries(x)) {
                                ok = true;
                                break;
                            }
                        }
                    }
                    if (!ok) {
                        e = 0.0;
                    }
                }
            }
        }
        res.eps[i] = e;
        fp[i] = newval(i, e);
    }

    // calcgrad with fm == null: forward difference over fx, with the
    // fp == fx zero guard (dcmpl).
    res.gradient.resize(static_cast<std::size_t>(dim));
    for (int i = 0; i < dim; ++i) {
        res.gradient[i] =
            fp[i] == res.fx ? 0.0 : (fp[i] - res.fx) / res.eps[i];
    }

    // calch: diagonal over hi = sqrt(|eps|), four-corner cross terms.
    std::vector<double> hi(static_cast<std::size_t>(dim));
    for (int i = 0; i < dim; ++i) {
        hi[i] = std::sqrt(std::fabs(res.eps[i]));
    }
    res.hessian.assign(static_cast<std::size_t>(dim) * dim, 0.0);
    for (int i = 0; i < dim; ++i) {
        const double h = hi[i];
        const double a = newval(i, h);
        const double b = newval(i, -h);
        const double r = (a - 2.0 * res.fx) + b;
        if (r != 0.0 && res.eps[i] != 0.0) {
            res.hessian[static_cast<std::size_t>(i) * dim + i] = r / (h * h);
        }
    }
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < i; ++j) {
            const double a = hi[i];
            const double b = hi[j];
            double s = newval2(i, j, a / 2.0, b / 2.0);
            s += newval2(i, j, (-a) / 2.0, (-b) / 2.0);
            s -= newval2(i, j, a / 2.0, (-b) / 2.0);
            s -= newval2(i, j, (-a) / 2.0, b / 2.0);
            if (s != 0.0 && a != 0.0 && b != 0.0) {
                const double val = s / (a * b);
                res.hessian[static_cast<std::size_t>(i) * dim + j] = val;
                res.hessian[static_cast<std::size_t>(j) * dim + i] = val;
            }
        }
    }

    // RegArmaProcessor.compute scalings (0.5 nobs / fx on the curvature,
    // (-0.5 nobs) / fx on the gradient).
    const double nobs_d = static_cast<double>(req.nobs);
    const double scale_i = (0.5 * nobs_d) / res.fx;
    const double scale_g = ((-0.5) * nobs_d) / res.fx;
    res.information.resize(static_cast<std::size_t>(dim) * dim);
    for (std::size_t t = 0; t < res.information.size(); ++t) {
        res.information[t] = res.hessian[t] * scale_i;
    }
    res.score.resize(static_cast<std::size_t>(dim));
    for (int i = 0; i < dim; ++i) {
        res.score[i] = res.gradient[i] * scale_g;
    }

    // SymmetricMatrix.inverse(information).
    try {
        res.asymptotic_covariance = symmetric_inverse(res.information, dim);
        res.asymptotic_covariance_diagonal.resize(static_cast<std::size_t>(dim));
        for (int i = 0; i < dim; ++i) {
            res.asymptotic_covariance_diagonal[i] =
                res.asymptotic_covariance[static_cast<std::size_t>(i) * dim + i];
        }
        res.status = final_information_status::available;
    } catch (const matrix_exception&) {
        res.status = final_information_status::inverse_failed;
    }
    return res;
}

} // namespace kalos
