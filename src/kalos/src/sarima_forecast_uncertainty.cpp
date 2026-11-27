#include "sarima_forecast_uncertainty.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/numerics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace kalos {
namespace {

// Column-major dense square-matrix helpers (element (i, j) of an n-row
// matrix lives at i + n * j). Products carry the (T P) T' association and
// the update order of the proven R reference.

std::vector<double> matrix_times_vector(const std::vector<double>& a, std::size_t n,
                                        const std::vector<double>& x) {
    std::vector<double> y(n, 0.);
    for (std::size_t i = 0; i < n; ++i) {
        double value = 0.;
        for (std::size_t k = 0; k < n; ++k) value += a[i + n * k] * x[k];
        y[i] = value;
    }
    return y;
}

// (A B) C' for n x n factors in column-major storage.
std::vector<double> sandwich(const std::vector<double>& a, const std::vector<double>& p,
                             std::size_t n) {
    std::vector<double> ap(n * n, 0.), out(n * n, 0.);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) {
            double value = 0.;
            for (std::size_t k = 0; k < n; ++k) value += a[i + n * k] * p[k + n * j];
            ap[i + n * j] = value;
        }
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) {
            double value = 0.;
            for (std::size_t k = 0; k < n; ++k) value += ap[i + n * k] * a[j + n * k];
            out[i + n * j] = value;
        }
    return out;
}

void symmetrize(std::vector<double>& p, std::size_t n) {
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < j; ++i) {
            const double mean = 0.5 * (p[i + n * j] + p[j + n * i]);
            p[i + n * j] = p[j + n * i] = mean;
        }
}

double maximum_abs(const std::vector<double>& p) {
    double best = 0.;
    for (const double v : p) best = std::fmax(best, std::fabs(v));
    return best;
}

// Gauss-Jordan solution of the (M^2) x (M^2) Lyapunov system with partial
// pivoting; b is consumed, the solution returned.
std::vector<double> solve_linear(std::vector<double> a, std::vector<double> b, std::size_t n) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        double best = std::fabs(a[col + n * col]);
        for (std::size_t row = col + 1; row < n; ++row)
            if (std::fabs(a[row + n * col]) > best) {
                best = std::fabs(a[row + n * col]);
                pivot = row;
            }
        if (!(best > 0.) || !std::isfinite(best))
            throw std::domain_error("stationary Lyapunov system is singular");
        if (pivot != col)
            for (std::size_t k = 0; k < n; ++k) std::swap(a[col + n * k], a[pivot + n * k]);
        std::swap(b[col], b[pivot]);
        const double diagonal = a[col + n * col];
        for (std::size_t k = 0; k < n; ++k) a[col + n * k] /= diagonal;
        b[col] /= diagonal;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = a[row + n * col];
            if (factor == 0.) continue;
            for (std::size_t k = 0; k < n; ++k) a[row + n * k] -= factor * a[col + n * k];
            b[row] -= factor * b[col];
        }
    }
    return b;
}

} // namespace

namespace detail {

std::optional<std::vector<double>> sarima_forecast_uncertainty(
    const std::vector<double>& response,
    const std::vector<std::vector<double>>& columns,
    const std::vector<double>& ar, const std::vector<double>& sar,
    const std::vector<double>& ma, const std::vector<double>& sma,
    int d, int bd, const std::vector<double>& coefficients,
    const std::vector<double>& covariance, double innovation_variance,
    bool log, int months) {
    if (months < 0) throw std::invalid_argument("sarima_forecast_uncertainty: negative months");
    const std::size_t n = response.size(), p = columns.size();
    if (n == 0) throw std::invalid_argument("sarima_forecast_uncertainty: empty response");
    if (p != coefficients.size())
        throw std::invalid_argument("sarima_forecast_uncertainty: columns and coefficients disagree");
    if (covariance.size() != p * p)
        throw std::invalid_argument("sarima_forecast_uncertainty: covariance shape");
    if (!(innovation_variance > 0.) || !std::isfinite(innovation_variance))
        throw std::invalid_argument("sarima_forecast_uncertainty: innovation variance");
    if (d < 0 || d > 2 || bd < 0 || bd > 2)
        throw std::invalid_argument("sarima_forecast_uncertainty: differencing order");
    if (ar.size() > 2 || sar.size() > 2 || ma.size() > 2 || sma.size() > 2)
        throw std::invalid_argument("sarima_forecast_uncertainty: ARMA order");
    for (const double value : response)
        if (!std::isfinite(value))
            throw std::invalid_argument("sarima_forecast_uncertainty: non-finite response");
    for (std::size_t j = 0; j < p; ++j) {
        if (columns[j].size() != n + static_cast<std::size_t>(months))
            throw std::invalid_argument("sarima_forecast_uncertainty: a column misses the extended domain");
        for (const double value : columns[j])
            if (!std::isfinite(value))
                throw std::invalid_argument("sarima_forecast_uncertainty: non-finite column");
    }
    for (const double value : covariance)
        if (!std::isfinite(value))
            throw std::invalid_argument("sarima_forecast_uncertainty: non-finite covariance");
    for (const double value : coefficients)
        if (!std::isfinite(value))
            throw std::invalid_argument("sarima_forecast_uncertainty: non-finite coefficient");
    if (!months) return std::vector<double>{};

    // Polynomials in ascending lag (JVM SarimaModel.parameters() convention):
    // A0(B) = (1 + sum ar_i B^i)(1 + sum sar_j B^(12 j)), likewise for the MA.
    const auto block_polynomial = [](const std::vector<double>& block, std::size_t stride) {
        if (block.empty()) return std::vector<double>{1.};
        std::vector<double> poly(stride * block.size() + 1, 0.);
        poly[0] = 1.;
        for (std::size_t i = 0; i < block.size(); ++i) poly[stride * (i + 1)] = block[i];
        return poly;
    };
    const std::vector<double> a0 = convolve(block_polynomial(ar, 1), block_polynomial(sar, 12));
    const std::vector<double> cc = convolve(block_polynomial(ma, 1), block_polynomial(sma, 12));
    const std::vector<double> dd = difference_polynomial(d, bd, 12);

    // Canonical ARMA block on the differenced response plus the exactly
    // diffuse differencing history (the reference build_ss state).
    const std::size_t m1 = std::fmax(a0.size() - 1, cc.size() - 1);
    if (m1 < 1) return std::nullopt;  // fully diffuse models are outside the surface
    const std::size_t mm = 2 * m1, k = dd.size() - 1, nn = mm + k;
    if (nn > 64) return std::nullopt;
    std::vector<double> t(nn * nn, 0.), loading(nn * nn, 0.), z(nn, 0.);
    for (std::size_t i = 0; i < m1; ++i)  // the -A0 lag row and the C row
        t[0 + nn * i] = i < a0.size() - 1 ? -a0[i + 1] : 0.;
    for (std::size_t i = 0; i < m1; ++i)
        t[0 + nn * (m1 + i)] = i < cc.size() - 1 ? cc[i + 1] : 0.;
    for (std::size_t i = 1; i < m1; ++i) t[i + nn * (i - 1)] = 1.;
    for (std::size_t i = 1; i < m1; ++i) t[(m1 + i) + nn * (m1 + i - 1)] = 1.;
    std::vector<double> r(nn, 0.);
    r[0] = 1.;
    r[m1] = 1.;  // G hits positions 1 and m1 + 1 (zero-indexed 0 and m1)
    z[0] = 1.;
    for (std::size_t i = 0; i < k; ++i) z[mm + i] = -dd[i + 1];
    if (k > 0) {
        for (std::size_t j = 0; j < nn; ++j) t[mm + nn * j] = z[j];  // the y_t row
        for (std::size_t j = 1; j < k; ++j) t[(mm + j) + nn * (mm + j - 1)] = 1.;
    }
    for (std::size_t j = 0; j < nn; ++j)
        for (std::size_t i = 0; i < nn; ++i) loading[i + nn * j] = r[i] * r[j];

    // Stationary initialization of the ARMA block: the exact discrete
    // Lyapunov solution vec(P) = (I - kron(F, F))^-1 vec(G G') in the
    // reference's Kronecker layout, then symmetrized.
    const std::size_t sys = mm * mm;
    std::vector<double> kk(sys * sys, 0.), rhs(sys, 0.);
    for (std::size_t row = 0; row < sys; ++row) {
        kk[row + sys * row] += 1.;
        const std::size_t i1 = row / mm, i2 = row % mm;
        for (std::size_t col = 0; col < sys; ++col)
            kk[row + sys * col] -= t[i1 + nn * (col / mm)] * t[i2 + nn * (col % mm)];
    }
    for (std::size_t j = 0; j < mm; ++j)
        for (std::size_t i = 0; i < mm; ++i) rhs[i + mm * j] = r[i] * r[j];
    std::vector<double> pstat;
    try {
        pstat = solve_linear(std::move(kk), std::move(rhs), sys);
    } catch (const std::domain_error&) {
        return std::nullopt;
    }
    if (!std::all_of(pstat.begin(), pstat.end(), [](double v) { return std::isfinite(v); }))
        return std::nullopt;
    std::vector<double> pm(nn * nn, 0.), pc(nn * nn, 0.);
    for (std::size_t i = 0; i < mm; ++i)
        for (std::size_t j = 0; j < mm; ++j)
            pc[i + nn * j] = 0.5 * (pstat[i + mm * j] + pstat[j + mm * i]);
    for (std::size_t i = 0; i < k; ++i) pm[(mm + i) + nn * (mm + i)] = 1.;

    std::vector<double> x(nn, 0.);
    std::vector<std::vector<double>> xs(p, std::vector<double>(nn, 0.));

    // Observation rows: the columns ride the response gains (the reference's
    // gain_of branch), then the response step (the ss_step dual diffuse
    // update with its branch thresholds, rank-annihilation kill, propagation
    // by T and symmetrization order).
    for (std::size_t step = 0; step < n; ++step) {
        const std::vector<double> zm = matrix_times_vector(pm, nn, z);
        const std::vector<double> zc = matrix_times_vector(pc, nn, z);
        double a = 0., b = 0.;
        for (std::size_t i = 0; i < nn; ++i) {
            a += z[i] * zm[i];
            b += z[i] * zc[i];
        }
        std::vector<double> gain0(nn, 0.), next_pm, next_pc;
        if (a > 1e-12) {
            for (std::size_t i = 0; i < nn; ++i) gain0[i] = zm[i] / a;
            next_pm = pm;
            next_pc = pc;
            const double inverse = 1. / a, scaled = b / (a * a);
            for (std::size_t j = 0; j < nn; ++j)
                for (std::size_t i = 0; i < nn; ++i) {
                    next_pm[i + nn * j] -= zm[i] * zm[j] * inverse;
                    next_pc[i + nn * j] -= (zm[i] * zc[j] + zc[i] * zm[j]) * inverse;
                    next_pc[i + nn * j] += scaled * zm[i] * zm[j];
                }
            if (maximum_abs(next_pm) < 1e-9) std::fill(next_pm.begin(), next_pm.end(), 0.);
            next_pm.swap(pm);
            next_pc.swap(pc);
        } else {
            if (maximum_abs(pm) < 1e-9) std::fill(pm.begin(), pm.end(), 0.);  // rank annihilated
            if (b > 1e-14) {
                for (std::size_t i = 0; i < nn; ++i) gain0[i] = zc[i] / b;
                for (std::size_t j = 0; j < nn; ++j)
                    for (std::size_t i = 0; i < nn; ++i)
                        pc[i + nn * j] -= zc[i] * zc[j] / b;
            }
        }
        for (std::size_t j = 0; j < p; ++j) {
            std::vector<double>& state = xs[j];
            double prediction = 0.;
            for (std::size_t i = 0; i < nn; ++i) prediction += z[i] * state[i];
            const double innovation = columns[j][step] - prediction;
            for (std::size_t i = 0; i < nn; ++i) state[i] += gain0[i] * innovation;
            state = matrix_times_vector(t, nn, state);
        }
        double prediction = 0.;
        for (std::size_t i = 0; i < nn; ++i) prediction += z[i] * x[i];
        const double innovation = response[step] - prediction;
        for (std::size_t i = 0; i < nn; ++i) x[i] += gain0[i] * innovation;
        x = matrix_times_vector(t, nn, x);
        pm = sandwich(t, pm, nn);
        pc = sandwich(t, pc, nn);
        for (std::size_t j = 0; j < nn * nn; ++j) pc[j] += loading[j];
        symmetrize(pm, nn);
        symmetrize(pc, nn);
    }

    // Future rows: predict-only response; the columns record their
    // sensitivities on the prior states and propagate under T alone (the
    // JVM FastDkFilter skips the gain where the response is missing).
    std::vector<double> errors(static_cast<std::size_t>(months));
    for (int h = 1; h <= months; ++h) {
        double f0 = 0., v = 0.;
        const std::vector<double> zc = matrix_times_vector(pc, nn, z);
        for (std::size_t i = 0; i < nn; ++i) {
            f0 += z[i] * x[i];
            v += z[i] * zc[i];
        }
        double forecast = f0, variance = innovation_variance * v;
        std::vector<double> s(p, 0.);
        for (std::size_t j = 0; j < p; ++j) {
            std::vector<double>& state = xs[j];
            double prediction = 0.;
            for (std::size_t i = 0; i < nn; ++i) prediction += z[i] * state[i];
            s[j] = columns[j][n + static_cast<std::size_t>(h) - 1] - prediction;
            forecast += s[j] * coefficients[j];
            state = matrix_times_vector(t, nn, state);
        }
        for (std::size_t i = 0; i < p; ++i)
            for (std::size_t j = 0; j < p; ++j) variance += s[i] * covariance[i * p + j] * s[j];
        if (!(variance > 0.) || !std::isfinite(variance) || !(v >= 0.)) return std::nullopt;
        double error = std::sqrt(variance);
        if (log) {
            // JD+ LogNormal.stdev on the raw forecast and its total variance.
            error = detail::fdlibm_exp(forecast + 0.5 * variance) *
                    std::sqrt(detail::fdlibm_exp(variance) - 1.);
        }
        if (!std::isfinite(error)) return std::nullopt;
        errors[static_cast<std::size_t>(h) - 1] = error;
        if (h < months) {  // the predict-only response step for the next row
            x = matrix_times_vector(t, nn, x);
            pm = sandwich(t, pm, nn);
            pc = sandwich(t, pc, nn);
            for (std::size_t j = 0; j < nn * nn; ++j) pc[j] += loading[j];
            symmetrize(pm, nn);
            symmetrize(pc, nn);
        }
    }
    return errors;
}

} // namespace detail
} // namespace kalos
