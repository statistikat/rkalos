#include "kalos/hannan_rissanen.hpp"

#include "kalos/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

// Pinned HannanRissanenInitializer.EPS and HannanRissanen bias guard.
constexpr double rms_guard = 1e-9;
constexpr double bias_guard = 1e16;

void checkpoint(const Interrupt& interrupt) {
    if (interrupt) interrupt();
}

// Pinned AbsMeanNormalizer: divide every entry by the mean absolute value of
// the finite entries; a zero sum is a no-op.
void normalize_abs_mean(std::vector<double>& values) {
    double sum = 0.0;
    std::size_t finite = 0;
    for (const double value : values) {
        if (std::isfinite(value)) {
            sum += std::abs(value);
            ++finite;
        }
    }
    if (sum == 0.0) return;
    const double scale = static_cast<double>(finite) / sum;
    for (double& value : values) value *= scale;
}

// Pinned LevinsonAlgorithm: Durbin recursion on the raw (non-demeaned)
// autocovariance sums; the returned coefficients are the negated prediction
// polynomial beyond its constant term.
std::vector<double> levinson_coefficients(const std::vector<double>& y, int order) {
    const std::size_t n = y.size();
    std::vector<double> x(static_cast<std::size_t>(order) + 1, 0.0);
    for (int k = 0; k <= order; ++k) {
        for (std::size_t t = 0; t + static_cast<std::size_t>(k) < n; ++t) {
            x[static_cast<std::size_t>(k)] += y[t] * y[t + static_cast<std::size_t>(k)];
        }
    }
    std::vector<double> phi(static_cast<std::size_t>(order) + 1, 0.0);
    phi[0] = 1.0;
    double delta = x[0];
    for (int m = 0; m < order; ++m) {
        double k = 0.0;
        for (int i = 0; i <= m; ++i) k -= phi[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(m + 1 - i)];
        k /= delta;
        for (int i = 0; i <= (m + 1) / 2; ++i) {
            const double previous = phi[static_cast<std::size_t>(m + 1 - i)];
            phi[static_cast<std::size_t>(m + 1 - i)] = previous + k * phi[static_cast<std::size_t>(i)];
            phi[static_cast<std::size_t>(i)] += k * previous;
        }
        delta *= 1.0 - k * k;
    }
    std::vector<double> a(static_cast<std::size_t>(order));
    for (std::size_t j = 0; j < a.size(); ++j) a[j] = -phi[j + 1];
    return a;
}

// Pinned AutoRegressiveEstimation.residuals(): top-down in-place filtering
// against the untouched history, truncating the AR window at the series start.
std::vector<double> ar_residuals(const std::vector<double>& y, const std::vector<double>& a) {
    std::vector<double> r = y;
    for (std::size_t i = y.size(); i-- > 0;) {
        const std::size_t window = std::min(a.size(), i);
        for (std::size_t j = 1; j <= window; ++j) r[i] -= a[j - 1] * y[i - j];
    }
    return r;
}

// Pinned HannanRissanen.npi(): AR approximation length of the prefilter.
int ar_window(std::size_t n, const HannanRissanenOrders& o) {
    const double logn = std::log(static_cast<double>(n));
    const double factor = static_cast<int>(logn * logn);
    const std::size_t dense_ar = static_cast<std::size_t>(o.p) + static_cast<std::size_t>(o.period) * static_cast<std::size_t>(o.bp);
    const std::size_t dense_ma = 2 * (static_cast<std::size_t>(o.q) + static_cast<std::size_t>(o.period) * static_cast<std::size_t>(o.bq));
    long long length = static_cast<long long>(factor) * static_cast<long long>(std::max(dense_ar, dense_ma));
    const long long nn = static_cast<long long>(n);
    if (length >= nn) length = nn - nn / 4;
    if (length > 50) length = 50;
    return static_cast<int>(length);
}

// Pinned m_pi layout [phi(p) | bp blocks of (p+1) | theta(q) | bq blocks of
// (q+1)] in the rjd3 sign convention (polynomials 1 + phi(B), 1 + theta(B)).
struct Pipeline {
    std::vector<double> data; // AbsMean-normalized response
    std::vector<double> a;    // Levinson prefilter residuals (m_a)
    std::vector<double> pi;
    HannanRissanenOrders o;
    bool biascorrection;
    bool finalcorrection;
    const Interrupt& interrupt;
};

std::size_t ar_count(const Pipeline& h) {
    return static_cast<std::size_t>(h.o.p) + static_cast<std::size_t>(h.o.bp) * (1 + static_cast<std::size_t>(h.o.p));
}

std::size_t ma_count(const Pipeline& h) {
    return static_cast<std::size_t>(h.o.q) + static_cast<std::size_t>(h.o.bq) * (1 + static_cast<std::size_t>(h.o.q));
}

// Pinned ls(): ordinary least squares on the dense design; the caller owns
// the sign conventions. DGELSY is the qualified Kalos solver adapter; it
// coincides with the pinned unpivoted Householder solution whenever the
// design has full column rank.
std::vector<double> least_squares(Pipeline& h, const DenseMatrix& x, const std::vector<double>& rhs) {
    const DenseMatrix b(rhs.size(), 1, rhs);
    const auto fit = qr_least_squares(x, b, 1e-10, h.interrupt);
    return fit.coefficients.values();
}

// Pinned minspq(): dense two-sided regression y[t] ~ [y-lags | seasonal
// y-lags | a-lags | seasonal a-lags] over every t beyond the warmup, then the
// AR section of m_pi is negated (the design estimates -phi).
void minspq(Pipeline& h) {
    const std::size_t p = static_cast<std::size_t>(h.o.p), bp = static_cast<std::size_t>(h.o.bp);
    const std::size_t q = static_cast<std::size_t>(h.o.q), bq = static_cast<std::size_t>(h.o.bq);
    const std::size_t period = static_cast<std::size_t>(h.o.period);
    const std::size_t nP = ar_count(h), nQ = ma_count(h);
    // Pinned m_i: the largest design lag, p + period*bp or q + period*bq.
    const std::size_t warmup = std::max(p + period * bp, q + period * bq);
    const std::size_t n = h.data.size();
    if (n <= warmup) {
        throw std::runtime_error("hannan_rissanen_start: response shorter than the design warmup");
    }
    const std::size_t rows = n - warmup, columns = nP + nQ;
    DenseMatrix x(rows, columns);
    std::vector<double> rhs(rows);
    for (std::size_t r = 0; r < rows; ++r) rhs[r] = h.data[warmup + r];
    std::size_t column = 0;
    const auto fill = [&](const std::vector<double>& source, std::size_t lag) {
        double* target = x.values().data() + column * rows;
        for (std::size_t r = 0; r < rows; ++r) target[r] = source[warmup + r - lag];
        ++column;
    };
    for (std::size_t j = 1; j <= p; ++j) fill(h.data, j);
    for (std::size_t block = 1; block <= bp; ++block) {
        for (std::size_t j = 0; j <= p; ++j) fill(h.data, block * period + j);
    }
    for (std::size_t j = 1; j <= q; ++j) fill(h.a, j);
    for (std::size_t block = 1; block <= bq; ++block) {
        for (std::size_t j = 0; j <= q; ++j) fill(h.a, block * period + j);
    }
    checkpoint(h.interrupt);
    h.pi = least_squares(h, x, rhs);
    for (std::size_t j = 0; j < nP; ++j) h.pi[j] = -h.pi[j];
}

// Pinned biascorrection(): ea/es sensitivity recursions against the current
// m_pi, one least-squares pass, m_pi += correction. A divergent innovation
// (|value| > 1e16) aborts with m_pi unchanged, as the pinned early return.
void biascorrection(Pipeline& h) {
    const std::size_t p = static_cast<std::size_t>(h.o.p), bp = static_cast<std::size_t>(h.o.bp);
    const std::size_t q = static_cast<std::size_t>(h.o.q), bq = static_cast<std::size_t>(h.o.bq);
    const std::size_t period = static_cast<std::size_t>(h.o.period);
    const std::size_t n = h.data.size();
    std::vector<double> ea(n, 0.0), es(n, 0.0), e(n, 0.0);
    DenseMatrix x(n, ar_count(h) + ma_count(h));
    double* storage = x.values().data();
    for (std::size_t t = 0; t < n; ++t) {
        if ((t & 1023U) == 0) checkpoint(h.interrupt);
        double value = h.data[t], acc_a = 0.0, acc_s = 0.0;
        std::size_t column = 0;
        for (std::size_t j = 1; j <= p; ++j, ++column) {
            if (t >= j) {
                value += h.pi[column] * h.data[t - j];
                acc_a -= h.pi[column] * ea[t - j];
                storage[column * n + t] = -ea[t - j];
            }
        }
        for (std::size_t block = 1; block <= bp; ++block) {
            for (std::size_t j = 0; j <= p; ++j, ++column) {
                const std::size_t lag = block * period + j;
                if (t >= lag) {
                    value += h.pi[column] * h.data[t - lag];
                    acc_a -= h.pi[column] * ea[t - lag];
                    storage[column * n + t] = -ea[t - lag];
                }
            }
        }
        for (std::size_t j = 1; j <= q; ++j, ++column) {
            if (t >= j) {
                value -= h.pi[column] * e[t - j];
                acc_s -= h.pi[column] * es[t - j];
                storage[column * n + t] = es[t - j];
            }
        }
        for (std::size_t block = 1; block <= bq; ++block) {
            for (std::size_t j = 0; j <= q; ++j, ++column) {
                const std::size_t lag = block * period + j;
                if (t >= lag) {
                    value -= h.pi[column] * e[t - lag];
                    acc_s -= h.pi[column] * es[t - lag];
                    storage[column * n + t] = es[t - lag];
                }
            }
        }
        if (std::abs(value) > bias_guard) return;
        e[t] = value;
        ea[t] = acc_a + value;
        es[t] = acc_s + value;
    }
    checkpoint(h.interrupt);
    const auto correction = least_squares(h, x, e);
    for (std::size_t j = 0; j < correction.size(); ++j) h.pi[j] += correction[j];
}

// Pinned HannanRissanen.calc().
void run(Pipeline& h);

// Pinned SarimaModel.getAr() as the full raw AR polynomial (product of the
// regular factor and one factor per seasonal block, no coefficient trimming).
std::vector<double> ar_polynomial(const Pipeline& h) {
    std::vector<double> product(1, 1.0);
    const auto multiply = [&](const std::vector<double>& factor) {
        std::vector<double> next(product.size() + factor.size() - 1, 0.0);
        for (std::size_t i = 0; i < product.size(); ++i) {
            for (std::size_t j = 0; j < factor.size(); ++j) next[i + j] += product[i] * factor[j];
        }
        product = std::move(next);
    };
    std::vector<double> regular(static_cast<std::size_t>(h.o.p) + 1, 0.0);
    regular[0] = 1.0;
    for (std::size_t j = 0; j < static_cast<std::size_t>(h.o.p); ++j) regular[j + 1] = h.pi[j];
    multiply(regular);
    const std::size_t period = static_cast<std::size_t>(h.o.period);
    for (std::size_t block = 1; block <= static_cast<std::size_t>(h.o.bp); ++block) {
        std::vector<double> factor(block * period + 1, 0.0);
        factor[0] = 1.0;
        factor[block * period] = h.pi[static_cast<std::size_t>(h.o.p) + (block - 1) * (1 + static_cast<std::size_t>(h.o.p))];
        multiply(factor);
    }
    return product;
}

// Pinned finalcorrection(): filter the normalized response through the
// estimated AR polynomial, refit the MA orders on the filtered series with a
// fresh pipeline (final correction disabled, bias inherited) and overwrite the
// MA tail of m_pi with the refit.
void finalcorrection(Pipeline& h) {
    const std::vector<double> filter = ar_polynomial(h);
    const std::size_t degree = filter.size() - 1;
    if (h.data.size() <= degree) {
        throw std::runtime_error("hannan_rissanen_start: filtered series is empty");
    }
    std::vector<double> filtered(h.data.size() - degree);
    for (std::size_t t = 0; t < filtered.size(); ++t) {
        double value = 0.0;
        for (std::size_t k = 0; k <= degree; ++k) value += filter[k] * h.data[t + k];
        filtered[t] = value;
    }
    checkpoint(h.interrupt);
    HannanRissanenOrders reduced{0, 0, h.o.q, h.o.bq, h.o.period};
    Pipeline second{std::move(filtered), {}, {}, reduced, h.biascorrection, false, h.interrupt};
    normalize_abs_mean(second.data);
    run(second);
    if (second.pi.size() > h.pi.size()) {
        throw std::logic_error("hannan_rissanen_start: correction tail exceeds m_pi");
    }
    std::copy(second.pi.begin(), second.pi.end(), h.pi.end() - static_cast<std::ptrdiff_t>(second.pi.size()));
}

// Pinned HannanRissanen.calc().
void run(Pipeline& h) {
    const std::size_t ar = static_cast<std::size_t>(h.o.p) + static_cast<std::size_t>(h.o.period) * static_cast<std::size_t>(h.o.bp);
    const std::size_t ma = static_cast<std::size_t>(h.o.q) + static_cast<std::size_t>(h.o.period) * static_cast<std::size_t>(h.o.bq);
    h.pi.assign(ar_count(h) + ma_count(h), 0.0);
    if (ar == 0 && ma == 0) return; // The pinned bic-only shortcut.
    if (ma > 0) {
        const int window = ar_window(h.data.size(), h.o);
        if (window < 0) throw std::runtime_error("hannan_rissanen_start: negative AR window");
        h.a = ar_residuals(h.data, levinson_coefficients(h.data, window));
    }
    minspq(h);
    if (ma > 0 && h.biascorrection) biascorrection(h);
    if (ma > 0 && ar > 0 && h.finalcorrection) finalcorrection(h);
}

// Pinned SarimaMapping.stabilize(DataBlock): trim trailing zeros, leave a
// stable block, invert roots inside the unit circle of 1 + b(B) by
// re-forming the polynomial from the inverted roots and normalizing the
// constant term to one. Returns whether the block was rewritten.
bool stabilize_block(std::vector<double>& block) {
    std::size_t length = block.size();
    while (length > 0 && block[length - 1] == 0.0) --length;
    if (length == 0) return false;
    if (length == 1) {
        if (std::abs(block[0]) > 1.0) {
            block[0] = 1.0 / block[0];
            return true;
        }
        return false;
    }
    std::vector<double> polynomial(length + 1, 0.0);
    polynomial[0] = 1.0;
    for (std::size_t i = 0; i < length; ++i) polynomial[i + 1] = block[i];
    const auto found = polynomial_roots(polynomial);
    if (found.empty()) return false;
    std::vector<std::complex<double>> roots;
    roots.reserve(found.size());
    for (const auto& root : found) roots.emplace_back(root.real, root.imag);
    bool changed = false;
    for (auto& root : roots) {
        const double modulus = std::abs(root);
        if (1.0 / modulus > 1.0) {
            if (modulus == 0.0) {
                throw std::runtime_error("hannan_rissanen_start: zero root cannot be inverted");
            }
            root = std::conj(root) / (modulus * modulus);
            changed = true;
        }
    }
    if (!changed) return false;
    // Pinned Polynomial.fromComplexRoots followed by division by the
    // constant term; only real parts cross back into the block.
    std::vector<std::complex<double>> z(roots.size() + 1, {0.0, 0.0});
    z[0] = {1.0, 0.0};
    for (std::size_t i = 0; i < roots.size(); ++i) {
        z[i + 1] = z[i];
        for (std::size_t j = i; j >= 1; --j) z[j] = z[j - 1] - z[j] * roots[i];
        z[0] = z[0] * (-roots[i]);
    }
    if (!(z[0].real() != 0.0)) {
        throw std::runtime_error("hannan_rissanen_start: roots must all be different from 0");
    }
    for (std::size_t i = 0; i < length; ++i) {
        if (!std::isfinite(z[i + 1].real())) {
            throw std::runtime_error("hannan_rissanen_start: unstable reconstruction");
        }
        block[i] = z[i + 1].real() / z[0].real();
    }
    return true;
}

// Pinned SarimaMapping.stabilize(boolean, SarimaOrders, DataBlock) over the
// model parameter layout: regular AR, one entry per seasonal AR block,
// regular MA, one entry per seasonal MA block, in that order.
bool stabilize_model(std::array<std::vector<double>, 4>& raw) {
    bool changed = false;
    changed |= stabilize_block(raw[0]);
    for (auto& entry : raw[1]) {
        std::vector<double> singleton{entry};
        if (stabilize_block(singleton)) {
            entry = singleton[0];
            changed = true;
        }
    }
    changed |= stabilize_block(raw[2]);
    for (auto& entry : raw[3]) {
        std::vector<double> singleton{entry};
        if (stabilize_block(singleton)) {
            entry = singleton[0];
            changed = true;
        }
    }
    return changed;
}

} // namespace

HannanRissanenStart hannan_rissanen_start(const std::vector<double>& response,
                                          const RegressionDesign& design,
                                          const HannanRissanenOrders& orders,
                                          const Interrupt& interrupt) {
    if (orders.p < 0 || orders.bp < 0 || orders.q < 0 || orders.bq < 0 || orders.period < 1) {
        throw std::invalid_argument("hannan_rissanen_start: orders need non-negative block sizes and a positive period");
    }
    const std::size_t n = response.size();
    if (n == 0) throw std::invalid_argument("hannan_rissanen_start: empty response");
    for (const double value : response) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("hannan_rissanen_start: non-finite response value");
        }
    }
    if (design.columns > 0) {
        if (design.values.size() != design.columns * n) {
            throw std::invalid_argument("hannan_rissanen_start: design does not match the response");
        }
        for (const double value : design.values) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("hannan_rissanen_start: non-finite design value");
            }
        }
    }

    std::exception_ptr interrupted;
    const Interrupt checked_interrupt = [&] {
        if (!interrupt) return;
        try { interrupt(); }
        catch (...) { interrupted = std::current_exception(); throw; }
    };

    HannanRissanenStart start{HannanRissanenStart::Source::Fitted, false,
        {std::vector<double>(static_cast<std::size_t>(orders.p), 0.0),
         std::vector<double>(static_cast<std::size_t>(orders.bp), 0.0),
         std::vector<double>(static_cast<std::size_t>(orders.q), 0.0),
         std::vector<double>(static_cast<std::size_t>(orders.bq), 0.0)}};

    try {
        // Pinned Ols.compute on the differenced linear model: the pipeline
        // runs on y - X b, or on y itself when the model has no regressors.
        std::vector<double> adjusted = response;
        if (design.columns > 0) {
            DenseMatrix x(n, design.columns, design.values);
            const DenseMatrix b(n, 1, response);
            const auto fit = qr_least_squares(x, b, 1e-10, checked_interrupt);
            adjusted = fit.residuals.values(); // y - X b, observation order.
        }
        // Pinned root-mean-square guard on those residuals.
        double ssq = 0.0;
        for (const double value : adjusted) ssq += value * value;
        if (std::sqrt(ssq / static_cast<double>(n)) < rms_guard) {
            start.source = HannanRissanenStart::Source::ZeroVariance;
            return start;
        }
        Pipeline h{std::move(adjusted), {}, {}, orders, true, true, checked_interrupt};
        normalize_abs_mean(h.data);
        run(h);
        for (const double value : h.pi) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("hannan_rissanen_start: pipeline produced non-finite coefficients");
            }
        }
        // Pinned updatemodel extraction (only the first entry of each
        // seasonal block reaches the model), then the stabilized SarimaModel.
        const std::size_t p_count = static_cast<std::size_t>(orders.p);
        const std::size_t q_count = static_cast<std::size_t>(orders.q);
        const std::size_t nP = p_count + static_cast<std::size_t>(orders.bp) * (1 + p_count);
        auto& raw = start.coefficients;
        raw[0].assign(h.pi.begin(), h.pi.begin() + static_cast<std::ptrdiff_t>(p_count));
        for (std::size_t block = 0; block < static_cast<std::size_t>(orders.bp); ++block) {
            raw[1][block] = h.pi[p_count + block * (1 + p_count)];
        }
        raw[2].assign(h.pi.begin() + static_cast<std::ptrdiff_t>(nP),
                      h.pi.begin() + static_cast<std::ptrdiff_t>(nP + q_count));
        for (std::size_t block = 0; block < static_cast<std::size_t>(orders.bq); ++block) {
            raw[3][block] = h.pi[nP + block * (1 + q_count)];
        }
        start.stabilized = stabilize_model(raw);
        // rjd3 stores raw 1 + phi(B); the Kalos convention is 1 - phi(B).
        for (double& value : raw[0]) value = -value;
        for (double& value : raw[1]) value = -value;
        checkpoint(checked_interrupt);
        return start;
    } catch (...) {
        if (interrupted) std::rethrow_exception(interrupted);
        // Pinned useDefaultIfFailed(true): the SarimaModel.Builder
        // setDefault() fill, raw phi = bphi = -0.1, th = bth = -0.2.
        for (double& value : start.coefficients[0]) value = 0.1;
        for (double& value : start.coefficients[1]) value = 0.1;
        for (double& value : start.coefficients[2]) value = -0.2;
        for (double& value : start.coefficients[3]) value = -0.2;
        start.source = HannanRissanenStart::Source::Defaulted;
        start.stabilized = false;
        return start;
    }
}

} // namespace kalos
