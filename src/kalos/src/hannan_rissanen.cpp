#include "kalos/hannan_rissanen.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/hannan_rissanen_pipeline.hpp"
#include "kalos/detail/sarima_stabilization.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/detail/ssq_levmar.hpp"
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
// Pinned AbsMeanNormalizer.normalize(DataBlock): scale = n / sum(|x|) over
// the finite entries (1.0 when the sum is 0), applied by multiplication to
// every entry. (getNormalizationFactor, which restricts to |x| > 1, is a
// different method not used on this path.)
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
            // Pinned temporary order (LevinsonAlgorithm 158-211): var11 =
            // a[m+1-i] + k*a[i] is computed from the OLD pair before either
            // slot is written, so the middle slot (i == m+1-i) receives
            // P*(1+k), not the double update P*(1+2k) an in-place sequence
            // without the temporary would produce.
            const double other = phi[static_cast<std::size_t>(m + 1 - i)];
            const double new_other = other + k * phi[static_cast<std::size_t>(i)];
            phi[static_cast<std::size_t>(i)] = phi[static_cast<std::size_t>(i)] + k * other;
            phi[static_cast<std::size_t>(m + 1 - i)] = new_other;
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

// Pinned HannanRissanen.npi() (bytecode 747-105): the Levinson prefilter
// window is max((int)(log n)^2, max(p + P*bp, 2*(q + P*bq))) (two Math.max
// calls, NOT a product), replaced by n - n/4 when it reaches n, then capped
// at 50.
int ar_window(std::size_t n, const HannanRissanenOrders& o) {
    const double logn = detail::fdlibm_log(static_cast<double>(n));
    const long long logfactor = static_cast<int>(logn * logn);
    const long long dense_ar = static_cast<long long>(o.p) + static_cast<long long>(o.period) * static_cast<long long>(o.bp);
    const long long dense_ma = 2 * (static_cast<long long>(o.q) + static_cast<long long>(o.period) * static_cast<long long>(o.bq));
    long long length = std::max(logfactor, std::max(dense_ar, dense_ma));
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
    double bic = 0.0; // Pinned HannanRissanen.bic, set by the minspq ls().
    const Interrupt& interrupt;
};

std::size_t ar_count(const Pipeline& h) {
    return static_cast<std::size_t>(h.o.p) + static_cast<std::size_t>(h.o.bp) * (1 + static_cast<std::size_t>(h.o.p));
}

std::size_t ma_count(const Pipeline& h) {
    return static_cast<std::size_t>(h.o.q) + static_cast<std::size_t>(h.o.bq) * (1 + static_cast<std::size_t>(h.o.q));
}

// Pinned ls(): QRLeastSquaresSolver.fastLeastSquares (unpivoted
// Householder2 chain, UpperTriangularMatrix.rank(R, 0) and solveUx) on the
// dense design, bit-reproduced by the shared pinned solver transcription;
// DGELSY only agrees at full rank and steered the Multiple-start chain into
// the wrong basins (21-D.2). The caller owns the sign conventions.
struct PinnedPipelineLs {
    std::vector<double> b; // getB(), zero padded to the column count
    double ssq_err = 0.0;
};

PinnedPipelineLs pinned_ls(const DenseMatrix& x, const std::vector<double>& rhs,
    const Interrupt& interrupt) {
    checkpoint(interrupt);
    const std::size_t rows = x.rows();
    const std::size_t columns = x.columns();
    auto solution = detail::pinned_least_squares(x.values(), rhs, rows, columns);
    if (solution.singular) {
        // The pinned solveUx throws before any residual is read, so the
        // (empty) tail must not be summed here.
        throw std::runtime_error("hannan_rissanen_start: singular design");
    }
    PinnedPipelineLs out;
    out.b = std::move(solution.x);
    // QRLeastSquaresSolution.getSsqErr: DataBlock.ssq over the residual
    // tail (the b coordinates from rank on, stored compactly by the
    // transcription), ascending, skipping the non-finite coordinates.
    double acc = 0.0;
    for (const double value : solution.residual_tail) {
        if (std::isfinite(value)) {
            acc += value * value;
        }
    }
    out.ssq_err = acc;
    return out;
}

std::vector<double> least_squares(Pipeline& h, const DenseMatrix& x, const std::vector<double>& rhs,
    bool bic) {
    const auto solution = pinned_ls(x, rhs, h.interrupt);
    if (bic) {
        // Pinned bic = log(ssqErr/n) + log(n) * count(b != 0) / n.
        const double n = static_cast<double>(rhs.size());
        std::size_t count = 0;
        for (const double value : solution.b) {
            if (value != 0.0) ++count;
        }
        h.bic = detail::fdlibm_log(solution.ssq_err / n) +
            detail::fdlibm_log(n) * static_cast<double>(count) / n;
    }
    return solution.b;
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
    if (columns > rows) {
        // Pinned short-window failure class: with more design columns than
        // rows the pinned JVM continues past the (well-defined) min-width
        // Householder stage into reflector and residual layouts that go out
        // of bounds (the corpus records ArrayIndexOutOfBounds/
        // NegativeArraySize for these cases). The qualified seam refuses the
        // degenerate design instead of reproducing a JVM memory failure.
        throw std::runtime_error("hannan_rissanen_start: design exceeds the short response window");
    }
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
    h.pi = least_squares(h, x, rhs, true);
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
    const auto correction = least_squares(h, x, e, false);
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
    // Pinned IFiniteFilter.apply(DataBlock, DataBlock) unit-increment fast
    // path: out[i] = sum_{k ascending} w[k] * data[i + k] with w the reversed
    // polynomial coefficients. The accumulation direction decides ulps, and
    // pass-2 sensitivity steers the Multiple-start chain (21-D.2).
    std::vector<double> filtered = detail::pinned_backfilter_apply(filter, h.data);
    if (filtered.empty()) {
        throw std::runtime_error("hannan_rissanen_start: filtered series is empty");
    }
    checkpoint(h.interrupt);
    HannanRissanenOrders reduced{0, 0, h.o.q, h.o.bq, h.o.period};
    Pipeline second{std::move(filtered), {}, {}, reduced, h.biascorrection, false, 0.0,
        h.interrupt};
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

// Pinned SarimaMapping.stabilize(boolean, SarimaOrders, DataBlock) over the
// initializer parameter layout: regular AR, one entry per seasonal AR block,
// regular MA, one entry per seasonal MA block, in that order.
bool stabilize_model(std::array<std::vector<double>, 4>& raw) {
    bool changed = false;
    changed |= detail::stabilize_coefficient_block(raw[0]);
    for (auto& entry : raw[1]) {
        std::vector<double> singleton{entry};
        if (detail::stabilize_coefficient_block(singleton)) {
            entry = singleton[0];
            changed = true;
        }
    }
    changed |= detail::stabilize_coefficient_block(raw[2]);
    for (auto& entry : raw[3]) {
        std::vector<double> singleton{entry};
        if (detail::stabilize_coefficient_block(singleton)) {
            entry = singleton[0];
            changed = true;
        }
    }
    return changed;
}

} // namespace

namespace detail {

std::array<std::vector<double>, 4> hannan_rissanen_pipeline(const std::vector<double>& response,
    const HannanRissanenOrders& orders, const Interrupt& interrupt) {
    if (orders.p < 0 || orders.bp < 0 || orders.q < 0 || orders.bq < 0 || orders.period < 1) {
        throw std::invalid_argument("hannan_rissanen_pipeline: negative orders");
    }
    if (response.empty()) throw std::invalid_argument("hannan_rissanen_pipeline: empty response");
    for (const double value : response) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("hannan_rissanen_pipeline: non-finite response value");
        }
    }
    std::array<std::vector<double>, 4> blocks{
        std::vector<double>(static_cast<std::size_t>(orders.p), 0.0),
        std::vector<double>(static_cast<std::size_t>(orders.bp), 0.0),
        std::vector<double>(static_cast<std::size_t>(orders.q), 0.0),
        std::vector<double>(static_cast<std::size_t>(orders.bq), 0.0)};
    Pipeline h{response, {}, {}, orders, true, true, 0.0, interrupt};
    normalize_abs_mean(h.data);
    run(h);
    for (const double value : h.pi) {
        if (!std::isfinite(value)) {
            throw std::runtime_error("hannan_rissanen_pipeline: pipeline produced non-finite coefficients");
        }
    }
    // Pinned updatemodel extraction: only the first entry of each seasonal
    // block reaches the model; the blocks stay in the rjd3 sign convention.
    const std::size_t p_count = static_cast<std::size_t>(orders.p);
    const std::size_t q_count = static_cast<std::size_t>(orders.q);
    const std::size_t nP = p_count + static_cast<std::size_t>(orders.bp) * (1 + p_count);
    blocks[0].assign(h.pi.begin(), h.pi.begin() + static_cast<std::ptrdiff_t>(p_count));
    for (std::size_t block = 0; block < static_cast<std::size_t>(orders.bp); ++block) {
        blocks[1][block] = h.pi[p_count + block * (1 + p_count)];
    }
    blocks[2].assign(h.pi.begin() + static_cast<std::ptrdiff_t>(nP),
                     h.pi.begin() + static_cast<std::ptrdiff_t>(nP + q_count));
    for (std::size_t block = 0; block < static_cast<std::size_t>(orders.bq); ++block) {
        // Pinned updatemodel: btheta takes m_pi[i1], then i1 += 1 + q, so
        // the seasonal MA section starts after the q regular theta entries.
        blocks[3][block] = h.pi[nP + q_count + block * (1 + q_count)];
    }
    return blocks;
}

} // namespace detail

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
        // Pinned OlsComputer guard: a response with norm2 below machine
        // epsilon makes Ols.compute return null, which the initializer
        // swallows into the defaulted start.
        double ynorm2 = 0.0;
        for (const double value : response) {
            ynorm2 += value * value;
        }
        if (std::sqrt(ynorm2) < 1.1102230246251565e-16) { // dcmpg: NaN passes
            start.source = HannanRissanenStart::Source::Defaulted;
            for (double& value : start.coefficients[0]) value = 0.1;
            for (double& value : start.coefficients[1]) value = 0.1;
            for (double& value : start.coefficients[2]) value = -0.2;
            for (double& value : start.coefficients[3]) value = -0.2;
            return start;
        }
        std::vector<double> adjusted = response;
        if (design.columns == 1) {
            // The pinned HouseholderWithPivot solve (bit-exact), not the
            // LAPACK QR: the ulp of this coefficient steers the whole
            // Multiple-start chain (21-D.2).
            const double b = detail::pinned_ols_1col(design.values, response);
            adjusted.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                adjusted[i] = response[i] - design.values[i] * b;
            }
        } else if (design.columns > 0) {
            std::vector<std::vector<double>> columns;
            for (std::size_t j = 0; j < design.columns; ++j) {
                checkpoint(checked_interrupt);
                auto first = design.values.begin() + static_cast<std::ptrdiff_t>(j * n);
                columns.emplace_back(first, first + static_cast<std::ptrdiff_t>(n));
            }
            const auto b = detail::pinned_ols_coefficients(response, columns);
            if (!b) throw std::runtime_error("Hannan-Rissanen initialization: OLS failed");
            for (std::size_t j = 0; j < columns.size(); ++j) {
                checkpoint(checked_interrupt);
                for (std::size_t i = 0; i < n; ++i)
                    adjusted[i] += -(*b)[j] * columns[j][i];
            }
        }
        // Pinned root-mean-square guard on those residuals.
        double ssq = 0.0;
        for (const double value : adjusted) ssq += value * value;
        if (std::sqrt(ssq / static_cast<double>(n)) < rms_guard) {
            start.source = HannanRissanenStart::Source::ZeroVariance;
            return start;
        }
        start.coefficients = detail::hannan_rissanen_pipeline(adjusted, orders, checked_interrupt);
        start.stabilized = stabilize_model(start.coefficients);
        // rjd3 stores raw 1 + phi(B); the Kalos convention is 1 - phi(B).
        for (double& value : start.coefficients[0]) value = -value;
        for (double& value : start.coefficients[1]) value = -value;
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
