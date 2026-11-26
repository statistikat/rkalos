#include "kalos/finite_state_factor.hpp"

// R supplies the LAPACK eigendecomposition and orthogonal QR factorizations.
// The finite-factor Joseph update and concatenated disturbance propagation are
// independently assembled from P=C C'; no numerical filter is translated.
#include "kalos/detail/lapack.hpp"
#include "kalos/detail/hostlineage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace detail {
namespace {

constexpr std::size_t max_states = 32;
constexpr double residual_tolerance = 512 * std::numeric_limits<double>::epsilon();

void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }

// long double is accumulation scratch only; the numeric contract of every
// guarded quantity is its correctly rounded binary64 value. A nonzero
// extended-precision intermediate that rounds to zero at binary64 granularity
// is zero on every platform, matching the double-width platforms on which
// this library is qualified; see
// docs/contracts/extended-precision-range-contract.md. Underflow contracts
// are enforced separately on named binary64 quantities below.
double checked(long double value) {
    const auto result = static_cast<double>(value);
    if (!std::isfinite(result)) throw std::overflow_error("finite state factor arithmetic produced a non-finite value");
    return result;
}

std::size_t initial_dimension(const DenseMatrix& covariance) {
    covariance.validate();
    if (covariance.rows() != covariance.columns() || covariance.rows() > max_states) {
        throw std::invalid_argument("finite state covariance must be square with at most 32 states");
    }
    return covariance.rows();
}

int workspace_size(double query, int minimum) {
    if (!std::isfinite(query) || query < 1 || query > std::numeric_limits<int>::max()) {
        throw std::length_error("finite state factor LAPACK workspace exceeds the integer range");
    }
    return std::max(minimum, static_cast<int>(std::ceil(query)));
}

void check_variance(double h) {
    if (!std::isfinite(h) || h < 0) throw std::invalid_argument("finite state observation variance must be finite and nonnegative");
}

} // namespace

FiniteStateFactor::FiniteStateFactor(const DenseMatrix& covariance,
                                     const Interrupt& interrupt)
    : factor_(initial_dimension(covariance), covariance.rows()) {
    const auto m = dimension();
    double scale = 0;
    for (const auto value : covariance.values()) scale = std::max(scale, std::abs(value));
    auto eigenvectors = covariance.values();
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            const auto first = covariance(i, j), second = covariance(j, i);
            if (std::abs(static_cast<long double>(first) - second) >
                64 * std::numeric_limits<double>::epsilon() * static_cast<long double>(scale)) {
                throw std::invalid_argument("finite state covariance must be symmetric");
            }
            const auto average = checked(.5L * first + .5L * second);
            eigenvectors[i + m * j] = eigenvectors[j + m * i] = average;
        }
    }
    const int n = static_cast<int>(m), query_size = -1;
    int info = 0;
    double query = 0;
    std::vector<double> eigenvalues(m);
    checkpoint(interrupt);
    lapack::dsyev("V", "U", &n, eigenvectors.data(), &n, eigenvalues.data(),
                    &query, &query_size, &info);
    checkpoint(interrupt);
    if (info != 0) throw std::runtime_error("finite state covariance eigen workspace query failed");
    const auto work_size = workspace_size(query, std::max(1, 3 * n - 1));
    std::vector<double> work(static_cast<std::size_t>(work_size));
    lapack::dsyev("V", "U", &n, eigenvectors.data(), &n, eigenvalues.data(),
                    work.data(), &work_size, &info);
    checkpoint(interrupt);
    if (info != 0) throw std::runtime_error("finite state covariance eigendecomposition failed");
    const long double negative_bound = static_cast<long double>(residual_tolerance) * static_cast<long double>(m) * scale;
    for (std::size_t j = 0; j < m; ++j) {
        auto value = eigenvalues[j];
        if (!std::isfinite(value) || static_cast<long double>(value) < -negative_bound) {
            throw std::invalid_argument("finite state covariance is not numerically positive semidefinite");
        }
        if (value < 0) { value = 0; ++clipped_negative_count_; }
        const auto root = std::sqrt(value);
        for (std::size_t i = 0; i < m; ++i) {
            factor_(i, j) = checked(static_cast<long double>(eigenvectors[i + m * j]) * root);
        }
    }
    const auto reconstructed = this->covariance(interrupt);
    long double error = 0;
    for (std::size_t i = 0; i < covariance.values().size(); ++i) {
        error = std::max(error, std::abs(static_cast<long double>(reconstructed.values()[i]) - covariance.values()[i]));
    }
    initialization_residual_ = scale == 0 ? checked(error) : checked(error / scale);
    if (initialization_residual_ > residual_tolerance * static_cast<double>(m)) {
        throw std::runtime_error("finite state covariance factor failed its reconstruction residual check");
    }
}

std::vector<double> FiniteStateFactor::observation_loading(const std::vector<double>& z,
                                                           const Interrupt& interrupt) const {
    const auto m = dimension();
    if (z.size() != m) throw std::invalid_argument("finite state observation dimensions must agree");
    for (const auto value : z) if (!std::isfinite(value)) throw std::invalid_argument("finite state observation coefficients must be finite");
    std::vector<double> result(m, 0);
    for (std::size_t j = 0; j < m; ++j) {
        checkpoint(interrupt);
        long double value = 0;
        for (std::size_t i = 0; i < m; ++i) value += static_cast<long double>(z[i]) * factor_(i, j);
        result[j] = checked(value);
    }
    return result;
}

FiniteObservation FiniteStateFactor::observation(const std::vector<double>& z, double h,
                                                 const Interrupt& interrupt) const {
    check_variance(h);
    const auto w = observation_loading(z, interrupt);
    const auto m = dimension();
    double norm = std::sqrt(h);
    for (const auto value : w) norm = checked(kalos::detail::glibc_hypot(norm, value));
    const auto variance = checked(static_cast<long double>(norm) * norm);
    if (norm != 0 && variance == 0) throw std::underflow_error("finite innovation variance underflows");
    std::vector<double> pz(m, 0);
    for (std::size_t i = 0; i < m; ++i) {
        long double value = 0;
        for (std::size_t j = 0; j < m; ++j) value += static_cast<long double>(factor_(i, j)) * w[j];
        pz[i] = checked(value);
    }
    return {variance, std::move(pz)};
}

void FiniteStateFactor::compress_transpose(std::size_t rows, const Interrupt& interrupt) {
    // loading' = Q R, hence loading loading' = R' R. All m rows of R
    // are retained, including zero diagonals: this is not rank truncation.
    const auto m = dimension();
    const int nrows = static_cast<int>(rows), cols = static_cast<int>(m);
    int info = 0;
    std::array<double, max_states> tau;
    checkpoint(interrupt);
    if (rows != qr_workspace_rows_) {
        const int query_size = -1;
        double query = 0;
        lapack::dgeqrf(&nrows, &cols, qr_input_.data(), &nrows, tau.data(), &query, &query_size, &info);
        checkpoint(interrupt);
        if (info != 0) throw std::runtime_error("finite state QR workspace query failed");
        qr_workspace_.resize(static_cast<std::size_t>(workspace_size(query, cols)));
        qr_workspace_rows_ = rows;
    }
    const int work_size = static_cast<int>(qr_workspace_.size());
    lapack::dgeqrf(&nrows, &cols, qr_input_.data(), &nrows, tau.data(), qr_workspace_.data(), &work_size, &info);
    checkpoint(interrupt);
    if (info != 0) throw std::runtime_error("finite state QR compression failed");
    DenseMatrix result(m, m);
    auto& result_values = result.values();
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = j; i < m; ++i) result_values[i + m * j] = checked(qr_input_[j + rows * i]);
    }
    factor_ = std::move(result);
    ++compression_count_;
}

void FiniteStateFactor::update(const std::vector<double>& gain,
                               const std::vector<double>& z, double h,
                               const Interrupt& interrupt) {
    check_variance(h);
    const auto m = dimension();
    if (gain.size() != m) throw std::invalid_argument("finite state gain dimensions must agree");
    for (const auto value : gain) if (!std::isfinite(value)) throw std::invalid_argument("finite state gain must be finite");
    const auto w = observation_loading(z, interrupt);
    const auto rows = m + 1;
    qr_input_.resize(rows * m);
    const auto& factor = factor_.values();
    const auto root_h = std::sqrt(h);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            qr_input_[j + rows * i] = checked(static_cast<long double>(factor[i + m * j]) -
                                              static_cast<long double>(gain[i]) * w[j]);
        }
        qr_input_[m + rows * i] = checked(static_cast<long double>(root_h) * gain[i]);
    }
    compress_transpose(rows, interrupt);
}

void FiniteStateFactor::predict(const DenseMatrix& transition,
                                const DenseMatrix& noise_loading,
                                const Interrupt& interrupt) {
    transition.validate(interrupt);
    noise_loading.validate(interrupt);
    const auto m = dimension();
    if (transition.rows() != m || transition.columns() != m ||
        noise_loading.rows() != m || noise_loading.columns() > max_states) {
        throw std::invalid_argument("finite state prediction dimensions must agree and have at most 32 disturbances");
    }
    // The monthly companion transition is mostly zeros. Store its nonzero
    // column indices for each row in their original accumulation order.
    const auto& transition_values = transition.values();
    const auto& factor_values = factor_.values();
    std::array<std::size_t, max_states + 1> row_offsets;
    std::array<std::size_t, max_states * max_states> nonzero_columns;
    std::size_t count = 0;
    for (std::size_t i = 0; i < m; ++i) {
        row_offsets[i] = count;
        for (std::size_t k = 0; k < m; ++k)
            if (transition_values[i + m * k] != 0) nonzero_columns[count++] = k;
    }
    row_offsets[m] = count;
    const auto rows = m + noise_loading.columns();
    qr_input_.resize(rows * m);
    for (std::size_t j = 0; j < m; ++j) {
        checkpoint(interrupt);
        for (std::size_t i = 0; i < m; ++i) {
            long double value = 0;
            for (std::size_t pos = row_offsets[i]; pos < row_offsets[i + 1]; ++pos) {
                const std::size_t k = nonzero_columns[pos];
                value += static_cast<long double>(transition_values[i + k * m]) * factor_values[k + j * m];
            }
            qr_input_[j + rows * i] = checked(value);
        }
    }
    const auto& disturbance = noise_loading.values();
    for (std::size_t j = 0; j < noise_loading.columns(); ++j) {
        for (std::size_t i = 0; i < m; ++i) qr_input_[m + j + rows * i] = disturbance[i + m * j];
    }
    compress_transpose(rows, interrupt);
}

DenseMatrix FiniteStateFactor::covariance(const Interrupt& interrupt) const {
    const auto m = dimension();
    DenseMatrix result(m, m);
    for (std::size_t j = 0; j < m; ++j) {
        checkpoint(interrupt);
        for (std::size_t i = 0; i <= j; ++i) {
            long double value = 0;
            for (std::size_t k = 0; k < m; ++k) value += static_cast<long double>(factor_(i, k)) * factor_(j, k);
            result(i, j) = result(j, i) = checked(value);
        }
    }
    return result;
}

// Item-29 design D4 helpers (docs/item29-diffuse-filter-probe-2026-11-26.md
// §6; qualification dossiers tmp_opencode/item29-d4/ and
// tmp_opencode/item29-d4b/). Every accumulation is long double scratch with a
// binary64 contract (docs/contracts/extended-precision-range-contract.md); the
// band is the house 512*eps*dimension residual tolerance used by the diffuse
// projection checks. The helpers only measure or replace this factor when the
// caller's premise holds — a rejected certificate changes nothing here, and
// the caller (the canonical SARIMA filter) records the route observably.
bool FiniteStateFactor::certify_observation_nullity(const std::vector<double>& z,
                                                    double& relative_residual) const {
    // P^+Z = 0 through the factor: with P^+ = C C', ||P^+ z|| <= ||C||_2 ||C'z||,
    // so a rounding-level relative ||C'z|| / (||z|| ||C||_F) certifies the
    // nullity premise of the post-pivot induction at every named binary64 step.
    const auto m = dimension();
    if (z.size() != m) throw std::invalid_argument("finite state nullity dimensions must agree");
    for (const auto value : z) if (!std::isfinite(value)) throw std::invalid_argument("finite state nullity coefficients must be finite");
    std::vector<double> projected(m, 0);
    for (std::size_t j = 0; j < m; ++j) {
        long double value = 0;
        for (std::size_t i = 0; i < m; ++i) value += static_cast<long double>(z[i]) * factor_(i, j);
        projected[j] = checked(value);
    }
    double projected_norm = 0, loading_norm = 0, factor_norm = 0;
    for (const auto value : projected) projected_norm = checked(kalos::detail::glibc_hypot(projected_norm, value));
    for (const auto value : z) loading_norm = checked(kalos::detail::glibc_hypot(loading_norm, value));
    for (const auto value : factor_.values()) factor_norm = checked(kalos::detail::glibc_hypot(factor_norm, value));
    const long double scale = static_cast<long double>(loading_norm) * factor_norm;
    relative_residual = scale == 0 ? 0 : checked(projected_norm / scale);
    return relative_residual <= residual_tolerance * static_cast<double>(m);
}

bool FiniteStateFactor::certify_history_block(std::size_t stationary_rows,
                                              double& relative_residual,
                                              const Interrupt& interrupt) const {
    // The post-pivot induction claims every factor row at or below index
    // stationary_rows - 1 is mathematically zero, so its entries must sit at
    // rounding level relative to the stationary rows. A genuine finite
    // direction there — however weak — keeps the entrywise maximum above the
    // band and fails the certificate closed (no weak-direction rank rule: the
    // premise is structural, not a magnitude decision).
    const auto m = dimension();
    if (stationary_rows == 0 || stationary_rows >= m) {
        relative_residual = 1;
        return false;
    }
    double stationary = 0, history = 0;
    for (std::size_t j = 0; j < m; ++j) {
        checkpoint(interrupt);
        for (std::size_t i = 0; i < m; ++i) {
            const auto magnitude = std::abs(factor_(i, j));
            if (i < stationary_rows) { if (magnitude > stationary) stationary = magnitude; }
            else if (magnitude > history) history = magnitude;
        }
    }
    if (stationary == 0) { relative_residual = history == 0 ? 0 : 1; return history == 0; }
    relative_residual = checked(history / stationary);
    return relative_residual <= residual_tolerance * static_cast<double>(m);
}

void FiniteStateFactor::compress_to_stationary_block(std::size_t stationary_rows,
                                                     const Interrupt& interrupt) {
    // The reduced route stores an r x r factor F with F F' equal to the
    // stationary block M M' of the full factor's Gram (M = the top r rows,
    // r x m). A Householder QR of the WIDE M would only reproduce the
    // column Gram M' M conjugated by Q, not the row Gram, so the sweep runs
    // on the tall transpose: M' = Q R gives M M' = R' R exactly up to the
    // factorization residual, and F = R' is adopted after the reconstruction
    // check below compares F F' against the long-double reference block.
    const auto m = dimension();
    if (stationary_rows == 0 || stationary_rows >= m) {
        throw std::invalid_argument("post-pivot stationary block reduction needs a strict history split");
    }
    const int rows = static_cast<int>(m), columns = static_cast<int>(stationary_rows);
    std::vector<double> block(static_cast<std::size_t>(m) * stationary_rows);
    const auto& factor = factor_.values();
    for (std::size_t i = 0; i < stationary_rows; ++i) {
        for (std::size_t j = 0; j < m; ++j) block[j + rows * i] = factor[i + m * j];
    }
    std::array<double, max_states> tau;
    int info = 0;
    checkpoint(interrupt);
    const int query_size = -1;
    double query = 0;
    lapack::dgeqrf(&rows, &columns, block.data(), &rows, tau.data(), &query, &query_size, &info);
    checkpoint(interrupt);
    if (info != 0) throw std::runtime_error("post-pivot stationary block workspace query failed");
    const auto work_size = workspace_size(query, columns);
    std::vector<double> work(static_cast<std::size_t>(work_size));
    const int work_extent = static_cast<int>(work_size);
    lapack::dgeqrf(&rows, &columns, block.data(), &rows, tau.data(), work.data(), &work_extent, &info);
    checkpoint(interrupt);
    if (info != 0) throw std::runtime_error("post-pivot stationary block compression failed");
    DenseMatrix factor_reduced(stationary_rows, stationary_rows);
    for (std::size_t a = 0; a < stationary_rows; ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            factor_reduced(a, b) = checked(block[b + rows * a]);
        }
    }
    // Reconstruction residual: F F' must reproduce M M' at the house band;
    // the reduced factor is adopted only on a passing check.
    long double error = 0, scale = 0;
    const auto& reduced = factor_reduced.values();
    for (std::size_t a = 0; a < stationary_rows; ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            long double previous = 0, current = 0;
            for (std::size_t k = 0; k < m; ++k) previous += static_cast<long double>(factor[a + m * k]) * factor[b + m * k];
            for (std::size_t k = 0; k < stationary_rows; ++k) current += static_cast<long double>(reduced[a + stationary_rows * k]) * reduced[b + stationary_rows * k];
            error = std::max(error, std::abs(previous - current));
            scale = std::max(scale, std::abs(previous));
        }
    }
    if (error > residual_tolerance * static_cast<long double>(m) * scale) {
        throw std::runtime_error("post-pivot stationary block reduction failed its reconstruction residual check");
    }
    post_pivot_rebuild_residual_ = scale == 0 ? 0 : checked(error / scale);
    factor_ = std::move(factor_reduced);
    ++post_pivot_rebuild_count_;
}

} // namespace detail
} // namespace kalos
