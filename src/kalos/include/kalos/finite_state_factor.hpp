#ifndef KALOS_FINITE_STATE_FACTOR_HPP
#define KALOS_FINITE_STATE_FACTOR_HPP

#include "linear_algebra.hpp"

#include <cstddef>
#include <vector>

namespace kalos {
namespace detail {

struct FiniteObservation {
    double variance;
    std::vector<double> pz;
};

// A unit-scale finite covariance P=C C'. Unlike a Cholesky covariance API,
// zero/rank-deficient factors are allowed. QR compression keeps all columns;
// it never decides the diffuse rank or drops weak positive finite directions.
class FiniteStateFactor {
public:
    explicit FiniteStateFactor(const DenseMatrix& covariance,
                               const Interrupt& interrupt = {});
    std::size_t dimension() const { return factor_.rows(); }
    std::size_t compression_count() const { return compression_count_; }
    std::size_t clipped_negative_count() const { return clipped_negative_count_; }
    double initialization_residual() const { return initialization_residual_; }

    FiniteObservation observation(const std::vector<double>& z, double h,
                                   const Interrupt& interrupt = {}) const;
    void update(const std::vector<double>& gain, const std::vector<double>& z,
                 double h, const Interrupt& interrupt = {});
    void predict(const DenseMatrix& transition, const DenseMatrix& noise_loading,
                  const Interrupt& interrupt = {});
    DenseMatrix covariance(const Interrupt& interrupt = {}) const;

    // Item-29 design D4 (docs/item29-diffuse-filter-probe-2026-11-26.md §4
    // E2, §6 D4). A canonical zero-observation-variance SARIMA filter proves
    // P^+Z = 0 after every observed update, so once its diffuse pivots are
    // exhausted the history rows of this factor carry only rounding-level
    // residuals of a mathematically zero block (the E2 induction of
    // docs/item29-diffuse-filter-probe-2026-11-26.md §4). The two certify_
    // helpers measure that premise at the current factor and return whether
    // it sits at rounding level against the house 512*eps*dimension band;
    // they never change the factor. compress_to_stationary_block replaces the
    // factor by the r-wide Householder compression of its stationary rows
    // (the one re-association of the certified-zero drop); it fails loudly on
    // a rejected precondition, never silently falls back. Selection between
    // the full-m and reduced-r recursions is owned by the filter and made
    // observable through DiffuseNumerics — both recursions stay exact.
    bool certify_observation_nullity(const std::vector<double>& z,
                                     double& relative_residual) const;
    bool certify_history_block(std::size_t stationary_rows,
                               double& relative_residual,
                               const Interrupt& interrupt = {}) const;
    void compress_to_stationary_block(std::size_t stationary_rows,
                                      const Interrupt& interrupt = {});
    std::size_t post_pivot_rebuild_count() const { return post_pivot_rebuild_count_; }
    double post_pivot_rebuild_residual() const { return post_pivot_rebuild_residual_; }

private:
    std::vector<double> observation_loading(const std::vector<double>& z,
                                            const Interrupt& interrupt) const;
    void compress_transpose(std::size_t rows, const Interrupt& interrupt);

    DenseMatrix factor_;
    std::size_t post_pivot_rebuild_count_ = 0;
    double post_pivot_rebuild_residual_ = 0;
    // QR uses the transposed loading directly. Both buffers belong to this filter.
    std::vector<double> qr_input_;
    std::vector<double> qr_workspace_;
    std::size_t qr_workspace_rows_ = 0;
    std::size_t compression_count_ = 0;
    std::size_t clipped_negative_count_ = 0;
    double initialization_residual_ = 0;
};

} // namespace detail
} // namespace kalos
#endif
