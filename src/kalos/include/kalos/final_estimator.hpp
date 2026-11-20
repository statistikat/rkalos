#ifndef KALOS_FINAL_ESTIMATOR_HPP
#define KALOS_FINAL_ESTIMATOR_HPP

// Item 21-H decision unit for jdplus.x13.base.core.x13.regarima
// .FinalEstimator.test: the post-estimation ARMA reduction cascade. Given
// the estimated SARIMA orders, the concatenated parameter vector
// [phi | bphi | theta | btheta] and the asymptotic covariance diagonal of
// the final (maximum-likelihood) estimation, it decides whether one
// coefficient block loses its last (highest-order) term. The bytecode pin
// lives in docs/contracts/final-estimator-contract.md.
//
// The root checks run the deterministic FilterUtility.checkRoots port
// (detail::check_roots); the pinned JVM companion solver is Random-seeded
// for degree >= 3, which the qualified corpus never reaches.

#include <cstddef>
#include <vector>

namespace kalos {

struct SarimaOrders {
    int p = 0;  // regular AR order
    int bp = 0; // seasonal AR order
    int q = 0;  // regular MA order
    int bq = 0; // seasonal MA order
};

// FinalEstimator.test(RegSarimaModelling). status == 0 keeps the
// specification; status == 1 drops one term (orders carries the reduced
// specification). status >= 2 makes the pinned estimate() abort the round
// loop as failure; with the pinned if-else resets it is unreachable, but
// the count is returned literally.
//
// parameters must hold p + bp + q + bq coefficients in the SarimaModel
// order; asy_covariance_diagonal indexes the same positions (the pinned
// SarimaModelling mapping places the ARMA coefficients first), and may be
// longer when regression coefficients follow. sample_length selects the
// sig2 plateau (<= 150 -> 0.15, else 0.1); tsig is the builder's
// significance limit (the ODM kernel leaves it at 1.0).
struct FinalEstimatorRound {
    int status;
    SarimaOrders orders;
};

FinalEstimatorRound final_estimator_test(const SarimaOrders& orders,
    const std::vector<double>& parameters,
    const std::vector<double>& asy_covariance_diagonal, int sample_length,
    double tsig);

} // namespace kalos

#endif
