#ifndef KALOS_FINAL_INFORMATION_HPP
#define KALOS_FINAL_INFORMATION_HPP

// Item 21-H.0: the final-information producer of the pinned jdplus 3.9.0
// regarima stack - the tail of RegArmaProcessor.compute that turns a final
// point into the score, the observed information and the asymptotic
// covariance (docs/contracts/final-information-contract.md).
//
// The stencil is NumericalDerivatives built as (point, false, false) by the
// IFunctionPoint.derivatives() default: per coordinate the domain epsilon,
// NumericalDerivatives.checkepsilon (boundary halving / flip / zeroing), the
// forward node newval(i, eps), then calcgrad ((fp - fx) / eps with the
// fp == fx zero guard) and calch (diagonal (f(+h) - 2 fx) + f(-h) over h =
// sqrt(|eps|) and the four-corner cross formula). RegArmaProcessor scales
// the raw arrays into the score and the information, and the asymptotic
// covariance is SymmetricMatrix.inverse: deepClone, Cholesky (lapack
// Cholesky.lcholesky, tol 0), LowerTriangularMatrix.inverse, LtL + fromLower,
// falling back to CroutDoolittle.decompose(Constants.getEpsilon()) plus
// LUDecomposition.solve(I) when that chain raises MatrixException. Every
// loop, dcmp/dcmpl/dcmpg choice and operation order is transcribed
// operation by operation from the pinned bytecode.
//
// The deviance itself is injected: the producer only ever calls the pinned
// RegArmaFunction, whose value here is the callback. A nullopt return means
// the JVM evaluation raised and NumericalDerivatives.newval caught it, which
// the stencil replaces by fx; a nullopt at the final point itself rejects the
// round (status final_point_rejected). Coordinates are the physical pinned
// coefficients [phi | bphi | theta | btheta] with the raw SarimaModel signs,
// exactly as the stencil sees them; any Kalos-convention sign flip and the
// differencing/nobs bookkeeping (nobs = len(y_d) - xcols - dim) belong to the
// caller (the 21-H.1 driver). Design columns are out of scope: the pinned
// tail reaches them only through the CompoundMapping coordinate extension,
// which is qualified with that driver slice.

#include <functional>
#include <optional>
#include <vector>

#include "kalos/final_estimator.hpp" // SarimaOrders

namespace kalos {

// SarimaMapping2's static epsilon: Math.sqrt(2.220446E-16).
double final_information_step();

enum class final_information_status {
    // the tail completed; score, information, asymptotic covariance filled.
    available,
    // both SymmetricMatrix.inverse branches raised; fx, gradient, score,
    // hessian and information stay valid, the covariance stays empty.
    inverse_failed,
    // the deviance at the final point itself was rejected (the JVM tail
    // throws out of the ssq evaluation before the stencil runs).
    final_point_rejected,
};

struct final_information_request {
    // stationary-model order counts; point must hold p + bp + q + bq
    // coordinates.
    SarimaOrders orders;
    std::vector<double> point;
    // the pinned compute() nobs argument, used literally (the JVM accepts
    // zero and negative values and scales with them).
    int nobs = 0;
    // SarimaMapping2 knobs.
    double ar_limit = 0.99999;
    double ma_limit = 1.0;
    bool check_all = true;
    double epsilon = final_information_step();
    // the pinned deviance machine at one coordinate vector; nullopt where
    // the JVM evaluation throws.
    std::function<std::optional<double>(const std::vector<double>&)> deviance;
};

struct final_information_result {
    final_information_status status = final_information_status::final_point_rejected;
    // the stencil epsilons after checkepsilon (always computed).
    std::vector<double> eps;
    double fx = 0.0;
    // NumericalDerivatives.calcgrad output.
    std::vector<double> gradient;
    // gradient * ((-0.5 * nobs) / fx).
    std::vector<double> score;
    // NumericalDerivatives.calch output, row-major flat.
    std::vector<double> hessian;
    // hessian * ((0.5 * nobs) / fx), row-major flat.
    std::vector<double> information;
    // SymmetricMatrix.inverse(information), row-major flat, and its
    // diagonal; empty unless status == available.
    std::vector<double> asymptotic_covariance;
    std::vector<double> asymptotic_covariance_diagonal;
};

final_information_result final_information(const final_information_request& req);

} // namespace kalos

#endif
