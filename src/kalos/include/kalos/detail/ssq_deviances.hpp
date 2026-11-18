#ifndef KALOS_DETAIL_SSQ_DEVIANCES_HPP
#define KALOS_DETAIL_SSQ_DEVIANCES_HPP

// Bytecode-pinned transcription of the jdplus 3.9.0 deviances machine that
// runs at the differencing-decision estimation seams (item 21-D.2). The JVM
// machine is RegArmaSsqFunction.Evaluation over ConcentratedLikelihood
// Computer.DEFAULT_FULL_COMPUTER: AutoCovarianceFunction (Yule-Walker with
// the LULinearSystemSolver2 column-pivoting kernel and its RationalFunction
// psi series) feeding the internal KalmanFilter (sfilter for a design-free
// model, mfilter plus the one-column HouseholderWithPivoting /
// QRLeastSquaresSolver pair once a constant regressor is present), the
// scaled DeterminantalTerm log-determinant, and the
// ConcentratedLikelihoodWithMissing assembly whose deviances vector
// residual * sqrt(exp(ldet / n)) is the Levenberg-Marquardt residual vector
// and whose ssqErr * factor is the sum-of-squares target. Every loop,
// comparison flavour (dcmp/dcmpl/dcmpg) and operation order is transcribed
// operation by operation from the pinned bytecode; the machine reproduces
// the JVM values bit for bit, including its failure semantics (an
// evaluation the JVM rejects returns std::nullopt).
//
// Coordinates are the physical pinned-model coefficients in the JVM
// convention, [phi(1..p), bphi(1..bp), theta(1..q), btheta(1..bq)]; the AR
// polynomials carry the raw JVM signs (1 + phi(B) storage with phi as
// stored by SarimaModel, the negation from the Kalos convention happening in
// the caller).

#include <cstddef>
#include <optional>
#include <vector>

namespace kalos::detail {

// (1 - B)^d (1 - B^period)^bd as an ascending-coefficient polynomial (the
// pinned UnitRoots.D product).
std::vector<double> pinned_differencing_polynomial(int period, int d, int bd);

// Pinned BackFilter.apply: the filter acts through the reversed weight array
// (weights(k) = f[deg - k]), accumulating over k ascending and skipping zero
// weights, so the highest power multiplies the oldest observation. Returns an
// empty vector when the filter is not shorter than the series (where the
// pinned machine fails on a negative-size array).
std::vector<double> pinned_backfilter_apply(
    const std::vector<double>& f, const std::vector<double>& x);

// The pinned OlsComputer coefficient for a single regressor column:
// HouseholderWithPivot.decompose plus QRLeastSquaresSolver.leastSquares at
// the robustLeastSquares rcond (1e-9), bit-identical to the JVM. A rank-0
// solve yields the zero coefficient.
double pinned_ols_1col(const std::vector<double>& col, const std::vector<double>& y);

// One evaluation point of the pinned deviances machine.
struct PinnedMachinePoint {
    // The LM residual vector: the ConcentratedLikelihoodWithMissing
    // deviances (full residuals scaled by sqrt(factor), non-finite
    // coordinates filtered the way DoubleSeq.select(isFinite) filters).
    std::vector<double> e;
    double ssq_err = 0.0;  // cl.ssq()
    double ldet = 0.0;     // cl.logDeterminant()
    double factor = 0.0;   // exp(ldet / n)
    double ssq_e = 0.0;    // ssq_err * factor, the deviance target
    double log_likelihood = 0.0; // cl.logLikelihood()
    int rank = 0;          // the whitened-design rank (0 or 1)
    double coefficient = 0.0;    // sol.b(0) (design present)
    double uc00 = 0.0;           // sol.unscaledCovariance()(0, 0)
    double tstat = 0.0;          // coefficient / sqrt((uc00 * ssq_err) / n)
};

// Evaluates the machine at one coordinate vector.
//   y          the already-differenced response (the machine never differences)
//   p, bp, q, bq  stationary model order counts; period the seasonal period
//   design     empty for the design-free (meanCorrection(false), sfilter)
//              machine, else the single whitened-model regressor column of
//              y.size() entries (the meanCorrection(true) machine)
//   jvm        the [phi, bphi, theta, btheta] coordinates (p+bp+q+bq)
// std::nullopt where the pinned machine throws (invalid ARIMA variances, a
// singular Yule-Walker system, a failed GLS solve) or where the JVM surface
// rejects the construction.
std::optional<PinnedMachinePoint> pinned_deviances_evaluate(
    const std::vector<double>& y, std::size_t p, std::size_t bp, std::size_t q, std::size_t bq,
    int period, const std::vector<double>& design, const std::vector<double>& jvm);

} // namespace kalos::detail

#endif
