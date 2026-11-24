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
    // The full least-squares solution vector (the pinned getB(), zero-padded
    // past the rank) = ConcentratedLikelihood.allCoefficients, filled by the
    // multi machine (item 23-F F.2); the 1x1 conveniences above mirror it.
    std::vector<double> b;
};

// Optional QR reporting at a fitted point. Optimizer evaluations can skip
// both the raw residual tail and the regression covariance calculation.
struct PinnedRegressionDiagnostics {
    std::vector<double> qr_residuals;
    std::vector<double> unscaled_covariance; // row major in input column order
};

// The two ConcentratedLikelihoodComputer shapes (item 23-F F.2): Full is
// DEFAULT_FULL_COMPUTER (xfixed = false: the whole whitened design enters
// HouseholderWithPivoting with column pivoting and there is no rank
// requirement); Fixed is DEFAULT_COMPUTER (xfixed = true: every design column
// is fixed and a rank below the column count is the "Gls failed" throw,
// mirrored as nullopt).
enum class PinnedComputerShape {
    Full,
    Fixed,
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

// The multi-regressor generalization of the same machine (item 23-F F.2):
// design carries the full whitened-model regressor block (each column
// y.size() entries, the differenced-domain view the pinned
// ConcentratedLikelihoodComputer sees: missing point columns first, the
// meanCorrection constant, then the differenced user columns). The QR runs
// the pivoting hwp_n decomposition of the DEFAULT_FULL_COMPUTER shape the
// RegArmaSsqFunction surface uses (xfixed = false, so nfixed = nmissing:
// full column pivoting) and the rcond-1e-12 solve (the same primitives the
// OLS and DEFAULT_COMPUTER seams replay), the residual assembly keeps the
// per-column dneg/addAY order and the e/ssqE/logLikelihood assembly is
// byte-identical to the single-column machine. A rank loss is not an error
// under the Full shape (the xfixed rank throw belongs to the landing seam,
// not to this computer); nullopt mirrors the pinned throws. pt.coefficient/
// uc00/tstat (the 1x1 conveniences) are not filled here; pt.b carries the
// full solution vector (allCoefficients). The shape parameter selects the
// computer: Fixed runs the DEFAULT_COMPUTER shape of the zero-parameter
// FinalEstimator round - every whitened column fixed, a rank below the
// column count mirroring the "Gls failed" throw as nullopt.
std::optional<PinnedMachinePoint> pinned_deviances_evaluate_multi(
    const std::vector<double>& y, std::size_t p, std::size_t bp, std::size_t q, std::size_t bq,
    int period, const std::vector<std::vector<double>>& design, const std::vector<double>& jvm,
    PinnedComputerShape shape = PinnedComputerShape::Full);

// At a fitted point the extra output includes the compact cl.e() and the
// full regression covariance. An empty design uses the scalar Kalman path.
std::optional<PinnedMachinePoint> pinned_deviances_evaluate_multi(
    const std::vector<double>& y, std::size_t p, std::size_t bp, std::size_t q, std::size_t bq,
    int period, const std::vector<std::vector<double>>& design, const std::vector<double>& jvm,
    PinnedComputerShape shape, PinnedRegressionDiagnostics* diagnostics);

// Pinned Polynomial.times (the full convolution behind the rational filter
// compositions of the item 23-B scan machine: composing the pi weights with
// an outlier filter representation).
std::vector<double> pinned_polynomial_times(
    const std::vector<double>& a, const std::vector<double>& b);

// LeastSquaresDivision.divide + isExact behind Polynomial.SimplifyingTool
// .simplifyExact (the exact-cancel step of RationalBackFilter.times): the
// divisor shifted across n rows, Householder2 least squares (no pivoting,
// rcond 0), quotient exact iff ssqErr / b.size() is below sqrt(epsilon).
// nullopt marks "divide refused" (a divisor longer than the dividend, a
// MatrixException) or "not exact", both of which send the simplifier to its
// roots pass.
std::optional<std::vector<double>> pinned_ls_divide(
    const std::vector<double>& a, const std::vector<double>& b);

// Pinned RationalFunction.coefficients(k): the first k power-series weights
// of num / den, the same literal recurrence the deviances machine runs for
// its psi series and the scan machine for the pi / (pi * rep) columns.
std::vector<double> pinned_rational_weights(
    const std::vector<double>& num, const std::vector<double>& den, std::size_t k);

// QRLeastSquaresSolver.robustLeastSquares behind OlsComputer.compute
// (HouseholderWithPivoting with nfixed = 0, column pivoting, rcond 1e-9):
// returns the getB() vector (zero-padded on rank deficiency). nullopt marks
// the null result for a response norm below Constants.epsilon or a
// MatrixException from the triangular tail, which the driver turns into
// "Ols failed".
struct PinnedOlsFit {
    std::vector<double> coefficients;
    std::size_t rank;
    double residual_ssq;
};

// QRLeastSquaresSolver.fastLeastSquares: fixed column order, zero rank tolerance.
std::optional<PinnedOlsFit> pinned_fast_ols_fit(const std::vector<double>& y,
    const std::vector<std::vector<double>>& xcols);

std::optional<PinnedOlsFit> pinned_ols_fit(const std::vector<double>& y,
    const std::vector<std::vector<double>>& xcols);
std::optional<PinnedOlsFit> pinned_ols_fit(const std::vector<double>& y,
    const std::vector<std::vector<double>>& xcols, std::vector<double>* unscaled_covariance);

std::optional<std::vector<double>> pinned_ols_coefficients(const std::vector<double>& y,
    const std::vector<std::vector<double>>& xcols);

// ConcentratedLikelihoodComputer.DEFAULT_COMPUTER over a multi-column design
// at one coordinate vector (the item 23-B initmodel seam): steady-state
// multi-use KalmanFilter whitening of the response and every design column,
// HouseholderWithPivoting with all columns fixed and
// QRLeastSquaresSolver.leastSquares at the machine rcond 1e-12. nullopt
// mirrors the "Gls failed" EcoException (rank below the column count) or a
// whitening/solve failure.
std::optional<std::vector<double>> pinned_gls_coefficients(const std::vector<double>& y,
    std::size_t p, std::size_t bp, std::size_t q, std::size_t bq, int period,
    const std::vector<std::vector<double>>& xcols, const std::vector<double>& jvm);

// The same DEFAULT_COMPUTER evaluation with the QRLeastSquaresSolution parts
// the outlier rounds' updateLikelihood reads (item 23-F F.2): the raw
// coefficient vector b (the ConcentratedLikelihoodWithMissing allCoefficients
// - missing-correction entries first), the unscaledCovariance() diagonal
// AFTER the missing-slice view (concentrated over the nx non-missing design
// columns, entry k = the diagonal of SymmetricMatrix.UUt applied to
// UpperTriangularMatrix.inverse of the leading rank block, the pivot applied
// when the decomposition carries one) and the raw ssqErr.
struct PinnedGlsPoint {
    std::vector<double> b;         // solution coefficients, design order
    std::vector<double> uc_diag;   // nx diagonal entries (missing slice)
    double ssq_err = 0.0;          // QRLeastSquaresSolution.getSsqErr()
};

std::optional<PinnedGlsPoint> pinned_gls_point(const std::vector<double>& y,
    std::size_t p, std::size_t bp, std::size_t q, std::size_t bq, int period,
    const std::vector<std::vector<double>>& xcols, const std::vector<double>& jvm,
    std::size_t nmissing);

// FastKalmanFilter.process(x).e(): the Durbin-Magenau innovations of the
// stationary model at the coordinates over x, with the pinned steady-state
// fast switch. nullopt mirrors the null Likelihood the pinned filter returns
// when a predicted variance is NaN or negative (where the caller
// dereferences it).
std::optional<std::vector<double>> pinned_fast_kalman_residuals(const std::vector<double>& x,
    std::size_t p, std::size_t bp, std::size_t q, std::size_t bq, int period,
    const std::vector<double>& jvm);

// FastKalmanFilter.fastFilter(x): the machine the differencing module's
// computeTMean runs (item 23-C). Distinct from the process(x) machine above:
// the warmup runs only max(arLen - 1, maLen - 1) steps, the stickiness fires
// when every c4 entry reaches |c4| <= eps, the h update runs on every warmup
// step (including after the switch), a NaN h after the update aborts the
// machine (the ArimaException), and every remaining step is the plain
// Durbin-Magenau recursion on raw series lags:
// e[t] = x[t] + sum_j ar[j] * x[t - j] - sum_j ma[j] * e[t - j]. A series
// shorter than the warmup trips the JVM array-index failure; that mirrors as
// nullopt.
std::optional<std::vector<double>> pinned_fast_filter_residuals(const std::vector<double>& x,
    std::size_t p, std::size_t bp, std::size_t q, std::size_t bq, int period,
    const std::vector<double>& jvm);

// The design-free ConcentratedLikelihoodComputer innovations (the item 23-B
// allCoefficients-empty branch: KalmanFilter.prepare plus sfilter, the
// meanCorrection(false) machine already bit-qualified through
// pinned_deviances_evaluate). nullopt mirrors the failed machine.
std::optional<std::vector<double>> pinned_kalman_innovations(const std::vector<double>& y,
    std::size_t p, std::size_t bp, std::size_t q, std::size_t bq, int period,
    const std::vector<double>& jvm);

} // namespace kalos::detail

#endif
