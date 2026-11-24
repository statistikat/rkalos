#ifndef KALOS_DETAIL_SSQ_LEVMAR_HPP
#define KALOS_DETAIL_SSQ_LEVMAR_HPP

// Internal seam for the pinned jdplus 3.9.0 sum-of-squares Levenberg-Marquardt
// machine (item 21-D.2): LevenbergMarquardtMinimizer (QR variant) over the
// SsqNumericalDerivatives one-sided stencil of the deviances residual vector
// of the RegArmaSsqFunction, with the SarimaMapping2 parameter domain. The
// coordinates are the physical pinned-model coefficients in JVM convention,
// [phi(1..p), bphi(1..bp), theta(1..q), btheta(1..bq)]; the caller owns the
// sign flip to the Kalos ARMA convention. The builder constants are the pinned
// LmBuilder defaults (eps1 = eps2 = 1e-15, tau = 1e-3, maxIter = 200,
// nu = 4, acceptance floor 1e-4, mu damping floor 0.3333333333333333, at most
// 100 damping retries per iteration and g_nsteps = 2 epsilon deflations).
// The least-squares solve reuses the pinned fastLeastSquares transcription
// (the same Householder2 chain as the deviances machine), keeping the LM
// trajectory on the pinned bits wherever the trajectory stays in the
// full-rank interior; near rank drops and domain stabilizations the landing
// can still sit in the same stationary basin at a ulp offset, which the
// corpus measures as a landing-drift tier (documented in
// docs/contracts/differencing-decision-contract.md).

#include "kalos/numerics.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace kalos::detail {

// The pinned QRLeastSquaresSolver.fastLeastSquares machine (unpivoted
// Householder2 chain, UpperTriangularMatrix.rank(R, 0), solveUx), shared by
// the sum-of-squares minimizer's QR step and the HannanRissanen pipeline's
// regressions: a LAPACK solver only agrees at full rank and drifts the
// pinned trajectories by ulps.
struct PinnedLsSolution {
    int rank = 0;
    std::vector<double> x; // dim entries, zero padded when rank < dim
    std::vector<double> residual_tail; // Q'b below the rank, for getSsqErr
    bool singular = false;
};

PinnedLsSolution pinned_least_squares(std::vector<double> qr, const std::vector<double>& e,
    std::size_t m, std::size_t dim);

// Pinned ParamValidation outcomes of SarimaMapping2.validate.
enum class SsqParamValidation { Invalid, Valid, Changed };

// Stationary domain over JVM-ordered coordinates. The default SarimaMapping2
// bounds AR roots away from one. GlsSarimaComputer uses SarimaMapping and
// inverts roots strictly inside the unit circle.
// The limits are the pinned of() defaults: arLimit 0.99999, maLimit 1.0,
// checkAll true, epsilon sqrt(2.220446E-16) directed toward zero.
class SarimaSsqDomain {
public:
    SarimaSsqDomain(int p, int bp, int q, int bq, bool gls_mapping = false);

    std::size_t dim() const { return offsets_.back(); }
    const std::vector<std::size_t>& offsets() const { return offsets_; }
    const std::vector<std::size_t>& sizes() const { return sizes_; }

    // Pinned SarimaMapping2.epsilon: +/- sqrt(2.220446E-16) * max(1, |x|),
    // directed toward zero.
    double epsilon(const std::vector<double>& x, std::size_t j) const;
    // Pinned checkBoundaries: FilterUtility.checkRoots per block.
    bool check_boundaries(const std::vector<double>& x) const;
    // Pinned validate: FilterUtility.stabilize per block in place; the length
    // mismatch (Invalid) cannot occur here because the coordinate vector is
    // built to the block layout.
    SsqParamValidation validate(std::vector<double>& x) const;

private:
    std::vector<std::size_t> sizes_;   // [p, bp, q, bq]
    std::vector<std::size_t> offsets_; // exclusive prefix, size 5
    bool gls_mapping_;
    std::vector<double> limits_;       // per block
};

// One evaluated point: e is the pinned deviances vector (cl.e() scaled by
// sqrt(cl.factor())) and ssq_e is the pinned objective cl.ssq() *
// cl.factor(). Providers return std::nullopt exactly where the pinned
// Evaluation constructor throws (non-stationary filter, failed Glsl solve)
// or reports non-finite values.
struct SsqPoint {
    std::vector<double> e;
    double ssq_e;
};

class SsqVectorFunction {
public:
    virtual ~SsqVectorFunction() = default;
    virtual std::optional<SsqPoint> ssq_evaluate(const std::vector<double>& x) const = 0;
    virtual const SarimaSsqDomain& domain() const = 0;
};

struct SsqLmResult {
    std::vector<double> parameters; // landing coordinates
    std::vector<double> e;          // landing residual vector
    double objective = 0.0;         // landing ssq_e
    int iterations = 0;             // pinned getIterationsCount()
    int evaluations = 0;            // probe evaluations (diagnostic)
    // Pinned stop code: 0 none, 1 gradient, 2 function/step precision,
    // 4 step overflow, 5 mu overflow, 6 initial scale, 7 non-finite.
    int stop = 0;
    bool converged = false;     // pinned calc() return
    bool has_converged = false; // pinned hasConverged(): stop == 2
};

// Pinned LevenbergMarquardtMinimizer.minimize + calc + iterateQR over the
// given function. fn_precision is the pinned builder functionPrecision (eps3);
// eps1, eps2, tau, maxIter keep the pinned defaults. A provider that fails at
// the start point throws std::runtime_error (the pinned Evaluation throws
// before the minimizer runs).
SsqLmResult minimize_ssq_levmar(const SsqVectorFunction& fn, std::vector<double> start,
    double fn_precision, const Interrupt& interrupt = {});

} // namespace kalos::detail

#endif
