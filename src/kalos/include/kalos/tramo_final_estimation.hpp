#ifndef KALOS_TRAMO_FINAL_ESTIMATION_HPP
#define KALOS_TRAMO_FINAL_ESTIMATION_HPP

// Item 23-F decision layer: the post-estimation spec cascade of the pinned
// jdplus 3.9.0 jdplus.tramoseats.base.core.tramo.FinalEstimator (the TRAMO
// side of the 21-H fork; docs/contracts/final-estimation-contract.md keeps
// the x13 side). This unit transcribes the three private decision methods -
// searchur + checkUnitRoots (regular/seasonal unit-root promotion),
// checkCommonRoots (seasonal AR/MA cancellation) and test (the last-
// coefficient reduction cascade) - exactly as the probe adjudicated them
// (tmp_opencode/kalos-m5-pinning.md SECTION 26, corpus
// validation/evidence/native/tramo-final-estimation/).
//
// The fork verdict (26.1): this is NOT the landed final_estimator unit.
// The tramo cascade has no checkRoots stability gate (the 0.95 slot of the
// bytecode is stored and never read), starts its takeover ratio at 99999.0
// (the x13 side uses 999), drives the degree-1/0 promotion pair
// p -> d and bp -> bd from its own searchur, and lets Bq fire without
// updating the ratio. The round-loop driver (estimate(), five rounds over
// the qualified estimation seams) and the multi-round outlier machine are
// the 23-F F.2 layer of this same unit; the process() assembly is tramo_fit
// (F.3).
//
// Coefficient convention: all parameter vectors are RAW jdplus SarimaModel
// coordinates, block order [phi | bphi | theta | btheta] with the AR
// polynomials read as 1 + phi(B) - the mirror().roots() scan therefore
// sees roots at the coefficients themselves and the promotion fires when
// an AR coefficient lands near +1 (root 1/phi near the unit root +1). The
// Kalos-convention [1 - phi(B)] blocks of the estimation layer must be
// flipped by the caller.

#include <string>
#include <vector>

namespace kalos {

// SarimaOrders of the description/module at the decision point; the period
// plays no role in any decision below (the seasonal scan reads the Bphi
// block polynomial directly, exactly as getSeasonalAR() does).
struct TramoFinalOrders {
    int p = 0, d = 0, q = 0;
    int bp = 0, bd = 0, bq = 0;
};

// SarimaOrders.toString(), byte-for-byte: "(p,d,q)(Bp,Bd,Bq)". This is the
// detail payload of the "arima spec changed" info event (26.A5).
std::string tramo_final_orders_to_string(const TramoFinalOrders& orders);

// One complex root (Complex.getRe()/getIm()).
struct TramoFinalComplex {
    double re = 0.0, im = 0.0;
};

// Polynomial.mirror().roots() on the fast solver: the coefficient array
// enters through Polynomial.of/ofInternal semantics (trailing - highest
// order - zeros shrink the effective degree; corpus mr06 pins this), the
// trimmed array is reversed (mirror() itself does not trim again) and the
// roots are computed with the bit-exact FastEigenValuesSolver mirror
// (degree <= 2 closed forms, degree >= 3 the 24-A.1 QR cascade). A QR
// cascade that raises in the JVM (null roots -> NPE) throws
// std::runtime_error here. The input is the RAW ascending coefficient
// array (no mirror applied yet): mirror_roots({1, phi1, ...}) returns the
// roots of reverse(trim({1, phi, ...})).
std::vector<TramoFinalComplex> tramo_mirror_roots(
    const std::vector<double>& ascending_coefficients);

// FinalEstimator.searchur over an explicit root array: counts roots with
// Complex.abs() > ur STRICT (dcmpl: NaN on either side excluded, infinity
// counted), |getIm()| <= 0.05 INCLUSIVE (dcmpg: NaN excluded, infinite
// imaginary part excluded) and getRe() > 0 STRICT (dcmpl: -0.0 excluded).
// Complex.abs() is the scaled modulus of the JVM, NOT std::hypot.
// A null array in the JVM returns 0; an empty array has no roots.
int tramo_searchur_roots(const std::vector<TramoFinalComplex>& roots, double ur);

// checkUnitRoots over one AR block: mirror().roots() of 1 + block(B) fed
// through the count above.
int tramo_searchur(const std::vector<double>& ar_block, double ur);

// FinalEstimator.checkUnitRoots. Gates: (p > 0 && d < 2) scans the regular
// block and (bp > 0 && bd < 1) scans the seasonal block; a nonzero count
// promotes p -> p - 1, d -> d + 1 (bp -> bp - 1, bd -> bd + 1) and the
// caller-side specification rewrite (setSpecification) RESETS every
// parameter to its default - the probe goldens pin the post-promotion
// landing at the -0.1/-0.2 defaults. unchanged == true means the pinned
// method returned true (no promotion; the estimate() driver then runs
// checkCommonRoots and the cascade; a false result skips them and re-
// estimates).
struct TramoUnitRootsRequest {
    TramoFinalOrders orders;
    std::vector<double> phi;  // regular AR block, SarimaModel coordinates
    std::vector<double> bphi; // seasonal AR block
    double ur = 0.96;         // builder default; the kernel passes options.ur
};

struct TramoUnitRootsResult {
    bool unchanged = true;
    bool regular_checked = false;   // gate (p > 0 && d < 2) opened
    bool seasonal_checked = false;  // gate (bp > 0 && bd < 1) opened
    int regular_count = 0;
    int seasonal_count = 0;
    bool regular_fired = false;
    bool seasonal_fired = false;
    TramoFinalOrders orders_after;
};

TramoUnitRootsResult tramo_check_unit_roots(const TramoUnitRootsRequest& request);

// FinalEstimator.checkCommonRoots. The live rule fires iff bp == 1 &&
// bq == 1 and |bphi(1) - btheta(1)| < cancel STRICT (dcmpg on the raw
// lag-1 coefficient difference: NaN never cancels); firing clears bp and
// bq and rewrites the specification (parameters reset). unchanged == true
// is the pinned return value (false feeds nnsig = 2 in the estimate()
// driver). The p/q root arrays are computed and inverted when p != 0 &&
// q != 0 and the results are never read - but the roots() call itself is
// observable through its failure modes (the pinned QR cascade exhausting
// raises where the JVM NPEs), so the calls are executed faithfully and
// only the arrays dropped (26.A3).
struct TramoCommonRootsRequest {
    TramoFinalOrders orders;
    std::vector<double> phi;   // regular AR block (dead-work roots call)
    std::vector<double> bphi;  // seasonal AR block; bphi[0] feeds the rule
    std::vector<double> theta; // regular MA block (dead-work roots call)
    std::vector<double> btheta; // seasonal MA block; btheta[0] feeds the rule
    double cancel = 0.044;    // builder default; the kernel never overrides it
};

struct TramoCommonRootsResult {
    bool unchanged = true;
    bool cancelled = false;
    TramoFinalOrders orders_after;
};

TramoCommonRootsResult tramo_check_common_roots(const TramoCommonRootsRequest& request);

// FinalEstimator.test: the last-coefficient reduction cascade. One block
// (P, Bp, Q, Bq, tested in that order on the running index
// p, p + bp, p + bp + q, p + bp + q + bq) fires when v > 0 STRICT (dcmpl),
// t = |a| / sqrt(v) satisfies dcmpg(t, tsig) < 0 and dcmpg(|a|, sig2) < 0
// (NaN never fires); P then always fires (setting the takeover ratio to t)
// while Bp/Q/Bq additionally require dcmpg(ratio, t) > 0 STRICT. Firing
// Bp clears P and updates the ratio; firing Q clears P and Bp and updates
// the ratio; firing Bq clears P, Bp and Q and does NOT update the ratio
// (26.A5, bytecode offsets 452-462). total = the sum of the four 0/1 flags;
// total == 0 returns with no specification change, otherwise each fired
// block loses exactly one order and the "arima spec changed" info event
// fires with the new orders as its detail. sig2 = sample_length <= 150 ?
// 0.15 : 0.1 on the ESTIMATION DOMAIN length; tsig is the field value
// (1.0 in every live build - the Builder has no setter).
struct TramoFinalEstimatorTestRequest {
    TramoFinalOrders orders;
    std::vector<double> parameters;            // [phi | bphi | theta | btheta]
    std::vector<double> asy_covariance_diagonal; // same positions; may be longer
    int sample_length = 0;
    double tsig = 1.0;
};

struct TramoFinalEstimatorTestResult {
    int total = 0;              // getChangedParametersCount contribution (0..4)
    bool fired_p = false;       // flags AFTER the takeover clears (decrement drivers)
    bool fired_bp = false;
    bool fired_q = false;
    bool fired_bq = false;
    double ratio_final = 0.0;   // audit: ratio after the cascade (0 when nothing set it)
    bool spec_changed = false;  // total > 0: setSpecification + info event
    std::string info_detail;    // orders_after string when spec_changed, else empty
    TramoFinalOrders orders_after;
};

TramoFinalEstimatorTestResult tramo_final_estimator_test(
    const TramoFinalEstimatorTestRequest& request);

// ---------------------------------------------------------------------------
// Item 23-F F.2: the live FinalEstimator.estimate(RegSarimaModelling) driver.
// One call replays the pinned round loop end to end (decodes in
// tmp_opencode/kalos-m5-pinning.md SECTION 29, per-round confirmation in the
// Scratch23n corpus validation/evidence/native/tramo-final-estimation/
// n23-capture-*.out): the single push("full estimation")/pop pair, at most
// five round bodies, and per round the RegSarimaModelling.estimate hand-over
// to RegSarimaComputer over the HannanRissanen start and the pinned
// Levenberg-Marquardt chain, the zero-parameter concentrated-likelihood
// shortcut (DEFAULT_COMPUTER, "Gls failed" refused), the exact-derivatives
// Point tail of RegArmaProcessor.compute (calcResiduals at the landing
// allCoefficients, the design-free ssq_e deviance re-run under the
// SarimaMapping2 stencil, scaled by 0.5 * nobs / fx on the cdf nobs;
// asymptoticCovariance = SymmetricMatrix.inverse), then the checkUnitRoots /
// checkCommonRoots / test cascade with the setSpecification resets (empty
// parameters, nulled estimation). The series is rebuilt bit-exactly from the
// family facts (SECTION 29.6: fdlibm sin/cos reproduce the probe's
// StrictMath series on
// all 1008 values); the design is the optional raw constant column and the
// response carries no missing values, exactly as the fe_est corpus.
// ---------------------------------------------------------------------------

// The SECTION-1 input family behind the fe_est cases: base = 50 + 0.05 t +
// 5 sin(2 pi t / 12) + 3 cos(4 pi t / 12) + 0.5 sin(14 pi t / 12) with
// t = 1..n; LevelShift adds 20 from t = 49; TwoOutliers adds 12 at index 10
// and subtracts 18 at index 60. Families reaching past the response throw
// std::invalid_argument (fail-closed; the JVM probe would fail the same
// case before estimate()).
enum class TramoEstFamily {
    None,
    LevelShift,
    TwoOutliers,
};

struct TramoFinalEstimatorRunRequest {
    TramoEstFamily family = TramoEstFamily::None;
    int n = 0;
    bool mean = true;                    // desc.setMean (concentrated constant)
    double eps = 1e-4;                   // FinalEstimator builder precision
    double ur = 0.96;                    // builder unitRootThreshold
    int pass = 0;                        // builder pass
    bool ami = true;                     // builder ami
    bool outliers = false;               // builder outliers
    bool airline = false;                // setAirline(true): (0,1,1)(0,1,1)
    TramoFinalOrders orders;             // the specification when not airline
};

// Per-round audit rows (the probe's m R<k> lines): the mapping dimension,
// the HannanRissanen start provenance (fitted/zerovariance/defaulted, empty
// on the zero-parameter shortcut), the landing in JVM coordinates, the
// ladder decisions and the round's exit.
struct TramoFinalEstimatorRoundRecord {
    int round = 0;
    int dim = 0;
    std::string hr_source;
    std::vector<double> landing; // JVM [phi | bphi | theta | btheta]; empty on dim == 0
    int checkur = -1;            // -1 not run, else unchanged flag 0/1
    int checkcr = -1;
    int test_total = -1;
    int nnsig_after = 0;
    bool shortcircuit = false;   // the dim == 0 or !ami early return
    std::string outcome;         // continue | return_true | return_false
    // Point-tail witnesses (empty when the round had no landing point): the
    // exact-derivatives score grad * (-0.5 * nobs / fx), the information
    // hessian * (0.5 * nobs / fx) row-major (both from the final_information
    // producer over the residual-machine deviance), and the tail fx =
    // Likelihood.ssq() * factor at the landing.
    std::vector<double> tail_score;
    std::vector<double> tail_information;
    double tail_ssq_e = 0.0;
};

struct TramoFinalEstimatorRunResult {
    bool returned = false;  // the estimate() return value
    int nnsig = 0;          // getChangedParametersCount after the call
    bool mean_out = false;
    TramoFinalOrders orders_after;
    std::vector<double> parameters_after; // desc.arima() parameters; empty
                                          // after a setSpecification
    // The final estimation state: 0 = null (a setSpecification was last),
    // 1 = full estimation with landing point, 2 = zero-parameter estimation
    // (no point: asymptoticCovariance would throw, the probe's NPE seam).
    int cl_state = 0;
    int cl_dim = 0;                 // concentratedLikelihood.dim()
    double cl_ssq = 0.0;            // concentratedLikelihood.ssq()
    std::vector<double> cl_coefs;   // allCoefficients
    std::vector<double> covariance_diagonal; // state 1 only
    // Point tail witness: the scaled information matrix (row-major dim x dim)
    // the last point-tail round inverted; empty when no point tail ran.
    std::vector<double> landing_information;
    // ProcessingLog mirror: (kind, depth) with 1 push, 2 pop, 3 remark,
    // 4 info; a normal call is exactly {(1, 0), (2, 0)} per witness.
    std::vector<std::pair<int, int>> events;
    std::vector<TramoFinalEstimatorRoundRecord> rounds;
};

// Static input violations (non-positive n, negative orders, a non-airline
// specification without parametersCount match, non-finite eps/ur) throw
// std::invalid_argument. A round whose pinned machine fails (the
// "estimation faied" remark path of estimate()) is NOT an exception here:
// it ends the call with returned == false, the remark event and the state
// left where the throw happened, exactly as the pinned catch (the corpus
// never takes this path).
TramoFinalEstimatorRunResult tramo_final_estimator_run(
    const TramoFinalEstimatorRunRequest& request);

// Run the same monthly estimator on complete observations in the model's
// level/log channel, with an optional mean. n may be zero (inferred) or
// equal the series length. The synthetic family setting is unused.
TramoFinalEstimatorRunResult tramo_final_estimator_run(
    const TramoFinalEstimatorRunRequest& request, const std::vector<double>& model_series);

} // namespace kalos

#endif
