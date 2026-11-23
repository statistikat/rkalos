#ifndef KALOS_TRAMO_SEATS_MODEL_VALIDATOR_HPP
#define KALOS_TRAMO_SEATS_MODEL_VALIDATOR_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

// Item 24-A: the SEATS model-validation and canonical-decomposition DECISION
// mirror. This unit reproduces, operation for operation, what the LIVE pinned
// jdplus 3.9.0 path decides between the SeatsKernel.validate seam and the
// SeatsKernel.decomposeModel seam:
//
//   * DefaultModelValidator - the XL parameter-cutoff machine (simplifyModel
//     trailing trim with cutOff, changeMA/stabilizeMA, changeAR, the
//     fixMaUnitRoots projection), bit-exact;
//   * the toolkit selector chain (TrendCycleSelector + SeasonalSelector as
//     SeatsToolkit.of wires them from the DecompositionSpec knobs) and the
//     RootDecomposer linear-system structure (the FastMatrix assembly and the
//     LULinearSystemSolver2 transcription give bit-exact symmetric filters);
//   * the decomposability decision of UcarimaModel.setVarianceMax(-1, noisy):
//     the Spectrum$Minimizer minimum-sum test (minimizeByRoots with the
//     FastEigenValuesSolver closed forms / minimizeByGrid with the
//     denominator-bracketed Newton iteration) whose sign decides whether the
//     toolkit decomposer yields a UcarimaModel or null, and the SeatsKernel
//     decomposeModel state machine flags/events that follow ("arima parameters
//     cut off", "non decomposable model", "canonical decomposition failed",
//     "noisy model used", modelChanged/parametersCutOff).
//
// It does NOT implement what happens AFTER a null decomposition: the
// DefaultModelApproximator approximation modes are decision register D6 and
// stay unimplemented here - the unit reports approximation_required and stops
// (the live record pins what the kernel then did). The spectral
// factorization of the split components (SymmetricFilter.factorize -> the
// one-sided MA / variance facts), the variance allocation itself and the
// UcarimaModel post-allocation values are item 24-B scope; the unit exposes
// the raw (pre-allocation) component structure only.
//
// Envelope (registers D5/D11, monthly-first, faithful-to-capture): the unit
// boundary accepts the fixed-parameter (all-Fixed) SarimaSpec shape the probe
// captured, monthly (period 12) and quarterly (period 4) series. The degree
// >= 3 root cascade (FastEigenValuesSolver.qr, ported bit-exact at 24-A.1)
// retired the old D11 refusals: stationary bphi, stationary AR degree > 2 and
// MA degree > 2 are all in scope now. The unit REFUSES (std::invalid_argument,
// never silently reinterpreting) only:
//   - period or frequency outside {4, 12} (the captured corpus shapes);
//   - any MA/AR input coefficient or knob non-finite; negative d/bd.
// Two cascade paths end in a JVM exception that DefaultModelDecomposer's
// catch(Exception) turns into a null result (it_max exhaustion -> null roots
// -> consumer NPE; a chase of the 2x2 tail block -> ArrayIndexOutOfBounds);
// the unit mirrors both as an internal abort caught at the decompose seam,
// reproducing the captured "non decomposable" decisions (q_airline/q_d2).
// Inside the cascade the one unpinnable JVM branch - the it_count % 25 random
// shift (new Random(); fires on v_simplify's degree-12 stationary AR in both
// the JVM and the unit, whose roots only feed gate-robust partition
// decisions) - is substituted deterministically, not refused: the JVM itself
// is not bit-reproducible there.
// JVM-observed crashes on the live path (the fixMa deg-1 quirk feeding a
// MatrixException into the Legacy approximator, the Noisy-mode NPE) are
// outside this unit and recorded as data in the corpus, never reproduced.

// Pinned DecompositionSpec.ModelApproximationMode. The live builder default
// (and the corpus default) is Legacy; None nulls the toolkit approximator,
// Noisy switches setVarianceMax to the noisy-model branch.
enum class seats_approx_mode { none, legacy, noisy };

// Pinned DecompositionSpec.ComponentsEstimationMethod / BiasCorrection
// (24-C). The live builder defaults (corpus default) are Burman and Legacy;
// the enum declaration orders are {Burman, KalmanSmoother} and {None, Legacy}
// — "New"/"Default" are NOT valid constant names (valueOf throws).
enum class seats_components_method { burman, kalman_smoother };
enum class seats_bias_correction { none, legacy };

// The frozen fixed-parameter SARIMA request (the probe's `spec`/`phi`/...
// lines). Orders ARE the parameter vector lengths (the SarimaSpec builder
// contract). innovation_variance is the SeatsModelSpec.innovationVariance
// echo the decomposition scales the symmetric MA by (convolutionOf(ma, var)).
struct tramo_seats_model_request {
    int frequency = 12;   // series annual frequency; period is the declared
    // SarimaSpec period; the corpus always pairs them (12/12, 4/4).
    int period = 12;
    int d = 0;
    int bd = 0;
    std::vector<double> phi, bphi, theta, btheta;
    // (orders p/bp/q/bq are the vector sizes; p is the phi length.)
    bool log = false;            // echoed in the result; no model-path effect
    bool mean_correction = true;  // echoed (SeatsModel.of carries it)
    double innovation_variance = 0.0;
    // DecompositionSpec knobs (builder defaults: xl 0.95, epsphi 2, rmod 0.5,
    // smod 0.8, smodAtPi 0.8, mode Legacy).
    double xl = 0.95;
    double epsphi = 2.0;
    double rmod = 0.5;
    double smod = 0.8;
    double smod_pi = 0.8;
    seats_approx_mode mode = seats_approx_mode::legacy;
    // ---- 24-C additions (defaults = DecompositionSpec builder defaults) ----
    int forecast_count = -2;  // extrapolationCount sign convention
    int backcast_count = 0;
    seats_components_method components_method = seats_components_method::burman;
    seats_bias_correction bias = seats_bias_correction::legacy;
    // The observation series (the probe's values row, RAW domain: the log
    // flag applies the transform inside the kernel). Only the 24-C
    // extraction unit consumes it; 24-A/24-B ignore it.
    std::vector<double> series;
    int start_year = 1990;  // observation start (TsPeriod month is 1-based
    int start_month = 1;    // here; TsPeriod.annualPosition() is 0-based)
};

// One raw (pre-variance-allocation) RootDecomposer split component, as the
// probe read it through the getters: stationary AR, unit roots and the
// symmetric difference filter (the linear-system solution). The one-sided
// MA and innovation variance of a symmetric model need the 24-B spectral
// factorization and are NOT part of this mirror.
struct seats_raw_component {
    std::vector<double> sar;  // getStationaryAr().asPolynomial()
    std::vector<double> ur;   // getNonStationaryAr().asPolynomial()
    std::vector<double> sfd;  // symmetricMa().coefficientsAsPolynomial()
};

// The unit's decision-and-event record for one request.
struct tramo_seats_model_result {
    // ---- DefaultModelValidator output (the v.* fact block) --------------
    bool parameters_cut_off = false;   // validate() returned false
    double xl = 0.95;                  // validator's xl (getXl())
    // The validated model (m1): orders = parameter lengths + d/bd.
    int orders_p = 0, orders_d = 0, orders_q = 0;
    int orders_bp = 0, orders_bd = 0, orders_bq = 0;
    std::vector<double> phi, bphi, theta, btheta;
    double innovation_variance = 0.0;  // SeatsModel echo (unchanged)
    int period = 12;

    // ---- decomposition decision layer ------------------------------------
    // The raw selector-chain split of the VALIDATED model (r.* facts;
    // r.count 3, structural lines only).
    std::vector<seats_raw_component> raw_components;
    // true when the chain threw (LU MatrixException / singular system):
    // r.decompose "null" + the DefaultModelDecomposer catch -> um null.
    bool decomposition_threw = false;
    // um != null after setVarianceMax(-1, mode == Noisy): the toolkit
    // decomposer's answer (d.decompose "ok"/"null").
    bool variance_max_ok = false;
    // sum of the raw component spectral minima (setVarianceMax's floor);
    // exposed for the sweep report, decisions use its sign only.
    double variance_floor = 0.0;
    // The decomposeModel state machine outcome:
    bool approximation_required = false;  // um null and mode != None
    bool seats_exception = false;         // um null and mode == None
    bool noisy_model_used = false;        // Noisy branch changed the model
    bool model_changed = false;           // SeatsModel.setModelChanged
                                          // (approximation or noisy)
    // The kernel validate+decompose frame as the probe's observer recorded
    // it (push/pop nesting with the remark/warning/error messages the unit
    // itself emits; the estimation-frame events are 23-F/G scope).
    struct event {
        int seq;
        std::string kind;  // push, pop, warning, remark, error
        int depth;
        std::string msg;
        std::string detail;  // "-" or "orders ..." (the approximation remark)
    };
    std::vector<event> events;
};

// Runs the mirror. Throws std::invalid_argument outside the envelope above.
tramo_seats_model_result validate_tramo_seats_model(
    const tramo_seats_model_request& request);

}  // namespace kalos

#endif
