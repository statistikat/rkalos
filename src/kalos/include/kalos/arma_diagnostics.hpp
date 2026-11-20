#ifndef KALOS_ARMA_DIAGNOSTICS_HPP
#define KALOS_ARMA_DIAGNOSTICS_HPP

// Item 21-H decision units for the automatic-model diagnostics: the pinned
// AutoCovariances.autoCorrelationFunction, the Ljung-Box test machine that
// ModelController.calcLb builds (LjungBox(seq, false) + TestsUtility.testOf
// with the Chi2 upper tail), the calcLbLength / calcResStat rules and the
// accept-model composition (modified Ljung-Box residuals -> rvr/rtval ->
// Ljung-Box acceptance). The bytecode pin lives in
// docs/contracts/arma-diagnostics-contract.md; the ML-residual filter is the
// shared src/mlj_engine.hpp transcription qualified by the 21-C corpus.
//
// Scope envelope: sequences must be non-empty and finite (the pinned data
// flow never carries NaN here; the non-finite branch of
// autoCorrelationFunction is not transcribed), lags within [0, n), the
// Ljung-Box sign filter stays at its default 0, and the ARMA blocks follow
// the kalos conventions (phi/theta/bphi/btheta as estimated, period >= 1).
// Violations raise std::invalid_argument before any JVM-visible state.

#include <cstddef>
#include <vector>

namespace kalos {

// AutoCovariances.autoCorrelationFunction(seq, mean).applyAsDouble(lag):
// varianceNoMissing(seq, mean) under 1e-38 (dcmpg) yields 0, otherwise
// autoCovarianceNoMissing(seq, mean, lag) / variance with the biased
// n - lag normalization.
double auto_correlation(const std::vector<double>& seq, double mean, int lag);

// ModelController.calcLbLength(annualFrequency, neff): 12 -> 24, 1 -> 8,
// else 4 * frequency with the (neff <= 22, frequency 4) -> 6 override, then
// the k >= neff halving (integer division).
int calc_lb_length(int annual_frequency, int neff);

// LjungBox(seq, false).autoCorrelationsCount(k).hyperParametersCount(nhp)
// .build() + TestsUtility.testOf(statistic, Chi2(df), Upper). The machine
// mean-corrects a copy of seq (the pinned calcLb runs tmp.sub(mean) on the
// caller's block; callers decide whether to keep the mutation). built ==
// false mirrors the Chi2 constructor rejection (df = k - nhp not > 0), which
// ModelController.calcLb catches as RuntimeException and maps onto false;
// p is NaN on the TestsUtility.pvalue catch surface.
struct LjungBoxTest {
    bool built;
    double value; // statistic = value() * n * (n + 2)
    double p;
    double df;
};

LjungBoxTest ljung_box_test(const std::vector<double>& seq, int k, int nhp);

// ModelController.calcResStat side effect: mean / sqrt(var / n) with
// mean = sum / n and var = ssq / n - mean * mean (var / n under the root).
double calc_res_stat_rtval(const std::vector<double>& seq);

// ModelController.accept(RegSarimaModelling) machine, given the already
// differenced linear-model residuals (LinearModel.calcResiduals(allCoefficients))
// and the estimation state it reads (effective observation count, annual
// frequency, free ARMA parameter count, ljungBoxLimit).
struct ModelControllerRequest {
    std::vector<double> residuals; // differenced model, length m >= 1
    std::vector<double> phi;
    std::vector<double> theta;
    std::vector<double> bphi;
    std::vector<double> btheta;
    int period = 1;
    int effective_observations = 0; // statistics().getEffectiveObservationsCount()
    int annual_frequency = 0;       // description.getAnnualFrequency()
    int free_parameters = 0;        // arimaSpec.freeParametersCount()
    double lb_limit = 0.95;         // AutoModelSpec.ljungBoxLimit
};

struct ModelControllerResult {
    bool accepted;         // accept() return value
    bool lb_built;         // calcLb survived the build
    int lb_length;         // calcLbLength result actually used
    double lb_value;       // statistic (0 when !lb_built)
    double lb_p;           // p-value (NaN when !lb_built)
    double rvr;            // sqrt(ssq(ml residuals) / neff), pre-drop
    double rtval;          // calcResStat on the kept tail
    std::vector<double> ml_residuals; // kept (last neff) ML residuals
};

ModelControllerResult model_controller_accept(const ModelControllerRequest& request);

} // namespace kalos

#endif
