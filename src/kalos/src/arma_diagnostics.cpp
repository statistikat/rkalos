// Item 21-H: transcription of the jdplus automatic-model diagnostics chain
// behind ModelController.accept: the AutoCovariances acf, the Ljung-Box
// statistic and Chi2 upper-tail p-value (via detail::chi2_upper_p), the
// calcLbLength / calcResStat rules and the accept composition on the shared
// ModifiedLjungBox machine (src/mlj_engine.hpp). Operand orders follow the
// pinned bytecode (docs/contracts/arma-diagnostics-contract.md); the JVM
// exception surfaces map onto built == false / std::invalid_argument per the
// header envelope.
#include "kalos/arma_diagnostics.hpp"

#include "kalos/detail/special_functions.hpp"
#include "mlj_engine.hpp"

#include <cmath>
#include <stdexcept>

namespace kalos {
namespace {

using mlj::Acvf;
using mlj::ModifiedLjungBox;
using mlj::expand_sarima_polynomials;
using mlj::rational_coefficients;

double sum_of(const std::vector<double>& x) {
    double s = 0.0;
    for (const double v : x) {
        s += v;
    }
    return s;
}

double ssq_of(const std::vector<double>& x) {
    double s = 0.0;
    for (const double v : x) {
        s += v * v;
    }
    return s;
}

} // namespace

double auto_correlation(const std::vector<double>& seq, double mean, int lag) {
    const int n = static_cast<int>(seq.size());
    if (lag < 0) {
        throw std::invalid_argument("auto_correlation needs a non-negative lag");
    }
    for (const double v : seq) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("auto_correlation needs finite values");
        }
    }
    if (static_cast<int>(seq.size()) == 0) {
        throw std::invalid_argument("auto_correlation needs a non-empty sequence");
    }
    // varianceNoMissing(seq, mean): all values count, divisor n.
    double var = 0.0;
    for (const double v : seq) {
        const double d = v - mean;
        var += d * d;
    }
    var = var / static_cast<double>(n);
    // Lambda guard: dcmpg(var, 1e-38) - a NaN variance compares greater and
    // takes the compute branch (yielding NaN), as in the JVM.
    if (var < 1e-38) {
        return 0.0;
    }
    // autoCovarianceNoMissing(seq, mean, lag): the pinned divisor is the
    // full length n (biased estimate), with +0.0 for a non-positive
    // n - lag count (out-of-range lag yields 0.0 / var).
    const int count = n - lag;
    if (!(count > 0)) {
        return 0.0;
    }
    double acc = 0.0;
    for (int t = 0; t < count; ++t) {
        acc += (seq[t] - mean) * (seq[t + lag] - mean);
    }
    return (acc / static_cast<double>(n)) / var;
}

int calc_lb_length(int annual_frequency, int neff) {
    int k;
    if (annual_frequency == 12) {
        k = 24;
    } else if (annual_frequency == 1) {
        k = 8;
    } else {
        k = 4 * annual_frequency;
        if (neff <= 22 && annual_frequency == 4) {
            k = 6;
        }
    }
    if (k >= neff) {
        k = neff / 2;
    }
    return k;
}

LjungBoxTest ljung_box_test(const std::vector<double>& seq, int k, int nhp) {
    const int n = static_cast<int>(seq.size());
    if (n == 0) {
        throw std::invalid_argument("ljung_box_test needs a non-empty sequence");
    }
    // calcLb: mean = tmp.sum() / n; tmp.sub(mean) (here on the local copy).
    const double mean = sum_of(seq) / static_cast<double>(n);
    std::vector<double> data(seq.size());
    for (std::size_t i = 0; i < seq.size(); ++i) {
        data[i] = seq[i] - mean;
    }
    // value(): sign filter is the default 0 (add every term); lag = 1.
    double acc = 0.0;
    for (int i = 1; i <= k; ++i) {
        const double r = auto_correlation(data, 0.0, i);
        acc += (r * r) / static_cast<double>(n - i);
    }
    // build(): statistic = value * n * (n + 2); df = (lag == 1) ? k - nhp : k.
    const double stat = acc * static_cast<double>(n) * static_cast<double>(n + 2);
    const double df = static_cast<double>(k - nhp);
    if (!(df > 0.0)) {
        // Chi2 constructor throws IllegalArgumentException; calcLb maps the
        // RuntimeException onto false.
        return LjungBoxTest{false, 0.0, std::nan(""), df};
    }
    const double p = detail::chi2_upper_p(df, stat);
    return LjungBoxTest{true, stat, p, df};
}

double calc_res_stat_rtval(const std::vector<double>& seq) {
    const int n = static_cast<int>(seq.size());
    if (n == 0) {
        throw std::invalid_argument("calc_res_stat_rtval needs a non-empty sequence");
    }
    const double s = sum_of(seq);
    const double ssq = ssq_of(seq);
    const double mean = s / static_cast<double>(n);
    const double var = ssq / static_cast<double>(n) - mean * mean;
    return mean / std::sqrt(var / static_cast<double>(n));
}

ModelControllerResult model_controller_accept(const ModelControllerRequest& request) {
    const int m = static_cast<int>(request.residuals.size());
    if (m == 0) {
        throw std::invalid_argument("accept machine needs differenced residuals");
    }
    for (const double v : request.residuals) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("accept machine needs finite residuals");
        }
    }
    if (request.period < 1 || request.effective_observations < 1
        || request.free_parameters < 0) {
        throw std::invalid_argument("accept machine needs a sane estimation state");
    }
    std::vector<double> ar_full;
    std::vector<double> ma_full;
    expand_sarima_polynomials(request.phi, request.theta, request.bphi, request.btheta,
        request.period, ar_full, ma_full);
    const int q = static_cast<int>(ma_full.size()) - 1;
    const int n = m + q; // filt.prepare(sarima, m) returns n + q
    if (request.effective_observations > n) {
        throw std::invalid_argument("effective observations exceed the ML residual count");
    }
    Acvf acvf(ma_full, ar_full);
    const auto psi = rational_coefficients(ma_full, ar_full, ma_full.size() - 1);
    const ModifiedLjungBox mlj(ar_full, ma_full, acvf, psi, static_cast<std::size_t>(m));
    std::vector<double> out(static_cast<std::size_t>(n), 0.0);
    mlj.apply(request.residuals, out.data());

    ModelControllerResult result{};
    // rvr reads the full ML-residual block before the neff drop.
    result.rvr = std::sqrt(ssq_of(out) / static_cast<double>(request.effective_observations));
    std::vector<double> tail(out.end() - request.effective_observations, out.end());
    result.rtval = calc_res_stat_rtval(tail);
    result.lb_length = calc_lb_length(request.annual_frequency, request.effective_observations);
    const LjungBoxTest lb = ljung_box_test(tail, result.lb_length, request.free_parameters);
    result.lb_built = lb.built;
    result.lb_value = lb.value;
    result.lb_p = lb.p;
    result.ml_residuals = std::move(tail);
    // accept(): if (!calcLb(...)) return false; return !((1 - p) > lbLimit);
    // A NaN p-value leaves the comparison false, i.e. the model is accepted.
    result.accepted = lb.built && !((1.0 - lb.p) > request.lb_limit);
    return result;
}

} // namespace kalos
