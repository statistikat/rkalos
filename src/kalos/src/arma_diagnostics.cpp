// Item 21-H: transcription of the jdplus automatic-model diagnostics chain
// behind ModelController.accept: the AutoCovariances acf, the Ljung-Box
// statistic and Chi2 upper-tail p-value (via detail::chi2_upper_p), the
// calcLbLength / calcResStat rules and the accept composition on the shared
// ModifiedLjungBox machine (src/mlj_engine.hpp). Operand orders follow the
// pinned bytecode (docs/contracts/arma-diagnostics-contract.md); the JVM
// exception surfaces map onto built == false / std::invalid_argument per the
// header envelope.
#include "kalos/arma_diagnostics.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/special_functions.hpp"
#include "kalos/outlier_scoring.hpp"
#include "mlj_engine.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

struct ResidualMoments {
    double skewness;
    double kurtosis;
};

ResidualMoments residual_moments(const std::vector<double>& residuals) {
    const double n = static_cast<double>(residuals.size());
    const double mean = sum_of(residuals) / n;
    double sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
    for (const double value : residuals) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("residual moments require finite residuals");
        }
        const double deviation = value - mean;
        const double square = deviation * deviation;
        const double cube = square * deviation;
        sum2 += square;
        sum3 += cube;
        sum4 += cube * deviation;
    }
    const double variance = sum2 / n;
    const double skewness = sum3 / (variance * std::sqrt(variance) * n);
    const double kurtosis = sum4 / (variance * variance * n);

    return {skewness, kurtosis};
}

struct ResidualSample {
    double mean;
    double variance;
};

ResidualSample residual_sample(const std::vector<double>& residuals, int start, int length) {
    double sum = 0.0, ssq = 0.0;
    for (int i = start; i < start + length; ++i) {
        const double x = residuals[i];
        sum += x;
        ssq += x * x;
    }
    return {sum / length, ssq / length};
}

double student_two_sided(double t, double df) {
    if (t == 0.0) return 1.0;
    const auto beta = detail::special_inc_beta(df / 2, 0.5, df / (df + t * t));
    if (!beta) return std::numeric_limits<double>::quiet_NaN();
    // Preserve rounding in SpecialFunctions.studentProbability followed by T.getProbability.
    const double delta = (1 - *beta) / 2;
    return 2 * (t < 0 ? 0.5 - delta : 1 - (0.5 + delta));
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

namespace {

LjungBoxTest ljung_box_at_zero(const std::vector<double>& data, int k, int nhp) {
    const int n = static_cast<int>(data.size());
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

} // namespace

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
    return ljung_box_at_zero(data, k, nhp);
}

LjungBoxTest rsa0_residual_ljung_box(const std::vector<double>& residuals) {
    if (residuals.empty() || residuals.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() - 2)) {
        throw std::invalid_argument("RSA0 residual count is outside the supported range");
    }
    if (residuals.size() <= 24) {
        for (double value : residuals) {
            if (!std::isfinite(value)) throw std::invalid_argument("RSA0 residuals must be finite");
        }
        return {true, std::nan(""), std::nan(""), 22.0};
    }
    return ljung_box_at_zero(residuals, 24, 2);
}

std::optional<ResidualNormalityTest> rsa0_residual_normality(const std::vector<double>& residuals) {
    if (residuals.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("rsa0_residual_normality requires a JVM-sized series");
    }
    if (residuals.size() <= 7) return std::nullopt;
    const double n = static_cast<double>(residuals.size());
    const auto [skewness, kurtosis] = residual_moments(residuals);

    // DoornikHansen.build() transforms sample skewness to a normal score.
    const double beta = 3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3)
                      / ((n - 2) * (n + 5) * (n + 7) * (n + 9));
    const double w2 = -1 + std::sqrt(2 * (beta - 1));
    const double delta = 1 / std::sqrt(0.5 * detail::fdlibm_log(w2));
    const double y = skewness * std::sqrt((w2 - 1) * (n + 1) * (n + 3) / (12 * (n - 2)));
    const double z1 = delta * detail::fdlibm_log(y + std::sqrt(y * y + 1));

    // The kurtosis transform is conditional on observed squared skewness.
    const double den = (n - 3) * (n + 1) * (n * n + 15 * n - 4);
    const double a = (n - 2) * (n + 5) * (n + 7) * (n * n + 27 * n - 70) / (den * 6);
    const double c = (n - 7) * (n + 5) * (n + 7) * (n * n + 2 * n - 5) / (den * 6);
    const double k = (n + 5) * (n + 7) * (n * n * n + n * n * 37 + 11 * n - 313) / (den * 12);
    const double alpha = a + c * skewness * skewness;
    const double chi = 2 * k * (kurtosis - 1 - skewness * skewness);
    const double z2 = std::sqrt(alpha * 9)
        * (std::pow(chi / (2 * alpha), 1.0 / 3.0) - 1 + 1 / (9 * alpha));
    const double value = z1 * z1 + z2 * z2;
    return ResidualNormalityTest{value, detail::chi2_upper_p(ResidualNormalityTest::df, value)};
}

TramoComparatorResidualTests tramo_comparator_residual_tests(const std::vector<double>& residuals) {
    if (residuals.size() < 4 || residuals.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("monthly TRAMO sample comparisons require at least four residuals");
    const int n = static_cast<int>(residuals.size());
    const double skew = residual_moments(residuals).skewness;
    const double lower = standard_normal_lower_cdf(skew / std::sqrt(6.0 / n));
    const double skew_p = 2 * (skew < 0 ? lower : 1 - lower);

    // ModelStatistics.of uses [0, n-drop) and [drop, n). These samples can
    // overlap on long series. Sample.ofResiduals measures variance around zero.
    const int drop = std::min(n / 2 + n % 2, 120);
    const int length = n - drop;
    const auto first = residual_sample(residuals, 0, length);
    const auto second = residual_sample(residuals, drop, length);
    const double df = length - 1;
    const double pooled = (first.variance * df + second.variance * df) / (2 * df);
    const double t = (second.mean - first.mean) / std::sqrt(pooled / length + pooled / length);
    const double ratio = second.variance / first.variance;
    return {std::abs(skew), skew_p, t, student_two_sided(t, 2 * df),
            ratio, detail::f_upper_p(ratio, df, df)};
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
