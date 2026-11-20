// Item 21-H: transcription of FinalEstimator.test (x13 regarima), the
// post-estimation ARMA reduction cascade (docs/contracts/
// final-estimator-contract.md). Operand order, the dcmpg/dcmpl comparison
// choices and the takeover resets follow the pinned bytecode literally.
#include "kalos/final_estimator.hpp"

#include "kalos/detail/filter_utility.hpp"

#include <cmath>
#include <stdexcept>

namespace kalos {

FinalEstimatorRound final_estimator_test(const SarimaOrders& orders,
    const std::vector<double>& parameters,
    const std::vector<double>& asy_covariance_diagonal, int sample_length,
    double tsig) {
    const int sum = orders.p + orders.bp + orders.q + orders.bq;
    if (orders.p < 0 || orders.bp < 0 || orders.q < 0 || orders.bq < 0) {
        throw std::invalid_argument("negative sarima orders");
    }
    if (static_cast<int>(parameters.size()) != sum) {
        throw std::invalid_argument("parameters must cover the declared orders");
    }
    if (static_cast<int>(asy_covariance_diagonal.size()) < sum) {
        throw std::invalid_argument("the asymptotic covariance diagonal is too short");
    }

    // sig2 plateau on the estimation domain length; sig1 roots limit.
    const double sig2 = sample_length <= 150 ? 0.15 : 0.1;
    const double sig1 = 0.95;
    const double limit = 1.0 / sig1;

    SarimaOrders reduced = orders;
    // FilterUtility.checkRoots per block (SarimaModel.parameters order
    // [phi | bphi | theta | btheta]); the cursor advances unconditionally.
    int cursor = 0;
    bool fP = false;
    if (orders.p > 0) {
        fP = detail::check_roots(&parameters[cursor], static_cast<std::size_t>(orders.p), limit);
    }
    cursor += orders.p;
    bool fBp = false;
    if (orders.bp > 0) {
        fBp = detail::check_roots(&parameters[cursor], static_cast<std::size_t>(orders.bp), limit);
    }
    cursor += orders.bp;
    bool fQ = false;
    if (orders.q > 0) {
        fQ = detail::check_roots(&parameters[cursor], static_cast<std::size_t>(orders.q), limit);
    }
    cursor += orders.q;
    bool fBq = false;
    if (orders.bq > 0) {
        fBq = detail::check_roots(&parameters[cursor], static_cast<std::size_t>(orders.bq), limit);
    }
    if (!fP && !fBp && !fQ && !fBq) {
        return FinalEstimatorRound{0, reduced};
    }

    int cP = 0, cBp = 0, cQ = 0, cBq = 0;
    double ratio = 999.0;
    int idx = -1;
    if (fP) {
        idx += orders.p;
        const double a = std::abs(parameters[idx]);
        const double v = asy_covariance_diagonal[idx];
        if (v > 0.0) {
            const double t = a / std::sqrt(v);
            if (t < tsig && a < sig2) {
                ++cP;
                ratio = t;
            }
        }
    }
    if (fBp) {
        idx += orders.bp;
        const double a = std::abs(parameters[idx]);
        const double v = asy_covariance_diagonal[idx];
        if (v > 0.0) {
            const double t = a / std::sqrt(v);
            if (t < tsig && a < sig2 && ratio > t) {
                ++cBp;
                ratio = t;
                cP = 0;
            }
        }
    }
    if (fQ) {
        idx += orders.q;
        const double a = std::abs(parameters[idx]);
        const double v = asy_covariance_diagonal[idx];
        if (v > 0.0) {
            const double t = a / std::sqrt(v);
            if (t < tsig && a < sig2 && ratio > t) {
                ++cQ;
                ratio = t;
                cP = 0;
                cBp = 0;
            }
        }
    }
    if (fBq) {
        idx += orders.bq;
        const double a = std::abs(parameters[idx]);
        const double v = asy_covariance_diagonal[idx];
        if (v > 0.0) {
            const double t = a / std::sqrt(v);
            if (t < tsig && a < sig2 && ratio > t) {
                ++cBq;
                // The pinned takeover clears the other counts but does not
                // update the ratio.
                cP = 0;
                cBp = 0;
                cQ = 0;
            }
        }
    }

    const int total = cP + cBp + cQ + cBq;
    if (total == 0) {
        return FinalEstimatorRound{0, reduced};
    }
    if (cP > 0) {
        reduced.p -= cP;
    } else if (cBp > 0) {
        reduced.bp -= cBp;
    } else if (cQ > 0) {
        reduced.q -= cQ;
    } else if (cBq > 0) {
        reduced.bq -= cBq;
    }
    return FinalEstimatorRound{total, reduced};
}

} // namespace kalos
