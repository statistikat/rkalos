#include "kalos/state_space.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

constexpr std::size_t max_states = 16;
constexpr double rank_tolerance = 1e-12;
constexpr double residual_tolerance = 1e-10;

void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }

double finite(double value) {
    if (!std::isfinite(value)) {
        throw std::overflow_error("state-space arithmetic produced a non-finite value");
    }
    return value;
}

double finite(long double value) { return finite(static_cast<double>(value)); }

void validate(const StationaryStateSpace& model, const Interrupt& interrupt) {
    model.transition.validate(interrupt);
    model.noise_loading.validate(interrupt);
    const auto m = model.transition.rows();
    if (m > max_states || model.transition.columns() != m ||
        model.observation.size() != m || model.noise_loading.rows() != m ||
        model.noise_loading.columns() > max_states) {
        throw std::invalid_argument("state-space dimensions must agree and have at most 16 states and loading columns");
    }
    for (const auto value : model.observation) {
        if (!std::isfinite(value)) throw std::invalid_argument("observation coefficients must be finite");
    }
    if (!std::isfinite(model.observation_variance) || model.observation_variance < 0) {
        throw std::invalid_argument("observation variance must be finite and nonnegative");
    }
}

DenseMatrix multiply(const DenseMatrix& a, const DenseMatrix& b) {
    if (a.columns() != b.rows()) throw std::invalid_argument("matrix multiplication dimensions do not agree");
    DenseMatrix result(a.rows(), b.columns());
    for (std::size_t j = 0; j < b.columns(); ++j) {
        for (std::size_t i = 0; i < a.rows(); ++i) {
            long double total = 0;
            for (std::size_t k = 0; k < a.columns(); ++k) total += static_cast<long double>(a(i, k)) * b(k, j);
            result(i, j) = finite(total);
        }
    }
    return result;
}

// A B A', with symmetric output. Loop sums use extended precision where the
// platform provides it, but correctness/tolerances do not depend on its width.
DenseMatrix sandwich(const DenseMatrix& a, const DenseMatrix& b) {
    const auto ab = multiply(a, b);
    DenseMatrix result(a.rows(), a.rows());
    for (std::size_t j = 0; j < a.rows(); ++j) {
        for (std::size_t i = 0; i <= j; ++i) {
            long double total = 0;
            for (std::size_t k = 0; k < a.columns(); ++k) total += static_cast<long double>(ab(i, k)) * a(j, k);
            result(i, j) = result(j, i) = finite(total);
        }
    }
    return result;
}

DenseMatrix symmetric_solution(const DenseMatrix& coefficients, std::size_t column,
                                std::size_t m) {
    DenseMatrix result(m, m);
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i <= j; ++i) {
            const long double left = coefficients(i + m * j, column);
            const long double right = coefficients(j + m * i, column);
            result(i, j) = result(j, i) = finite((left + right) / 2);
        }
    }
    return result;
}

#ifdef KALOS_WITH_SLICOT
extern "C" {
// SLICOT SB03MD from the vendored closure (SLICOT-Reference @ 5eeb845d,
// BSD-3-Clause; provenance, license and the closure/symbol argument in
// vendor/slicot/UPSTREAM.md). For DICO='D' it solves
// op(A) X op(A)' - X + SCALE*C = 0, so with TRANA='T' and C=-R the returned X
// satisfies X/SCALE - A'(X/SCALE)A = R; for symmetric R the unique symmetric
// solution of that transposed equation is also the solution of
// P - T P T' = R (verified by the preserved residual check against the
// original equation). gfortran passes hidden string lengths.
void sb03md_(const char* dico, const char* job, const char* fact, const char* transa,
             const int* n, double* a, const int* lda, double* u, const int* ldu,
             double* c, const int* ldc, double* scale, double* sep, double* ferr,
             double* wr, double* wi, int* iwork, double* dwork, const int* ldwork,
             int* info, size_t l_dico, size_t l_job, size_t l_fact, size_t l_transa);
}

struct SlicotSolve {
    DenseMatrix solution{1, 1}; // Symmetrized X / SCALE (the true solution).
    int info = -999;
    double scale = 0.0;
};

// Optimal workspace size for JOB='B' at state order m (LDWORK=-1 query,
// floored at the documented lower bound).
std::size_t slicot_workspace(std::size_t m) {
    const int n = static_cast<int>(m);
    DenseMatrix a(m, m), u(m, m), c(m, m);
    for (std::size_t j = 0; j < m; ++j) u(j, j) = 1;
    std::vector<int> iwork(m * m + 1, 0);
    std::vector<double> wr(m), wi(m), dwork(1, 0.0);
    int ldwork = -1, info = 0;
    double scale = 0, sep = 0, ferr = 0;
    sb03md_("D", "B", "N", "T", &n, a.values().data(), &n, u.values().data(), &n,
            c.values().data(), &n, &scale, &sep, &ferr, wr.data(), wi.data(),
            iwork.data(), dwork.data(), &ldwork, &info, 1, 1, 1, 1);
    if (info != 0) throw std::runtime_error("SB03MD workspace query failed with INFO=" + std::to_string(info));
    const std::size_t lower_bound = 2 * m * m + 2 * m + 16;
    return std::max(lower_bound, static_cast<std::size_t>(std::max(1.0, dwork[0])));
}

// One solve of Y - T Y T' = right_hand_side. `factored` reuses the Schur
// pair (A, U) produced by the previous call (FACT='F'); otherwise the
// transition matrix is factored here (FACT='N') and (A, U) are written out
// for the next call. JOB='B' also computes SEP/FERR; those estimates are
// informational only and deliberately not retained: the certificate measure
// remains the Cholesky reciprocal_condition computed by the caller.
SlicotSolve slicot_solve(const DenseMatrix& t_or_schur, const DenseMatrix& u_in,
                         const DenseMatrix& right_hand_side, bool factored,
                         std::vector<double>& dwork, DenseMatrix* schur_out,
                         DenseMatrix* transformation_out) {
    const std::size_t m = t_or_schur.rows();
    const int n = static_cast<int>(m);
    SlicotSolve out;
    DenseMatrix a(t_or_schur), u(u_in), c(m, m);
    for (std::size_t i = 0; i < m * m; ++i) c.values()[i] = -right_hand_side.values()[i];
    std::vector<int> iwork(m * m + 1, 0);
    std::vector<double> wr(m), wi(m);
    int ldwork = static_cast<int>(dwork.size());
    double sep = 0, ferr = 0;
    sb03md_("D", "B", factored ? "F" : "N", "T", &n, a.values().data(), &n,
            u.values().data(), &n, c.values().data(), &n, &out.scale, &sep, &ferr,
            wr.data(), wi.data(), iwork.data(), dwork.data(), &ldwork, &out.info,
            1, 1, 1, 1);
    // Explicit status mapping; no silent fallback and no flattened status.
    if (out.info == n + 1) {
        // Almost-reciprocal eigenvalue products: the Lyapunov operator is
        // numerically singular. SLICOT continues on perturbed values; that
        // solution is NOT accepted. Maps onto the existing singular/rank-
        // deficient rejection of the QR route (std::invalid_argument), with
        // the SB03MD cause kept visible in the message.
        throw std::invalid_argument(
            "stationary initialization is singular or numerically rank deficient; "
            "diffuse initialization is unsupported (SB03MD INFO=N+1: almost-reciprocal "
            "eigenvalue products; the perturbed-value solution is not accepted)");
    }
    if (out.info != 0) {
        // QR algorithm failure (INFO=1..N) or rejected arguments (INFO<0):
        // hard numerical failure, surfaced, never retried on the QR path.
        throw std::runtime_error("stationary initialization SB03MD failure: INFO=" +
                                 std::to_string(out.info));
    }
    if (!std::isfinite(out.scale) || out.scale <= 0 || out.scale > 1) {
        // The documented overflow guard returns 0 < SCALE <= 1; anything else
        // (including the undocumented SCALE=0) is an explicit failure.
        throw std::runtime_error("stationary initialization SB03MD returned an unusable SCALE=" +
                                 std::to_string(out.scale));
    }
    out.solution = DenseMatrix(m, m);
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i <= j; ++i) {
            const long double value = (static_cast<long double>(c(i, j)) + c(j, i)) / 2 / out.scale;
            out.solution(i, j) = out.solution(j, i) = finite(value);
        }
    }
    if (schur_out) *schur_out = std::move(a);
    if (transformation_out) *transformation_out = std::move(u);
    return out;
}

// Size gate for the SLICOT route (M9 measurement dossier,
// tmp_opencode/m9-perf/): the microbench (SLICOT-scoping-2026-11-25.md)
// measured speedups of 2.9-140x over the Kronecker-QR route at state order
// m >= 8 and 0.8-2x (timer noise) below it; the seam qualification
// (seam/SEAM-QUALIFICATION-2026-11-25.md) measured per-init division by 11-12
// and -26..-36% end-to-end at m = 14 against ~0/+microseconds at m <= 4 (the
// workspace query + opaque Fortran call costs ~1 us/evaluation there). The
// gain window therefore starts at order 8; below it the QR route is retained
// bit-identically.
constexpr std::size_t slicot_minimum_state_order = 8;
#endif // KALOS_WITH_SLICOT

double lyapunov_residual(const DenseMatrix& transition, const DenseMatrix& p,
                         const DenseMatrix& rhs) {
    const auto propagated = sandwich(transition, p);
    long double error = 0;
    long double scale = 0;
    for (std::size_t i = 0; i < p.values().size(); ++i) {
        const long double a = p.values()[i], b = propagated.values()[i], c = rhs.values()[i];
        error = std::max(error, std::abs(a - b - c));
        scale = std::max(scale, std::abs(a) + std::abs(b) + std::abs(c));
    }
    return scale == 0 ? 0 : finite(error / scale);
}

std::vector<double> lag_polynomial(const std::vector<double>& coefficients,
                                    int step, double sign) {
    if (step < 1 || coefficients.size() > max_states ||
        (!coefficients.empty() && static_cast<std::size_t>(step) > max_states / coefficients.size())) {
        throw std::invalid_argument("ARMA polynomial order exceeds the 16-state pilot bound");
    }
    std::vector<double> result(coefficients.size() * static_cast<std::size_t>(step) + 1, 0);
    result[0] = 1;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (!std::isfinite(coefficients[i])) throw std::invalid_argument("ARMA coefficients must be finite");
        result[(i + 1) * static_cast<std::size_t>(step)] = sign * coefficients[i];
    }
    return result;
}

} // namespace

StationaryInitialization stationary_initialization(const StationaryStateSpace& model,
                                                    const Interrupt& interrupt) {
    validate(model, interrupt);
    const auto m = model.transition.rows();
    DenseMatrix process(m, m);
    for (std::size_t j = 0; j < m; ++j) {
        checkpoint(interrupt);
        for (std::size_t i = 0; i <= j; ++i) {
            long double value = 0;
            for (std::size_t k = 0; k < model.noise_loading.columns(); ++k) {
                value += static_cast<long double>(model.noise_loading(i, k)) * model.noise_loading(j, k);
            }
            process(i, j) = process(j, i) = finite(value);
        }
    }
    DenseMatrix covariance(m, m), certificate(m, m), identity(m, m);
#ifdef KALOS_WITH_SLICOT
    if (m >= slicot_minimum_state_order) {
        // Two SB03MD calls, JOB='B': the covariance equation factors T
        // (FACT='N') and the certificate equation reuses the Schur (FACT='F').
        // Interruption checkpoints sit before, between and after the opaque
        // Fortran calls (no checks are possible inside them); a cancellation
        // between the calls discards the Schur pair with the stack frame, so a
        // later call re-factors from the transition matrix (no stale reuse).
        for (std::size_t j = 0; j < m; ++j) identity(j, j) = 1;
        const auto workspace = slicot_workspace(m);
        std::vector<double> dwork(workspace, 0.0);
        DenseMatrix unitary(m, m);
        for (std::size_t j = 0; j < m; ++j) unitary(j, j) = 1;
        checkpoint(interrupt);
        DenseMatrix schur(m, m), schur_transformation(m, m);
        auto covariance_solve = slicot_solve(model.transition, unitary, process, false, dwork,
                                             &schur, &schur_transformation);
        checkpoint(interrupt);
        const auto certificate_solve = slicot_solve(schur, schur_transformation, identity, true,
                                                    dwork, nullptr, nullptr);
        checkpoint(interrupt);
        covariance = std::move(covariance_solve.solution);
        certificate = std::move(certificate_solve.solution);
    } else
#endif // KALOS_WITH_SLICOT
    {
        // The pivoted-QR route: the m^2 x m^2 Kronecker assembly with the
        // two-RHS solve, verbatim in operation order (statements, checkpoints
        // and exception texts as in the pre-seam core). With
        // KALOS_WITH_SLICOT=OFF this is the only route, for every state order.
        const auto dimension = m * m;
        DenseMatrix system(dimension, dimension), rhs(dimension, 2);
        for (std::size_t j = 0; j < m; ++j) {
            checkpoint(interrupt);
            identity(j, j) = 1;
            for (std::size_t i = 0; i < m; ++i) {
                const auto row = i + m * j;
                rhs(row, 0) = process(i, j);
                rhs(row, 1) = i == j ? 1 : 0;
                for (std::size_t l = 0; l < m; ++l) {
                    for (std::size_t k = 0; k < m; ++k) {
                        const auto column = k + m * l;
                        system(row, column) = finite(static_cast<long double>(row == column ? 1 : 0) -
                            static_cast<long double>(model.transition(i, k)) * model.transition(j, l));
                    }
                }
            }
        }
        const auto solved = qr_least_squares(system, rhs, rank_tolerance, interrupt);
        if (solved.rank != static_cast<int>(dimension)) {
            throw std::invalid_argument("stationary initialization is singular or numerically rank deficient; diffuse initialization is unsupported");
        }
        covariance = symmetric_solution(solved.coefficients, 0, m);
        certificate = symmetric_solution(solved.coefficients, 1, m);
    }
    const double error = std::max(lyapunov_residual(model.transition, covariance, process),
                                  lyapunov_residual(model.transition, certificate, identity));
    if (!std::isfinite(error) || error > residual_tolerance) {
        throw std::runtime_error("stationary Lyapunov solve failed its scaled residual check");
    }
    // For real T, existence of positive definite S satisfying S-TST'=I is
    // equivalent to Schur stability. This tests every mode, including modes
    // with zero process noise or zero observation loading.
    CholeskyFactor factor = [&]() {
        try { return cholesky_factor(certificate, interrupt); }
        catch (const std::invalid_argument&) {
            throw std::invalid_argument("transition is not certified stationary: the Lyapunov stability certificate is not positive definite");
        }
    }();
    if (factor.reciprocal_condition <= rank_tolerance) {
        throw std::invalid_argument("stationary stability certificate is numerically ill-conditioned at the pilot tolerance");
    }
    return {std::move(covariance), std::move(process), error, factor.reciprocal_condition};
}

InnovationsLikelihood stationary_innovations(const MaskedSeries& series,
                                              const StationaryStateSpace& model,
                                              double mean, double variance,
                                              const Interrupt& interrupt) {
    series.validate();
    if (!std::isfinite(mean)) throw std::invalid_argument("observation mean must be finite");
    if (!std::isfinite(variance) || variance <= 0) throw std::invalid_argument("innovation scale variance must be finite and positive");
    auto initialization = stationary_initialization(model, interrupt);
    const auto m = model.transition.rows(), n = series.values.size();
    std::vector<double> state(m, 0), predictions(n, 0), variances(n, 0);
    MaskedSeries innovations{std::vector<double>(n, 0), series.missing};
    MaskedSeries standardized{std::vector<double>(n, 0), series.missing};
    DenseMatrix covariance = initialization.covariance;
    std::size_t observed = 0;
    long double squares = 0, log_variances = 0;
    for (std::size_t t = 0; t < n; ++t) {
        if ((t & 63U) == 0) checkpoint(interrupt);
        long double prediction = mean;
        std::vector<double> pz(m, 0);
        for (std::size_t i = 0; i < m; ++i) {
            prediction += static_cast<long double>(model.observation[i]) * state[i];
            long double total = 0;
            for (std::size_t j = 0; j < m; ++j) total += static_cast<long double>(covariance(i, j)) * model.observation[j];
            pz[i] = finite(total);
        }
        predictions[t] = finite(prediction);
        long double f = model.observation_variance;
        for (std::size_t i = 0; i < m; ++i) f += static_cast<long double>(model.observation[i]) * pz[i];
        variances[t] = finite(f);
        if (variances[t] < 0) throw std::runtime_error("negative innovation variance: stationary covariance is numerically unresolved");
        if (!series.missing[t]) {
            if (variances[t] == 0) throw std::invalid_argument("observed data require a positive innovation variance; deterministic observation likelihood is unsupported");
            const double innovation = finite(static_cast<long double>(series.values[t]) - prediction);
            const double standardized_value = finite(innovation / std::sqrt(variances[t]));
            // The guard and the accumulation both observe the binary64 square,
            // so a nonzero standardized innovation with a sub-binary64 exact
            // square is rejected identically on every platform.
            const double squared = static_cast<double>(static_cast<long double>(standardized_value) * standardized_value);
            if ((innovation != 0 && standardized_value == 0) ||
                (standardized_value != 0 && squared == 0)) {
                throw std::underflow_error("standardized innovation square underflows the finite numeric range");
            }
            innovations.values[t] = innovation;
            standardized.values[t] = standardized_value;
            squares += squared;
            log_variances += std::log(static_cast<long double>(variances[t]));
            ++observed;
            std::vector<double> gain(m, 0);
            DenseMatrix update(m, m);
            for (std::size_t i = 0; i < m; ++i) {
                gain[i] = finite(pz[i] / variances[t]);
                state[i] = finite(static_cast<long double>(state[i]) + static_cast<long double>(gain[i]) * innovation);
                for (std::size_t j = 0; j < m; ++j) {
                    update(i, j) = finite(static_cast<long double>(i == j ? 1 : 0) -
                                           static_cast<long double>(gain[i]) * model.observation[j]);
                }
            }
            // Joseph form: (I-KZ)P(I-KZ)' + h K K'.
            auto posterior = sandwich(update, covariance);
            for (std::size_t j = 0; j < m; ++j) {
                for (std::size_t i = 0; i <= j; ++i) {
                    posterior(i, j) = posterior(j, i) = finite(static_cast<long double>(posterior(i, j)) +
                        static_cast<long double>(model.observation_variance) * gain[i] * gain[j]);
                }
            }
            covariance = std::move(posterior);
        }
        if (t + 1 < n) {
            std::vector<double> next(m, 0);
            for (std::size_t i = 0; i < m; ++i) {
                long double value = 0;
                for (std::size_t j = 0; j < m; ++j) value += static_cast<long double>(model.transition(i, j)) * state[j];
                next[i] = finite(value);
            }
            state = std::move(next);
            covariance = sandwich(model.transition, covariance);
            for (std::size_t i = 0; i < covariance.values().size(); ++i) {
                covariance.values()[i] = finite(static_cast<long double>(covariance.values()[i]) +
                                                 initialization.process_covariance.values()[i]);
            }
        }
    }
    const double sum_squares = finite(squares), sum_logs = finite(log_variances);
    if (squares > 0 && sum_squares == 0) {
        throw std::underflow_error("standardized sum of squares underflows the finite numeric range");
    }
    const long double log_two_pi = std::log(2 * std::acos(-1.0L));
    const auto n_observed = static_cast<long double>(observed);
    const double fixed_loglik = finite(-0.5L * (n_observed * (log_two_pi + std::log(static_cast<long double>(variance))) +
                                               log_variances + squares / variance));
    double profile_variance = std::numeric_limits<double>::quiet_NaN();
    double profile_loglik = std::numeric_limits<double>::quiet_NaN();
    std::string status = "no_observations";
    if (observed > 0) {
        profile_variance = finite(squares / n_observed);
        if (squares > 0 && profile_variance == 0) {
            throw std::underflow_error("profiled variance underflows the finite numeric range");
        }
        if (profile_variance == 0) {
            status = "zero_scale_boundary";
            profile_loglik = std::numeric_limits<double>::infinity();
        } else {
            status = "available";
            profile_loglik = finite(-0.5L * (n_observed * (log_two_pi + 1 + std::log(static_cast<long double>(profile_variance))) + log_variances));
        }
    }
    checkpoint(interrupt);
    return {std::move(predictions), std::move(innovations), std::move(standardized),
            std::move(variances), std::move(state), std::move(covariance),
            std::move(initialization), observed, sum_squares, sum_logs, fixed_loglik,
            profile_variance, profile_loglik, std::move(status)};
}

StationaryArma stationary_arma(const std::vector<double>& ar, const std::vector<double>& ma,
                                const std::vector<double>& seasonal_ar,
                                const std::vector<double>& seasonal_ma,
                                int period, const Interrupt& interrupt) {
    if (period < 1) throw std::invalid_argument("seasonal period must be a positive integer");
    auto ar_poly = convolve(lag_polynomial(ar, 1, -1), lag_polynomial(seasonal_ar, period, -1), interrupt);
    auto ma_poly = convolve(lag_polynomial(ma, 1, 1), lag_polynomial(seasonal_ma, period, 1), interrupt);
    const auto m = std::max(ar_poly.size() - 1, ma_poly.size());
    if (m > max_states) throw std::invalid_argument("expanded ARMA model exceeds the 16-state pilot bound");
    std::vector<double> expanded_ar(ar_poly.size() - 1), expanded_ma(ma_poly.begin() + 1, ma_poly.end());
    DenseMatrix transition(m, m), loading(m, 1);
    std::vector<double> observation(m, 0);
    observation[0] = 1;
    for (std::size_t i = 0; i < expanded_ar.size(); ++i) transition(i, 0) = expanded_ar[i] = -ar_poly[i + 1];
    for (std::size_t i = 0; i + 1 < m; ++i) transition(i, i + 1) = 1;
    for (std::size_t i = 0; i < ma_poly.size(); ++i) loading(i, 0) = ma_poly[i];
    return {{std::move(transition), std::move(observation), std::move(loading), 0},
            std::move(expanded_ar), std::move(expanded_ma)};
}

} // namespace kalos
