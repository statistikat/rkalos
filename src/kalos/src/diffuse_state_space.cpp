#include "kalos/diffuse_state_space.hpp"
#include "kalos/exact_history.hpp"
#include "kalos/sparse_history.hpp"
#include "kalos/finite_state_factor.hpp"

// DLARFG is called through R's declared LAPACK ABI. Only its orthogonal
// reflector is used; no third-party filtering implementation is incorporated.
#include "kalos/detail/lapack.hpp"
#include "kalos/detail/hostlineage.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

constexpr std::size_t max_states = 32;
constexpr double rank_tolerance = 1e-12;
constexpr double projection_tolerance = 512 * std::numeric_limits<double>::epsilon();

void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }

// long double is accumulation scratch only; the numeric contract of every
// guarded quantity is its correctly rounded binary64 value, so a nonzero
// extended-precision intermediate that rounds to zero at binary64 granularity
// is zero on every platform. Underflow contracts are enforced separately on
// named binary64 quantities; see
// docs/contracts/extended-precision-range-contract.md.
double checked(long double value) {
    const auto result = static_cast<double>(value);
    if (!std::isfinite(result)) throw std::overflow_error("diffuse state-space arithmetic produced a non-finite value");
    return result;
}

double norm(const std::vector<double>& values) {
    double result = 0;
    for (const auto value : values) result = checked(kalos::detail::glibc_hypot(result, value));
    return result;
}

double square(double value) {
    // The guard observes the binary64 product, so a nonzero binary64 value
    // with a sub-binary64 exact square is rejected identically on every
    // platform regardless of the scratch precision used for the product.
    const double product = static_cast<double>(static_cast<long double>(value) * value);
    if (value != 0 && product == 0) throw std::underflow_error("diffuse state-space square underflows the finite numeric range");
    return product;
}

struct Basis {
    std::size_t rows;
    std::size_t columns;
    std::vector<double> values;
    double operator()(std::size_t i, std::size_t j) const { return values.at(i + rows * j); }
};

void validate_basis(const Basis& basis, const Interrupt& interrupt, const char* context) {
    if (basis.columns > basis.rows || basis.values.size() != basis.rows * basis.columns) {
        throw std::invalid_argument("diffuse loading dimensions do not match its stated rank");
    }
    for (const auto value : basis.values) {
        if (!std::isfinite(value)) throw std::invalid_argument("diffuse loading values must be finite");
    }
    if (basis.columns == 0) return;
    const DenseMatrix loading(basis.rows, basis.columns, basis.values);
    const DenseMatrix zero_rhs(basis.rows, 1);
    const auto solved = qr_least_squares(loading, zero_rhs, rank_tolerance, interrupt);
    if (solved.rank != static_cast<int>(basis.columns)) throw std::invalid_argument(context);
}

DenseMatrix gram(const Basis& loading) {
    DenseMatrix result(loading.rows, loading.rows);
    for (std::size_t i = 0; i < loading.rows; ++i) {
        double row_norm = 0;
        for (std::size_t k = 0; k < loading.columns; ++k) row_norm = checked(kalos::detail::glibc_hypot(row_norm, loading(i, k)));
        result(i, i) = square(row_norm);
        for (std::size_t j = i + 1; j < loading.rows; ++j) {
            long double total = 0;
            for (std::size_t k = 0; k < loading.columns; ++k) total += static_cast<long double>(loading(i, k)) * loading(j, k);
            result(i, j) = result(j, i) = checked(total);
        }
    }
    return result;
}

DenseMatrix gram(const DenseMatrix& loading) {
    return gram(Basis{loading.rows(), loading.columns(), loading.values()});
}

DenseMatrix sandwich(const DenseMatrix& a, const DenseMatrix& p) {
    DenseMatrix ap(a.rows(), p.columns()), result(a.rows(), a.rows());
    for (std::size_t j = 0; j < p.columns(); ++j) {
        for (std::size_t i = 0; i < a.rows(); ++i) {
            long double value = 0;
            for (std::size_t k = 0; k < a.columns(); ++k) value += static_cast<long double>(a(i, k)) * p(k, j);
            ap(i, j) = checked(value);
        }
    }
    for (std::size_t j = 0; j < a.rows(); ++j) {
        for (std::size_t i = 0; i <= j; ++i) {
            long double value = 0;
            for (std::size_t k = 0; k < a.columns(); ++k) value += static_cast<long double>(ap(i, k)) * a(j, k);
            result(i, j) = result(j, i) = checked(value);
        }
    }
    return result;
}

std::vector<double> product(const DenseMatrix& a, const std::vector<double>& x) {
    std::vector<double> result(a.rows(), 0);
    for (std::size_t i = 0; i < a.rows(); ++i) {
        long double value = 0;
        for (std::size_t j = 0; j < a.columns(); ++j) value += static_cast<long double>(a(i, j)) * x[j];
        result[i] = checked(value);
    }
    return result;
}

std::vector<double> observation_loading(const std::vector<double>& z, const Basis& basis) {
    std::vector<double> result(basis.columns, 0);
    for (std::size_t j = 0; j < basis.columns; ++j) {
        long double value = 0;
        for (std::size_t i = 0; i < basis.rows; ++i) value += static_cast<long double>(z[i]) * basis(i, j);
        result[j] = checked(value);
    }
    return result;
}

Basis propagate(const DenseMatrix& transition, const Basis& basis,
                 const std::vector<double>& z, bool observation_processed) {
    Basis result{basis.rows, basis.columns, std::vector<double>(basis.values.size(), 0)};
    for (std::size_t j = 0; j < basis.columns; ++j) {
        for (std::size_t i = 0; i < basis.rows; ++i) {
            // After either observed update, Z B_post is mathematically zero.
            // An exactly identical transition row therefore has an exactly
            // zero propagated loading. This structural identity avoids
            // cancellation when SARIMA stores the newly observed history;
            // no approximate row match or small-value threshold is used.
            bool observed_history_row = observation_processed;
            for (std::size_t k = 0; observed_history_row && k < basis.rows; ++k) {
                observed_history_row = transition(i, k) == z[k];
            }
            if (observed_history_row) continue;
            long double value = 0;
            for (std::size_t k = 0; k < basis.rows; ++k) value += static_cast<long double>(transition(i, k)) * basis(k, j);
            result.values[i + basis.rows * j] = checked(value);
        }
    }
    return result;
}

// A canonical history frame is normalized in the ORIGINAL initial-history
// coordinates. Only its orientation changes after observations; this is not
// a renormalization of the propagated state loading or the diffuse measure.
Basis history_loading(const detail::SparseHistoryBasis& history, std::size_t states,
                       std::uint64_t time, const Interrupt& interrupt) {
    const auto k = history.dimension(), rank = history.remaining_rank();
    if (k > states) throw std::logic_error("history frame exceeds the state dimension");
    const auto values = history.loading(time, interrupt);
    if (values.size() != k * rank) throw std::logic_error("history frame has inconsistent dimensions");
    Basis result{states, rank, std::vector<double>(states * rank, 0)};
    for (std::size_t j = 0; j < rank; ++j) {
        for (std::size_t i = 0; i < k; ++i) result.values[states - k + i + states * j] = values[i + k * j];
    }
    return result;
}

// w H = [beta, 0,...,0] for an orthogonal LAPACK reflector H. Retaining the
// remaining columns of B H removes one identified coordinate without a
// subtractive P_inf update or an arbitrary covariance threshold.
Basis remove_observed_direction(const Basis& basis, std::vector<double> w,
                                 const std::vector<double>& z, double w_norm,
                                 double& maximum_residual, const Interrupt& interrupt) {
    if (basis.columns == 1) return Basis{basis.rows, 0, {}};
    for (auto& value : w) value /= w_norm;
    const int n = static_cast<int>(w.size()), increment = 1;
    double tau = 0;
    checkpoint(interrupt);
    lapack::dlarfg(&n, &w[0], w.data() + 1, &increment, &tau);
    checkpoint(interrupt);
    if (!std::isfinite(tau) || tau < 0 || tau > 2) throw std::runtime_error("LAPACK returned an invalid diffuse reflector");
    w[0] = 1;
    Basis result{basis.rows, basis.columns - 1,
                 std::vector<double>(basis.rows * (basis.columns - 1), 0)};
    for (std::size_t i = 0; i < basis.rows; ++i) {
        long double dot = 0;
        for (std::size_t k = 0; k < basis.columns; ++k) dot += static_cast<long double>(basis(i, k)) * w[k];
        for (std::size_t j = 0; j < result.columns; ++j) {
            result.values[i + result.rows * j] = checked(static_cast<long double>(basis(i, j + 1)) - tau * dot * w[j + 1]);
        }
    }
    validate_basis(result, interrupt, "diffuse projection is numerically rank deficient");
    const auto residual = norm(observation_loading(z, result));
    const long double scale = static_cast<long double>(norm(z)) * norm(result.values);
    const double relative = scale == 0 ? 0 : checked(residual / scale);
    maximum_residual = std::max(maximum_residual, relative);
    if (relative > projection_tolerance * static_cast<double>(basis.rows)) {
        throw std::runtime_error("diffuse orthogonal projection failed its residual check");
    }
    return result;
}

DenseMatrix joseph(const DenseMatrix& p, const std::vector<double>& gain,
                    const std::vector<double>& z, double h) {
    const auto m = p.rows();
    DenseMatrix update(m, m);
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i < m; ++i) {
            update(i, j) = checked(static_cast<long double>(i == j ? 1 : 0) - static_cast<long double>(gain[i]) * z[j]);
        }
    }
    auto posterior = sandwich(update, p);
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i <= j; ++i) {
            posterior(i, j) = posterior(j, i) = checked(static_cast<long double>(posterior(i, j)) +
                static_cast<long double>(h) * gain[i] * gain[j]);
        }
    }
    return posterior;
}

// ---- Item-30 option (ii): structure-aware dense k = 0 steps ----
// (docs/item30-k0-routing-probe-2026-11-28.md §6-§7; maintainer adoption
// ruling of 2026-11-28: "yes go with option (ii)"). Scope is the k = 0
// dense regime ONLY: the pointer-local accessor idiom and the exact-zero
// structural skips must not widen beyond the k = 0 gate below, and every
// k > 0 or generic-diffuse shape keeps the dense/factorized code above
// bit-for-bit. The witnesses are shape-only and evaluated once per call;
// a shape failure keeps today's dense path. The k > 0 canonical route is
// structurally excluded (its FiniteStateFactor is armed, and the structural
// branches sit behind the factor in every branch chain).
//
// Exactness argument (the same one the shipped D4 case E1 uses): under the
// witnesses the skipped products are EXACT +-0.0 terms in accumulators that
// start at +0; +0 + (+-0) = +0 and a nonzero accumulator is unchanged by a
// +-0 add, so every stored long double accumulator — and hence every named
// binary64 value — equals the generic sandwich/joseph/product result bit for
// bit, in the same ascending-k term order. The claim was additionally tested
// byte-for-byte on both frozen k = 0 shapes (probe corpus, OPTION-II-BIT-
// IDENTICAL) and is re-pinned by tests/test_diffuse_k0_goldens.cpp.
bool witness_e1_loading(const std::vector<double>& z) {
    for (std::size_t i = 0; i < z.size(); ++i)
        if (z[i] != (i == 0 ? 1.0 : 0.0)) return false;
    return true;
}

bool witness_companion_columns(const DenseMatrix& t) {
    // A(i,0) arbitrary; A(i,i+1) exactly 1.0 for i < m-1; every other entry
    // exactly 0.0 (the stationary_arma companion: coefficient column 0 and a
    // superdiagonal shift).
    const auto m = t.rows();
    if (t.columns() != m) return false;
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 1; j < m; ++j)
            if (t(i, j) != (j == i + 1 ? 1.0 : 0.0)) return false;
    return true;
}

struct StructuralScratch {
    DenseMatrix ap{1, 1};         // persistent (TP / UP) product buffer
    std::vector<double> u_col0;   // (delta_i0 - gain_i) column of U = I - K Z
};

// pz = P z and f = h + z' P z for z = e_1 (witness-checked by the caller):
// pz[i] = checked(0 + P(i,0)*1 + sum_{k>0} P(i,k)*0), the generic product's
// accumulator reduced to its only nonzero term, and f starts from h and adds
// z[0]*pz[0] = pz[0] before its +-0 tail. pz is resized, never reallocated
// mid-run; the structural path performs no per-step allocation.
void struct_observation(const DenseMatrix& p, double h, std::vector<double>& pz, long double& f) {
    const auto m = p.rows();
    pz.resize(m);
    const double* pp = p.values().data();
    for (std::size_t i = 0; i < m; ++i) {
        // Accumulator semantics match the generic product term for term,
        // including the sign of zero: +0 + (-0) = +0, so a -0.0 entry of the
        // first covariance column reads as the +0.0 the generic accumulation
        // stores (a raw read would not; the skeleton's raw read is hardened
        // here, a value-identical deviation on every captured shape).
        long double value = 0;
        value += static_cast<long double>(pp[i]);
        pz[i] = checked(value);
    }
    f = static_cast<long double>(h);
    f += static_cast<long double>(pz[0]);
}

// joseph with z = e_1: U = I - K e_1' has column 0 = (1-K_0, -K_1, ...) and
// identity columns elsewhere; the posterior is written into p in place after
// the UP product completed into scratch (the generic path's return-and-assign
// has the same read/write history). Raw values() views replace the
// out-of-line bounds-checked DenseMatrix::operator() — pointer-local
// accumulation is a k = 0 idiom per the 2026-11-28 ruling and the dense
// k > 0 fallback above is deliberately left on the accessor path.
void struct_joseph(DenseMatrix& p, const std::vector<double>& gain,
                   double h, StructuralScratch& scratch) {
    const auto m = p.rows();
    scratch.u_col0.resize(m);
    for (std::size_t i = 0; i < m; ++i)
        scratch.u_col0[i] = checked((i == 0 ? 1.0L : 0.0L) - static_cast<long double>(gain[i]));
    DenseMatrix& ap = scratch.ap;  // persistent buffer: no per-step allocation
    if (ap.rows() != m || ap.columns() != m) ap = DenseMatrix(m, m);
    double* pp = p.values().data();
    double* pap = ap.values().data();
    // ap = U * P: k = 0 term (u_col0[i] * P(0,j)) then, for i > 0, the k = i
    // term +1.0 * P(i,j); the generic loop accumulates k in ascending order
    // from a +0 accumulator — the folds below keep that start, so a product
    // that is exactly +-0 folds to +0 exactly as the generic loop's does.
    for (std::size_t j = 0; j < m; ++j) {
        const double p0j = pp[m * j];
        for (std::size_t i = 0; i < m; ++i) {
            long double value = 0;
            value += static_cast<long double>(scratch.u_col0[i]) * p0j;
            if (i > 0) value += static_cast<long double>(pp[i + m * j]);
            pap[i + m * j] = checked(value);
        }
    }
    // posterior(i,j) = sum_k ap(i,k) U(j,k): j = 0 keeps the k = 0 term only;
    // j > 0 accumulates k = 0 (ap(i,0) * u_col0[j]) then k = j (+ap(i,j)).
    for (std::size_t j = 0; j < m; ++j) {
        const double u0j = scratch.u_col0[j];
        for (std::size_t i = 0; i <= j; ++i) {
            long double value = 0;
            value += static_cast<long double>(pap[i]) * u0j;
            if (j > 0) value += static_cast<long double>(pap[i + m * j]);
            const double v = checked(value);
            pp[i + m * j] = v;
            pp[j + m * i] = v;
        }
    }
    for (std::size_t j = 0; j < m; ++j) {
        const double gj = gain[j];
        for (std::size_t i = 0; i <= j; ++i) {
            const double v = checked(static_cast<long double>(pp[i + m * j]) +
                static_cast<long double>(h) * gain[i] * gj);
            pp[i + m * j] = v;
            pp[j + m * i] = v;
        }
    }
}

// sandwich(t, p) + process for the witness_companion_columns transition:
// (T P)(i,j) accumulates k = 0 (phi_i * P(0,j)) then k = i+1 (+P(i+1,j));
// (T P T')(i,j) accumulates k = 0 (ap(i,0) * phi_j) then k = j+1 (+ap(i,j+1)).
void struct_predict_add(DenseMatrix& p, const DenseMatrix& t, const DenseMatrix& process,
                        StructuralScratch& scratch) {
    const auto m = p.rows();
    DenseMatrix& ap = scratch.ap;
    if (ap.rows() != m || ap.columns() != m) ap = DenseMatrix(m, m);
    double* pp = p.values().data();
    double* pap = ap.values().data();
    const double* pt = t.values().data();
    // The folds below start from the same +0 accumulator as the generic
    // sandwich loops (see the note in struct_joseph), so stored covariance
    // entries carry the same sign of zero as the generic path's do.
    for (std::size_t j = 0; j < m; ++j) {
        const double p0j = pp[m * j];
        for (std::size_t i = 0; i < m; ++i) {
            long double value = 0;
            value += static_cast<long double>(pt[i]) * p0j;
            if (i + 1 < m) value += static_cast<long double>(pp[i + 1 + m * j]);
            pap[i + m * j] = checked(value);
        }
    }
    for (std::size_t j = 0; j < m; ++j) {
        const double phi_j = pt[j];
        const double* apj1 = j + 1 < m ? pap + m * (j + 1) : nullptr;
        for (std::size_t i = 0; i <= j; ++i) {
            long double value = 0;
            value += static_cast<long double>(pap[i]) * phi_j;
            if (apj1) value += static_cast<long double>(apj1[i]);
            const double v = checked(value);
            pp[i + m * j] = v;
            pp[j + m * i] = v;
        }
    }
    const double* pr = process.values().data();
    for (std::size_t i = 0; i < p.values().size(); ++i)
        pp[i] = checked(static_cast<long double>(pp[i]) + pr[i]);
}

void validate_model(const DiffuseStateSpace& model, const std::vector<double>& state,
                     double mean, double variance, const Interrupt& interrupt) {
    model.transition.validate(interrupt);
    model.noise_loading.validate(interrupt);
    const auto m = model.transition.rows();
    if (m > max_states || model.transition.columns() != m || model.observation.size() != m ||
        model.noise_loading.rows() != m || model.noise_loading.columns() > max_states || state.size() != m) {
        throw std::invalid_argument("diffuse dimensions must agree and have at most 32 states and loading columns");
    }
    for (const auto value : state) if (!std::isfinite(value)) throw std::invalid_argument("initial state must be finite");
    for (const auto value : model.observation) if (!std::isfinite(value)) throw std::invalid_argument("observation coefficients must be finite");
    if (!std::isfinite(model.observation_variance) || model.observation_variance < 0) throw std::invalid_argument("observation variance must be finite and nonnegative");
    if (!std::isfinite(mean)) throw std::invalid_argument("mean must be finite");
    if (!std::isfinite(variance) || variance <= 0) throw std::invalid_argument("variance must be finite and positive");
}

DiffuseLikelihood filter(const MaskedSeries& series, const DiffuseStateSpace& model,
                          std::vector<double> state, DenseMatrix p, Basis basis,
                          double mean, double variance, const Interrupt& interrupt,
                          detail::CountsHistoryCertificate* history = nullptr,
                          detail::SparseHistoryBasis* stable_history = nullptr,
                          detail::FiniteStateFactor* finite_factor = nullptr,
                          bool save_gains = false, bool k0_canonical = false) {
    series.validate();
    validate_model(model, state, mean, variance, interrupt);
    validate_basis(basis, interrupt, "initial diffuse loading is numerically rank deficient");
    const auto m = model.transition.rows(), n = series.values.size(), initial_rank = basis.columns;
    if (save_gains && n > static_cast<std::size_t>(std::numeric_limits<int>::max()) / m) {
        throw std::length_error("regression gain history exceeds the supported matrix range");
    }
    std::vector<double> observation_gains(save_gains ? m * n : 0, 0);
    if (history) {
        if (history->dimension() != initial_rank || history->rank() != 0) {
            throw std::logic_error("exact history certificate disagrees with initial diffuse rank");
        }
        static_assert(std::numeric_limits<std::size_t>::digits <=
                      std::numeric_limits<std::uint64_t>::digits,
                      "exact history certificate requires size_t to fit uint64");
        history->validate_length(static_cast<std::uint64_t>(n));
    }
    if (stable_history && (!history || stable_history->dimension() != initial_rank ||
                           stable_history->rank() != 0)) {
        throw std::logic_error("stable history frame disagrees with the initial diffuse certificate");
    }
    if (finite_factor && (!stable_history || finite_factor->dimension() != m)) {
        throw std::logic_error("finite factor disagrees with the canonical state dimension");
    }
    // Item-29 design D4 (design dossier docs/item29-diffuse-filter-probe-
    // 2026-11-26.md §4 E2, §6 D4): after the diffuse pivots are exhausted, a
    // canonical zero-observation-variance SARIMA filter carries the k history
    // coordinates of its finite covariance as a certified-zero block; the
    // post-pivot recursion may run on the r-wide stationary factor while the
    // state keeps the observed-history ring (the m-wide state bookkeeping
    // below is exact and cheap). Eligibility is the canonical SARIMA route —
    // exact-history certificate, stable history propagation, square-root
    // finite factor, a nonzero history block, and zero observation variance
    // (the P^+Z = 0 nullity is a zero-variance identity; generic h > 0 models
    // and the generic diffuse_innovations entry keep the full path, route 0).
    // Route states: -1 ineligible, 0 awaiting pivot exhaustion, 1 nullity
    // certified (armed: the pivot-exhaust predict still runs full-m, then the
    // rounding-level history block is certified and the factor compressed),
    // 2 reduced r-wide recursion, 3 certificate rejected (full path kept).
    // Nothing ever switches on a weak direction and nothing falls back
    // silently: the route, its step count and its residual land in
    // DiffuseNumerics.
    const bool d4_eligible = history != nullptr && stable_history != nullptr &&
                             finite_factor != nullptr && initial_rank > 0 &&
                             model.observation_variance == 0;
    int d4_state = d4_eligible ? 0 : -1;
    std::size_t d4_certified_steps = 0, d4_route = 0;
    double d4_certificate_residual = 0;
    const std::size_t d4_r = m - initial_rank;
    std::vector<double> d4_observation;                        // stationary Z block (r)
    std::unique_ptr<DenseMatrix> d4_transition, d4_noise;      // r-blocks, built at the switch
    // Trim stages of the qualification dossier (tmp_opencode/item29-d4b/
    // D4B-TRIM-QUALIFICATION-2026-11-26.md, exactness cases E1/E3/E4):
    // d4_pz keeps the r-wide projected loading directly (case E4),
    // d4_state_scratch backs the block-structured state advance (case E1),
    // and d4_block_state records the fail-closed companion-shape witness
    // verified once at the switch — the E1 fast path is a bit-identical
    // statement of product(model.transition, state) only under that shape.
    // The measured-zero stages 3/4 (witness silence in the window, r-long
    // state adds) are NOT shipped: the ring/frame/rank witnesses run at every
    // step and the state-add loop keeps its m-long form.
    std::vector<double> d4_pz, d4_state_scratch, d4_gain_r;
    bool d4_block_state = false;
    // Item-30 option (ii) gate (structure-aware dense k = 0 steps; helper
    // comments above carry the scope note and the exactness argument). The
    // witness is shape-only and evaluated once per call, fail-closed: any
    // failure — non-canonical route, an armed factor (every canonical k > 0
    // call arms one), a nonzero initial diffuse rank, nonzero observation
    // variance, or a shape departure from the e_1 loading / companion
    // transition — keeps today's dense fallback bit-for-bit. Under the armed
    // witness the diffuse pivot branch is unreachable (a rank-zero basis
    // keeps w empty, so w_norm == 0 at every step) and the post-pivot route
    // stays 0, exactly as today at k = 0.
    StructuralScratch k0_scratch;
    const bool k0_structural = k0_canonical && finite_factor == nullptr &&
                               initial_rank == 0 && model.observation_variance == 0 &&
                               witness_e1_loading(model.observation) &&
                               witness_companion_columns(model.transition);
    // The projected-loading buffer lives outside the loop: the structural
    // observation branch resizes it (capacity retained, no mid-run
    // allocation) and every other branch assigns it in full before reading,
    // exactly as the per-step local did (one buffer acquisition per
    // evaluation for the factor/dense branches, none for the structural
    // branch after the first step).
    std::vector<double> pz;
    const auto process = gram(model.noise_loading);
    const auto initial_finite = p;
    const auto initial_diffuse = gram(basis);
    std::vector<double> predictions(n, 0), finite_variances(n, 0), diffuse_variances(n, 0);
    std::vector<int> update_kind(n, 0);
    std::vector<std::size_t> rank_before(n, 0), rank_after(n, 0);
    MaskedSeries innovations{std::vector<double>(n, 0), series.missing};
    MaskedSeries standardized{std::vector<double>(n, 0), std::vector<std::uint8_t>(n, 1)};
    std::size_t observed = 0, diffuse_count = 0, finite_count = 0, last_diffuse = 0, rank_checks = 0;
    std::size_t exact_zero_count = 0, exact_pivot_count = 0, history_reconstructions = 0;
    long double squares = 0, finite_logs = 0, diffuse_logs = 0;
    double maximum_projection_residual = 0;
    const double z_norm = norm(model.observation);
    bool scalar_transition = true;
    for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i < m; ++i) {
            const double required = i == j ? model.transition(0, 0) : 0;
            if (model.transition(i, j) != required) scalar_transition = false;
        }
    }
    bool observation_null_certified = false;
    for (std::size_t t = 0; t < n; ++t) {
        if ((t & 31U) == 0) checkpoint(interrupt);
        rank_before[t] = basis.columns;
        bool exact_independent = false;
        bool exact_zero = false;
        if (history && basis.columns > 0) {
            if (basis.columns != initial_rank - history->rank()) {
                throw std::runtime_error("floating diffuse rank disagrees with exact history certificate");
            }
            exact_independent = history->classify_current(!series.missing[t]);
            exact_zero = !exact_independent;
        }
        long double predicted = mean;
        for (std::size_t i = 0; i < m; ++i) predicted += static_cast<long double>(model.observation[i]) * state[i];
        predictions[t] = checked(predicted);
        // Case E4: the reduced route keeps the r-wide projected loading
        // (d4_pz) without the m-wide zero padding; every full-path branch
        // below assigns the m-wide vector in full, as before.
        long double finite_f = model.observation_variance;
        if (d4_state == 2) {
            // Reduced observation on the certified r-wide factor: every
            // dropped term multiplies a certified-zero history row, so the
            // loading sum is the full-path sum with its certified-zero tail
            // terms omitted in place. pz is padded with exact zeros, keeping
            // the m-wide gain, state and gain-history sums in their original
            // term order (zeros appended at the certified positions).
            auto projected = finite_factor->observation(d4_observation, model.observation_variance, interrupt);
            // Case E4: the r-wide loading moves out directly; the gain
            // division below reads it (case E3), so the m-wide zero
            // padding and its copy — which only re-created exact +0
            // entries at the certified positions — are not computed.
            d4_pz = std::move(projected.pz);
            finite_f = projected.variance;
            ++d4_certified_steps;
        } else if (finite_factor) {
            auto projected = finite_factor->observation(model.observation, model.observation_variance, interrupt);
            pz = std::move(projected.pz);
            finite_f = projected.variance;
        } else if (k0_structural) {
            struct_observation(p, model.observation_variance, pz, finite_f);
        } else {
            pz = product(p, model.observation);
            for (std::size_t i = 0; i < m; ++i) finite_f += static_cast<long double>(model.observation[i]) * pz[i];
        }
        finite_variances[t] = checked(finite_f);
        if (finite_variances[t] < 0) throw std::runtime_error("negative finite innovation variance: diffuse arithmetic is numerically unresolved");
        auto w = observation_loading(model.observation, basis);
        double w_norm = norm(w);
        // An observed update makes Z B_post = 0 exactly. T=cI preserves this
        // identity through any number of missing periods. Only the exact
        // scalar-identity pattern enables the certificate; approximate
        // proportionality is never used to discard weak observability.
        // Canonical SARIMA has an additional exact certificate from the
        // observed homogeneous-history design. It establishes dependence
        // without numerical rank thresholds or changing B's coordinates.
        if ((observation_null_certified || exact_zero) && basis.columns > 0) {
            const long double scale = static_cast<long double>(z_norm) * norm(basis.values);
            const double relative = scale == 0 ? 0 : checked(w_norm / scale);
            maximum_projection_residual = std::max(maximum_projection_residual, relative);
            if (relative > projection_tolerance * static_cast<double>(m)) {
                throw std::runtime_error("certified diffuse observation-null identity failed its residual check");
            }
            std::fill(w.begin(), w.end(), 0);
            w_norm = 0;
        }
        if (exact_independent && w_norm == 0) {
            throw std::runtime_error("exact history pivot is numerically unresolved in the diffuse loading");
        }
        diffuse_variances[t] = square(w_norm);
        if (!series.missing[t]) {
            ++observed;
            const double innovation = checked(static_cast<long double>(series.values[t]) - predicted);
            innovations.values[t] = innovation;
            std::vector<double> gain(m, 0);
            if (w_norm > 0) {
                const long double scale = static_cast<long double>(z_norm) * norm(basis.values);
                if (w_norm <= projection_tolerance * static_cast<long double>(m) * scale) {
                    throw std::invalid_argument("numerically ambiguous diffuse innovation; rescale the diffuse coordinates or model");
                }
                for (std::size_t i = 0; i < m; ++i) {
                    long double total = 0;
                    for (std::size_t j = 0; j < basis.columns; ++j) total += static_cast<long double>(basis(i, j)) * (w[j] / w_norm);
                    gain[i] = checked(total / w_norm);
                }
                if (finite_factor) finite_factor->update(gain, model.observation, model.observation_variance, interrupt);
                else p = joseph(p, gain, model.observation, model.observation_variance);
                basis = remove_observed_direction(basis, w, model.observation, w_norm,
                                                    maximum_projection_residual, interrupt);
                diffuse_logs += 2 * std::log(static_cast<long double>(w_norm));
                ++diffuse_count;
                if (history) ++exact_pivot_count;
                last_diffuse = t + 1;
                update_kind[t] = 1;
                if (d4_state == 0 && basis.columns == 0) {
                    // Pivot exhaustion: certify the P^+Z = 0 nullity on the
                    // just-updated factor — the induction premise of the
                    // certified-zero history block. On rejection the full-m
                    // path runs to the end (route 3); nothing changes.
                    double relative = 0;
                    const bool nullity = finite_factor->certify_observation_nullity(model.observation, relative);
                    if (relative > d4_certificate_residual) d4_certificate_residual = relative;
                    // The induction re-pins the history block through every
                    // observation update: a step without an observation
                    // leaves the process noise loaded on the history
                    // coordinates unpinned, so the window must be fully
                    // observed (qualification finding: the monitor census on
                    // real-airline-gaps-x reaches 0.83 relative after a gap
                    // inside the window versus 4e-16 with none).
                    bool window_observed = true;
                    for (std::size_t s = t + 1; s < n; ++s) {
                        if (series.missing[s]) { window_observed = false; break; }
                    }
                    d4_state = (nullity && window_observed) ? 1 : 3;
                }
            } else {
                if (finite_variances[t] == 0) throw std::invalid_argument("finite observations require positive innovation variance; deterministic likelihood is unsupported");
                const double residual = checked(innovation / std::sqrt(finite_variances[t]));
                if (innovation != 0 && residual == 0) throw std::underflow_error("finite standardized innovation underflows");
                const double residual_square = square(residual);
                squares += residual_square;
                finite_logs += std::log(static_cast<long double>(finite_variances[t]));
                standardized.values[t] = residual;
                standardized.missing[t] = 0;
                if (d4_state == 2) {
                    // Case E3: divide the r-wide loading directly; the
                    // certified-zero gain tail is written as the exact +0.0
                    // the padded division produced (checked(+0/F) == +0 for
                    // F > 0), so the m-wide gain buffer — read by the state
                    // update and the regression whitening lane — keeps its
                    // bit content.
                    for (std::size_t i = 0; i < d4_r; ++i) gain[i] = checked(d4_pz[i] / finite_variances[t]);
                    for (std::size_t i = d4_r; i < m; ++i) gain[i] = 0;
                } else
                for (std::size_t i = 0; i < m; ++i) gain[i] = checked(pz[i] / finite_variances[t]);
                if (d4_state == 2) {
                    // Joseph compression on the reduced factor: the history
                    // gain entries are the padded exact zeros, so restricting
                    // to r stationary entries drops only certified-zero terms.
                    // Case E3 (buffer reuse): the persistent r-wide copy
                    // holds the same entries the per-step temporary did
                    // (same bit content, same call), no per-step allocation.
                    for (std::size_t i = 0; i < d4_r; ++i) d4_gain_r[i] = gain[i];
                    finite_factor->update(d4_gain_r, d4_observation, model.observation_variance, interrupt);
                } else if (finite_factor) finite_factor->update(gain, model.observation, model.observation_variance, interrupt);
                else if (k0_structural) struct_joseph(p, gain, model.observation_variance, k0_scratch);
                else p = joseph(p, gain, model.observation, model.observation_variance);
                ++finite_count;
                if (exact_zero) ++exact_zero_count;
                update_kind[t] = 2;
            }
            if (save_gains) std::copy(gain.begin(), gain.end(), observation_gains.begin() + static_cast<std::ptrdiff_t>(m * t));
            // Case E3b: the state additions stay m-long in every route. The
            // r-long form (dossier stage-4 probe) measured zero gain and
            // would leave a sign-of-zero residual on a -0.0 coordinate; the
            // m-long loop is bit-identical to the full path BY CONSTRUCTION
            // for O(m) long-double FMAs.
            for (std::size_t i = 0; i < m; ++i) state[i] = checked(static_cast<long double>(state[i]) + static_cast<long double>(gain[i]) * innovation);
            observation_null_certified = true;
            if (stable_history) stable_history->observe(static_cast<std::uint64_t>(t) + 1, interrupt);
        }
        rank_after[t] = basis.columns;
        if (history && basis.columns != initial_rank - history->rank()) {
            throw std::runtime_error("diffuse update disagrees with exact history rank");
        }
        if (stable_history && stable_history->remaining_rank() != basis.columns) {
            throw std::runtime_error("stable history frame disagrees with the exact diffuse rank");
        }
        if (t + 1 < n) {
            if (history) history->advance();
            observation_null_certified = observation_null_certified && scalar_transition;
            if (d4_state == 2 && d4_block_state) {
                // Case E1 (block-structured state advance): bit-identical to
                // product(model.transition, state) under the companion-shape
                // witness verified at the switch. Stationary rows of the
                // canonical companion carry exact zeros on the history
                // columns, and adding an exact ±0 term never changes a long
                // double accumulator (it starts +0; +0+(±0)=+0; a nonzero
                // value is unchanged by ±0), so their m-term sums equal the
                // r-term sums entry by entry. The pivot row keeps the full
                // m-term accumulation in the original ascending order — this
                // is the observed-history ring entry, a genuine m-wide
                // quantity. The shift rows have a single 1.0 coefficient:
                // exact copies.
                for (std::size_t i = 0; i < d4_r; ++i) {
                    long double value = 0;
                    for (std::size_t j = 0; j < d4_r; ++j)
                        value += static_cast<long double>((*d4_transition)(i, j)) * state[j];
                    d4_state_scratch[i] = checked(value);
                }
                {
                    long double value = 0;
                    for (std::size_t j = 0; j < m; ++j)
                        value += static_cast<long double>(model.transition(d4_r, j)) * state[j];
                    d4_state_scratch[d4_r] = checked(value);
                }
                for (std::size_t i = d4_r + 1; i < m; ++i) d4_state_scratch[i] = state[i - 1];
                state.swap(d4_state_scratch);
            } else {
                state = product(model.transition, state);
            }
            if (d4_state == 2) {
                finite_factor->predict(*d4_transition, *d4_noise, interrupt);
            } else if (d4_state == 1) {
                // Pivot-exhaust step: this predict still runs full-m. The last
                // unidentified diffuse coordinate leaves the state with the
                // shift, so EVERY history row of the PREDICTED factor is a
                // rounding residual of a mathematically zero row; certify that
                // block, then compress the stationary rows once. On rejection
                // the full-m path continues (route 3).
                finite_factor->predict(model.transition, model.noise_loading, interrupt);
                double relative = 0;
                const bool block = finite_factor->certify_history_block(d4_r, relative, interrupt);
                if (relative > d4_certificate_residual) d4_certificate_residual = relative;
                if (block) {
                    d4_observation.assign(model.observation.begin(),
                                          model.observation.begin() + static_cast<std::ptrdiff_t>(d4_r));
                    d4_transition = std::make_unique<DenseMatrix>(d4_r, d4_r);
                    d4_noise = std::make_unique<DenseMatrix>(d4_r, 1);
                    for (std::size_t j = 0; j < d4_r; ++j) {
                        for (std::size_t i = 0; i < d4_r; ++i) (*d4_transition)(i, j) = model.transition(i, j);
                        (*d4_noise)(j, 0) = model.noise_loading(j, 0);
                    }
                    finite_factor->compress_to_stationary_block(d4_r, interrupt);
                    // Case E1 premise (fail-closed): the case-E1 state
                    // advance is used only under the canonical companion
                    // shape of model.transition — stationary rows exactly
                    // zero on the history columns, shift rows a single 1.0
                    // coefficient; the pivot row is summed m-wide either
                    // way. Checked once per fit; a shape failure simply
                    // keeps the general product() (the route and every
                    // output are untouched either way).
                    bool companion = true;
                    for (std::size_t i = 0; i < d4_r && companion; ++i)
                        for (std::size_t j = d4_r; j < m; ++j)
                            if (model.transition(i, j) != 0) { companion = false; break; }
                    for (std::size_t i = d4_r + 1; i < m && companion; ++i)
                        for (std::size_t j = 0; j < m; ++j)
                            if (model.transition(i, j) != (j + 1 == i ? 1.0 : 0.0)) { companion = false; break; }
                    d4_block_state = companion;
                    d4_state_scratch.assign(m, 0);
                    d4_gain_r.assign(d4_r, 0);
                    d4_state = 2;
                } else {
                    d4_state = 3;
                }
            }
            else if (finite_factor) finite_factor->predict(model.transition, model.noise_loading, interrupt);
            else if (k0_structural) struct_predict_add(p, model.transition, process, k0_scratch);
            else {
                p = sandwich(model.transition, p);
                for (std::size_t i = 0; i < p.values().size(); ++i) p.values()[i] = checked(static_cast<long double>(p.values()[i]) + process.values()[i]);
            }
            if (basis.columns > 0) {
                if (stable_history) {
                    basis = history_loading(*stable_history, m, static_cast<std::uint64_t>(t) + 2, interrupt);
                    ++history_reconstructions;
                } else basis = propagate(model.transition, basis, model.observation, !series.missing[t]);
                validate_basis(basis, interrupt, "transition loses diffuse rank or is numerically rank deficient; initial coordinates remain unidentified");
                ++rank_checks;
            }
        }
    }
    if (d4_state == 2) {
        // Reassemble the m-wide output covariance from the certified block
        // structure blockdiag(P_r, 0): the history block of the full path
        // holds the same zeros at rounding level (certified above), so the
        // reduced block plus exact zeros is the same named quantity within
        // the certificate residual.
        const auto reduced_covariance = finite_factor->covariance(interrupt);
        DenseMatrix full(m, m);
        for (std::size_t j = 0; j < d4_r; ++j)
            for (std::size_t i = 0; i < d4_r; ++i) full(i, j) = reduced_covariance(i, j);
        p = std::move(full);
    } else if (finite_factor) p = finite_factor->covariance(interrupt);
    d4_route = d4_state == 2 ? 2 : d4_state == 3 ? 3 : d4_state < 0 ? 0 : 1;
    const double sum_squares = checked(squares), sum_finite_logs = checked(finite_logs), sum_diffuse_logs = checked(diffuse_logs);
    const double unavailable = std::numeric_limits<double>::quiet_NaN();
    double fixed_loglik = unavailable, profiled_variance = unavailable, profiled_loglik = unavailable;
    std::string likelihood_status = "unresolved_diffuse", profile_status = "unresolved_diffuse";
    if (basis.columns == 0) {
        likelihood_status = observed == 0 ? "no_observations" : "available";
        profile_status = "no_finite_innovations";
        const auto n_finite = static_cast<long double>(finite_count);
        const long double log_two_pi = std::log(2 * std::acos(-1.0L));
        fixed_loglik = checked(-.5L * (diffuse_logs + finite_logs +
            n_finite * (log_two_pi + std::log(static_cast<long double>(variance))) + squares / variance));
        if (finite_count > 0) {
            // The named profiled variance must be bit-identical to the binary64
            // quotient consumers recompute from the named sum (forecast,
            // interpolation and arma information divide the stored binary64
            // sum_squared_standardized in binary64), so it is a binary64
            // division of binary64 operands, not a once-rounded scratch
            // quotient (extended-precision-range-contract rule 3). The
            // quotient cannot overflow: it is bounded by the finite sum.
            profiled_variance = sum_squares / static_cast<double>(finite_count);
            if (sum_squares > 0 && profiled_variance == 0) throw std::underflow_error("profiled diffuse variance underflows");
            if (profiled_variance == 0) {
                profile_status = "zero_scale_boundary";
                profiled_loglik = std::numeric_limits<double>::infinity();
            } else {
                profile_status = "available";
                profiled_loglik = checked(-.5L * (diffuse_logs + finite_logs +
                    n_finite * (log_two_pi + 1 + std::log(static_cast<long double>(profiled_variance)))));
            }
        }
    }
    checkpoint(interrupt);
    return {std::move(predictions), std::move(innovations), std::move(standardized),
            std::move(finite_variances), std::move(diffuse_variances), std::move(update_kind),
            std::move(rank_before), std::move(rank_after), observed, diffuse_count, finite_count,
            initial_rank, basis.columns, last_diffuse, sum_squares, sum_finite_logs, sum_diffuse_logs,
            fixed_loglik, profiled_variance, profiled_loglik, std::move(likelihood_status),
            std::move(profile_status), std::move(state), initial_finite, initial_diffuse,
            std::move(p), gram(basis),
            {rank_tolerance, projection_tolerance, rank_checks, maximum_projection_residual,
             history != nullptr, exact_zero_count, exact_pivot_count,
             stable_history != nullptr, history_reconstructions,
             stable_history ? stable_history->rebuild_count() : 0,
             stable_history ? stable_history->normalization_residual() : 0,
             finite_factor != nullptr,
             finite_factor ? finite_factor->compression_count() : 0,
             finite_factor ? finite_factor->clipped_negative_count() : 0,
             finite_factor ? finite_factor->initialization_residual() : 0,
             d4_route, d4_certified_steps, d4_certificate_residual,
             finite_factor ? finite_factor->post_pivot_rebuild_count() : 0,
             finite_factor ? finite_factor->post_pivot_rebuild_residual() : 0
             },
            std::move(observation_gains)};
}

} // namespace

DiffuseLikelihood diffuse_innovations(const MaskedSeries& series, const DiffuseStateSpace& model,
                                      const DiffuseInitialState& initial, double mean,
                                      double variance, const Interrupt& interrupt) {
    validate_model(model, initial.state, mean, variance, interrupt);
    initial.finite_loading.validate(interrupt);
    const auto m = model.transition.rows();
    if (initial.finite_loading.rows() != m || initial.finite_loading.columns() > max_states) {
        throw std::invalid_argument("initial finite loading dimensions must agree and have at most 32 columns");
    }
    return filter(series, model, initial.state, gram(initial.finite_loading),
                  Basis{m, initial.diffuse_rank, initial.diffuse_loading}, mean, variance, interrupt);
}

SarimaLikelihood sarima_innovations(const MaskedSeries& series, const std::vector<double>& ar,
                                   const std::vector<double>& ma, const std::vector<double>& seasonal_ar,
                                   const std::vector<double>& seasonal_ma, int period, int d,
                                   int seasonal_d, double mean, double variance, const Interrupt& interrupt,
                                   bool save_gains) {
    if (period != 4 && period != 12) throw std::invalid_argument("SARIMA period must be 4 or 12");
    if (d < 0 || d > 2 || seasonal_d < 0 || seasonal_d > 2) throw std::invalid_argument("d and seasonal_d must be between 0 and 2");
    if (!std::isfinite(mean)) throw std::invalid_argument("mean must be finite");
    if ((d > 0 || seasonal_d > 0) && mean != 0) throw std::invalid_argument("integrated SARIMA requires mean zero; drift is not implemented");
    auto arma = stationary_arma(ar, ma, seasonal_ar, seasonal_ma, period, interrupt);
    auto difference = difference_polynomial(d, seasonal_d, period, interrupt);
    const auto r = arma.state_space.transition.rows(), k = difference.size() - 1, m = r + k;
    if (m > max_states) throw std::invalid_argument("combined SARIMA model exceeds the 32-state diffuse pilot bound");
    const auto stationary = stationary_initialization(arma.state_space, interrupt);
    DenseMatrix transition(m, m), loading(m, 1), initial_finite(m, m);
    std::vector<double> z(m, 0);
    for (std::size_t i = 0; i < r; ++i) {
        z[i] = arma.state_space.observation[i];
        loading(i, 0) = arma.state_space.noise_loading(i, 0);
        for (std::size_t j = 0; j < r; ++j) {
            transition(i, j) = arma.state_space.transition(i, j);
            initial_finite(i, j) = stationary.covariance(i, j);
        }
    }
    Basis basis{m, k, std::vector<double>(m * k, 0)};
    for (std::size_t i = 0; i < k; ++i) {
        z[r + i] = -difference[i + 1];
        basis.values[r + i + m * i] = 1;
    }
    if (k > 0) {
        for (std::size_t j = 0; j < m; ++j) transition(r, j) = z[j];
        for (std::size_t i = 1; i < k; ++i) transition(r + i, r + i - 1) = 1;
    }
    DiffuseStateSpace model{std::move(transition), std::move(z), std::move(loading), 0};
    detail::CountsHistoryCertificate history(d, seasonal_d, period);
    detail::SparseHistoryBasis stable_history(d, seasonal_d, period, interrupt);
    std::unique_ptr<detail::FiniteStateFactor> finite_factor;
    if (k > 0) finite_factor = std::make_unique<detail::FiniteStateFactor>(initial_finite, interrupt);
    auto likelihood = filter(series, model, std::vector<double>(m, 0), std::move(initial_finite),
                              std::move(basis), mean, variance, interrupt,
                              k > 0 ? &history : nullptr, k > 0 ? &stable_history : nullptr,
                              finite_factor.get(), save_gains, k == 0);
    return {std::move(model), std::move(arma.ar), std::move(arma.ma),
            std::move(difference), std::move(likelihood)};
}

} // namespace kalos
