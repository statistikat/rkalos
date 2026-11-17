#ifndef KALOS_OUTLIER_DETAIL_HPP
#define KALOS_OUTLIER_DETAIL_HPP

// Internal helpers shared by the two outlier-scoring translation units.
// Everything here is a literal transcription of the pinned jdplus 3.9.0
// bytecode qualified in docs/contracts/outlier-scoring-contract.md; the
// semantics notes name the pinned source. Not part of the public API.

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace kalos {
namespace outlier_detail {

// JVM double comparisons with explicit NaN ordering, per the class-file
// semantics (jvmse 6.5): dcmpl pushes -1 when either operand is NaN, dcmpg
// pushes +1. javac emits dcmpl for a ">" test and dcmpg for a "<" test, so
// both are false on NaN; the pinned call sites fix which one a given
// comparison used and the NaN side effects follow from it.
inline int dcmpl(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return -1;
    return a > b ? 1 : (a < b ? -1 : 0);
}

inline int dcmpg(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return 1;
    return a > b ? 1 : (a < b ? -1 : 0);
}

// java.lang.Math.max/min: NaN on either side yields NaN (std::fmax does not).
inline double java_max(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::nan("");
    return a >= b ? a : b;
}

inline double java_min(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::nan("");
    return a <= b ? a : b;
}

// Pinned LULinearSystemSolver2(1e-13) (LinearSystemSolver.robustSolver()),
// transcribed operation for operation from the bytecode: a column-major
// n x n matrix a and a length-n right-hand side x are solved in place.
// Row i of the augmented system takes the largest |entry| over columns not
// yet pivoted (dcmpl: a NaN entry loses, ties keep the earlier column); the
// Gauss-Jordan elimination runs over the unpivoted columns of the FULL
// n-column augmented storage (the right-hand side travels in the trailing
// column, mirroring the pinned columnOf wrapper).
//
// Pinned failure mapping: no unpivoted non-NaN-losing candidate (ip < 0) or
// |pivot| < 1e-13 (dcmpg) throws the pinned MatrixException; both surface as
// std::runtime_error here. The pinned update subtracts (pivotInverse * v) *
// row_i before row i is scaled, which is bit-distinct from scaling first;
// the order below is the pinned one.
inline void robust_lu_solve(std::vector<double>& a, std::size_t n, std::vector<double>& x) {
    constexpr double eps = 1.0e-13; // LS_Factory ROBUST_FACTORY: LULinearSystemSolver2(1e-13)
    if (a.size() < n * n || x.size() < n) {
        throw std::invalid_argument("robust_lu_solve: dimensions");
    }
    // Augmented storage: n right-hand side columns appended after A (the
    // pinned solve builds a n x (n + m) linear block; m == 1 here).
    const std::size_t m = n + 1;
    std::vector<double> storage(n * m);
    for (std::size_t c = 0; c < n; ++c) {
        for (std::size_t r = 0; r < n; ++r) storage[c * n + r] = a[c * n + r];
    }
    for (std::size_t r = 0; r < n; ++r) storage[n * n + r] = x[r];

    std::vector<int> piv(m, -1);
    for (std::size_t i = 0; i < n; ++i) {
        double big = std::numeric_limits<double>::denorm_min(); // 4.9E-324
        int ip = -1;
        for (std::size_t k = 0; k < n; ++k) {
            if (piv[k] != -1) continue;
            const double xk = std::fabs(storage[k * n + i]);
            if (dcmpl(xk, big) > 0) {
                ip = static_cast<int>(k);
                big = xk;
            }
        }
        if (ip < 0) throw std::runtime_error("matrix is singular (m_err_dim)");
        piv[static_cast<std::size_t>(ip)] = static_cast<int>(i);
        const double pivot = storage[static_cast<std::size_t>(ip) * n + i];
        if (dcmpg(std::fabs(pivot), eps) < 0) {
            throw std::runtime_error("matrix is singular (m_err_sing)");
        }
        const double pivot_inverse = 1.0 / pivot;
        const std::size_t ip_base = static_cast<std::size_t>(ip) * n;
        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const double v = storage[ip_base + j];
            if (dcmpl(v, 0.0) == 0) continue;
            const double mult = pivot_inverse * v;
            for (std::size_t k = 0; k < m; ++k) {
                if (piv[k] != -1) continue;
                storage[k * n + j] -= mult * storage[k * n + i];
            }
        }
        for (std::size_t k = 0; k < m; ++k) {
            if (piv[k] != -1) continue;
            storage[k * n + i] *= pivot_inverse;
        }
    }
    // Pinned read-back: over B rows k (each matrix column got a pivot row),
    // B[k] = storage[piv[k] + rhsColumnOffset]; row piv[k] holds the reduced
    // solution of variable k and was already scaled by the pivot inverse.
    for (std::size_t k = 0; k < n; ++k) {
        if (piv[k] == -1) continue;
        x[k] = storage[n * n + static_cast<std::size_t>(piv[k])];
    }
}

} // namespace outlier_detail
} // namespace kalos
#endif
