#ifndef KALOS_DETAIL_FILTER_UTILITY_HPP
#define KALOS_DETAIL_FILTER_UTILITY_HPP

// Transcription of jdplus.toolkit.base.core.math.linearfilters.FilterUtility
// .checkRoots(DoubleSeq, double) (pinned javap dump item21h/
// v390_FilterUtility.txt), qualified against the JVM by the 21-B ssq-domain
// goldens. The effective degree is the position after the last non-zero
// coefficient; degrees 0..2 are closed form and exact; degree >= 3 runs the
// deterministic Durand-Kerner port of polynomial_roots (the JVM companion
// solver is seeded by java.util.Random, so its >= 3 verdicts are the
// documented non-reproducible surface - the corpus never reaches it).

#include <cstddef>

namespace kalos::detail {

// coeffs[0..degree) are the polynomial coefficients c1..ck of
// 1 + c1 z + ... + ck z^k (the raw AR/MA block, leading 1 implicit).
// Returns true when every inverse root modulus stays below limit.
bool check_roots(const double* coeffs, std::size_t length, double limit);

} // namespace kalos::detail

#endif
