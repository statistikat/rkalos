#ifndef KALOS_DETAIL_SARIMA_STABILIZATION_HPP
#define KALOS_DETAIL_SARIMA_STABILIZATION_HPP

// Internal seam for the pinned SarimaMapping.stabilize(DataBlock) block
// routine (item 21-A.2 transcription, shared by the Hannan-Rissanen
// initializer and the differencing decision unit of item 21-D). The
// routine is convention-equivariant under the AR sign flip: applying it to
// a block of 1 + b(B) coefficients or to the negated coefficients of
// 1 - b(B) yields sign-mirrored results with identical magnitudes.
#include "kalos/linear_algebra.hpp"
#include <vector>

namespace kalos::detail {

std::vector<ComplexRoot> jvm_closed_roots(const std::vector<double>& ascending_coefficients);
double jvm_complex_abs(double real, double imaginary);

// Pinned SarimaMapping.stabilize(DataBlock): trim trailing zeros, leave a
// stable block, invert roots inside the unit circle of 1 + b(B) by
// re-forming the polynomial from the inverted roots and normalizing the
// constant term to one. Returns whether the block was rewritten.
// Numerical dead ends (zero root, non-finite reconstruction) throw;
// callers map them onto their pinned failure policies.
bool stabilize_coefficient_block(std::vector<double>& block);

} // namespace kalos::detail
#endif
