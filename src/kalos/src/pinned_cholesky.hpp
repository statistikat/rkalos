#ifndef KALOS_PINNED_CHOLESKY_HPP
#define KALOS_PINNED_CHOLESKY_HPP

#include <vector>

namespace kalos::detail {
// Factor a column-major symmetric matrix in place using the JDemetra
// lapack.Cholesky arithmetic. A non-positive pivot makes the factor unusable.
bool pinned_lcholesky(std::vector<double>& matrix, int size);
}

#endif
