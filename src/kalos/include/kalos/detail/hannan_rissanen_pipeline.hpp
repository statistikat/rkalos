#ifndef KALOS_DETAIL_HANNAN_RISSANEN_PIPELINE_HPP
#define KALOS_DETAIL_HANNAN_RISSANEN_PIPELINE_HPP

// Internal seam for the pinned raw HannanRissanen pipeline
// (jdplus.toolkit.base.core.sarima.estimation.HannanRissanen, builder
// defaults biasCorrection = finalCorrection = true, Levinson
// initialization). The pipeline itself has neither the
// HannanRissanenInitializer root-mean-square guard nor the initializer
// stabilization and default fill, and it never reports a soft failure:
// sarima.estimation.HannanRissanen.calc() returns true on every path, so
// a numerical problem always surfaces as an exception. The pinned public
// entry point kalos::hannan_rissanen_start (item 21-A.2) wraps this seam
// with the initializer guard, stabilization, and default fill; the
// differencing decision unit of item 21-D calls the raw seam directly.
#include "kalos/regression_likelihood.hpp" // Interrupt, Response-only seam
#include "kalos/hannan_rissanen.hpp"       // HannanRissanenOrders

#include <array>
#include <vector>

namespace kalos::detail {

// Runs the pinned pipeline on a finite response with no regression stage
// (the differencing call site hands the filter- and mean-processed series
// directly). Returns the updatemodel-extracted blocks [ar(p), sar(bp),
// ma(q), sma(bq)] in the JVM raw convention (polynomials 1 + phi(B),
// 1 + theta(B)), un-stabilized. Pipeline numerical failures throw;
// interrupt exceptions always propagate.
std::array<std::vector<double>, 4> hannan_rissanen_pipeline(const std::vector<double>& response,
    const HannanRissanenOrders& orders, const Interrupt& interrupt = {});

} // namespace kalos::detail
#endif
