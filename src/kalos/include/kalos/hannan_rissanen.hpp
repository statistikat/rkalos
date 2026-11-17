#ifndef KALOS_HANNAN_RISSANEN_HPP
#define KALOS_HANNAN_RISSANEN_HPP

// Transcription of the pinned jdplus HannanRissanenInitializer and its
// HannanRissanen pipeline (Hannan-Rissanen three-stage ARMA initialization on
// a differenced response), built on the qualified linear-algebra adapters.
// This header owns no solver and no likelihood; it produces the starting
// coefficients the estimation-control layer feeds to the optimizer.
#include "regression_likelihood.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace kalos {

// Stationary SarimaMapping orders of the pinned initializer: p regular AR
// coefficients, bp seasonal AR blocks, q regular MA, bq seasonal MA blocks.
// The pinned initializer runs this pipeline on the differenced model, so the
// integrated orders d and D never enter here.
struct HannanRissanenOrders {
    int p;
    int bp;
    int q;
    int bq;
    int period;
};

// Starting coefficients with their provenance. Blocks use the shared Kalos
// sign convention (AR 1 - phi(B), MA 1 + theta(B)) and are ordered
// [ar, sar, ma, sma] with sizes p, bp, q, bq.
struct HannanRissanenStart {
    enum class Source {
        // The pipeline completed; stabilized records whether the pinned
        // SarimaMapping.stabilize step rewrote at least one block.
        Fitted,
        // Pinned root-mean-square guard on the regression residuals
        // (rms < 1e-9): the start is zero on every block.
        ZeroVariance,
        // Pinned failure branch (process false or any exception): the
        // builder setDefault fill, Kalos ar/sar +0.1, ma/sma -0.2.
        Defaulted,
    };
    Source source;
    bool stabilized;
    std::array<std::vector<double>, 4> coefficients;
};

// response and design are the differenced (stationary) model quantities:
// complete finite response of length n and its column-major n-by-k design
// (design.columns == 0 for none). The transcription follows the pinned
// pipeline: OLS residuals (or the response itself without design), the
// rms variance guard, AbsMean normalization, the Levinson AR prefilter of the
// pinned npi() length, the dense two-sided regression on lagged response and
// residual columns, the bias correction, and the final correction by
// refiltering through the estimated AR polynomial and re-running the MA part.
// Builder defaults biascorrection = finalcorrection = true and the
// initializer settings stabilize = true, useDefaultIfFailed = true,
// failifunstable = false are the pinned RegSarimaComputer call-site values
// and are not exposed.
//
// Static input violations (negative orders, empty or non-finite response,
// design/response mismatch) throw std::invalid_argument. Numerical failures
// anywhere in the pipeline mirror the pinned catch-all: the result is the
// Defaulted fill, never an exception; interrupt exceptions always propagate.
HannanRissanenStart hannan_rissanen_start(const std::vector<double>& response,
                                          const RegressionDesign& design,
                                          const HannanRissanenOrders& orders,
                                          const Interrupt& interrupt = {});

} // namespace kalos
#endif
