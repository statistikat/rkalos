#ifndef KALOS_RCPP_X11_TABLE_SUPPORT_HPP
#define KALOS_RCPP_X11_TABLE_SUPPORT_HPP
// Binding-only serialization shared by standalone primitives and composition.
#include <Rcpp.h>
#include "kalos/x11_seasonal_processing.hpp"
#include "kalos/x11_normalization.hpp"
#include "kalos/x11_decomposition.hpp"
namespace kalos_x11_rcpp {
// Decomposition-mode parsing shared by the four pass boundaries: NULL keeps
// the additive default (the pre-mode positional calls stay valid), the four
// pinned mode names are accepted verbatim.
kalos::X11DecompositionMode parse_mode(SEXP mode,const char* label);
const char* mode_name(kalos::X11DecompositionMode mode);
Rcpp::NumericVector values(const kalos::MaskedSeries& series);
Rcpp::List seasonal_result(const kalos::X11SeasonalProcessing& result,int frequency,
    const std::vector<int>& choices);
Rcpp::List normalization_result(const kalos::X11SeasonalNormalization& result,int frequency,
    const std::vector<int>& choices,int extension);
}
#endif
