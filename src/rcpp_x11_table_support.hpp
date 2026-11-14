#ifndef KALOS_RCPP_X11_TABLE_SUPPORT_HPP
#define KALOS_RCPP_X11_TABLE_SUPPORT_HPP
// Binding-only serialization shared by standalone primitives and composition.
#include <Rcpp.h>
#include "kalos/x11_seasonal_processing.hpp"
#include "kalos/x11_normalization.hpp"
namespace kalos_x11_rcpp {
Rcpp::NumericVector values(const kalos::MaskedSeries& series);
Rcpp::List seasonal_result(const kalos::X11SeasonalProcessing& result,int frequency,
    const std::vector<int>& choices);
Rcpp::List normalization_result(const kalos::X11SeasonalNormalization& result,int frequency,
    const std::vector<int>& choices,int extension);
}
#endif
