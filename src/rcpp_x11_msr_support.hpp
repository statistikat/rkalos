#ifndef KALOS_RCPP_X11_MSR_SUPPORT_HPP
#define KALOS_RCPP_X11_MSR_SUPPORT_HPP
// Binding helpers only. No core header includes this Rcpp-dependent adapter.
#include <Rcpp.h>
#include "kalos/x11_msr.hpp"
namespace kalos_msr_rcpp {
int integer(SEXP value,const char* label);
kalos::MaskedSeries series(SEXP value);
kalos::X11MsrMode mode(SEXP value);
Rcpp::List table_result(const kalos::X11MsrTable& result,int frequency,kalos::X11MsrMode mode);
}
#endif
