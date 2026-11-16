#ifndef KALOS_RCPP_X13_SHARED_HPP
#define KALOS_RCPP_X13_SHARED_HPP
// X-13 boundary marshalling shared by the allocation-seam and restricted-fit
// binding TUs: the allocation vocabulary (pinned variable classes, exact
// regeffect buckets, coefficient/column shapes) and the stored-table records
// keep one shape wherever the seams are crossed.
#include <Rcpp.h>
#include "kalos/x13_allocation.hpp"
#include <string>
#include <vector>
namespace kalos_x13_rcpp {
std::vector<double> reals_in(SEXP value,const char* label);
Rcpp::NumericVector reals_out(const std::vector<double>& values);
// Fetch one named element of a pairlist-carrying VECSXP (pure R API, so a
// partially hand-built object fails with a clear error instead of an
// extractor surprise).
SEXP named(SEXP object,const char* name,const char* label);
// The pinned model-order variable list: one list per variable with class,
// effect, one coefficient per regression column and the columns themselves
// (the pinned Regression.matrix columns over the input calendar).
std::vector<kalos::X13AllocationVariable> variables(SEXP value);
Rcpp::List preadjustment_list(const kalos::X13Preadjustment& pre,bool log);
Rcpp::List finals_list(const kalos::X13Finals& finals);
// The restore consumes the caller's finished preadjustment; the pinned seam
// enforces the coupling (available status, non-empty calendar, shared table
// lengths) on the rebuilt record, so the boundary only rebuilds faithfully.
kalos::X13Preadjustment rebuild_preadjustment(SEXP value);
}
#endif
