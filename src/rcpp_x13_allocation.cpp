#include "rcpp_x13_shared.hpp"
#include "rcpp_x11_table_support.hpp"
#include "kalos/x13_allocation.hpp"
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
// Item 22c M3c: the two deterministic-effect allocation seams, exposed one
// entry point each. The preadjust takes the complete input plus the pinned
// model-order variables (class, exact regeffect bucket, one coefficient per
// regression column, the columns themselves); the restore consumes the
// caller's finished preadjustment object with the same faithful-rebuild/
// passthrough pattern as the C-on-B and D-on-C boundaries: the boundary only
// rebuilds the pinned record and the pinned seam enforces the coupling
// (available preadjustment, shared calendar). The pinned shape violations
// throw verbatim through Rcpp.
// [[Rcpp::export(name = ".native_x13_preadjust")]]
Rcpp::List native_x13_preadjust(SEXP values,SEXP log,SEXP variables_list) {
    if(TYPEOF(log)!=LGLSXP || Rf_xlength(log)!=1 || LOGICAL(log)[0]==NA_LOGICAL ||
        !Rf_isNull(Rf_getAttrib(log,R_DimSymbol)))
        throw std::invalid_argument("x13 preadjust log flag must be TRUE or FALSE");
    kalos::X13PreadjustInput input;
    input.input=kalos_x13_rcpp::reals_in(values,"x13 preadjust input");
    input.log=LOGICAL(log)[0]==TRUE;
    input.variables=kalos_x13_rcpp::variables(variables_list);
    if(input.input.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("x13 preadjust input within the supported length");
    const auto result=kalos::x13_preadjust(input);
    return kalos_x13_rcpp::preadjustment_list(result,input.log);
}
// [[Rcpp::export(name = ".native_x13_restore")]]
Rcpp::List native_x13_restore(SEXP preadjustment,SEXP d10,SEXP d11,SEXP d12,SEXP d13,
    SEXP c17,SEXP mode=R_NilValue) {
    kalos::X13RestoreInput input;
    input.pre=kalos_x13_rcpp::rebuild_preadjustment(preadjustment); // never recomputed here
    input.d10=kalos_x13_rcpp::reals_in(d10,"x13 restore D10");
    input.d11=kalos_x13_rcpp::reals_in(d11,"x13 restore D11");
    input.d12=kalos_x13_rcpp::reals_in(d12,"x13 restore D12");
    input.d13=kalos_x13_rcpp::reals_in(d13,"x13 restore D13");
    input.c17=kalos_x13_rcpp::reals_in(c17,"x13 restore C17");
    input.mode=kalos_x11_rcpp::parse_mode(mode,"x13 restore decomposition mode");
    const auto result=kalos::x13_restore(input);
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x13_finals",
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(input.mode),
        Rcpp::_["mode_note"]="the effective mode: non-log runs additively; "
            "log runs multiplicatively unless log_additive or pseudo_additive was requested",
        Rcpp::_["tables"]=kalos_x13_rcpp::finals_list(result),Rcpp::_["status"]=result.status,
        Rcpp::_["pre_coupling"]="checked by the pinned seam against the supplied preadjustment result",
        Rcpp::_["scope"]="deterministic_effect_allocation_only",Rcpp::_["backend"]="native");
}
