#include "rcpp_x13_shared.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
namespace kalos_x13_rcpp {
namespace {
kalos::X13VariableClass variable_class(SEXP value) {
    using C=kalos::X13VariableClass;
    if(TYPEOF(value)!=STRSXP || Rf_xlength(value)!=1 || STRING_ELT(value,0)==NA_STRING)
        throw std::invalid_argument("x13 variable class must be outlier, days_related, moving_holiday or user");
    const auto name=Rcpp::as<std::string>(value);
    if(name=="outlier")return C::outlier;
    if(name=="days_related")return C::days_related;
    if(name=="moving_holiday")return C::moving_holiday;
    if(name=="user")return C::user;
    throw std::invalid_argument("x13 variable class must be outlier, days_related, moving_holiday or user");
}
kalos::X13Component variable_effect(SEXP value) {
    using E=kalos::X13Component;
    if(TYPEOF(value)!=STRSXP || Rf_xlength(value)!=1 || STRING_ELT(value,0)==NA_STRING)
        throw std::invalid_argument(
            "x13 variable effect must be undefined, series, trend, seasonal, "
            "seasonally_adjusted, irregular or calendar_effect");
    const auto name=Rcpp::as<std::string>(value);
    if(name=="undefined")return E::undefined;
    if(name=="series")return E::series;
    if(name=="trend")return E::trend;
    if(name=="seasonal")return E::seasonal;
    if(name=="seasonally_adjusted")return E::seasonally_adjusted;
    if(name=="irregular")return E::irregular;
    if(name=="calendar_effect")return E::calendar_effect;
    throw std::invalid_argument(
        "x13 variable effect must be undefined, series, trend, seasonal, "
        "seasonally_adjusted, irregular or calendar_effect");
}
}
std::vector<double> reals_in(SEXP value,const char* label) {
    if((TYPEOF(value)!=REALSXP && TYPEOF(value)!=INTSXP) || Rf_isFactor(value) ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(std::string(label)+" must be a numeric vector");
    const auto y=Rcpp::as<Rcpp::NumericVector>(value);
    std::vector<double> out;out.reserve(static_cast<std::size_t>(y.size()));
    for(R_xlen_t i=0;i<y.size();++i)out.push_back(y[i]); // NA/NaN pass through; the pinned seam judges
    return out;
}
Rcpp::NumericVector reals_out(const std::vector<double>& values) {
    Rcpp::NumericVector out(values.size());
    for(std::size_t i=0;i<values.size();++i)
        out[static_cast<R_xlen_t>(i)]=std::isnan(values[i]) ? NA_REAL : values[i];
    return out;
}
// Fetch one named element of a pairlist-carrying VECSXP (pure R API, so a
// partially hand-built object fails with a clear error instead of an
// extractor surprise).
SEXP named(SEXP object,const char* name,const char* label) {
    const SEXP names=Rf_getAttrib(object,R_NamesSymbol);
    if(TYPEOF(object)!=VECSXP || TYPEOF(names)!=STRSXP)
        throw std::invalid_argument(label);
    for(R_xlen_t i=0;i<Rf_xlength(names);++i)
        if(STRING_ELT(names,i)!=NA_STRING && std::string(name)==CHAR(STRING_ELT(names,i)))
            return VECTOR_ELT(object,i);
    throw std::invalid_argument(std::string(label)+": missing '"+name+"'");
}
// The pinned model-order variable list: one list per variable with class,
// effect, one coefficient per regression column and the columns themselves
// (the pinned Regression.matrix columns over the input calendar).
std::vector<kalos::X13AllocationVariable> variables(SEXP value) {
    if(TYPEOF(value)!=VECSXP)
        throw std::invalid_argument("x13 model variables must be a list of variable specifications");
    const char* shape="each x13 variable needs class, effect, coefficient and columns";
    std::vector<kalos::X13AllocationVariable> out;out.reserve(static_cast<std::size_t>(Rf_xlength(value)));
    for(R_xlen_t k=0;k<Rf_xlength(value);++k) {
        const SEXP row=VECTOR_ELT(value,k);
        if(TYPEOF(row)!=VECSXP)throw std::invalid_argument(shape);
        kalos::X13AllocationVariable v;
        v.cls=variable_class(named(row,"class",shape));
        v.effect=variable_effect(named(row,"effect",shape));
        v.coefficient=reals_in(named(row,"coefficient",shape),"x13 variable coefficient");
        const SEXP columns=named(row,"columns",shape);
        if(TYPEOF(columns)!=VECSXP)
            throw std::invalid_argument("x13 variable columns must be a list of numeric column vectors");
        for(R_xlen_t j=0;j<Rf_xlength(columns);++j)
            v.columns.push_back(reals_in(VECTOR_ELT(columns,j),"x13 variable column"));
        out.push_back(std::move(v));
    }
    return out;
}
Rcpp::List preadjustment_list(const kalos::X13Preadjustment& pre,bool log) {
    Rcpp::List tables=Rcpp::List::create(
        Rcpp::_["a1"]=reals_out(pre.a1),Rcpp::_["a6"]=reals_out(pre.a6),
        Rcpp::_["a7"]=reals_out(pre.a7),Rcpp::_["a8"]=reals_out(pre.a8),
        Rcpp::_["a8t"]=reals_out(pre.a8t),Rcpp::_["a8s"]=reals_out(pre.a8s),
        Rcpp::_["a8i"]=reals_out(pre.a8i),Rcpp::_["a9"]=reals_out(pre.a9),
        Rcpp::_["a9u"]=reals_out(pre.a9u),Rcpp::_["a9cal"]=reals_out(pre.a9cal),
        Rcpp::_["a9sa"]=reals_out(pre.a9sa),Rcpp::_["a9ser"]=reals_out(pre.a9ser),
        Rcpp::_["out"]=reals_out(pre.out));
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x13_preadjustment",
        Rcpp::_["log"]=log,Rcpp::_["length"]=static_cast<int>(pre.a1.size()),
        Rcpp::_["tables"]=tables,Rcpp::_["status"]=pre.status,
        Rcpp::_["scope"]="deterministic_effect_allocation_only",Rcpp::_["backend"]="native");
}
kalos::X13Preadjustment rebuild_preadjustment(SEXP value) {
    if(TYPEOF(value)!=VECSXP)
        throw std::invalid_argument("x13 restore requires a preadjustment result object");
    const char* shape="x13 restore requires a preadjustment result object";
    if(Rcpp::as<std::string>(named(value,"type",shape))!="x13_preadjustment")
        throw std::invalid_argument(shape);
    const SEXP status=named(value,"status",shape);
    if(TYPEOF(status)!=STRSXP || Rf_xlength(status)!=1 || STRING_ELT(status,0)==NA_STRING)
        throw std::invalid_argument(shape);
    kalos::X13Preadjustment pre;
    pre.status=Rcpp::as<std::string>(status);
    const SEXP tables=named(value,"tables",shape);
    const std::array<const char*,13> names={"a1","a6","a7","a8","a8t","a8s","a8i","a9",
        "a9u","a9cal","a9sa","a9ser","out"};
    using Member=std::vector<double> kalos::X13Preadjustment::*;
    static const Member members[13]=
        {&kalos::X13Preadjustment::a1,&kalos::X13Preadjustment::a6,&kalos::X13Preadjustment::a7,
        &kalos::X13Preadjustment::a8,&kalos::X13Preadjustment::a8t,&kalos::X13Preadjustment::a8s,
        &kalos::X13Preadjustment::a8i,&kalos::X13Preadjustment::a9,&kalos::X13Preadjustment::a9u,
        &kalos::X13Preadjustment::a9cal,&kalos::X13Preadjustment::a9sa,
        &kalos::X13Preadjustment::a9ser,&kalos::X13Preadjustment::out};
    for(std::size_t k=0;k<names.size();++k)
        pre.*(members[k])=reals_in(named(tables,names[k],shape),"x13 preadjustment table");
    return pre;
}
Rcpp::List finals_list(const kalos::X13Finals& result) {
    return Rcpp::List::create(
        Rcpp::_["d11final"]=reals_out(result.d11final),
        Rcpp::_["d12final"]=reals_out(result.d12final),
        Rcpp::_["d13final"]=reals_out(result.d13final),
        Rcpp::_["d16"]=reals_out(result.d16),Rcpp::_["d18"]=reals_out(result.d18),
        Rcpp::_["e1"]=reals_out(result.e1),Rcpp::_["e2"]=reals_out(result.e2),
        Rcpp::_["e3"]=reals_out(result.e3),Rcpp::_["e11"]=reals_out(result.e11));
}
} // namespace kalos_x13_rcpp
