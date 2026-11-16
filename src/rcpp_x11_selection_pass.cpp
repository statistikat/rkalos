#include "rcpp_x11_selection_records.hpp"
#include "kalos/x11_selection_pass.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
double scalar(SEXP value,const char* label) {
    if((TYPEOF(value)!=REALSXP && TYPEOF(value)!=INTSXP) || Rf_isFactor(value) ||
        Rf_xlength(value)!=1 || !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(std::string(label)+" must be one real number");
    const auto out=Rcpp::as<double>(value);
    if(!std::isfinite(out))throw std::invalid_argument(std::string(label)+" must be finite");
    return out;
}
int integer(SEXP value,const char* label) {
    const auto out=scalar(value,label);
    if(out!=std::floor(out) || out<-2147483647.0 || out>2147483647.0)
        throw std::invalid_argument(std::string(label)+" must be a supported integer");
    return static_cast<int>(out);
}
kalos::MaskedSeries masked(SEXP values) {
    if((TYPEOF(values)!=REALSXP && TYPEOF(values)!=INTSXP) || Rf_isFactor(values) ||
        !Rf_isNull(Rf_getAttrib(values,R_DimSymbol)) ||
        Rf_xlength(values)>std::numeric_limits<int>::max())
        throw std::invalid_argument("selection-pass input must be a real numeric vector within the supported length");
    const auto y=Rcpp::as<Rcpp::NumericVector>(values);const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=y[static_cast<R_xlen_t>(i)];
        if(std::isnan(value))out.missing[i]=1;
        else if(!std::isfinite(value))throw std::invalid_argument("selection-pass input must be finite or missing");
        else out.values[i]=value;
    }
    return out;
}
}
// [[Rcpp::export(name = ".native_x11_selection_pass")]]
Rcpp::List native_x11_selection_pass(SEXP values,SEXP frequency,SEXP start_period,
    SEXP requests_vector,SEXP trend_length,SEXP policy,SEXP msr_policy,
    SEXP lower_sigma,SEXP upper_sigma,SEXP backcast,SEXP forecast,SEXP mode=R_NilValue) {
    if(TYPEOF(policy)!=STRSXP || Rf_xlength(policy)!=1 || STRING_ELT(policy,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(policy,R_DimSymbol)))
        throw std::invalid_argument("selection-pass sigma weight policy must be legacy or accurate");
    const auto name=Rcpp::as<std::string>(policy);
    kalos_x11_rcpp::parse_policy(name,"selection-pass sigma weight policy");
    if(TYPEOF(msr_policy)!=STRSXP || Rf_xlength(msr_policy)!=1 ||
        STRING_ELT(msr_policy,0)==NA_STRING || !Rf_isNull(Rf_getAttrib(msr_policy,R_DimSymbol)))
        throw std::invalid_argument("MSR iteration policy must be jdplus or x11_five_year_cap");
    const auto msr_name=Rcpp::as<std::string>(msr_policy);
    const auto period=integer(frequency,"frequency"),first=integer(start_period,"start period");
    const auto trend=integer(trend_length,"trend filter length");
    const auto lower=scalar(lower_sigma,"lower sigma limit"),upper=scalar(upper_sigma,"upper sigma limit");
    const auto back=integer(backcast,"backcast horizon"),fc=integer(forecast,"forecast horizon");
    const auto series=masked(values);
    const auto decomp=kalos_x11_rcpp::parse_mode(mode,"selection-pass decomposition mode");
    const auto result=kalos::x11_selection_pass(series,period,first,
        kalos_x11_rcpp::requests(requests_vector),
        trend,kalos_x11_rcpp::parse_policy(name,"selection-pass sigma weight policy"),
        kalos_x11_rcpp::parse_msr_policy(msr_name),lower,upper,back,fc,decomp,
        [] {Rcpp::checkUserInterrupt();});
    // The serialization lives in the shared records TU so the restricted
    // X-13 fit delivers the identical nested record for its selection run.
    return kalos_x11_rcpp::selection_result(result,name,msr_name,
        kalos_x11_rcpp::mode_name(decomp),"pass_level_selection_only",
        lower,upper,period,first,trend,back,fc,Rcpp::CharacterVector(requests_vector));
}
