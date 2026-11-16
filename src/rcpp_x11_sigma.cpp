#include <Rcpp.h>
#include "kalos/x11_sigma.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
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
bool flag(SEXP value,const char* label) {
    if(TYPEOF(value)!=LGLSXP || Rf_xlength(value)!=1 || LOGICAL(value)[0]==NA_LOGICAL ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(std::string(label)+" must be TRUE or FALSE");
    return LOGICAL(value)[0]!=0;
}
kalos::MaskedSeries masked(SEXP values) {
    if((TYPEOF(values)!=REALSXP && TYPEOF(values)!=INTSXP) || Rf_isFactor(values) ||
        !Rf_isNull(Rf_getAttrib(values,R_DimSymbol)) ||
        Rf_xlength(values)>std::numeric_limits<int>::max())
        throw std::invalid_argument("Sigma input must be a real numeric vector within the supported length");
    const auto y=Rcpp::as<Rcpp::NumericVector>(values);const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=y[static_cast<R_xlen_t>(i)];
        if(std::isnan(value))out.missing[i]=1;
        else if(!std::isfinite(value))throw std::invalid_argument("Sigma input must be finite or missing");
        else out.values[i]=value;
    }
    return out;
}
constexpr std::array<const char*,5> kWindowKinds={
    "global_single","leading_extension","central","trailing_extension","boundary_reuse"};
constexpr std::array<const char*,6> kWeightStates={
    "full","ramp","zero","full_below_dispersion_cutoff","full_without_dispersion","missing_input"};
Rcpp::List window_list(const std::vector<kalos::X11SigmaWindow>& windows) {
    Rcpp::NumericVector dispersion;Rcpp::IntegerVector start,length,usable;
    Rcpp::CharacterVector kind;
    for(const auto& window:windows) {
        dispersion.push_back(window.dispersion);
        start.push_back(static_cast<int>(window.start));
        length.push_back(static_cast<int>(window.length));
        usable.push_back(static_cast<int>(window.usable));
        kind.push_back(kWindowKinds.at(static_cast<std::size_t>(window.kind)));
    }
    return Rcpp::List::create(Rcpp::_["dispersions"]=dispersion,Rcpp::_["starts"]=start,
        Rcpp::_["lengths"]=length,Rcpp::_["usable"]=usable,Rcpp::_["kinds"]=kind);
}
}
// [[Rcpp::export(name = ".native_x11_sigma_weights")]]
Rcpp::List native_x11_sigma_weights(SEXP values,SEXP frequency,SEXP start_period,SEXP policy,
    SEXP backcasts,SEXP forecasts,SEXP exclude,SEXP lower_sigma,SEXP upper_sigma,
    SEXP multiplicative) {
    if(TYPEOF(policy)!=STRSXP || Rf_xlength(policy)!=1 || STRING_ELT(policy,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(policy,R_DimSymbol)))
        throw std::invalid_argument("Sigma weight policy must be legacy or accurate");
    const auto name=Rcpp::as<std::string>(policy);
    if(name!="legacy" && name!="accurate")
        throw std::invalid_argument("Sigma weight policy must be legacy or accurate");
    const auto period=integer(frequency,"frequency"),first=integer(start_period,"start period");
    const auto head=integer(backcasts,"backcast count"),tail=integer(forecasts,"forecast count");
    if(head<0||tail<0)throw std::invalid_argument("projection counts must be nonnegative");
    const auto lower=scalar(lower_sigma,"lower sigma limit"),upper=scalar(upper_sigma,"upper sigma limit");
    const auto series=masked(values);
    const auto result=kalos::x11_sigma_weights(series,period,first,
        name=="legacy"?kalos::X11SigmaPolicy::legacy:kalos::X11SigmaPolicy::accurate,
        head,tail,flag(exclude,"projection exclusion"),lower,upper,
        flag(multiplicative,"multiplicative mode"),[] {Rcpp::checkUserInterrupt();});
    Rcpp::NumericVector weights(result.weights.begin(),result.weights.end());
    Rcpp::NumericVector first_pass(result.first_pass_weights.begin(),result.first_pass_weights.end());
    Rcpp::CharacterVector states;
    for(const auto state:result.states)states.push_back(kWeightStates.at(static_cast<std::size_t>(state)));
    Rcpp::IntegerVector zeroed;
    for(const auto row:result.first_pass_zeroed)
        zeroed.push_back(static_cast<int>(row)+1); // 1-based observation rows
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_sigma_weights",
        Rcpp::_["policy"]=name,Rcpp::_["frequency"]=period,Rcpp::_["start_period"]=first,
        Rcpp::_["backcasts"]=head,Rcpp::_["forecasts"]=tail,
        Rcpp::_["exclude_projections"]=flag(exclude,"projection exclusion"),
        Rcpp::_["lower_sigma"]=lower,Rcpp::_["upper_sigma"]=upper,
        Rcpp::_["multiplicative"]=flag(multiplicative,"multiplicative mode"),
        Rcpp::_["weights"]=weights,Rcpp::_["states"]=states,
        Rcpp::_["windows"]=window_list(result.windows),
        Rcpp::_["first_pass_weights"]=first_pass,
        Rcpp::_["first_pass_windows"]=window_list(result.first_pass_windows),
        Rcpp::_["first_pass_zeroed"]=zeroed,Rcpp::_["passes"]=result.passes,
        Rcpp::_["status"]=result.status,Rcpp::_["scope"]="sigma_weights_only",
        Rcpp::_["backend"]="native");
}
