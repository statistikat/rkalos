#include <Rcpp.h>
#include "kalos/x11_replacement.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
namespace {
int integer(SEXP value,const char* label) {
    if((TYPEOF(value)!=REALSXP && TYPEOF(value)!=INTSXP) || Rf_isFactor(value) ||
        Rf_xlength(value)!=1 || !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(std::string(label)+" must be one real number");
    const auto out=Rcpp::as<double>(value);
    if(!std::isfinite(out) || out!=std::floor(out) || out<-2147483647.0 || out>2147483647.0)
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
        throw std::invalid_argument("Replacement input must be a real numeric vector within the supported length");
    const auto y=Rcpp::as<Rcpp::NumericVector>(values);const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=y[static_cast<R_xlen_t>(i)];
        if(std::isnan(value))out.missing[i]=1;
        else if(!std::isfinite(value))throw std::invalid_argument("Replacement input must be finite or missing");
        else out.values[i]=value;
    }
    return out;
}
std::vector<double> weights_of(SEXP values,std::size_t n) {
    if((TYPEOF(values)!=REALSXP && TYPEOF(values)!=INTSXP) || Rf_isFactor(values) ||
        !Rf_isNull(Rf_getAttrib(values,R_DimSymbol)) ||
        static_cast<std::size_t>(Rf_xlength(values))!=n)
        throw std::invalid_argument("Observation weights must be a finite real vector of the series length");
    const auto w=Rcpp::as<Rcpp::NumericVector>(values);
    std::vector<double> out(n,0);
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=w[static_cast<R_xlen_t>(i)];
        if(!std::isfinite(value))throw std::invalid_argument("Observation weights must be finite");
        out[i]=value;
    }
    return out;
}
constexpr std::array<const char*,4> kKinds={"none","four_neighbors","period_average","outside_span"};
}
// [[Rcpp::export(name = ".native_x11_replacements")]]
Rcpp::List native_x11_replacements(SEXP values,SEXP frequency,SEXP weights,SEXP backcasts,
    SEXP forecasts,SEXP exclude,SEXP multiplicative) {
    const auto period=integer(frequency,"frequency");
    const auto head=integer(backcasts,"backcast count"),tail=integer(forecasts,"forecast count");
    if(head<0||tail<0)throw std::invalid_argument("projection counts must be nonnegative");
    const auto series=masked(values);
    const auto w=weights_of(weights,series.values.size());
    const auto plan=kalos::x11_replacements(series,period,w,head,tail,
        flag(exclude,"projection exclusion"),flag(multiplicative,"multiplicative mode"),
        [] {Rcpp::checkUserInterrupt();});
    const auto n=static_cast<std::size_t>(plan.replacements.size());
    Rcpp::NumericVector replacements(n),corrected(n),factors(n),averages(
        plan.period_averages.begin(),plan.period_averages.end());
    Rcpp::CharacterVector kinds;
    Rcpp::IntegerVector counts;
    Rcpp::NumericMatrix neighbors(n,4);
    auto put=[](double value) {return std::isnan(value)?NA_REAL:value;};
    for(std::size_t row=0;row<n;++row) {
        if((row&4095U)==0)Rcpp::checkUserInterrupt();
        replacements[static_cast<R_xlen_t>(row)]=put(plan.replacements[row]);
        corrected[static_cast<R_xlen_t>(row)]=put(plan.corrected[row]);
        factors[static_cast<R_xlen_t>(row)]=put(plan.correction_factors[row]);
        kinds.push_back(kKinds.at(static_cast<std::size_t>(plan.kinds[row])));
        counts.push_back(static_cast<int>(plan.neighbor_counts[row]));
        for(std::size_t slot=0;slot<4;++slot)
            neighbors(static_cast<R_xlen_t>(row),static_cast<R_xlen_t>(slot))=
                static_cast<double>(plan.neighbors[row][slot]); // pinned 0-based order, -1 padded
    }
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_replacements",
        Rcpp::_["frequency"]=period,Rcpp::_["backcasts"]=head,Rcpp::_["forecasts"]=tail,
        Rcpp::_["exclude_projections"]=flag(exclude,"projection exclusion"),
        Rcpp::_["multiplicative"]=flag(multiplicative,"multiplicative mode"),
        Rcpp::_["replacements"]=replacements,Rcpp::_["kinds"]=kinds,
        Rcpp::_["neighbors"]=neighbors,Rcpp::_["neighbor_counts"]=counts,
        Rcpp::_["period_averages"]=averages,Rcpp::_["correction_factors"]=factors,
        Rcpp::_["corrected"]=corrected,Rcpp::_["fallbacks"]=static_cast<int>(plan.fallbacks),
        Rcpp::_["status"]=plan.status,Rcpp::_["scope"]="replacement_pass_only",
        Rcpp::_["backend"]="native");
}
