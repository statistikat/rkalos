#include "rcpp_x11_table_support.hpp"
#include "kalos/x11_b_pass.hpp"
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
        throw std::invalid_argument("B-pass input must be a real numeric vector within the supported length");
    const auto y=Rcpp::as<Rcpp::NumericVector>(values);const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=y[static_cast<R_xlen_t>(i)];
        if(std::isnan(value))out.missing[i]=1;
        else if(!std::isfinite(value))throw std::invalid_argument("B-pass input must be finite or missing");
        else out.values[i]=value;
    }
    return out;
}
// One fixed seasonal request per calendar period (0 stable, N for 3xN),
// matching the pinned per-phase arrays; value support is validated by the core.
std::vector<int> requests(SEXP value,const char* label) {
    if((TYPEOF(value)!=REALSXP && TYPEOF(value)!=INTSXP) || Rf_isFactor(value) ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(std::string(label)+" must be a numeric vector of fixed seasonal requests");
    const auto y=Rcpp::as<Rcpp::NumericVector>(value);
    std::vector<int> out;out.reserve(static_cast<std::size_t>(y.size()));
    for(R_xlen_t i=0;i<y.size();++i) {
        const auto row=y[i];
        if(!std::isfinite(row) || row!=std::floor(row) || row<0.0 || row>2147483647.0)
            throw std::invalid_argument(std::string(label)+" must contain nonnegative integer requests");
        out.push_back(static_cast<int>(row));
    }
    return out;
}
Rcpp::NumericVector replacements(const std::vector<double>& values) {
    Rcpp::NumericVector out(values.size());
    for(std::size_t i=0;i<values.size();++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        out[static_cast<R_xlen_t>(i)]=std::isnan(values[i]) ? NA_REAL : values[i];
    }
    return out;
}
const char* policy_name(kalos::X11SigmaPolicy policy) {
    return policy==kalos::X11SigmaPolicy::legacy ? "legacy" : "accurate";
}
Rcpp::List correction_round(const kalos::X11BCorrectionRound& round) {
    const auto plan=replacements(round.replacements.replacements);
    int count=0;
    for(R_xlen_t i=0;i<plan.size();++i) if(!ISNAN(plan[i])) ++count;
    Rcpp::NumericVector weights(round.sigma.weights.begin(),round.sigma.weights.end());
    const auto corrections=replacements(round.replacements.correction_factors);
    return Rcpp::List::create(Rcpp::_["detection"]=kalos_x11_rcpp::values(round.detection),
        Rcpp::_["analysis_position"]=static_cast<int>(round.analysis_position)+1, // 1-based rows
        Rcpp::_["policy"]=policy_name(round.policy),
        Rcpp::_["sigma_weights"]=weights,Rcpp::_["replacements"]=plan,
        Rcpp::_["corrections"]=corrections,Rcpp::_["replacement_count"]=count);
}
}
// [[Rcpp::export(name = ".native_x11_b_pass")]]
Rcpp::List native_x11_b_pass(SEXP values,SEXP frequency,SEXP start_period,SEXP policy,
    SEXP initial_filters,SEXP final_filters,SEXP trend_length,SEXP lower_sigma,SEXP upper_sigma) {
    if(TYPEOF(policy)!=STRSXP || Rf_xlength(policy)!=1 || STRING_ELT(policy,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(policy,R_DimSymbol)))
        throw std::invalid_argument("B-pass sigma weight policy must be legacy or accurate");
    const auto name=Rcpp::as<std::string>(policy);
    if(name!="legacy" && name!="accurate")
        throw std::invalid_argument("B-pass sigma weight policy must be legacy or accurate");
    const auto period=integer(frequency,"frequency"),first=integer(start_period,"start period");
    const auto trend=integer(trend_length,"trend filter length");
    const auto lower=scalar(lower_sigma,"lower sigma limit"),upper=scalar(upper_sigma,"upper sigma limit");
    const auto initial=requests(initial_filters,"Initial seasonal filters");
    const auto finals=requests(final_filters,"Final seasonal filters");
    const auto series=masked(values);
    const auto result=kalos::x11_b_pass(series,period,first,initial,finals,trend,
        name=="legacy"?kalos::X11SigmaPolicy::legacy:kalos::X11SigmaPolicy::accurate,
        lower,upper,[] {Rcpp::checkUserInterrupt();});
    using kalos_x11_rcpp::values;
    Rcpp::IntegerVector resolved_initial(result.initial_filters.begin(),result.initial_filters.end());
    Rcpp::IntegerVector resolved_final(result.final_filters.begin(),result.final_filters.end());
    // The 19 pinned member tables in label order; the B9d detection series
    // lives in the correction history (round 2 detection) below.
    Rcpp::List tables=Rcpp::List::create(
        Rcpp::_["B1"]=values(result.b1),Rcpp::_["B2"]=values(result.b2),
        Rcpp::_["B3"]=values(result.b3),Rcpp::_["B4"]=values(result.b4),
        Rcpp::_["B4a"]=values(result.b4a),Rcpp::_["B4anorm"]=values(result.b4anorm),
        Rcpp::_["B4d"]=values(result.b4d),Rcpp::_["B4g"]=values(result.b4g),
        Rcpp::_["B5"]=values(result.b5),Rcpp::_["B6"]=values(result.b6),
        Rcpp::_["B7"]=values(result.b7),Rcpp::_["B8"]=values(result.b8),
        Rcpp::_["B9"]=values(result.b9),Rcpp::_["B9g"]=values(result.b9g),
        Rcpp::_["B10"]=values(result.b10),Rcpp::_["B11"]=values(result.b11),
        Rcpp::_["B13"]=values(result.b13),Rcpp::_["B17"]=values(result.b17),
        Rcpp::_["B20"]=values(result.b20));
    Rcpp::List rounds=Rcpp::List::create(correction_round(result.rounds[0]),
        correction_round(result.rounds[1]),correction_round(result.rounds[2]));
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_b_pass",
        Rcpp::_["policy"]=name,Rcpp::_["frequency"]=period,Rcpp::_["start_period"]=first,
        Rcpp::_["trend_length"]=result.trend_length,Rcpp::_["b2drop"]=static_cast<int>(result.b2drop),
        Rcpp::_["lower_sigma"]=lower,Rcpp::_["upper_sigma"]=upper,
        Rcpp::_["initial_filters"]=resolved_initial,Rcpp::_["final_filters"]=resolved_final,
        Rcpp::_["tables"]=tables,Rcpp::_["rounds"]=rounds,
        Rcpp::_["status"]=result.status,Rcpp::_["scope"]="fixed_additive_b_pass_only",
        Rcpp::_["backend"]="native");
}
