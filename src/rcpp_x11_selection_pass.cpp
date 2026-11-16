#include "rcpp_x11_table_support.hpp"
#include "kalos/x11_selection_pass.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
using Request=kalos::X11SeasonalRequest;
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
// One pinned SeasonalFilterOption request per calendar period, named exactly
// as in the qualified seasonal-plan boundary.
std::vector<Request> requests(SEXP value) {
    if(TYPEOF(value)!=STRSXP || Rf_xlength(value)<1)
        throw std::invalid_argument("selection requests must be a character vector of seasonal requests");
    std::vector<Request> out;out.reserve(static_cast<std::size_t>(Rf_xlength(value)));
    for(R_xlen_t i=0;i<Rf_xlength(value);++i) {
        if(STRING_ELT(value,i)==NA_STRING)
            throw std::invalid_argument("seasonal requests cannot be missing");
        const auto name=Rcpp::as<std::string>(STRING_ELT(value,i));
        if(name=="stable")out.push_back(Request::stable);
        else if(name=="s3x1")out.push_back(Request::s3x1);
        else if(name=="s3x3")out.push_back(Request::s3x3);
        else if(name=="s3x5")out.push_back(Request::s3x5);
        else if(name=="s3x9")out.push_back(Request::s3x9);
        else if(name=="s3x15")out.push_back(Request::s3x15);
        else if(name=="x11default")out.push_back(Request::x11default);
        else if(name=="msr")out.push_back(Request::msr);
        else throw std::invalid_argument("unsupported seasonal request");
    }
    return out;
}
const char* policy_name(kalos::X11SigmaPolicy policy) {
    return policy==kalos::X11SigmaPolicy::legacy ? "legacy" : "accurate";
}
kalos::X11SigmaPolicy parse_policy(const std::string& name,const char* label) {
    if(name=="legacy")return kalos::X11SigmaPolicy::legacy;
    if(name=="accurate")return kalos::X11SigmaPolicy::accurate;
    throw std::invalid_argument(std::string(label)+" must be legacy or accurate");
}
kalos::X11MsrIterationPolicy parse_msr_policy(const std::string& name) {
    if(name=="jdplus")return kalos::X11MsrIterationPolicy::jdplus;
    if(name=="x11_five_year_cap")return kalos::X11MsrIterationPolicy::x11_five_year_cap;
    throw std::invalid_argument("MSR iteration policy must be jdplus or x11_five_year_cap");
}
const char* stage_status_name(kalos::X11TrendStageStatus status) {
    using S=kalos::X11TrendStageStatus;
    switch(status) {
    case S::fixed:return "fixed";
    case S::selected:return "selected";
    case S::unsupported:return "unsupported";
    }
    throw std::logic_error("unknown trend stage status");
}
const char* band_name(kalos::X11MsrBand band) {
    using B=kalos::X11MsrBand;
    switch(band) {
    case B::unavailable:return "unavailable";
    case B::s3x3:return "s3x3";
    case B::lower_gap:return "lower_gap";
    case B::s3x5:return "s3x5";
    case B::upper_gap:return "upper_gap";
    case B::s3x9:return "s3x9";
    }
    throw std::logic_error("unknown MSR band");
}
const char* msr_status_name(kalos::X11MsrSelectionStatus status) {
    using S=kalos::X11MsrSelectionStatus;
    switch(status) {
    case S::available:return "available";
    case S::missing_observations:return "missing_observations";
    case S::nonpositive_relative_input:return "nonpositive_relative_input";
    case S::numerical_range:return "numerical_range";
    case S::incomplete_table:return "incomplete_table";
    case S::budget_exhausted:return "budget_exhausted";
    case S::undefined_msr:return "undefined_msr";
    }
    throw std::logic_error("unknown MSR selection status");
}
const char* msr_reason_name(kalos::X11MsrSelectionReason reason) {
    using R=kalos::X11MsrSelectionReason;
    switch(reason) {
    case R::unavailable:return "unavailable";
    case R::decisive_ratio:return "decisive_ratio";
    case R::insufficient_initial_history:return "insufficient_initial_history";
    case R::insufficient_remaining_history:return "insufficient_remaining_history";
    case R::five_year_cap:return "five_year_cap";
    case R::iteration_budget:return "iteration_budget";
    }
    throw std::logic_error("unknown MSR selection reason");
}
const char* resolution_name(kalos::X11SeasonalResolution resolution) {
    using R=kalos::X11SeasonalResolution;
    switch(resolution) {
    case R::fixed:return "fixed_choice";
    case R::x11default_initial:return "x11default_initial";
    case R::x11default_final:return "x11default_final";
    case R::msr_initial:return "msr_initial";
    case R::msr_preliminary:return "msr_preliminary_final";
    case R::msr_mixed:return "msr_mixed_request";
    case R::requires_msr:return "requires_msr";
    }
    throw std::logic_error("unknown seasonal resolution");
}
double real(double value) {return std::isnan(value) ? NA_REAL : value;}
Rcpp::NumericVector reals(const std::vector<double>& values) {
    Rcpp::NumericVector out(values.size());
    for(std::size_t i=0;i<values.size();++i)
        out[static_cast<R_xlen_t>(i)]=std::isnan(values[i]) ? NA_REAL : values[i];
    return out;
}
Rcpp::IntegerVector integers(const std::vector<int>& values) {
    return Rcpp::IntegerVector(values.begin(),values.end());
}
// The B and C correction-round records share their fields, so one converter
// serves both (the 1-based row convention matches the fixed-pass boundaries).
template<class Round>
Rcpp::List correction_round(const Round& round) {
    auto plan=reals(round.replacements.replacements);
    int count=0;
    for(R_xlen_t i=0;i<plan.size();++i) if(!ISNAN(plan[i])) ++count;
    return Rcpp::List::create(Rcpp::_["detection"]=kalos_x11_rcpp::values(round.detection),
        Rcpp::_["analysis_position"]=static_cast<int>(round.analysis_position)+1, // 1-based rows
        Rcpp::_["policy"]=policy_name(round.policy),
        Rcpp::_["sigma_weights"]=reals(round.sigma.weights),Rcpp::_["replacements"]=plan,
        Rcpp::_["corrections"]=reals(round.replacements.correction_factors),
        Rcpp::_["replacement_count"]=count);
}
Rcpp::List msr_table(const kalos::X11MsrDiagnostics& msr) {
    return Rcpp::List::create(
        Rcpp::_["mean_irregular_evolutions"]=reals(msr.mean_irregular_evolutions),
        Rcpp::_["mean_seasonal_evolutions"]=reals(msr.mean_seasonal_evolutions),
        Rcpp::_["counts"]=Rcpp::IntegerVector(msr.counts.begin(),msr.counts.end()),
        Rcpp::_["msr"]=reals(msr.msr),Rcpp::_["global_msr"]=real(msr.global_msr));
}
Rcpp::List stage(const kalos::X11TrendStage& stage) {
    return Rcpp::List::create(Rcpp::_["status"]=stage_status_name(stage.status),
        Rcpp::_["length"]=stage.length,Rcpp::_["ic_ratio"]=real(stage.ic_ratio));
}
Rcpp::List plan(const kalos::X11SeasonalPlan& plan) {
    Rcpp::CharacterVector resolutions(plan.resolutions.size());
    for(std::size_t i=0;i<plan.resolutions.size();++i)
        resolutions[static_cast<R_xlen_t>(i)]=resolution_name(plan.resolutions[i]);
    return Rcpp::List::create(Rcpp::_["base_filters"]=integers(plan.base_filters),
        Rcpp::_["resolutions"]=resolutions,Rcpp::_["global_msr"]=plan.global_msr,
        Rcpp::_["requires_msr"]=plan.requires_msr);
}
// The nested selected-state passes keep the exact shapes of the fixed-pass
// boundaries, so a selection result's $b chains into .native_x11_c_pass and
// its $c chains into .native_x11_d_pass unchanged.
Rcpp::List b_result(const kalos::X11BPass& b,const std::string& policy,double lower,
        double upper,int frequency,int first) {
    using kalos_x11_rcpp::values;
    Rcpp::List tables=Rcpp::List::create(Rcpp::_["B1"]=values(b.b1),Rcpp::_["B2"]=values(b.b2),
        Rcpp::_["B3"]=values(b.b3),Rcpp::_["B4"]=values(b.b4),Rcpp::_["B4a"]=values(b.b4a),
        Rcpp::_["B4anorm"]=values(b.b4anorm),Rcpp::_["B4d"]=values(b.b4d),
        Rcpp::_["B4g"]=values(b.b4g),Rcpp::_["B5"]=values(b.b5),Rcpp::_["B6"]=values(b.b6),
        Rcpp::_["B7"]=values(b.b7),Rcpp::_["B8"]=values(b.b8),Rcpp::_["B9"]=values(b.b9),
        Rcpp::_["B9g"]=values(b.b9g),Rcpp::_["B10"]=values(b.b10),Rcpp::_["B11"]=values(b.b11),
        Rcpp::_["B13"]=values(b.b13),Rcpp::_["B17"]=values(b.b17),Rcpp::_["B20"]=values(b.b20));
    Rcpp::List rounds=Rcpp::List::create(correction_round(b.rounds[0]),
        correction_round(b.rounds[1]),correction_round(b.rounds[2]));
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_b_pass",
        Rcpp::_["policy"]=policy,Rcpp::_["frequency"]=frequency,Rcpp::_["start_period"]=first,
        Rcpp::_["trend_length"]=b.trend_length,Rcpp::_["b2drop"]=static_cast<int>(b.b2drop),
        Rcpp::_["lower_sigma"]=lower,Rcpp::_["upper_sigma"]=upper,
        Rcpp::_["initial_filters"]=integers(b.initial_filters),
        Rcpp::_["final_filters"]=integers(b.final_filters),
        Rcpp::_["tables"]=tables,Rcpp::_["rounds"]=rounds,Rcpp::_["status"]=b.status,
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(b.mode),
        Rcpp::_["scope"]=std::string("selected_state_")
            +kalos_x11_rcpp::mode_name(b.mode)+"_b_pass",Rcpp::_["backend"]="native");
}
Rcpp::List c_result(const kalos::X11CPass& c,const std::string& policy,double lower,
        double upper,int frequency) {
    using kalos_x11_rcpp::values;
    Rcpp::List tables=Rcpp::List::create(Rcpp::_["C1"]=values(c.c1),Rcpp::_["C2"]=values(c.c2),
        Rcpp::_["C4"]=values(c.c4),Rcpp::_["C5"]=values(c.c5),Rcpp::_["C5a"]=values(c.c5a),
        Rcpp::_["C6"]=values(c.c6),Rcpp::_["C7"]=values(c.c7),Rcpp::_["C9"]=values(c.c9),
        Rcpp::_["C10"]=values(c.c10),Rcpp::_["C10a"]=values(c.c10a),
        Rcpp::_["C11"]=values(c.c11),Rcpp::_["C13"]=values(c.c13),
        Rcpp::_["C17"]=values(c.c17),Rcpp::_["C20"]=values(c.c20));
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_c_pass",
        Rcpp::_["policy"]=policy,Rcpp::_["frequency"]=frequency,
        Rcpp::_["start_period"]=c.first_period,Rcpp::_["trend_length"]=c.trend_length,
        Rcpp::_["c2drop"]=static_cast<int>(c.c2drop),Rcpp::_["lower_sigma"]=lower,
        Rcpp::_["upper_sigma"]=upper,Rcpp::_["initial_filters"]=integers(c.initial_filters),
        Rcpp::_["final_filters"]=integers(c.final_filters),Rcpp::_["tables"]=tables,
        Rcpp::_["round"]=correction_round(c.round),
        Rcpp::_["reference_series"]=values(c.reference_series),Rcpp::_["status"]=c.status,
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(c.mode),
        Rcpp::_["scope"]=std::string("selected_state_")
            +kalos_x11_rcpp::mode_name(c.mode)+"_c_pass",Rcpp::_["backend"]="native");
}
Rcpp::List d_result(const kalos::X11DPass& d,const std::string& policy,int frequency,int first) {
    using kalos_x11_rcpp::values;
    Rcpp::List tables=Rcpp::List::create(Rcpp::_["D1"]=values(d.d1),Rcpp::_["D2"]=values(d.d2),
        Rcpp::_["D4"]=values(d.d4),Rcpp::_["D5"]=values(d.d5),Rcpp::_["D6"]=values(d.d6),
        Rcpp::_["D7"]=values(d.d7),Rcpp::_["D8"]=values(d.d8),Rcpp::_["D9"]=values(d.d9),
        Rcpp::_["D9bis"]=values(d.d9bis),Rcpp::_["D9_g_bis"]=values(d.d9bis),
        Rcpp::_["D10"]=values(d.d10),Rcpp::_["D10bis"]=values(d.d10bis),
        Rcpp::_["D11"]=values(d.d11),Rcpp::_["D11bis"]=values(d.d11bis),
        Rcpp::_["D12"]=values(d.d12),Rcpp::_["D13"]=values(d.d13));
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_d_pass",
        Rcpp::_["policy"]=policy,Rcpp::_["frequency"]=frequency,Rcpp::_["start_period"]=first,
        Rcpp::_["trend_length"]=d.trend_length,Rcpp::_["d2drop"]=static_cast<int>(d.d2drop),
        Rcpp::_["backcast"]=d.backcast,Rcpp::_["forecast"]=d.forecast,
        Rcpp::_["d7_henderson_length"]=d.d7_henderson_length,
        Rcpp::_["final_henderson_length"]=d.final_henderson_length,
        Rcpp::_["msr_iterations"]=d.msr_iterations,
        Rcpp::_["initial_filters"]=integers(d.initial_filters),
        Rcpp::_["final_filters"]=integers(d.final_filters),
        Rcpp::_["seas_filter"]=integers(d.seas_filter),Rcpp::_["tables"]=tables,
        Rcpp::_["msr_table"]=msr_table(d.msr_table),Rcpp::_["status"]=d.status,
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(d.mode),
        Rcpp::_["scope"]=std::string("selected_state_")
            +kalos_x11_rcpp::mode_name(d.mode)+"_d_pass",Rcpp::_["backend"]="native");
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
    parse_policy(name,"selection-pass sigma weight policy");
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
    const auto result=kalos::x11_selection_pass(series,period,first,requests(requests_vector),
        trend,parse_policy(name,"selection-pass sigma weight policy"),
        parse_msr_policy(msr_name),lower,upper,back,fc,decomp,[] {Rcpp::checkUserInterrupt();});
    Rcpp::List stages=Rcpp::List::create(Rcpp::_["b7"]=stage(result.b7),
        Rcpp::_["c7"]=stage(result.c7),Rcpp::_["d7"]=stage(result.d7),
        Rcpp::_["final"]=stage(result.final));
    std::vector<int> att_count,att_years;std::vector<double> att_msr;
    std::vector<std::string> att_band,att_status;
    for(const auto& row:result.msr.attempts) {
        att_count.push_back(static_cast<int>(row.count));
        att_years.push_back(static_cast<int>(row.years_removed));
        att_msr.push_back(real(row.global_msr));
        att_band.push_back(band_name(row.band));
        att_status.push_back(msr_status_name(row.status));
    }
    Rcpp::List attempts=Rcpp::List::create(Rcpp::_["count"]=integers(att_count),
        Rcpp::_["years_removed"]=integers(att_years),Rcpp::_["global_msr"]=reals(att_msr),
        Rcpp::_["band"]=Rcpp::wrap(att_band),Rcpp::_["status"]=Rcpp::wrap(att_status));
    Rcpp::List msr=Rcpp::List::create(Rcpp::_["engaged"]=result.msr.engaged,
        Rcpp::_["policy"]=msr_name,Rcpp::_["status"]=msr_status_name(result.msr.status),
        Rcpp::_["reason"]=msr_reason_name(result.msr.reason),
        Rcpp::_["selected_filter"]=result.msr.selected_filter,
        Rcpp::_["using_default"]=result.msr.using_default,
        Rcpp::_["iterations"]=static_cast<int>(result.msr.iterations),
        Rcpp::_["global_msr"]=real(result.msr.global_msr),Rcpp::_["attempts"]=attempts);
    Rcpp::List out=Rcpp::List::create(Rcpp::_["schema_version"]=1,
        Rcpp::_["type"]="x11_selection_pass",Rcpp::_["status"]=result.status,
        Rcpp::_["frequency"]=period,Rcpp::_["start_period"]=first,
        Rcpp::_["requested_trend_length"]=trend,Rcpp::_["automatic"]=trend==0,
        Rcpp::_["policy"]=name,Rcpp::_["msr_policy"]=msr_name,
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(decomp),
        Rcpp::_["lower_sigma"]=lower,Rcpp::_["upper_sigma"]=upper,
        Rcpp::_["backcast"]=back,Rcpp::_["forecast"]=fc,
        Rcpp::_["requests"]=Rcpp::CharacterVector(requests_vector),
        Rcpp::_["initial_filters"]=integers(result.initial_filters),
        Rcpp::_["final_filters"]=integers(result.final_filters),
        Rcpp::_["stages"]=stages,Rcpp::_["msr"]=msr,
        Rcpp::_["initial_plan"]=plan(result.initial_plan),
        Rcpp::_["final_plan"]=plan(result.final_plan),
        Rcpp::_["scope"]="pass_level_selection_only",Rcpp::_["backend"]="native");
    if(result.status=="available") { // unsupported runs keep records only
        out["b"]=b_result(result.b,name,lower,upper,period,first);
        out["c"]=c_result(result.c,name,lower,upper,period);
        out["d"]=d_result(result.d,name,period,first);
    }
    return out;
}
