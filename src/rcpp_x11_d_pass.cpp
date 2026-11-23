#include "rcpp_x11_table_support.hpp"
#include "kalos/x11_d_pass.hpp"
#include <array>
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
        throw std::invalid_argument("D-pass input must be a real numeric vector within the supported length");
    const auto y=Rcpp::as<Rcpp::NumericVector>(values);const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=y[static_cast<R_xlen_t>(i)];
        if(std::isnan(value))out.missing[i]=1;
        else if(!std::isfinite(value))throw std::invalid_argument("D-pass input must be finite or missing");
        else out.values[i]=value;
    }
    return out;
}
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
kalos::X11SigmaPolicy parse_policy(const std::string& name,const char* label) {
    if(name=="legacy")return kalos::X11SigmaPolicy::legacy;
    if(name=="accurate")return kalos::X11SigmaPolicy::accurate;
    throw std::invalid_argument(std::string(label)+" must be legacy or accurate");
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
// Optional variant: an absent field yields R_NilValue (used for the mode
// echo, which pre-mode companion objects legitimately lack — those are
// additive by construction).
SEXP named_opt(SEXP object,const char* name) {
    const SEXP names=Rf_getAttrib(object,R_NamesSymbol);
    if(TYPEOF(object)!=VECSXP || TYPEOF(names)!=STRSXP)return R_NilValue;
    for(R_xlen_t i=0;i<Rf_xlength(names);++i)
        if(STRING_ELT(names,i)!=NA_STRING && std::string(name)==CHAR(STRING_ELT(names,i)))
            return VECTOR_ELT(object,i);
    return R_NilValue;
}
// Rebuild the pinned MaskedSeries from a .native_x11_c_pass table (NA marks
// the pinned missing flags; ts decorations on the public wrapper's tables
// are irrelevant here).
kalos::MaskedSeries untable(SEXP value,const char* label) {
    if((TYPEOF(value)!=REALSXP && TYPEOF(value)!=INTSXP) || Rf_isFactor(value) ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(std::string(label)+" must be a numeric vector");
    const auto y=Rcpp::as<Rcpp::NumericVector>(value);
    const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto row=y[static_cast<R_xlen_t>(i)];
        if(ISNAN(row))out.missing[i]=1;
        else if(!std::isfinite(row))throw std::invalid_argument(std::string(label)+" must be finite or missing");
        else out.values[i]=row;
    }
    return out;
}
// The D pass consumes the caller's finished C-pass result; the pinned core
// enforces the coupling (status, setting echoes, round policy, reference
// series equality, C20 coverage) on the reconstructed object, so the boundary
// only rebuilds it faithfully and never recomputes the C pass here.
kalos::X11CPass rebuild_c(SEXP value) {
    if(TYPEOF(value)!=VECSXP)
        throw std::invalid_argument("D pass requires a C-pass result object");
    const char* shape="D pass requires a C-pass result object";
    if(Rcpp::as<std::string>(named(value,"type",shape))!="x11_c_pass")
        throw std::invalid_argument(shape);
    const SEXP tables=named(value,"tables",shape);
    const std::array<const char*,14> names={"C1","C2","C4","C5","C5a","C6","C7","C9",
        "C10","C10a","C11","C13","C17","C20"};
    using Member=kalos::MaskedSeries kalos::X11CPass::*;
    static const Member members[14]=
        {&kalos::X11CPass::c1,&kalos::X11CPass::c2,&kalos::X11CPass::c4,&kalos::X11CPass::c5,
        &kalos::X11CPass::c5a,&kalos::X11CPass::c6,&kalos::X11CPass::c7,&kalos::X11CPass::c9,
        &kalos::X11CPass::c10,&kalos::X11CPass::c10a,&kalos::X11CPass::c11,&kalos::X11CPass::c13,
        &kalos::X11CPass::c17,&kalos::X11CPass::c20};
    kalos::X11CPass c;
    c.c2drop=static_cast<std::size_t>(integer(named(value,"c2drop",shape),"C-pass c2drop"));
    c.first_period=integer(named(value,"start_period",shape),"C-pass start period");
    const SEXP status=named(value,"status",shape);
    if(TYPEOF(status)!=STRSXP || Rf_xlength(status)!=1 || STRING_ELT(status,0)==NA_STRING)
        throw std::invalid_argument(shape);
    c.status=Rcpp::as<std::string>(status);
    // The decomposition mode is carried by the C-pass result (itself
    // inherited from its B pass); the D pass runs on that same context. An
    // object predating the mode echo is additive by construction.
    c.mode=kalos_x11_rcpp::parse_mode(named_opt(value,"mode"),"C-pass decomposition mode");
    c.trend_length=integer(named(value,"trend_length",shape),"C-pass trend filter length");
    c.initial_filters=requests(named(value,"initial_filters",shape),"C-pass initial filters");
    c.final_filters=requests(named(value,"final_filters",shape),"C-pass final filters");
    c.reference_series=untable(named(value,"reference_series",shape),"C-pass reference series");
    for(std::size_t k=0;k<names.size();++k)
        c.*(members[k])=untable(named(tables,names[k],shape),"C-pass table");
    const SEXP round=named(value,"round",shape);
    if(TYPEOF(round)!=VECSXP)throw std::invalid_argument("the C-pass correction round is incomplete");
    const char* rshape="the C-pass correction round is incomplete";
    auto& target=c.round;
    target.detection=untable(named(round,"detection",rshape),"C-pass round detection");
    target.analysis_position=static_cast<std::int64_t>(
        integer(named(round,"analysis_position",rshape),"C-pass round analysis position"))-1; // back to 0-based
    target.policy=parse_policy(Rcpp::as<std::string>(named(round,"policy",rshape)),
        "C-pass round policy");
    const auto weights=Rcpp::as<Rcpp::NumericVector>(named(round,"sigma_weights",rshape));
    target.sigma.weights.assign(weights.begin(),weights.end());
    const auto vals=Rcpp::as<Rcpp::NumericVector>(named(round,"replacements",rshape));
    target.replacements.replacements.assign(vals.begin(),vals.end());
    const auto factors=Rcpp::as<Rcpp::NumericVector>(named(round,"corrections",rshape));
    target.replacements.correction_factors.assign(factors.begin(),factors.end());
    return c;
}
Rcpp::NumericVector reals(const std::vector<double>& values) {
    Rcpp::NumericVector out(values.size());
    for(std::size_t i=0;i<values.size();++i)
        out[static_cast<R_xlen_t>(i)]=std::isnan(values[i]) ? NA_REAL : values[i];
    return out;
}
}
// [[Rcpp::export(name = ".native_x11_d_pass")]]
Rcpp::List native_x11_d_pass(SEXP values,SEXP frequency,SEXP start_period,SEXP policy,
    SEXP initial_filters,SEXP final_filters,SEXP trend_length,SEXP cpass,
    SEXP backcast,SEXP forecast) {
    if(TYPEOF(policy)!=STRSXP || Rf_xlength(policy)!=1 || STRING_ELT(policy,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(policy,R_DimSymbol)))
        throw std::invalid_argument("D-pass sigma weight policy must be legacy or accurate");
    const auto name=Rcpp::as<std::string>(policy);
    if(name!="legacy" && name!="accurate")
        throw std::invalid_argument("D-pass sigma weight policy must be legacy or accurate");
    const auto period=integer(frequency,"frequency"),first=integer(start_period,"start period");
    const auto trend=integer(trend_length,"trend filter length");
    const auto back=integer(backcast,"backcast horizon"),fc=integer(forecast,"forecast horizon");
    const auto initial=requests(initial_filters,"Initial seasonal filters");
    const auto finals=requests(final_filters,"Final seasonal filters");
    const auto series=masked(values);
    const auto c=rebuild_c(cpass); // caller-supplied C pass; never recomputed here
    const auto result=kalos::x11_d_pass(series,period,first,initial,finals,trend,
        name=="legacy"?kalos::X11SigmaPolicy::legacy:kalos::X11SigmaPolicy::accurate,
        c,back,fc,[] {Rcpp::checkUserInterrupt();});
    Rcpp::IntegerVector resolved_initial(result.initial_filters.begin(),result.initial_filters.end());
    Rcpp::IntegerVector resolved_final(result.final_filters.begin(),result.final_filters.end());
    Rcpp::IntegerVector seas(result.seas_filter.begin(),result.seas_filter.end());
    // The 16 public pinned D tables in label order; D9_g_bis is the pinned
    // alias of D9bis in the fixed corrector mode (the null D9g is not
    // delivered, exactly like the core).
    Rcpp::List tables=Rcpp::List::create(
        Rcpp::_["D1"]=kalos_x11_rcpp::values(result.d1),Rcpp::_["D2"]=kalos_x11_rcpp::values(result.d2),
        Rcpp::_["D4"]=kalos_x11_rcpp::values(result.d4),Rcpp::_["D5"]=kalos_x11_rcpp::values(result.d5),
        Rcpp::_["D6"]=kalos_x11_rcpp::values(result.d6),Rcpp::_["D7"]=kalos_x11_rcpp::values(result.d7),
        Rcpp::_["D8"]=kalos_x11_rcpp::values(result.d8),Rcpp::_["D9"]=kalos_x11_rcpp::values(result.d9),
        Rcpp::_["D9bis"]=kalos_x11_rcpp::values(result.d9bis),Rcpp::_["D9_g_bis"]=kalos_x11_rcpp::values(result.d9bis),
        Rcpp::_["D10"]=kalos_x11_rcpp::values(result.d10),Rcpp::_["D10bis"]=kalos_x11_rcpp::values(result.d10bis),
        Rcpp::_["D11"]=kalos_x11_rcpp::values(result.d11),Rcpp::_["D11bis"]=kalos_x11_rcpp::values(result.d11bis),
        Rcpp::_["D12"]=kalos_x11_rcpp::values(result.d12),Rcpp::_["D13"]=kalos_x11_rcpp::values(result.d13));
    // The D9 MSRE diagnostic: the pinned per-phase table on the dropped
    // original calendar (phase indices 0..frequency-1 of the pinned table
    // at getPosition(backcast)) plus the count-weighted global ratio.
    Rcpp::IntegerVector counts(result.msr_table.counts.begin(),result.msr_table.counts.end());
    Rcpp::List msr_table=Rcpp::List::create(
        Rcpp::_["mean_irregular_evolutions"]=reals(result.msr_table.mean_irregular_evolutions),
        Rcpp::_["mean_seasonal_evolutions"]=reals(result.msr_table.mean_seasonal_evolutions),
        Rcpp::_["counts"]=counts,Rcpp::_["msr"]=reals(result.msr_table.msr),
        Rcpp::_["global_msr"]=std::isnan(result.msr_table.global_msr) ? NA_REAL
            : result.msr_table.global_msr);
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_d_pass",
        Rcpp::_["policy"]=name,Rcpp::_["frequency"]=period,Rcpp::_["start_period"]=first,
        Rcpp::_["trend_length"]=result.trend_length,
        Rcpp::_["d2drop"]=static_cast<int>(result.d2drop),
        Rcpp::_["backcast"]=result.backcast,Rcpp::_["forecast"]=result.forecast,
        Rcpp::_["d7_henderson_length"]=result.d7_henderson_length,
        Rcpp::_["final_henderson_length"]=result.final_henderson_length,
        Rcpp::_["msr_iterations"]=result.msr_iterations,
        Rcpp::_["initial_filters"]=resolved_initial,Rcpp::_["final_filters"]=resolved_final,
        Rcpp::_["seas_filter"]=seas,
        Rcpp::_["tables"]=tables,Rcpp::_["msr_table"]=msr_table,
        Rcpp::_["c_coupling"]="checked by the pinned core against the supplied C-pass result",
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(result.mode),
        Rcpp::_["mode_source"]="inherited from the C-pass result",
        Rcpp::_["status"]=result.status,
        Rcpp::_["scope"]=std::string("fixed_")+kalos_x11_rcpp::mode_name(result.mode)
            +"_d_pass_only",Rcpp::_["backend"]="native");
}
