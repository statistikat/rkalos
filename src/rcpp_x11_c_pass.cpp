#include "rcpp_x11_table_support.hpp"
#include "kalos/x11_c_pass.hpp"
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
        throw std::invalid_argument("C-pass input must be a real numeric vector within the supported length");
    const auto y=Rcpp::as<Rcpp::NumericVector>(values);const auto n=static_cast<std::size_t>(y.size());
    kalos::MaskedSeries out{std::vector<double>(n,0),std::vector<std::uint8_t>(n,0)};
    for(std::size_t i=0;i<n;++i) {
        if((i&4095U)==0)Rcpp::checkUserInterrupt();
        const auto value=y[static_cast<R_xlen_t>(i)];
        if(std::isnan(value))out.missing[i]=1;
        else if(!std::isfinite(value))throw std::invalid_argument("C-pass input must be finite or missing");
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
kalos::X11SigmaPolicy parse_policy(const std::string& name,const char* label) {
    if(name=="legacy")return kalos::X11SigmaPolicy::legacy;
    if(name=="accurate")return kalos::X11SigmaPolicy::accurate;
    throw std::invalid_argument(std::string(label)+" must be legacy or accurate");
}
// Rebuild the pinned MaskedSeries from a .native_x11_b_pass table (NA marks
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
// The C pass consumes the caller's finished B-pass result; the pinned core
// enforces the coupling (status, setting echoes, round policies, series
// equality) on the reconstructed object, so the boundary only rebuilds it
// faithfully and never recomputes the B pass here.
kalos::X11BPass rebuild_b(SEXP value) {
    if(TYPEOF(value)!=VECSXP)
        throw std::invalid_argument("C pass requires a B-pass result object");
    const char* shape="C pass requires a B-pass result object";
    if(Rcpp::as<std::string>(named(value,"type",shape))!="x11_b_pass")
        throw std::invalid_argument(shape);
    const SEXP tables=named(value,"tables",shape);
    const std::array<const char*,19> names={"B1","B2","B3","B4","B4a","B4anorm","B4d","B4g",
        "B5","B6","B7","B8","B9","B9g","B10","B11","B13","B17","B20"};
    using Member=kalos::MaskedSeries kalos::X11BPass::*;
    static const Member members[19]=
        {&kalos::X11BPass::b1,&kalos::X11BPass::b2,&kalos::X11BPass::b3,&kalos::X11BPass::b4,
        &kalos::X11BPass::b4a,&kalos::X11BPass::b4anorm,&kalos::X11BPass::b4d,&kalos::X11BPass::b4g,
        &kalos::X11BPass::b5,&kalos::X11BPass::b6,&kalos::X11BPass::b7,&kalos::X11BPass::b8,
        &kalos::X11BPass::b9,&kalos::X11BPass::b9g,&kalos::X11BPass::b10,&kalos::X11BPass::b11,
        &kalos::X11BPass::b13,&kalos::X11BPass::b17,&kalos::X11BPass::b20};
    kalos::X11BPass b;
    b.b2drop=static_cast<std::size_t>(integer(named(value,"b2drop",shape),"B-pass b2drop"));
    b.first_period=integer(named(value,"start_period",shape),"B-pass start period");
    const SEXP status=named(value,"status",shape);
    if(TYPEOF(status)!=STRSXP || Rf_xlength(status)!=1 || STRING_ELT(status,0)==NA_STRING)
        throw std::invalid_argument(shape);
    b.status=Rcpp::as<std::string>(status);
    // The decomposition mode is carried by the B-pass result: the pinned
    // kernel runs the C step on the same context, so the C mode is exactly
    // the mode the caller's B pass reported (never restated here; an object
    // predating the mode echo is additive by construction).
    b.mode=kalos_x11_rcpp::parse_mode(named_opt(value,"mode"),"B-pass decomposition mode");
    b.trend_length=integer(named(value,"trend_length",shape),"B-pass trend filter length");
    b.initial_filters=requests(named(value,"initial_filters",shape),"B-pass initial filters");
    b.final_filters=requests(named(value,"final_filters",shape),"B-pass final filters");
    for(std::size_t k=0;k<names.size();++k)
        b.*(members[k])=untable(named(tables,names[k],shape),"B-pass table");
    const SEXP rounds=named(value,"rounds",shape);
    if(TYPEOF(rounds)!=VECSXP || Rf_xlength(rounds)!=3)
        throw std::invalid_argument("the B-pass object must carry three correction rounds");
    for(std::size_t k=0;k<3;++k) {
        const SEXP round=VECTOR_ELT(rounds,k);
        const char* rshape="the B-pass correction rounds are incomplete";
        if(TYPEOF(round)!=VECSXP)throw std::invalid_argument(rshape);
        auto& target=b.rounds[k];
        target.detection=untable(named(round,"detection",rshape),"B-pass round detection");
        target.analysis_position=static_cast<std::int64_t>(
            integer(named(round,"analysis_position",rshape),"B-pass round analysis position"))-1; // back to 0-based
        target.policy=parse_policy(Rcpp::as<std::string>(named(round,"policy",rshape)),
            "B-pass round policy");
        const auto weights=Rcpp::as<Rcpp::NumericVector>(named(round,"sigma_weights",rshape));
        target.sigma.weights.assign(weights.begin(),weights.end());
        const auto vals=Rcpp::as<Rcpp::NumericVector>(named(round,"replacements",rshape));
        target.replacements.replacements.assign(vals.begin(),vals.end());
        const auto factors=Rcpp::as<Rcpp::NumericVector>(named(round,"corrections",rshape));
        target.replacements.correction_factors.assign(factors.begin(),factors.end());
    }
    return b;
}
Rcpp::List correction_round(const kalos::X11CCorrectionRound& round) {
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
// [[Rcpp::export(name = ".native_x11_c_pass")]]
Rcpp::List native_x11_c_pass(SEXP values,SEXP frequency,SEXP start_period,SEXP policy,
    SEXP initial_filters,SEXP final_filters,SEXP trend_length,SEXP bpass,
    SEXP lower_sigma,SEXP upper_sigma) {
    if(TYPEOF(policy)!=STRSXP || Rf_xlength(policy)!=1 || STRING_ELT(policy,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(policy,R_DimSymbol)))
        throw std::invalid_argument("C-pass sigma weight policy must be legacy or accurate");
    const auto name=Rcpp::as<std::string>(policy);
    if(name!="legacy" && name!="accurate")
        throw std::invalid_argument("C-pass sigma weight policy must be legacy or accurate");
    const auto period=integer(frequency,"frequency"),first=integer(start_period,"start period");
    const auto trend=integer(trend_length,"trend filter length");
    const auto lower=scalar(lower_sigma,"lower sigma limit"),upper=scalar(upper_sigma,"upper sigma limit");
    const auto initial=requests(initial_filters,"Initial seasonal filters");
    const auto finals=requests(final_filters,"Final seasonal filters");
    const auto series=masked(values);
    const auto b=rebuild_b(bpass); // caller-supplied B pass; never recomputed here
    const auto result=kalos::x11_c_pass(series,period,first,initial,finals,trend,
        name=="legacy"?kalos::X11SigmaPolicy::legacy:kalos::X11SigmaPolicy::accurate,
        b,lower,upper,[] {Rcpp::checkUserInterrupt();});
    Rcpp::IntegerVector resolved_initial(result.initial_filters.begin(),result.initial_filters.end());
    Rcpp::IntegerVector resolved_final(result.final_filters.begin(),result.final_filters.end());
    // The 14 pinned C tables in label order.
    Rcpp::List tables=Rcpp::List::create(
        Rcpp::_["C1"]=kalos_x11_rcpp::values(result.c1),Rcpp::_["C2"]=kalos_x11_rcpp::values(result.c2),
        Rcpp::_["C4"]=kalos_x11_rcpp::values(result.c4),Rcpp::_["C5"]=kalos_x11_rcpp::values(result.c5),
        Rcpp::_["C5a"]=kalos_x11_rcpp::values(result.c5a),Rcpp::_["C6"]=kalos_x11_rcpp::values(result.c6),
        Rcpp::_["C7"]=kalos_x11_rcpp::values(result.c7),Rcpp::_["C9"]=kalos_x11_rcpp::values(result.c9),
        Rcpp::_["C10"]=kalos_x11_rcpp::values(result.c10),Rcpp::_["C10a"]=kalos_x11_rcpp::values(result.c10a),
        Rcpp::_["C11"]=kalos_x11_rcpp::values(result.c11),Rcpp::_["C13"]=kalos_x11_rcpp::values(result.c13),
        Rcpp::_["C17"]=kalos_x11_rcpp::values(result.c17),Rcpp::_["C20"]=kalos_x11_rcpp::values(result.c20));
    return Rcpp::List::create(Rcpp::_["schema_version"]=1,Rcpp::_["type"]="x11_c_pass",
        Rcpp::_["policy"]=name,Rcpp::_["frequency"]=period,Rcpp::_["start_period"]=first,
        Rcpp::_["trend_length"]=result.trend_length,Rcpp::_["c2drop"]=static_cast<int>(result.c2drop),
        Rcpp::_["lower_sigma"]=lower,Rcpp::_["upper_sigma"]=upper,
        Rcpp::_["initial_filters"]=resolved_initial,Rcpp::_["final_filters"]=resolved_final,
        Rcpp::_["tables"]=tables,Rcpp::_["round"]=correction_round(result.round),
        // The pinned reference-series echo (S3): the exact series this C pass
        // consumed, which the D pass coupling is checked against.
        Rcpp::_["reference_series"]=kalos_x11_rcpp::values(result.reference_series),
        Rcpp::_["mode"]=kalos_x11_rcpp::mode_name(result.mode),
        Rcpp::_["mode_source"]="inherited from the B-pass result",
        Rcpp::_["status"]=result.status,
        Rcpp::_["scope"]=std::string("fixed_")+kalos_x11_rcpp::mode_name(result.mode)
            +"_c_pass_only",Rcpp::_["backend"]="native");
}
