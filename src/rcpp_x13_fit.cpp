#include "rcpp_x13_shared.hpp"
#include "rcpp_combined_seasonality.hpp"
#include "rcpp_x11_selection_records.hpp"
#include "rcpp_x11_table_support.hpp"
#include "kalos/x13_fit.hpp"
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
// The REQUIRED sigma policy arrives as NULL (unset — the fit then answers the
// pinned fail-closed status, never a default) or the policy name string.
std::optional<kalos::X11SigmaPolicy> sigma_policy(SEXP value) {
    if(Rf_isNull(value))return std::nullopt;
    if(TYPEOF(value)!=STRSXP || Rf_xlength(value)!=1 || STRING_ELT(value,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument("x13 fit sigma weight policy must be legacy or accurate");
    return kalos_x11_rcpp::parse_policy(Rcpp::as<std::string>(value),
        "x13 fit sigma weight policy");
}
}
// [[Rcpp::export(name = ".native_x13_fit")]]
Rcpp::List native_x13_fit(SEXP values,SEXP variables_list,SEXP log,SEXP frequency,
    SEXP start_period,SEXP seasonal,SEXP trend_length,SEXP sigma,SEXP lower_sigma,
    SEXP upper_sigma,SEXP msr_policy,SEXP mode,SEXP backcast,SEXP forecast) {
    kalos::X13FitRequest request;
    request.input=kalos_x13_rcpp::reals_in(values,"x13 fit input");
    request.variables=kalos_x13_rcpp::variables(variables_list);
    if(TYPEOF(log)!=LGLSXP || Rf_xlength(log)!=1 || LOGICAL(log)[0]==NA_LOGICAL ||
        !Rf_isNull(Rf_getAttrib(log,R_DimSymbol)))
        throw std::invalid_argument("x13 fit log flag must be TRUE or FALSE");
    request.log=LOGICAL(log)[0]==TRUE;
    // Scalars are validated by the pinned composition and its seams (the
    // frequency/seasonal pairing, the Henderson parity, the horizons); the
    // boundary only marshals. Rcpp converts the integer/real vectors itself.
    request.frequency=Rcpp::as<int>(frequency);
    request.first_period=Rcpp::as<int>(start_period);
    request.seasonal=kalos_x11_rcpp::requests(seasonal);
    request.trend_length=Rcpp::as<int>(trend_length);
    request.sigma=sigma_policy(sigma);
    request.lower_sigma=Rcpp::as<double>(lower_sigma);
    request.upper_sigma=Rcpp::as<double>(upper_sigma);
    if(TYPEOF(msr_policy)!=STRSXP || Rf_xlength(msr_policy)!=1 ||
        STRING_ELT(msr_policy,0)==NA_STRING || !Rf_isNull(Rf_getAttrib(msr_policy,R_DimSymbol)))
        throw std::invalid_argument("MSR iteration policy must be jdplus or x11_five_year_cap");
    const auto msr_name=Rcpp::as<std::string>(msr_policy);
    request.msr_policy=kalos_x11_rcpp::parse_msr_policy(msr_name);
    request.mode=kalos_x11_rcpp::parse_mode(mode,"x13 fit decomposition mode");
    request.backcast=Rcpp::as<int>(backcast);
    request.forecast=Rcpp::as<int>(forecast);
    const auto result=kalos::x13_fit(request);
    const Rcpp::Nullable<Rcpp::CharacterVector> sigma_field=result.sigma_set
        ? Rcpp::Nullable<Rcpp::CharacterVector>(
            Rcpp::wrap(std::string(kalos_x11_rcpp::policy_name(result.sigma))))
        : Rcpp::Nullable<Rcpp::CharacterVector>(R_NilValue);
    Rcpp::List out=Rcpp::List::create(Rcpp::_["schema_version"]=1,
        Rcpp::_["type"]="x13_fit",Rcpp::_["status"]=result.status,
        // The provenance echo is complete whether or not the run completed,
        // exactly as the pinned contract requires for a refused request.
        Rcpp::_["frequency"]=result.frequency,Rcpp::_["start_period"]=result.first_period,
        Rcpp::_["log"]=result.log,Rcpp::_["seasonal"]=Rcpp::CharacterVector(seasonal),
        Rcpp::_["trend_length"]=result.trend_length,
        Rcpp::_["sigma_set"]=result.sigma_set,Rcpp::_["sigma"]=sigma_field,
        Rcpp::_["lower_sigma"]=result.lower_sigma,Rcpp::_["upper_sigma"]=result.upper_sigma,
        Rcpp::_["msr_policy"]=msr_name,
        Rcpp::_["mode_requested"]=kalos_x11_rcpp::mode_name(result.mode_requested),
        Rcpp::_["mode_effective"]=kalos_x11_rcpp::mode_name(result.mode_effective),
        Rcpp::_["mode_forcing"]="pinned X13Kernel.updateSpec: preserve an explicit pseudo_additive mode; "
            "otherwise non-log -> additive and log -> multiplicative unless log_additive was requested",
        Rcpp::_["backcast"]=result.backcast,Rcpp::_["forecast"]=result.forecast,
        Rcpp::_["scope"]="restricted_fixed_model_x13_fit_only",Rcpp::_["backend"]="native");
    if(!result.preadjustment.status.empty()) // past the fail-closed guards
        out["preadjustment"]=kalos_x13_rcpp::preadjustment_list(result.preadjustment,result.log);
    if(!result.selection.status.empty()) { // the run over out, same record shape
        out["selection"]=kalos_x11_rcpp::selection_result(result.selection,
            result.sigma_set ? kalos_x11_rcpp::policy_name(result.sigma) : "",
            msr_name,kalos_x11_rcpp::mode_name(result.mode_effective),
            "restricted_fit_selection_run",result.lower_sigma,result.upper_sigma,
            result.frequency,result.first_period,result.trend_length,0,0,
            Rcpp::CharacterVector(seasonal));
    }
    if(result.status=="available") { // tables only on a completed run
        out["diagnostics"] = wrap_si_diagnostics(kalos::x13_si_seasonality(result));
        out["results"]=Rcpp::List::create(
            Rcpp::_["d10"]=kalos_x13_rcpp::reals_out(result.d10),
            Rcpp::_["d11"]=kalos_x13_rcpp::reals_out(result.d11),
            Rcpp::_["d12"]=kalos_x13_rcpp::reals_out(result.d12),
            Rcpp::_["d13"]=kalos_x13_rcpp::reals_out(result.d13),
            Rcpp::_["c17"]=kalos_x13_rcpp::reals_out(result.c17),
            Rcpp::_["note"]=result.mode_effective==kalos::X11DecompositionMode::log_additive
                ? "the X11Results tables as finals consumed them (log-additive "
                    "results-builder stage, Ratio bias default)"
                : "the X11Results tables as finals consumed them (raw step tables)");
        out["finals"]=Rcpp::List::create(Rcpp::_["schema_version"]=1,
            Rcpp::_["type"]="x13_finals",Rcpp::_["scope"]="restricted_fit_finals",
            Rcpp::_["tables"]=kalos_x13_rcpp::finals_list(result.finals));
        // The pinned final-components triplet (X13Finals exports): there is
        // no separate pinned final-trend/final-irregular table; these are
        // bit-identical to finals.d11final/d12final/d13final.
        out["finalSA"]=kalos_x13_rcpp::reals_out(result.finalSA);
        out["finalTrend"]=kalos_x13_rcpp::reals_out(result.finalTrend);
        out["finalIrregular"]=kalos_x13_rcpp::reals_out(result.finalIrregular);
    }
    return out;
}
