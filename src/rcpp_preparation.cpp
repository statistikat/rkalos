#include <Rcpp.h>
#include "kalos/preparation.hpp"
#include "kalos/version.hpp"
#include "kalos/detail/lapack.hpp"
#include <cmath>
#include <limits>

namespace {
std::int64_t integer(double x) {
    if (!std::isfinite(x) || x != std::floor(x) || std::abs(x) > 9007199254740991.0)
        Rcpp::stop("Domain values must be exact integers.");
    return static_cast<std::int64_t>(x);
}
kalos::TimeDomain domain(Rcpp::List x) {
    Rcpp::NumericVector start=x["start"];
    if(start.size()!=2 || !std::isfinite(start[1]) || start[1]<1 || start[1]>12)Rcpp::stop("Invalid domain start.");
    return {integer(start[0]), static_cast<int>(integer(start[1])), Rcpp::as<int>(x["frequency"]),
            integer(Rcpp::as<double>(x["length"]))};
}
Rcpp::NumericVector period(const std::pair<std::int64_t,int>& x) {
    return Rcpp::NumericVector::create(static_cast<double>(x.first),x.second);
}
Rcpp::List output_domain(const kalos::TimeDomain& x) {
    // Preserve the exact-double gate used by the original R ts conversion.
    const auto year=static_cast<double>(x.at(0).first);
    const auto n=static_cast<double>(x.length());
    if(std::abs(year)*x.frequency()+std::abs(n)+x.frequency()>9007199254740991.0)
        Rcpp::stop("time-domain ordinals exceed R's exact integer range");
    return Rcpp::List::create(Rcpp::_["start"]=period(x.at(0)),
        Rcpp::_["end"]=period(x.at(x.length()-1)),Rcpp::_["next"]=period(x.at(x.length())),
        Rcpp::_["frequency"]=x.frequency(),Rcpp::_["length"]=n,
        Rcpp::_["start_ordinal"]=static_cast<double>(x.start_ordinal()));
}
kalos::Interrupt interrupt() {return []{Rcpp::checkUserInterrupt();};}
Rcpp::IntegerVector indices(const std::vector<std::size_t>& x) {
    Rcpp::IntegerVector result(x.size());
    for(std::size_t i=0;i<x.size();++i) {
        if(x[i]>=static_cast<std::size_t>(std::numeric_limits<int>::max()))Rcpp::stop("Row exceeds R integer range.");
        result[i]=static_cast<int>(x[i]+1);
    }
    return result;
}
}

// [[Rcpp::export(name = ".native_estimation_window")]]
Rcpp::List native_estimation_window(Rcpp::List source, Rcpp::NumericVector endpoints) {
    if(endpoints.size()!=2)Rcpp::stop("Expected two window endpoints.");
    const auto result=kalos::estimation_window(domain(source),integer(endpoints[0]),integer(endpoints[1]));
    return Rcpp::List::create(Rcpp::_["first"]=static_cast<int>(result.first+1),
        Rcpp::_["length"]=static_cast<int>(result.domain.length()),Rcpp::_["domain"]=output_domain(result.domain));
}

// [[Rcpp::export(name = ".native_align_regressors")]]
Rcpp::List native_align_regressors(Rcpp::List observed, Rcpp::List columns,
    Rcpp::List domains, Rcpp::CharacterVector labels, double backcasts, double forecasts) {
    if(columns.size()!=domains.size() || columns.size()!=labels.size())Rcpp::stop("Inconsistent regressor metadata.");
    std::vector<kalos::NamedRegressor> regressors;
    for(R_xlen_t i=0;i<columns.size();++i) {
        if(labels[i]==NA_STRING)Rcpp::stop("Regressor names must not be missing.");
        regressors.push_back({Rcpp::as<std::string>(labels[i]),domain(Rcpp::as<Rcpp::List>(domains[i])),
            Rcpp::as<std::vector<double>>(columns[i])});
    }
    const auto x=kalos::align_regressors(domain(observed),regressors,integer(backcasts),integer(forecasts),interrupt());
    Rcpp::NumericMatrix values(static_cast<int>(x.required.length()),static_cast<int>(x.names.size()));
    std::copy(x.values.begin(),x.values.end(),values.begin());
    if(labels.size())Rcpp::colnames(values)=labels;
    return Rcpp::List::create(Rcpp::_["required_domain"]=output_domain(x.required),Rcpp::_["regressors"]=values);
}

// [[Rcpp::export(name = ".native_prepare_regression")]]
Rcpp::List native_prepare_regression(Rcpp::NumericMatrix values, Rcpp::CharacterVector labels,
    Rcpp::List observed, Rcpp::List required, int backcasts, Rcpp::NumericVector coefficients,
    Rcpp::CharacterVector coefficient_names, int d, int seasonal_d, double mean,
    bool include_mean, bool include_drift, int source_first) {
    if(values.ncol()!=labels.size() || coefficients.size()!=coefficient_names.size() || backcasts<0 || source_first<1)
        Rcpp::stop("Inconsistent preparation metadata.");
    std::vector<kalos::NamedCoefficient> fixed;
    for(R_xlen_t i=0;i<coefficients.size();++i) {
        if(coefficient_names[i]==NA_STRING || R_IsNaN(coefficients[i]))Rcpp::stop("Invalid coefficients.");
        const bool is_fixed=!R_IsNA(coefficients[i]);
        fixed.push_back({Rcpp::as<std::string>(coefficient_names[i]),is_fixed?coefficients[i]:0,is_fixed});
    }
    kalos::PreparedRegressors prepared{domain(observed),domain(required),static_cast<std::size_t>(backcasts),
        Rcpp::as<std::vector<std::string>>(labels),Rcpp::as<std::vector<double>>(values)};
    const auto x=kalos::prepare_regression(std::move(prepared),fixed,d,seasonal_d,mean,include_mean,include_drift,
        static_cast<std::size_t>(source_first-1),interrupt());
    const auto names=Rcpp::wrap(x.regressors.names);
    Rcpp::NumericMatrix all(static_cast<int>(x.regressors.required.length()),static_cast<int>(x.design.columns));
    Rcpp::NumericMatrix design(static_cast<int>(x.regressors.observed.length()),static_cast<int>(x.design.columns));
    std::copy(x.regressors.values.begin(),x.regressors.values.end(),all.begin());
    std::copy(x.design.values.begin(),x.design.values.end(),design.begin());
    if(x.design.columns){Rcpp::colnames(all)=names;Rcpp::colnames(design)=names;}
    Rcpp::NumericVector states(x.design.columns);
    for(std::size_t i=0;i<x.design.columns;++i)states[i]=x.design.fixed[i]?x.design.coefficients[i]:NA_REAL;
    states.attr("names")=names;
    return Rcpp::List::create(Rcpp::_["regressors"]=all,Rcpp::_["design"]=design,
        Rcpp::_["labels"]=names,Rcpp::_["fixed"]=states,Rcpp::_["generated"]=x.generated);
}

// [[Rcpp::export(name = ".native_likelihood_adjustment")]]
Rcpp::List native_likelihood_adjustment(Rcpp::NumericVector values, bool logarithm,
    int degree, int source_first) {
    if(degree<0 || source_first<1)Rcpp::stop("Invalid likelihood adjustment range.");
    kalos::MaskedSeries input{std::vector<double>(values.size()),std::vector<std::uint8_t>(values.size())};
    for(R_xlen_t i=0;i<values.size();++i) {
        if((i&4095)==0)Rcpp::checkUserInterrupt();
        if(std::isnan(values[i]))input.missing[i]=1;else input.values[i]=values[i];
    }
    const auto x=kalos::likelihood_adjustment(input,logarithm?kalos::Transform::logarithm:kalos::Transform::none,
        static_cast<std::size_t>(degree),static_cast<std::size_t>(source_first-1),interrupt());
    return Rcpp::List::create(Rcpp::_["local_rows"]=indices(x.local_rows),Rcpp::_["source_rows"]=indices(x.source_rows),
        Rcpp::_["terms"]=x.terms,Rcpp::_["log_jacobian"]=x.log_jacobian);
}

// [[Rcpp::export(name = ".native_core_version")]]
Rcpp::List native_core_version() {
    // Build state of the vendored SLICOT SB03MD seam, mirroring the core's
    // KALOS_WITH_SLICOT option (configure's Fortran probe): true when the
    // closure is compiled into this shared object (stationary initialization
    // uses SB03MD for state order >= 8), false in the QR-everywhere build.
#ifdef KALOS_WITH_SLICOT
    const bool slicot_stationary_initialization = true;
#else
    const bool slicot_stationary_initialization = false;
#endif
    return Rcpp::List::create(Rcpp::_["version"]=kalos::version,Rcpp::_["commit"]=kalos::core_revision,
        Rcpp::_["provider"]=kalos::lapack::provider_name(),Rcpp::_["provider_version"]=kalos::lapack::provider_version(),
        Rcpp::_["blas_integer_bits"]=32L,Rcpp::_["cpp_standard"]=17L,
        Rcpp::_["slicot_stationary_initialization"]=slicot_stationary_initialization);
}
