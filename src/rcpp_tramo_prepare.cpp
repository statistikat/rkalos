#include <Rcpp.h>
#include "kalos/tramo_model_builder.hpp"

// [[Rcpp::export(name = ".native_tramo_prepare_model", rng = false)]]
Rcpp::List native_tramo_prepare_model(std::vector<double> values, std::string name,
        int year, int month, int frequency) {
    const auto* spec = kalos::tramo_template(name);
    if (!spec) Rcpp::stop("Unknown TRAMO template: %s", name);
    kalos::tramo_build_request request;
    request.values.swap(values);
    request.start_year = year;
    request.start_month = month;
    request.annual_frequency = frequency;
    request.using_auto_model = spec->using_auto_model;
    request.transform_function = spec->transform_function;
    request.transform_adjust = spec->transform_adjust;
    request.arima = spec->arima;
    request.mean = spec->mean;
    request.trading_days = spec->trading_days;
    request.easter = spec->easter;
    request.delta_tc = spec->delta_tc;
    const auto result = kalos::build_tramo_model(request);
    Rcpp::IntegerVector missing(result.missing.size());
    for (std::size_t i = 0; i < result.missing.size(); ++i) missing[i] = result.missing[i] + 1;
    const auto& a = result.arima;
    Rcpp::IntegerVector orders = Rcpp::IntegerVector::create(
        Rcpp::_["p"] = static_cast<int>(a.phi.size()), Rcpp::_["d"] = a.d,
        Rcpp::_["q"] = static_cast<int>(a.theta.size()), Rcpp::_["bp"] = static_cast<int>(a.bphi.size()),
        Rcpp::_["bd"] = a.bd, Rcpp::_["bq"] = static_cast<int>(a.btheta.size()));
    Rcpp::CharacterVector variables(result.variables.size());
    for (std::size_t i = 0; i < result.variables.size(); ++i) variables[i] = result.variables[i].name;
    return Rcpp::List::create(Rcpp::_["template"] = spec->name,
        Rcpp::_["automatic_model"] = spec->using_auto_model, Rcpp::_["log"] = result.log_transformation,
        Rcpp::_["mean"] = result.mean, Rcpp::_["orders"] = orders, Rcpp::_["period"] = a.period,
        Rcpp::_["missing"] = missing, Rcpp::_["interpolated"] = result.interpolated,
        Rcpp::_["variables"] = variables, Rcpp::_["warnings"] = result.warnings);
}
