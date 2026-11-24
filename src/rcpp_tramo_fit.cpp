#include "rcpp_tramo_fit.hpp"
#include "rcpp_monthly_residuals.hpp"

Rcpp::List wrap_monthly_tramo_residuals(const std::vector<double>& values, double standard_error, int arma_parameters) {
    const auto lb = kalos::monthly_tramo_residual_ljung_box(values, arma_parameters);
    const auto normality = kalos::tramo_residual_normality(values);
    return wrap_monthly_residuals(values, standard_error, lb, normality);
}


Rcpp::List wrap_tramo_preprocessing(const kalos::TramoFitResult& fit) {
    const auto& o = fit.orders;
    auto residuals = wrap_monthly_tramo_residuals(fit.full_residuals,
        kalos::tramo_residual_standard_error(fit.ls_ssq, fit.ls_eff, fit.ls_npar),
        o.p + o.q + o.bp + o.bq);
    const auto n = fit.variables.size();
    Rcpp::CharacterVector names(n), types(n);
    Rcpp::IntegerVector positions(n);
    Rcpp::NumericVector estimates(n);
    const char* outlier_types[] = {"AO", "LS", "TC", "SO"};
    for (std::size_t i = 0; i < n; ++i) {
        const auto& variable = fit.variables[i];
        names[i] = variable.name;
        types[i] = outlier_types[static_cast<int>(variable.type)];
        positions[i] = variable.position + 1;
        estimates[i] = variable.coefficient;
    }
    return Rcpp::List::create(
        Rcpp::_["status"] = fit.status,
        Rcpp::_["model"] = Rcpp::List::create(
            Rcpp::_["orders"] = Rcpp::IntegerVector::create(Rcpp::_["p"] = o.p, Rcpp::_["d"] = o.d,
                Rcpp::_["q"] = o.q, Rcpp::_["bp"] = o.bp, Rcpp::_["bd"] = o.bd, Rcpp::_["bq"] = o.bq),
            Rcpp::_["parameters"] = fit.parameters, Rcpp::_["log"] = fit.log,
            Rcpp::_["mean_correction"] = fit.mean,
            Rcpp::_["mean_coefficient"] = fit.mean ? fit.regression_coefficients.at(0) : NA_REAL,
            Rcpp::_["outliers"] = Rcpp::DataFrame::create(Rcpp::_["name"] = names, Rcpp::_["type"] = types,
                Rcpp::_["position"] = positions, Rcpp::_["coefficient"] = estimates)),
        Rcpp::_["likelihood"] = Rcpp::List::create(Rcpp::_["log_likelihood"] = fit.ls_ll,
            Rcpp::_["transformation_adjustment"] = fit.ls_tadj, Rcpp::_["adjusted_log_likelihood"] = fit.ls_all,
            Rcpp::_["ssq"] = fit.ls_ssq, Rcpp::_["effective"] = fit.ls_eff, Rcpp::_["estimated"] = fit.ls_npar),
        Rcpp::_["residuals"] = residuals);
}

// [[Rcpp::export(name = ".native_tramo_tr1", rng = false)]]
Rcpp::List native_tramo_tr1(const std::vector<double>& series, int year, int month) {
    return wrap_tramo_preprocessing(kalos::estimate_tramo_tr1(series, year, month));
}

// [[Rcpp::export(name = ".native_tramo_tr3", rng = false)]]
Rcpp::List native_tramo_tr3(const std::vector<double>& series, int year, int month) {
    return wrap_tramo_preprocessing(kalos::estimate_tramo_tr3(series, year, month));
}

void attach_tramo_calendar(Rcpp::List& result, const kalos::TramoFittedCalendar& calendar) {
    const auto& selection = calendar.selection;
    std::string days = "none";
    std::vector<std::string> labels;
    if (selection.trading_days == kalos::TramoTradingDays::WorkingDays) {
        days = "working_days";
        labels.push_back("td");
    } else if (selection.trading_days == kalos::TramoTradingDays::TradingDays) {
        days = "trading_days";
        for (int i = 0; i < 6; ++i) labels.push_back("td[" + std::to_string(i) + "]");
    }
    if (selection.leap_year) labels.push_back("lp");
    if (selection.easter) labels.push_back("easter");
    Rcpp::NumericVector coefficients = Rcpp::wrap(calendar.coefficients);
    coefficients.attr("names") = labels;
    Rcpp::List model = result["model"];
    model["calendar"] = Rcpp::List::create(Rcpp::_["trading_days"] = days,
        Rcpp::_["leap_year"] = selection.leap_year, Rcpp::_["easter"] = selection.easter,
        Rcpp::_["coefficients"] = coefficients);
    result["model"] = model;
}

// [[Rcpp::export(name = ".native_tramo_trfull", rng = false)]]
Rcpp::List native_tramo_trfull(const std::vector<double>& series, int year, int month) {
    const auto fit = kalos::estimate_tramo_trfull(series, year, month);
    auto result = wrap_tramo_preprocessing(fit);
    attach_tramo_calendar(result, fit.calendar);
    return result;
}
