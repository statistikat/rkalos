#include "rcpp_tramo_fit.hpp"
#include <kalos/arma_diagnostics.hpp>

Rcpp::List wrap_monthly_airline_residuals(const std::vector<double>& values, double standard_error) {
    const auto lb = kalos::monthly_airline_residual_ljung_box(values);
    const auto normality = kalos::tramo_residual_normality(values);
    Rcpp::RObject dh = R_NilValue;
    if (normality) {
        dh = Rcpp::List::create(Rcpp::_["value"] = normality->value,
            Rcpp::_["pvalue"] = normality->p, Rcpp::_["df"] = normality->df);
    }
    return Rcpp::List::create(Rcpp::_["tsres"] = values, Rcpp::_["ser"] = standard_error,
        Rcpp::_["lb"] = Rcpp::List::create(Rcpp::_["value"] = lb.value, Rcpp::_["pvalue"] = lb.p,
            Rcpp::_["df"] = lb.df, Rcpp::_["lags"] = 24),
        Rcpp::_["doornikhansen"] = dh);
}


Rcpp::List wrap_tramo_preprocessing(const kalos::TramoFitResult& fit) {
    const auto& o = fit.orders;
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
        Rcpp::_["residuals"] = wrap_monthly_airline_residuals(fit.full_residuals,
            kalos::tramo_residual_standard_error(fit.ls_ssq, fit.ls_eff, fit.ls_npar)));
}

// [[Rcpp::export(name = ".native_tramo_tr1", rng = false)]]
Rcpp::List native_tramo_tr1(const std::vector<double>& series, int year, int month) {
    return wrap_tramo_preprocessing(kalos::estimate_tramo_tr1(series, year, month));
}
