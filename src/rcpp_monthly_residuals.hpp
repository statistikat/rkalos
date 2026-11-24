#pragma once

#include <Rcpp.h>
#include <kalos/arma_diagnostics.hpp>

inline Rcpp::List wrap_monthly_residuals(const std::vector<double>& values, double standard_error,
        const kalos::LjungBoxTest& lb, const std::optional<kalos::ResidualNormalityTest>& normality) {
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
