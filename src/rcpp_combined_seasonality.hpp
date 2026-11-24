#pragma once

#include <Rcpp.h>
#include <kalos/combined_seasonality.hpp>

inline Rcpp::List wrap_si_diagnostics(const std::optional<kalos::CombinedSiSeasonality>& test) {
    Rcpp::RObject combined = R_NilValue;
    if (test) {
        const auto statistic = [](double value, double p) {
            return Rcpp::List::create(Rcpp::_["value"] = value, Rcpp::_["pvalue"] = p);
        };
        combined = Rcpp::List::create(Rcpp::_["summary"] = test->summary,
            Rcpp::_["stable"] = statistic(test->stable_f, test->stable_p),
            Rcpp::_["moving"] = statistic(test->moving_f, test->moving_p),
            Rcpp::_["kruskal_wallis"] = statistic(test->kruskal_wallis, test->kruskal_p));
    }
    return Rcpp::List::create(Rcpp::_["combined_seasonality"] = combined);
}
