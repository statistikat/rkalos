#include <Rcpp.h>
#include <kalos/combined_seasonality.hpp>

namespace {
Rcpp::RObject adjusted_test_report(const std::optional<kalos::AdjustedSeasonalityTest>& test) {
    if (!test) return R_NilValue;
    return Rcpp::List::create(Rcpp::_["value"] = test->value,
                              Rcpp::_["p_value"] = test->p_value);
}
}

// [[Rcpp::export(name = ".native_adjusted_seasonality")]]
Rcpp::List native_adjusted_seasonality(const std::vector<double>& values,
        int differencing_order, bool mean) {
    const auto test = kalos::monthly_adjusted_seasonality(values, differencing_order, mean);
    return Rcpp::List::create(Rcpp::_["qs"] = adjusted_test_report(test.qs),
                              Rcpp::_["seasonal_f"] = adjusted_test_report(test.seasonal_f));
}
