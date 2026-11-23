#include <Rcpp.h>
#include <kalos/tramo_transform_decision.hpp>

namespace {
Rcpp::List test_report(const kalos::TramoSeasonalityTest& test) {
    return Rcpp::List::create(Rcpp::_["value"] = test.value,
                              Rcpp::_["p_value"] = test.p_value,
                              Rcpp::_["valid"] = test.valid);
}

Rcpp::List seasonality_report(const kalos::TramoSeasonalityReport& report) {
    return Rcpp::List::create(Rcpp::_["verdict"] = report.verdict,
                              Rcpp::_["to_int"] = report.to_int,
                              Rcpp::_["qs"] = test_report(report.qs),
                              Rcpp::_["friedman"] = test_report(report.friedman),
                              Rcpp::_["differenced"] = report.differenced);
}
}

// [[Rcpp::export(name = ".native_tramo_critical_value")]]
double native_tramo_critical_value(int observations) {
    if (observations < 0) Rcpp::stop("observations must be nonnegative");
    return kalos::tramo_critical_value(static_cast<std::size_t>(observations));
}

// [[Rcpp::export(name = ".native_tramo_seasonality")]]
Rcpp::List native_tramo_seasonality(const std::vector<double>& values) {
    return seasonality_report(kalos::tramo_seasonality(values));
}

// [[Rcpp::export(name = ".native_tramo_transform_decision")]]
Rcpp::List native_tramo_transform_decision(std::vector<double> values, bool mean,
        bool detect_seasonality, bool spec_seasonality, double log_preference,
        double precision) {
    kalos::TramoTransformDecisionRequest request;
    request.values.swap(values);
    request.period = 12;
    request.mean = mean;
    request.detect_seasonality = detect_seasonality;
    request.spec_seasonality = spec_seasonality;
    request.log_preference = log_preference;
    request.precision = precision;
    const auto decision = kalos::decide_tramo_transform(request);
    Rcpp::RObject detector = R_NilValue;
    if (detect_seasonality) detector = seasonality_report(decision.seasonality);
    return Rcpp::List::create(
        Rcpp::_["seasonality"] = detector, Rcpp::_["spec_pre"] = decision.spec_pre,
        Rcpp::_["seasonal"] = decision.seasonal, Rcpp::_["nonpositive"] = decision.nonpositive,
        Rcpp::_["loglevel_status"] = decision.loglevel_status,
        Rcpp::_["choosing_log"] = decision.choosing_log,
        Rcpp::_["log_ll"] = decision.log_ll, Rcpp::_["level_ll"] = decision.level_ll,
        Rcpp::_["log_correction"] = decision.log_correction,
        Rcpp::_["log_preference"] = decision.log_preference,
        Rcpp::_["log_transform"] = decision.log_transform, Rcpp::_["remarks"] = decision.remarks);
}
