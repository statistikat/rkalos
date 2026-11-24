#include <Rcpp.h>
#include "kalos/tramo_differencing_decision.hpp"

namespace {

const char* result_name(kalos::TramoDifferencingDecision::Result result) {
    using Result = kalos::TramoDifferencingDecision::Result;
    switch (result) {
    case Result::Changed: return "Changed";
    case Result::Unchanged: return "Unchanged";
    case Result::Unprocessed: return "Unprocessed";
    }
    Rcpp::stop("Unknown TRAMO differencing decision.");
}

} // namespace

// [[Rcpp::export(name = ".native_tramo_differencing", rng = false)]]
Rcpp::List native_tramo_differencing(std::vector<double> values, int period,
        Rcpp::IntegerVector orders, bool mean, bool log, bool seasonal, bool initial) {
    if (orders.size() != 6) Rcpp::stop("Expected six TRAMO SARIMA orders.");
    kalos::TramoDifferencingRequest request;
    request.series.swap(values);
    request.period = period;
    request.p0 = orders[0];
    request.d0 = orders[1];
    request.q0 = orders[2];
    request.bp0 = orders[3];
    request.bd0 = orders[4];
    request.bq0 = orders[5];
    request.mean0 = mean;
    request.log_input = log;
    request.seasonal = seasonal;
    request.initial_round = initial;
    const auto decision = kalos::decide_tramo_differencing(request);
    return Rcpp::List::create(
        Rcpp::_["decision"] = result_name(decision.result),
        Rcpp::_["orders"] = Rcpp::IntegerVector::create(
            Rcpp::_["p"] = decision.spec_p, Rcpp::_["d"] = decision.spec_d,
            Rcpp::_["q"] = decision.spec_q, Rcpp::_["bp"] = decision.spec_bp,
            Rcpp::_["bd"] = decision.spec_bd, Rcpp::_["bq"] = decision.spec_bq),
        Rcpp::_["mean"] = decision.mean, Rcpp::_["mean_correction"] = decision.mean_correction,
        Rcpp::_["tmean"] = decision.tmean, Rcpp::_["default_model"] = decision.default_model,
        Rcpp::_["selection_failed"] = decision.selection_failed,
        Rcpp::_["airline_applied"] = decision.airline_applied,
        Rcpp::_["events"] = decision.events);
}
