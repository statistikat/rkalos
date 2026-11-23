#include <Rcpp.h>
#include <kalos/seats_extract_wk.hpp>

namespace {
Rcpp::DataFrame component_table(const std::vector<kalos::seats_component_series>& rows) {
    const char* components[] = {"Series", "Trend", "Seasonal", "SeasonallyAdjusted", "Irregular", "CalendarEffect"};
    const char* information[] = {"Value", "Stdev", "Forecast", "StdevForecast", "Backcast", "StdevBackcast"};
    Rcpp::CharacterVector type(rows.size()), info(rows.size());
    Rcpp::NumericVector year(rows.size());
    Rcpp::IntegerVector month(rows.size());
    Rcpp::List values(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        type[i] = components[static_cast<int>(row.component)];
        info[i] = information[static_cast<int>(row.information)];
        year[i] = static_cast<double>(row.start_year);
        month[i] = row.start_month;
        values[i] = row.values;
    }
    values.attr("class") = "AsIs";
    return Rcpp::DataFrame::create(Rcpp::_["component"] = type, Rcpp::_["information"] = info,
        Rcpp::_["start_year"] = year, Rcpp::_["start_month"] = month, Rcpp::_["values"] = values);
}
}

// [[Rcpp::export(name = ".native_seats_extract_wk", rng = false)]]
Rcpp::List native_seats_extract_wk(const Rcpp::List& spec, const Rcpp::List& calendar,
        std::vector<double> series) {
    kalos::tramo_seats_model_request r;
    r.series.swap(series);
    const Rcpp::NumericVector start = calendar["start"];
    r.frequency = r.period = Rcpp::as<int>(calendar["frequency"]);
    r.start_year = start[0];
    r.start_month = start[1];
    r.phi = Rcpp::as<std::vector<double>>(spec["phi"]);
    r.theta = Rcpp::as<std::vector<double>>(spec["theta"]);
    r.bphi = Rcpp::as<std::vector<double>>(spec["bphi"]);
    r.btheta = Rcpp::as<std::vector<double>>(spec["btheta"]);
    r.d = spec["d"];
    r.bd = spec["bd"];
    r.log = spec["log"];
    r.mean_correction = spec["mean_correction"];
    r.innovation_variance = spec["innovation_variance"];
    r.forecast_count = spec["forecast"];
    r.backcast_count = spec["backcast"];
    r.xl = spec["xl"];
    r.epsphi = spec["epsphi"];
    r.rmod = spec["rmod"];
    r.smod = spec["smod"];
    r.smod_pi = spec["smod_pi"];
    const std::string mode = spec["approximation"];
    r.mode = mode == "legacy" ? kalos::seats_approx_mode::legacy :
             mode == "none" ? kalos::seats_approx_mode::none : kalos::seats_approx_mode::noisy;
    r.components_method = Rcpp::as<std::string>(spec["method"]) == "burman" ?
        kalos::seats_components_method::burman : kalos::seats_components_method::kalman_smoother;
    r.bias = Rcpp::as<std::string>(spec["bias"]) == "legacy" ?
        kalos::seats_bias_correction::legacy : kalos::seats_bias_correction::none;
    const auto result = kalos::seats_extract_wk(r);
    return Rcpp::List::create(Rcpp::_["status"] = result.fact("status"),
        Rcpp::_["initial"] = component_table(result.initial_components()),
        Rcpp::_["final"] = component_table(result.final_components()));
}
