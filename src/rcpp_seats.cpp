#include "rcpp_tramo_fit.hpp"
#include <kalos/arma_diagnostics.hpp>
#include <kalos/tramo_seats_pipeline.hpp>

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

// [[Rcpp::export(name = ".native_seats_from_tramo", rng = false)]]
Rcpp::List native_seats_from_tramo(const Rcpp::List& spec, const Rcpp::List& calendar,
        std::vector<double> series, Rcpp::DataFrame outliers) {
    kalos::TramoFitResult fit;
    fit.status = "available";
    auto& builder = fit.builder;
    const Rcpp::NumericVector start = calendar["start"];
    builder.annual_frequency = Rcpp::as<int>(calendar["frequency"]);
    builder.domain_year = builder.estimation_year = static_cast<int>(start[0]);
    builder.domain_position = builder.estimation_position = static_cast<int>(start[1]) - 1;
    builder.domain_length = builder.estimation_length = static_cast<int>(series.size());
    builder.interpolated.swap(series);
    fit.orders.d = spec["d"];
    fit.orders.bd = spec["bd"];
    auto append = [&](const char* name, int& order) {
        const auto values = Rcpp::as<std::vector<double>>(spec[name]);
        order = static_cast<int>(values.size());
        fit.parameters.insert(fit.parameters.end(), values.begin(), values.end());
    };
    append("phi", fit.orders.p);
    append("bphi", fit.orders.bp);
    append("theta", fit.orders.q);
    append("btheta", fit.orders.bq);
    fit.log = spec["log"];
    fit.mean = spec["mean_correction"];
    fit.ls_ssq = spec["ssq"];
    fit.ls_eff = spec["effective"];
    fit.ls_npar = spec["estimated"];
    const Rcpp::CharacterVector types = outliers["type"];
    const Rcpp::IntegerVector positions = outliers["position"];
    const Rcpp::NumericVector coefficients = outliers["coefficient"];
    for (R_xlen_t i = 0; i < types.size(); ++i) {
        kalos::TramoFitVariable variable;
        if (types[i] == "AO") variable.type = kalos::TramoOutlierKind::Ao;
        else if (types[i] == "LS") variable.type = kalos::TramoOutlierKind::Ls;
        else if (types[i] == "TC") variable.type = kalos::TramoOutlierKind::Tc;
        else Rcpp::stop("Unknown fitted outlier type.");
        variable.position = positions[i] - 1;
        variable.coefficient = coefficients[i];
        fit.variables.push_back(variable);
    }
    const auto result = kalos::run_tramo_seats_pipeline(fit, spec["backcast"], spec["forecast"]);
    return Rcpp::List::create(Rcpp::_["status"] = result.extraction.fact("status"),
        Rcpp::_["initial"] = component_table(result.extraction.initial_components()),
        Rcpp::_["final"] = component_table(result.extraction.final_components()),
        Rcpp::_["merged"] = component_table(result.finals));
}

// [[Rcpp::export(name = ".native_tramoseats_rsa0", rng = false)]]
Rcpp::List native_tramoseats_rsa0(const std::vector<double>& series, int year, int month,
        int backcast, int forecast) {
    const auto result = kalos::run_tramo_seats_rsa0(series, year, month, backcast, forecast);
    const auto& fitted = result.preprocessing;
    const auto& model = fitted.model;
    const auto& decomposition = result.decomposition;
    // Both projections receive the month counts normalized at the R boundary.
    const auto errors = kalos::tramo_tr0_forecast_errors(fitted, forecast);
    const auto lb = kalos::rsa0_residual_ljung_box(fitted.residuals.values);
    const auto normality = kalos::rsa0_residual_normality(fitted.residuals.values);
    Rcpp::RObject dh = R_NilValue;
    if (normality) {
        dh = Rcpp::List::create(Rcpp::_["value"] = normality->value,
            Rcpp::_["pvalue"] = normality->p, Rcpp::_["df"] = normality->df);
    }
    return Rcpp::List::create(
        Rcpp::_["status"] = decomposition.extraction.fact("status"),
        Rcpp::_["residuals"] = Rcpp::List::create(
            Rcpp::_["start"] = Rcpp::NumericVector::create(
                static_cast<double>(fitted.residuals.start_year), fitted.residuals.start_month),
            Rcpp::_["tsres"] = fitted.residuals.values, Rcpp::_["ser"] = fitted.residuals.standard_error,
            Rcpp::_["lb"] = Rcpp::List::create(Rcpp::_["value"] = lb.value, Rcpp::_["pvalue"] = lb.p,
                Rcpp::_["df"] = lb.df, Rcpp::_["lags"] = 24),
            Rcpp::_["doornikhansen"] = dh),
        Rcpp::_["series_forecast_stdev_status"] = errors.status,
        Rcpp::_["series_forecast_stdev"] = errors.standard_errors,
        Rcpp::_["model"] = Rcpp::List::create(
            Rcpp::_["d"] = model.d, Rcpp::_["bd"] = model.bd,
            Rcpp::_["orders"] = Rcpp::IntegerVector::create(
                Rcpp::_["p"] = static_cast<int>(model.phi.size()), Rcpp::_["d"] = model.d,
                Rcpp::_["q"] = static_cast<int>(model.theta.size()), Rcpp::_["bp"] = static_cast<int>(model.bphi.size()),
                Rcpp::_["bd"] = model.bd, Rcpp::_["bq"] = static_cast<int>(model.btheta.size())),
            Rcpp::_["phi"] = model.phi, Rcpp::_["bphi"] = model.bphi,
            Rcpp::_["theta"] = model.theta, Rcpp::_["btheta"] = model.btheta,
            Rcpp::_["log"] = model.log, Rcpp::_["mean_correction"] = model.mean_correction,
            Rcpp::_["mean_coefficient"] = fitted.mean_coefficient,
            Rcpp::_["innovation_variance"] = model.innovation_variance),
        Rcpp::_["likelihood"] = Rcpp::List::create(
            Rcpp::_["log_likelihood"] = fitted.log_likelihood, Rcpp::_["ssq"] = fitted.residual_ssq,
            Rcpp::_["effective"] = fitted.effective_observations, Rcpp::_["estimated"] = fitted.estimated_parameters),
        Rcpp::_["initial"] = component_table(decomposition.extraction.initial_components()),
        Rcpp::_["final"] = component_table(decomposition.extraction.final_components()),
        Rcpp::_["merged"] = component_table(decomposition.finals));
}

// [[Rcpp::export(name = ".native_tramoseats_rsa1", rng = false)]]
Rcpp::List native_tramoseats_rsa1(const std::vector<double>& series, int year, int month,
        int backcast, int forecast) {
    const auto fit = kalos::run_tramo_seats_rsa1(series, year, month, backcast, forecast);
    const auto& decomposition = fit.decomposition;
    auto result = wrap_tramo_preprocessing(fit.preprocessing);
    result["status"] = decomposition.extraction.fact("status");
    result["model_changed"] = decomposition.extraction.model_changed;
    result["initial"] = component_table(decomposition.extraction.initial_components());
    result["final"] = component_table(decomposition.extraction.final_components());
    result["merged"] = component_table(decomposition.finals);
    return result;
}
