#include "kalos/tramo_seats_bridge.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/outlier_scoring.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace kalos {
namespace {

void require(bool ok, const char* detail) {
    if (!ok) throw std::invalid_argument(std::string("make_tramo_seats_model: ") + detail);
}

void copy_parameters(const TramoFitResult& fit, tramo_seats_model_request& model) {
    const auto& o = fit.orders;
    require(o.d >= 0 && o.bd >= 0, "negative differencing order");
    model.d = o.d;
    model.bd = o.bd;
    auto cursor = fit.parameters.cbegin();
    const auto last = fit.parameters.cend();
    auto take = [&](int count, std::vector<double>& dst) {
        require(count >= 0 && count <= last - cursor, "fitted parameter count disagrees with orders");
        dst.assign(cursor, cursor + count);
        cursor += count;
        require(std::all_of(dst.begin(), dst.end(), [](double v) { return std::isfinite(v); }),
                "non-finite fitted parameter");
    };
    take(o.p, model.phi);
    take(o.bp, model.bphi);
    take(o.q, model.theta);
    take(o.bq, model.btheta);
    require(cursor == last, "extra fitted parameters");
}

}  // namespace

tramo_seats_model_request make_tramo_seats_model(const TramoFitResult& fit) {
    require(fit.status == "available", "requires an available TRAMO fit");
    const auto& build = fit.builder;
    require(build.annual_frequency == 12, "monthly fits only");
    const auto& input = build.interpolated;
    require(!input.empty() && build.domain_length >= 0 &&
            static_cast<std::size_t>(build.domain_length) == input.size() &&
            build.domain_position >= 0 && build.domain_position < 12, "invalid input domain");
    require(build.estimation_year == build.domain_year &&
            build.estimation_position == build.domain_position &&
            build.estimation_length == build.domain_length, "requires a full estimation domain");
    require(build.missing.empty() && fit.cl_nmissing == 0, "missing observations are outside the bridge gate");
    require(build.preadjustment == tramo_lop::none, "length-of-period correction is outside the bridge gate");
    for (const auto& variable : build.variables) {
        require(variable.kind == tramo_variable_kind::trend_constant,
                "pre-registered regressors are outside the bridge gate");
    }
    require(fit.ls_npar >= 0 && fit.ls_eff > fit.ls_npar &&
            fit.ls_eff <= build.domain_length && std::isfinite(fit.ls_ssq) && fit.ls_ssq > 0,
            "invalid fitted likelihood for the positive-variance SEATS envelope");
    tramo_seats_model_request model;
    model.start_year = build.domain_year;
    model.start_month = build.domain_position + 1;
    model.log = fit.log;
    model.mean_correction = fit.mean;
    model.innovation_variance = fit.ls_ssq / (fit.ls_eff - fit.ls_npar);
    require(model.innovation_variance > 0, "transmitted variance underflow");
    copy_parameters(fit, model);

    std::vector<double> effect(input.size(), 0.0);
    for (const auto& var : fit.variables) {
        require(var.type == TramoOutlierKind::Ao || var.type == TramoOutlierKind::Ls ||
                (var.type == TramoOutlierKind::Tc && var.tc_rate == 0.7),
                "qualified effects are AO, LS and TC with monthly decay 0.7");
        require(var.position >= 0 && static_cast<std::size_t>(var.position) < input.size() &&
                std::isfinite(var.coefficient), "invalid fitted outlier");
        std::vector<double> column(input.size(), 0.0);
        fill_outlier_column(static_cast<OutlierType>(var.type), static_cast<std::size_t>(var.position),
                            12, var.tc_rate, column);
        require(var.column.empty() || var.column == column, "outlier column disagrees with its factory");
        for (std::size_t i = 0; i < input.size(); ++i) effect[i] += var.coefficient * column[i];
    }
    model.series.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        require(std::isfinite(input[i]) && (!fit.log || input[i] > 0), "invalid interpolated observation");
        const double value = fit.log ? input[i] / detail::fdlibm_exp(effect[i]) : input[i] - effect[i];
        require(std::isfinite(value) && (!fit.log || value > 0), "adjusted observation outside SEATS domain");
        model.series.push_back(value);
    }
    return model;
}

}  // namespace kalos
