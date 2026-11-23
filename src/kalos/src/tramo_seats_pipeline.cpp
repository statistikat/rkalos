#include "kalos/tramo_seats_pipeline.hpp"

#include "kalos/detail/fdlibm.hpp"
#include "kalos/outlier_scoring.hpp"
#include "kalos/tramo_seats_bridge.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace kalos {
namespace {

using Row = seats_component_series;
using Component = seats_component;
using Info = seats_information;
constexpr std::array<Info, 3> kinds{Info::value, Info::forecast, Info::backcast};
constexpr std::array<Component, 5> components{Component::series, Component::trend,
    Component::seasonal, Component::seasonally_adjusted, Component::irregular};
using Grid = std::array<std::array<const Row*, 3>, 5>;

void require(bool ok, const char* reason) {
    if (!ok) throw std::invalid_argument(std::string("merge_tramo_seats_components: ") + reason);
}

Grid component_grid(const std::vector<Row>& wk) {
    Grid rows{};
    for (const auto& row : wk) {
        if (row.information == Info::stdev || row.information == Info::stdev_forecast ||
            row.information == Info::stdev_backcast) continue;
        const auto t = static_cast<std::size_t>(row.component);
        require(t < components.size(), "unsupported WK component");
        const auto span = std::find(kinds.begin(), kinds.end(), row.information);
        require(span != kinds.end(), "unsupported WK information");
        auto& slot = rows[t][static_cast<std::size_t>(span - kinds.begin())];
        require(!slot, "duplicate WK component");
        slot = &row;
    }
    return rows;
}

// Both series and deterministic effects are monthly throughout the bridge envelope.
std::vector<Row> merge(const std::vector<double>& original, const std::vector<TramoFitVariable>& variables,
                       const tramo_seats_model_request& model,
                       const std::vector<Row>& wk) {
    const auto rows = component_grid(wk);
    const auto& series = rows.front();
    require(series[0] && series[0]->values == model.series,
            "observed WK series disagrees with the fitted bridge input");
    const std::size_t n = model.series.size();
    const std::size_t nf = series[1] ? series[1]->values.size() : 0;
    const std::size_t nb = series[2] ? series[2]->values.size() : 0;
    const std::size_t limit = std::numeric_limits<int>::max();
    require(n <= limit && nf <= limit - n && nb <= limit - n - nf, "extended domain exceeds index range");
    const auto first = static_cast<std::int64_t>(model.start_year) * 12 + model.start_month - 1;
    const std::array<std::int64_t, 3> starts{first, first + static_cast<std::int64_t>(n),
                                           first - static_cast<std::int64_t>(nb)};
    const std::array<std::size_t, 3> lengths{n, nf, nb};
    const std::array<std::size_t, 3> offsets{nb, nb + n, 0};
    std::vector<double> trend(n + nf + nb), irregular(trend.size());
    for (const auto& var : variables) {
        const bool ls = var.type == TramoOutlierKind::Ls;
        auto& target = ls ? trend : irregular;
        std::vector<double> column(target.size());
        fill_outlier_column(static_cast<OutlierType>(var.type),
                            nb + static_cast<std::size_t>(var.position), 12, var.tc_rate, column);
        for (std::size_t i = 0; i < column.size(); ++i) target[i] += column[i] * var.coefficient;
    }
    if (model.log) {
        for (double& x : trend) x = detail::fdlibm_exp(x);
        for (double& x : irregular) x = detail::fdlibm_exp(x);
    }
    auto combine = [log = model.log](double a, double b) { return log ? a * b : a + b; };
    std::vector<Row> output;
    for (std::size_t t = 0; t < components.size(); ++t) {
        for (std::size_t s = 0; s < kinds.size(); ++s) {
            const Row* row = rows[t][s];
            if (!lengths[s]) {
                require(!row || row->values.empty(), "nonempty WK component outside its series domain");
                continue;
            }
            require(row && row->values.size() == lengths[s], "incomplete WK domain");
            const int month = static_cast<int>((starts[s] % 12 + 12) % 12) + 1;
            require(row->start_month == month && row->start_year == (starts[s] - (month - 1)) / 12,
                    "misaligned WK component");
            output.push_back(*row);
            auto& values = output.back().values;
            if (components[t] == Component::series && kinds[s] == Info::value) {
                values = original;
                continue;
            }
            for (std::size_t j = 0; j < values.size(); ++j) {
                const std::size_t i = offsets[s] + j;
                switch (components[t]) {
                case Component::trend: values[j] = combine(trend[i], values[j]); break;
                case Component::irregular: values[j] = combine(irregular[i], values[j]); break;
                case Component::series:
                case Component::seasonally_adjusted:
                    values[j] = combine(combine(trend[i], irregular[i]), values[j]); break;
                default: break;  // the qualified fit has no seasonal/calendar regressors
                }
            }
        }
    }
    return output;
}

tramo_seats_pipeline_result extract_and_merge(const std::vector<double>& original,
        const std::vector<TramoFitVariable>& variables, const tramo_seats_model_request& model) {
    tramo_seats_pipeline_result result;
    result.extraction = seats_extract_wk(model);
    if (result.extraction.fact("status") == "ok")
        result.finals = merge(original, variables, model, result.extraction.final_components());
    return result;
}

}  // namespace

std::vector<seats_component_series> merge_tramo_seats_components(
        const TramoFitResult& fit, const std::vector<seats_component_series>& wk_final) {
    return merge(fit.builder.interpolated, fit.variables, make_tramo_seats_model(fit), wk_final);
}

tramo_seats_pipeline_result run_tramo_seats_pipeline(const TramoFitResult& fit,
        int backcast_count, int forecast_count) {
    auto model = make_tramo_seats_model(fit);
    model.backcast_count = backcast_count;
    model.forecast_count = forecast_count;
    return extract_and_merge(fit.builder.interpolated, fit.variables, model);
}

tramo_seats_rsa0_result run_tramo_seats_rsa0(const std::vector<double>& series, int start_year,
        int start_month, int backcast_count, int forecast_count) {
    tramo_seats_rsa0_result result;
    result.preprocessing = estimate_tramo_tr0(series, start_year, start_month);
    auto& model = result.preprocessing.model;
    model.backcast_count = backcast_count;
    model.forecast_count = forecast_count;
    result.decomposition = extract_and_merge(series, {}, model);
    return result;
}

tramo_seats_automatic_result run_tramo_seats_rsa1(const std::vector<double>& series,
    int start_year, int start_month, int backcast_count, int forecast_count) {
    auto fit = estimate_tramo_tr1(series, start_year, start_month);
    auto decomposition = run_tramo_seats_pipeline(fit, backcast_count, forecast_count);
    return {std::move(fit), std::move(decomposition)};
}

tramo_seats_automatic_result run_tramo_seats_rsa3(const std::vector<double>& series,
    int start_year, int start_month, int backcast_count, int forecast_count) {
    auto fit = estimate_tramo_tr3(series, start_year, start_month);
    auto decomposition = run_tramo_seats_pipeline(fit, backcast_count, forecast_count);
    return {std::move(fit), std::move(decomposition)};
}

}  // namespace kalos
