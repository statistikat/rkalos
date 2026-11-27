#pragma once

#include <optional>
#include <vector>

namespace kalos::detail {

// Original-scale forecast uncertainty (the JVM y_ef surface) of a fixed-coefficient
// regression SARIMA at any landed stationary orders and differencing (M13-P5).
// This mirrors the pinned RegArimaForecasts.calcForecast semantics established
// by the dual-JVM capture proof (tmp_opencode/m13-p5/PROOF.md, agreement at
// most 9e-15 on the air/ao/seas landings):
//   * response is the transformed series INCLUDING the in-sample regression
//     effects (the JVM regY() route removes nothing while no variable is
//     flagged preadjustment);
//   * columns span the observation plus forecast domain in fitted model order,
//     paired with coefficients; covariance pairs with the columns and enters
//     UNSCALED (the innovation variance multiplies only the model term);
//   * variance_h = innovation_variance * V_h + s_h' Cov s_h, where V_h is the
//     unit-variance diffuse-filter prediction variance and s_h the future-row
//     innovations of the columns filtered on the response's per-time gains
//     (the JVM FastDkFilter skips the gain update where the response is
//     missing, so future column states propagate under T alone);
//   * raw forecast f_h = f0_h + s_h' coefficients; a log fit back-transforms
//     with LogNormal.stdev: exp(f + v/2) * sqrt(exp(v) - 1).
// The filter itself is the JVM kappa-limit dual: the canonical ARMA block on
// the differenced response (stationarily initialized) plus the exactly diffuse
// differencing history, carried as P = Pc + Pm/eps (JVM kappa = 1e12), with
// the pinned branch thresholds (diffuse pivot a > 1e-12, finite innovation
// b > 1e-14) and the 1e-9 rank-annihilation kill of the 1/eps component. The
// polynomials follow the JVM SarimaModel.parameters() convention directly:
// ar, sar, ma, sma carry the polynomial coefficients of A0(B) = (1+sum
// ar_i B^i)(1+sum sar_j B^(12j)) and likewise for MA; D(B) = (1-B)^d (1-B^12)^bd.
// A missing-time-free, complete response is required; shape violations throw
// std::invalid_argument and a non-finite or non-positive variance answers
// nullopt.
std::optional<std::vector<double>> sarima_forecast_uncertainty(
    const std::vector<double>& response,
    const std::vector<std::vector<double>>& columns,
    const std::vector<double>& ar, const std::vector<double>& sar,
    const std::vector<double>& ma, const std::vector<double>& sma,
    int d, int bd, const std::vector<double>& coefficients,
    const std::vector<double>& covariance, double innovation_variance,
    bool log, int months);

} // namespace kalos::detail
