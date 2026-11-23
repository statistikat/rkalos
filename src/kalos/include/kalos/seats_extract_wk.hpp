#pragma once

// kalos M5 slice 24-C.1: the SEATS Wiener-Kolmogorov extraction machine,
// the live SeatsKernel path from the validated model to the initial and
// final component decompositions (WienerKolmogorovEstimator.decompose +
// DefaultBiasCorrector.correctBias), as a bit-exact native mirror.
//
// Composes on the public 24-A/24-B surfaces: the request is re-validated
// (validate_tramo_seats_model), the canonical UCARIMA split is re-derived
// (the 24-B selector chain + UcarimaModel.setVarianceMax(-1, false)), then
// the kernel continues: UcarimaModel.compact(2, 2), the Burman/WK estimates
// (seats_burman_estimates.hpp), the initial decomposition assembly (the
// WKE builder surface), the bias corrector (fillLevels / correctLogs plus
// the fill{Fore,Back}casts passes) and the result domains. The emitted
// fact lines reproduce validation/reference/SeatsWkExtractionProbe.java:
// the "status" row, the three "dom" rows (SeatsResults reads the INITIAL
// components), then the fixed (ComponentType x ComponentInformation)
// enumeration of the initial and final decompositions at %.17g.
//
// Envelope (fail-closed, pinned docs/contracts/seats-wk-extraction-contract
// .md): monthly series (D11), the Burman/WK estimator only - KalmanSmoother
// is D1, a non-zero innovation variance (var 0 takes the JVM
// ConcentratedLikelihood variance path), not Noisy (24-B parity) and no
// approximation_required decision (D6: the DefaultModelApproximator retry
// loop is NOT ported; the validator already classifies those requests -
// the classical frozen airline and the fixMa degenerate case live there).
// A decomposition that came back null under mode None is a DATA outcome,
// not a refusal: the JVM SeatsException status row is mirrored. Mid-run
// JVM aborts (EcoException / DemetraException from the prediction-error
// sink, NullPointerException, MatrixException) are likewise mirrored as
// status rows - exactly what the probe records - and no fact rows follow.
//
// Error contract: std::invalid_argument (envelope refusal, pinned reason;
// validator refusals propagate unchanged). Otherwise the unit returns the
// status/dom/init/fact grammar and no exception escapes for corpus-shaped
// inputs; an unexpected C++ exception is a mirror divergence, not a JVM
// behaviour.

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "kalos/tramo_seats_model_validator.hpp"

namespace kalos {

enum class seats_component {
  series, trend, seasonal, seasonally_adjusted, irregular, calendar_effect
};

enum class seats_information {
  value, stdev, forecast, stdev_forecast, backcast, stdev_backcast
};

struct seats_component_series {
  seats_component component;
  seats_information information;
  std::int64_t start_year;
  int start_month;  // 1 = January
  std::vector<double> values;
};

struct seats_extract_wk_result {
  // JVM-faithful fact lines in probe order: "status ..." first, then
  // "dom actual|backcast|forecast n N start <year> <annualPosition>"
  // (only on status ok), then the init rows and finally the final rows,
  // each "<init|final> <Type> <Info> n N start <year> <annualPosition>
  // vals <v ...>" with the fixed enum order Series,Trend,Seasonal,
  // SeasonallyAdjusted,Irregular,CalendarEffect x Value,Stdev,Forecast,
  // StdevForecast,Backcast,StdevBackcast and absent items skipped.
  // Record-only corpus rows (toolkit, v.*, model, r/d/cmp, k.*) belong
  // to the 24-A/24-B units and are NOT emitted.
  const std::vector<std::string>& facts() const { return facts_; }

  // Convenience: the fact value for `key' (rest of the line after the
  // first whitespace run) or "" when the key is absent ("status" and the
  // dotted-free keys "dom actual" etc. match like the 24-B accessor).
  std::string fact(const std::string& key) const;

  // Present rows follow the JVM ComponentType x ComponentInformation order.
  // A present, zero-length series retains its start. Error results have no rows.
  const std::vector<seats_component_series>& initial_components() const { return initial_; }
  const std::vector<seats_component_series>& final_components() const { return final_; }

  seats_extract_wk_result() = default;
  explicit seats_extract_wk_result(std::vector<std::string> facts)
      : facts_(std::move(facts)) {}

  seats_extract_wk_result(std::vector<std::string> facts,
                          std::vector<seats_component_series> initial,
                          std::vector<seats_component_series> final)
      : facts_(std::move(facts)), initial_(std::move(initial)), final_(std::move(final)) {}

 private:
  std::vector<std::string> facts_;
  std::vector<seats_component_series> initial_, final_;
};

// Runs the extraction mirror over `request' (the frozen fixed-parameter
// grammar of the 24-A/24-B probes plus the series itself). See the header
// comment for the envelope and the failure-surface contract.
seats_extract_wk_result seats_extract_wk(const tramo_seats_model_request& request);

}  // namespace kalos
