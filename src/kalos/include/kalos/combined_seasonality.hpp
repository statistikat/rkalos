#ifndef KALOS_COMBINED_SEASONALITY_HPP
#define KALOS_COMBINED_SEASONALITY_HPP

#include <optional>
#include <string>
#include <vector>

namespace kalos {

struct CombinedSiSeasonality {
    std::string summary;
    double stable_f, stable_p;
    double moving_f, moving_p;
    double kruskal_wallis, kruskal_p;
};

// All-years public SI test for complete monthly or quarterly data. first_period
// is one-based. Fewer than three years of observations have no public result.
std::optional<CombinedSiSeasonality> combined_si_seasonality(const std::vector<double>& si,
    int frequency, int first_period, bool multiplicative);

struct seats_component_series;
struct X13FitResult;
std::optional<CombinedSiSeasonality> seats_si_seasonality(const std::vector<seats_component_series>& wk_final,
    int frequency, bool multiplicative);
std::optional<CombinedSiSeasonality> x13_si_seasonality(const X13FitResult& fit);

}  // namespace kalos
#endif
