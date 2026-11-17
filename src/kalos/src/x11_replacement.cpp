#include "kalos/x11_replacement.hpp"
#include "kalos/parallel.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace kalos {
namespace {
// Independent reimplementation of the pinned default corrector replacement
// stage (item 22b.3b); no upstream source is translated. The rules, walk
// order and fallbacks are documented in docs/contracts/x11-replacement-
// contract.md and qualified by validation/evidence/native/x11-replacement.
void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }
// The pinned walk compares against exactly one; masked rows keep weight one
// and are selectable neighbors carrying effective NaN values.
double effective(const MaskedSeries& series, std::size_t row) {
    return series.missing[row] ? std::numeric_limits<double>::quiet_NaN() : series.values[row];
}
} // namespace

X11ReplacementPlan x11_replacements(const MaskedSeries& series,int frequency,
        const std::vector<double>& weights,int backcasts,int forecasts,
        bool exclude_projections,bool multiplicative,const Interrupt& interrupt) {
    series.validate();
    if(frequency != 4 && frequency != 12)
        throw std::invalid_argument("X-11 replacements support monthly and quarterly series");
    const auto n = series.values.size();
    if(n == 0) throw std::invalid_argument("X-11 replacements require at least one observation");
    if(weights.size() != n)
        throw std::invalid_argument("X-11 replacement weights cover every observation");
    for(const double weight:weights) if(!std::isfinite(weight))
        throw std::invalid_argument("X-11 replacement weights must be finite");
    if(backcasts < 0 || forecasts < 0)
        throw std::invalid_argument("X-11 replacement projection horizons must be nonnegative");
    if(exclude_projections &&
        static_cast<std::size_t>(backcasts) + static_cast<std::size_t>(forecasts) > n)
        // The pinned projection masking writes outside the output array once
        // the horizons exceed the series; undefined behavior is not a result.
        throw std::invalid_argument("X-11 replacement projection horizons exceed the series");
    checkpoint(interrupt);
    const auto f = static_cast<std::size_t>(frequency);
    X11ReplacementPlan plan;
    plan.replacements.assign(n, std::nan(""));
    plan.kinds.assign(n, X11ReplacementKind::none);
    plan.neighbors.assign(n, { {-1,-1,-1,-1} });
    plan.neighbor_counts.assign(n, 0);
    plan.correction_factors.assign(n, multiplicative ? 1.0 : 0.0);
    plan.corrected.resize(n);
    // Pinned period averages: ascending same-phase sums including the masked
    // (effective NaN) values, so a poisoned phase stays poisoned.
    plan.period_averages.assign(f, std::nan(""));
    for(std::size_t phase=0;phase<f;++phase) {
        double sum=0; std::size_t count=0;
        for(std::size_t row=phase;row<n;row+=f) { sum+=effective(series,row); ++count; }
        plan.period_averages[phase] = sum / static_cast<double>(count);
    }
    const auto lo = exclude_projections ? static_cast<std::size_t>(backcasts) : std::size_t{0};
    const auto hi = exclude_projections ? n - static_cast<std::size_t>(forecasts) : n;
    // Rows are independent; the neighbor walk and the blend accumulation are
    // serial per row, so OpenMP and serial runs agree bit for bit.
    parallel_for(n,1,interrupt,[&](std::size_t row) {
        if(exclude_projections && (row < lo || row >= hi)) {
            plan.kinds[row]=X11ReplacementKind::outside_span; return;
        }
        const auto weight=weights[row];
        if(weight == 1) return; // Kind none: no replacement, corrections cancel.
        // Pinned discovery walk: up to two fully weighted same-phase rows per
        // side, then the opposite side extended to four in total.
        auto l=static_cast<std::int64_t>(row), r=static_cast<std::int64_t>(row);
        std::int64_t left=0, right=0, count=0;
        auto& positions = plan.neighbors[row];
        while(l >= static_cast<std::int64_t>(f) && left != 2) {
            l -= static_cast<std::int64_t>(f);
            if(weights[static_cast<std::size_t>(l)] == 1) { positions[count++]=l; ++left; }
        }
        while(r < static_cast<std::int64_t>(n - f) && right != 2) {
            r += static_cast<std::int64_t>(f);
            if(weights[static_cast<std::size_t>(r)] == 1) { positions[count++]=r; ++right; }
        }
        if(left < 2) while(r < static_cast<std::int64_t>(n - f) && count < 4) {
            r += static_cast<std::int64_t>(f);
            if(weights[static_cast<std::size_t>(r)] == 1) { positions[count++]=r; ++left; }
        }
        if(right < 2) while(l >= static_cast<std::int64_t>(f) && count < 4) {
            l -= static_cast<std::int64_t>(f);
            if(weights[static_cast<std::size_t>(l)] == 1) { positions[count++]=l; ++right; }
        }
        if(left + right != 4) { // Insufficient neighbors: pinned fallback.
            positions={ {-1,-1,-1,-1} };
            plan.kinds[row]=X11ReplacementKind::period_average;
            plan.replacements[row]=plan.period_averages[row % f];
            return;
        }
        plan.neighbor_counts[row]=4;
        plan.kinds[row]=X11ReplacementKind::four_neighbors;
        // Pinned arithmetic: v = w*x, add the four neighbors in discovery
        // order, then multiply by the reciprocal of four plus the weight.
        auto value=weight*series.values[row];
        for(std::size_t slot=0;slot<4;++slot) value+=effective(series,
            static_cast<std::size_t>(positions[slot]));
        plan.replacements[row]=value*(1.0/(4.0+weight));
    });
    for(std::size_t row=0;row<n;++row) {
        if((row & 1023U) == 0) checkpoint(interrupt);
        const auto weight=weights[row];
        if(weight < 1) { // Exactly the pinned getter branch condition.
            const auto x=effective(series,row);
            plan.correction_factors[row] = multiplicative ?
                x/(1.0+weight*(x-1.0)) : x*(1.0-weight);
        }
        plan.corrected[row] = std::isnan(plan.replacements[row]) ?
            effective(series,row) : plan.replacements[row];
    }
    plan.fallbacks = static_cast<std::int64_t>(std::count(plan.kinds.begin(),
        plan.kinds.end(), X11ReplacementKind::period_average));
    plan.status="available";
    checkpoint(interrupt);
    return plan;
}
} // namespace kalos
