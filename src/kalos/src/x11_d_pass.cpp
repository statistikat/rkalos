#include "kalos/x11_d_pass.hpp"
#include "kalos/detail/x11_pinned.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace kalos {
namespace {
using namespace detail;
} // namespace

X11DPass x11_d_pass(const MaskedSeries& series_in,int frequency,int first_period,
        const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
        int trend_length,X11SigmaPolicy policy,const X11CPass& c,
        int backcast,int forecast,const Interrupt& interrupt) {
    series_in.validate();
    // The decomposition mode is inherited from the C pass, re-gated here on
    // the raw values and re-transformed by the same kernel input transform
    // (log(y) under LogAdditive, item 22b.5 M3b-i) so the C-pass input echo
    // check below compares against the very sequence the C pass consumed.
    const auto mode=c.mode;
    check_decomposition("D pass",mode,series_in);
    const MaskedSeries series=mode==X11DecompositionMode::log_additive
        ?log_domain(series_in):series_in;
    if(frequency!=4 && frequency!=12)
        throw std::invalid_argument("D pass supports monthly and quarterly series");
    if(first_period<1 || first_period>frequency)
        throw std::invalid_argument("D pass first period is within the annual frequency");
    if(trend_length<5 || trend_length%2==0)
        throw std::invalid_argument("D pass trend filter is an odd Henderson length >= 5");
    const auto f=static_cast<std::size_t>(frequency),n=series.values.size(),t=f/2;
    check_filters(initial_filters,f,"initial seasonal filters");
    check_filters(final_filters,f,"final seasonal filters");
    if(n<static_cast<std::size_t>(trend_length))
        throw std::invalid_argument("D pass requires at least the trend filter length");
    if(n<3*t+f) throw std::invalid_argument("D pass requires three annual radii plus one cycle");
    if(backcast<0 || forecast<0)
        throw std::invalid_argument("D pass horizons are non-negative");
    // The MSRE stride blocks of the dropped original calendar need at least
    // the 7-MA endpoint radius (four points) per phase, which the pinned
    // endpoint ranges silently truncate; a full four cycles keeps them whole.
    if(static_cast<std::int64_t>(n)-backcast-forecast<static_cast<std::int64_t>(4*f))
        throw std::invalid_argument("D pass requires four cycles of original calendar");
    // The pinned X11Kernel composes the D pass only onto the finished C pass
    // of the same context and series; reproduce that coupling as validation.
    if(c.status!="available")
        throw std::invalid_argument("D pass requires a completed C pass");
    if(c.first_period!=first_period || c.trend_length!=trend_length
            || c.initial_filters!=initial_filters || c.final_filters!=final_filters)
        throw std::invalid_argument("D pass settings must match the C pass");
    if(c.round.policy!=policy)
        throw std::invalid_argument("D pass sigma policy must match the C pass");
    if(c.c20.values.size()!=n)
        throw std::invalid_argument("D pass C result does not cover the series");
    for(std::size_t row=0;row<n;++row) { // the C pass consumed this exact series
        if(c.reference_series.missing[row]!=static_cast<std::uint8_t>(series.missing[row]))
            throw std::invalid_argument("D pass series differs from the C pass input");
        if(!series.missing[row] && c.reference_series.values[row]!=series.values[row])
            throw std::invalid_argument("D pass series differs from the C pass input");
    }
    const auto mult=x11_mode_is_multiplicative(mode);
    const auto pseudo=mode==X11DecompositionMode::pseudo_additive;
    checkpoint(interrupt);
    X11DPass pass;
    pass.d2drop=t;
    pass.mode=mode;
    pass.backcast=backcast; pass.forecast=forecast;
    pass.initial_filters=initial_filters; pass.final_filters=final_filters;
    pass.seas_filter=final_filters; // the fixed requests; the pinned
        // seasFilter equals the final filter because MSRE never engages
    pass.trend_length=trend_length;
    pass.d7_henderson_length=trend_length; pass.final_henderson_length=trend_length;
    pass.msr_iterations=0; // pinned constant: no MSRE filter selection
    // X11DStep.process(refSeries, C20, ctx, cvFilters): d1 = remove(refSeries, C20).
    const auto y=to_pinned(series);
    // X11DStepPseudoAdd.d1 = adjustRefSeries(refSeries, C20, C7, C13)
    // = ref + C7*(C13/C20 - C13) (the constructor carries C7/C13 from the
    // kernel's finished C step; C20 is the reference correction).
    pass.d1=pseudo?pseudo_adjust_ref_series(y,c.c20,c.c7,c.c13,interrupt)
        :remove_pinned(y,c.c20,mult,interrupt);
    // d2Step: X11TrendCycleFilterFactory.makeTrendFilter(period) over the
    // central span only (the pinned table has length n - 2*d2drop).
    const auto annual=x11_symmetric_filter(repo_view(pass.d1),X11FilterKind::trend,
        frequency,0,interrupt);
    pass.d2=to_pinned(MaskedSeries{std::vector<double>(annual.values.begin()+t,annual.values.begin()+n-t),
        std::vector<std::uint8_t>(annual.missing.begin()+t,annual.missing.begin()+n-t)});
    MaskedSeries cropped;
    cropped.values.assign(pass.d1.values.begin()+t,pass.d1.values.begin()+n-t);
    cropped.missing.assign(pass.d1.missing.begin()+t,pass.d1.missing.begin()+n-t);
    pass.d4=remove_pinned(cropped,pass.d2,mult,interrupt);
    const auto first=static_cast<std::size_t>(first_period-1);
    const auto center_position=(first+t)%f; // pinned getPosition(d2drop)
    checkpoint(interrupt);
    // d5Step: the INITIAL requests over the cropped span at getPosition(
    // d2drop), normalized with iStart = d2drop (extension back to the span).
    const auto d5raw=process_pinned(pass.d4,frequency,center_position,initial_filters,interrupt);
    pass.d5=normalize_pinned(d5raw,frequency,center_position,0,t,final_filters,mult,interrupt);
    // X11DStepPseudoAdd.d6: ratio D1/D5 at the borders, interior
    // D1.drop - D2*(D5.drop - 1).
    pass.d6=pseudo?pseudo_with_borders(pass.d1,pass.d5,pass.d2,t,interrupt)
        :remove_pinned(pass.d1,pass.d5,mult,interrupt);
    checkpoint(interrupt);
    // d7Step: the fixed Henderson filter with mmsre endpoints over D6.
    pass.d7=pinned_trend_filter(pass.d6,frequency,trend_length,interrupt);
    // The pinned d7Step positivity repair under a multiplicative context.
    if(mult) make_positivity(pass.d7.values);
    // X11DStepPseudoAdd.d8: (D1/D7)*C20 — the ratio ctx.remove runs on D1
    // (not the reference) and the result is scaled by the C20 correction.
    pass.d8=pseudo?product_pinned(remove_pinned(pass.d1,pass.d7,true,interrupt),
        c.c20,interrupt)
        :remove_pinned(y,pass.d7,mult,interrupt); // remove(refSeries, d7)
    pass.d9bis=remove_pinned(pass.d1,pass.d7,mult,interrupt);
    pass.d9=d9_mask_pinned(pass.d9bis,pass.d8,interrupt); // pinned dcmpg mask
    checkpoint(interrupt);
    // dFinalStep: the final requests over D9bis at getPosition(0), normalized
    // with iStart = 0; the final Henderson trend over D11bis.
    pass.d10bis=process_pinned(pass.d9bis,frequency,first,final_filters,interrupt);
    pass.d10=normalize_pinned(pass.d10bis,frequency,first,0,0,final_filters,mult,interrupt);
    // X11DStepPseudoAdd.d11bis: removeSeasonalAdjusted(D1, D10, D7)
    // = D1 - D7*(D10 - 1).
    pass.d11bis=pseudo?pseudo_remove_seasonal_adjusted(pass.d1,pass.d10,pass.d7,interrupt)
        :remove_pinned(pass.d1,pass.d10,mult,interrupt);
    pass.d12=pinned_trend_filter(pass.d11bis,frequency,trend_length,interrupt);
    // The pinned dFinalStep positivity repair of the final Henderson trend,
    // before D13 composes it away (the second X11DStep makePositivity site).
    if(mult) make_positivity(pass.d12.values);
    // X11DStepPseudoAdd.d11: ref - D12*(D10 - 1)
    // (removeSeasonalAdjustedFinal(refSeries, D12*(D10 - 1))).
    pass.d11=pseudo?pseudo_remove_seasonal_adjusted(y,pass.d10,pass.d12,interrupt)
        :remove_pinned(y,pass.d10,mult,interrupt); // remove(refSeries, d10)
    pass.d13=remove_pinned(pass.d11,pass.d12,mult,interrupt);
    checkpoint(interrupt);
    // The pinned MSRE diagnostic runs on the dropped original calendar at
    // the original-calendar phase getPosition(backcast).
    std::vector<double> original(pass.d9bis.values.begin()+static_cast<std::ptrdiff_t>(backcast),
        pass.d9bis.values.end()-static_cast<std::ptrdiff_t>(forecast));
    std::vector<double> ri,rs;
    std::vector<std::int64_t> counts;
    pinned_msr_table(original,frequency,(first+static_cast<std::size_t>(backcast))%f,
        mode,ri,rs,counts);
    auto& msr=pass.msr_table;
    msr.mean_irregular_evolutions=std::move(ri);
    msr.mean_seasonal_evolutions=std::move(rs);
    msr.counts=std::move(counts);
    msr.msr.assign(f,nan());
    double numerator=0,denominator=0;
    for(std::size_t phase=0;phase<f;++phase) {
        msr.msr[phase]=msr.mean_irregular_evolutions[phase]
            /msr.mean_seasonal_evolutions[phase];
        numerator=numerator+msr.mean_irregular_evolutions[phase]
            *static_cast<double>(msr.counts[phase]);
        denominator=denominator+msr.mean_seasonal_evolutions[phase]
            *static_cast<double>(msr.counts[phase]);
    }
    msr.global_msr=numerator/denominator;
    pass.status="available";
    checkpoint(interrupt);
    return pass;
}
} // namespace kalos
