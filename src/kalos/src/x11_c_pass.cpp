#include "kalos/x11_c_pass.hpp"
#include "kalos/detail/x11_pinned.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace kalos {
namespace {
using namespace detail;
} // namespace

X11CPass x11_c_pass(const MaskedSeries& series_in,int frequency,int first_period,
        const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
        int trend_length,X11SigmaPolicy policy,const X11BPass& b,
        double lower_sigma,double upper_sigma,const Interrupt& interrupt) {
    series_in.validate();
    // The decomposition mode is inherited from the B pass (the pinned
    // kernel runs the C step on the same context), re-gated here on the raw
    // values and re-transformed by the same kernel input transform
    // (log(y) under LogAdditive, item 22b.5 M3b-i) so the B-pass input echo
    // check below compares against the very sequence the B pass consumed.
    const auto mode=b.mode;
    check_decomposition("C pass",mode,series_in);
    const MaskedSeries series=mode==X11DecompositionMode::log_additive
        ?log_domain(series_in):series_in;
    if(frequency!=4 && frequency!=12)
        throw std::invalid_argument("C pass supports monthly and quarterly series");
    if(first_period<1 || first_period>frequency)
        throw std::invalid_argument("C pass first period is within the annual frequency");
    if(trend_length<5 || trend_length%2==0)
        throw std::invalid_argument("C pass trend filter is an odd Henderson length >= 5");
    const auto f=static_cast<std::size_t>(frequency),n=series.values.size(),t=f/2;
    check_filters(initial_filters,f,"initial seasonal filters");
    check_filters(final_filters,f,"final seasonal filters");
    if(n<static_cast<std::size_t>(trend_length))
        throw std::invalid_argument("C pass requires at least the trend filter length");
    if(n<3*t+f) throw std::invalid_argument("C pass requires three annual radii plus one cycle");
    // The pinned X11Kernel composes the C pass only onto the finished B pass
    // of the same context and series; reproduce that coupling as validation.
    if(b.status!="available")
        throw std::invalid_argument("C pass requires a completed B pass");
    if(b.first_period!=first_period || b.trend_length!=trend_length
            || b.initial_filters!=initial_filters || b.final_filters!=final_filters)
        throw std::invalid_argument("C pass settings must match the B pass");
    for(const auto& round : b.rounds)
        if(round.policy!=policy)
            throw std::invalid_argument("C pass sigma policy must match the B pass");
    if(b.b20.values.size()!=n)
        throw std::invalid_argument("C pass B result does not cover the series");
    for(std::size_t row=0;row<n;++row) { // the B pass consumed this exact series
        if(b.b1.missing[row]!=static_cast<std::uint8_t>(series.missing[row]))
            throw std::invalid_argument("C pass series differs from the B pass input");
        if(!series.missing[row] && b.b1.values[row]!=series.values[row])
            throw std::invalid_argument("C pass series differs from the B pass input");
    }
    const auto mult=x11_mode_is_multiplicative(mode);
    const auto pseudo=mode==X11DecompositionMode::pseudo_additive;
    checkpoint(interrupt);
    X11CPass pass;
    pass.c2drop=t;
    pass.mode=mode;
    pass.first_period=first_period;
    pass.initial_filters=initial_filters; pass.final_filters=final_filters;
    pass.trend_length=trend_length;
    // X11CStep.process(refSeries, B20, ctx): c1 = remove(refSeries, B20).
    const auto y=to_pinned(series);
    pass.reference_series=y;
    // X11CStepPseudoAdd.c1 = adjustRefSeries(refSeries, B20, B7, B13)
    // = ref + B7*(B13/B20 - B13) (the constructor carries B7/B13 from the
    // kernel's finished B step).
    pass.c1=pseudo?pseudo_adjust_ref_series(y,b.b20,b.b7,b.b13,interrupt)
        :remove_pinned(y,b.b20,mult,interrupt);
    // c2Step: X11TrendCycleFilterFactory.makeTrendFilter(period) over the
    // central span only (the pinned table has length n - 2*c2drop).
    const auto annual=x11_symmetric_filter(repo_view(pass.c1),X11FilterKind::trend,
        frequency,0,interrupt);
    pass.c2=to_pinned(MaskedSeries{std::vector<double>(annual.values.begin()+t,annual.values.begin()+n-t),
        std::vector<std::uint8_t>(annual.missing.begin()+t,annual.missing.begin()+n-t)});
    MaskedSeries cropped;
    cropped.values.assign(pass.c1.values.begin()+t,pass.c1.values.begin()+n-t);
    cropped.missing.assign(pass.c1.missing.begin()+t,pass.c1.missing.begin()+n-t);
    pass.c4=remove_pinned(cropped,pass.c2,mult,interrupt);
    const auto first=static_cast<std::size_t>(first_period-1);
    const auto center_position=(first+t)%f; // pinned getPosition(c2drop)
    checkpoint(interrupt);
    // c5Step: the initial requests over the cropped span, normalized with
    // iStart = c2drop (the pinned extension back to the full span).
    pass.c5a=process_pinned(pass.c4,frequency,center_position,initial_filters,interrupt);
    pass.c5=normalize_pinned(pass.c5a,frequency,center_position,0,t,final_filters,mult,interrupt);
    // X11CStepPseudoAdd.c6: ratio C1/C5 at the borders, interior
    // C1.drop - C2*(C5.drop - 1).
    pass.c6=pseudo?pseudo_with_borders(pass.c1,pass.c5,pass.c2,t,interrupt)
        :remove_pinned(pass.c1,pass.c5,mult,interrupt);
    checkpoint(interrupt);
    // c7Step: the fixed Henderson filter with mmsre endpoints over C6.
    pass.c7=pinned_trend_filter(pass.c6,frequency,trend_length,interrupt);
    // The pinned c7Step positivity repair of the trend table under a
    // multiplicative context (the one X11CStep makePositivity site).
    if(mult) make_positivity(pass.c7.values);
    pass.c9=remove_pinned(pass.c1,pass.c7,mult,interrupt);
    checkpoint(interrupt);
    // cFinalStep: the final requests over C9, normalized with iStart = 0.
    pass.c10a=process_pinned(pass.c9,frequency,first,final_filters,interrupt);
    pass.c10=normalize_pinned(pass.c10a,frequency,first,0,0,final_filters,mult,interrupt);
    // X11CStepPseudoAdd.c11: removeSeasonalAdjusted(refSeries, C10, C7)
    // = ref - C7*(C10 - 1).
    pass.c11=pseudo?pseudo_remove_seasonal_adjusted(y,pass.c10,pass.c7,interrupt)
        :remove_pinned(y,pass.c10,mult,interrupt); // remove(refSeries, c10)
    pass.c13=remove_pinned(pass.c11,pass.c7,mult,interrupt);
    checkpoint(interrupt);
    auto& round=pass.round;
    round.detection=pass.c13;
    round.analysis_position=static_cast<std::int64_t>(first);
    round.policy=policy;
    round.sigma=x11_sigma_weights(repo_view(pass.c13),frequency,
        static_cast<int>(first)+1,policy,0,0,false,lower_sigma,upper_sigma,mult,interrupt);
    round.replacements=x11_replacements(repo_view(pass.c13),frequency,round.sigma.weights,
        0,0,false,mult,interrupt);
    pass.c17=MaskedSeries{round.sigma.weights,missing_of(round.sigma.weights)};
    pass.c20=MaskedSeries{round.replacements.correction_factors,
        missing_of(round.replacements.correction_factors)};
    pass.status="available";
    checkpoint(interrupt);
    return pass;
}
} // namespace kalos
