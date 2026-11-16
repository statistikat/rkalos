#include "kalos/x11_c_pass.hpp"
#include "kalos/detail/x11_pinned.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace kalos {
namespace {
using namespace detail;
} // namespace

X11CPass x11_c_pass(const MaskedSeries& series,int frequency,int first_period,
        const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
        int trend_length,X11SigmaPolicy policy,const X11BPass& b,
        double lower_sigma,double upper_sigma,const Interrupt& interrupt) {
    series.validate();
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
    checkpoint(interrupt);
    X11CPass pass;
    pass.c2drop=t;
    pass.initial_filters=initial_filters; pass.final_filters=final_filters;
    pass.trend_length=trend_length;
    // X11CStep.process(refSeries, B20, ctx): c1 = remove(refSeries, B20).
    const auto y=to_pinned(series);
    pass.c1=minus_pinned(y,b.b20,interrupt);
    // c2Step: X11TrendCycleFilterFactory.makeTrendFilter(period) over the
    // central span only (the pinned table has length n - 2*c2drop).
    const auto annual=x11_symmetric_filter(repo_view(pass.c1),X11FilterKind::trend,
        frequency,0,interrupt);
    pass.c2=to_pinned(MaskedSeries{std::vector<double>(annual.values.begin()+t,annual.values.begin()+n-t),
        std::vector<std::uint8_t>(annual.missing.begin()+t,annual.missing.begin()+n-t)});
    MaskedSeries cropped;
    cropped.values.assign(pass.c1.values.begin()+t,pass.c1.values.begin()+n-t);
    cropped.missing.assign(pass.c1.missing.begin()+t,pass.c1.missing.begin()+n-t);
    pass.c4=minus_pinned(cropped,pass.c2,interrupt);
    const auto first=static_cast<std::size_t>(first_period-1);
    const auto center_position=(first+t)%f; // pinned getPosition(c2drop)
    checkpoint(interrupt);
    // c5Step: the initial requests over the cropped span, normalized with
    // iStart = c2drop (the pinned extension back to the full span).
    pass.c5a=process_pinned(pass.c4,frequency,center_position,initial_filters,interrupt);
    pass.c5=normalize_pinned(pass.c5a,frequency,center_position,0,t,final_filters,interrupt);
    pass.c6=minus_pinned(pass.c1,pass.c5,interrupt);
    checkpoint(interrupt);
    // c7Step: the fixed Henderson filter with mmsre endpoints over C6.
    pass.c7=pinned_trend_filter(pass.c6,frequency,trend_length,interrupt);
    pass.c9=minus_pinned(pass.c1,pass.c7,interrupt);
    checkpoint(interrupt);
    // cFinalStep: the final requests over C9, normalized with iStart = 0.
    pass.c10a=process_pinned(pass.c9,frequency,first,final_filters,interrupt);
    pass.c10=normalize_pinned(pass.c10a,frequency,first,0,0,final_filters,interrupt);
    pass.c11=minus_pinned(y,pass.c10,interrupt); // remove(refSeries, c10)
    pass.c13=minus_pinned(pass.c11,pass.c7,interrupt);
    checkpoint(interrupt);
    auto& round=pass.round;
    round.detection=pass.c13;
    round.analysis_position=static_cast<std::int64_t>(first);
    round.policy=policy;
    round.sigma=x11_sigma_weights(repo_view(pass.c13),frequency,
        static_cast<int>(first)+1,policy,0,0,false,lower_sigma,upper_sigma,false,interrupt);
    round.replacements=x11_replacements(repo_view(pass.c13),frequency,round.sigma.weights,
        0,0,false,false,interrupt);
    pass.c17=MaskedSeries{round.sigma.weights,missing_of(round.sigma.weights)};
    pass.c20=MaskedSeries{round.replacements.correction_factors,
        missing_of(round.replacements.correction_factors)};
    pass.status="available";
    checkpoint(interrupt);
    return pass;
}
} // namespace kalos
