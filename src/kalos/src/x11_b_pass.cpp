#include "kalos/x11_b_pass.hpp"
#include "kalos/detail/x11_pinned.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace kalos {
namespace {
using namespace detail;
} // namespace

X11BPass x11_b_pass(const MaskedSeries& series_in,int frequency,int first_period,
        const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
        int trend_length,X11SigmaPolicy policy,
        double lower_sigma,double upper_sigma,X11DecompositionMode mode,
        const Interrupt& interrupt) {
    series_in.validate();
    check_decomposition("B pass",mode,series_in);
    // The pinned X11Kernel.process input transform: the plain additive steps
    // run on DoubleSeq.log() under LogAdditive (KERNEL.process, item 22b.5
    // M3b-i), so every table below is a step-getter log-domain table and the
    // gated positivity check above ran on the raw values.
    const MaskedSeries series=mode==X11DecompositionMode::log_additive
        ?log_domain(series_in):series_in;
    const auto mult=x11_mode_is_multiplicative(mode);
    // X11{B,C,D}StepPseudoAdd override the composition sites (item 22b.5
    // M3b-ii); every other site keeps the multiplicative context they
    // inherit (isMultiplicative() is true for PseudoAdditive).
    const auto pseudo=mode==X11DecompositionMode::pseudo_additive;
    if(frequency!=4 && frequency!=12)
        throw std::invalid_argument("B pass supports monthly and quarterly series");
    if(first_period<1 || first_period>frequency)
        throw std::invalid_argument("B pass first period is within the annual frequency");
    if(trend_length<5 || trend_length%2==0)
        throw std::invalid_argument("B pass trend filter is an odd Henderson length >= 5");
    const auto f=static_cast<std::size_t>(frequency),n=series.values.size(),t=f/2;
    check_filters(initial_filters,f,"initial seasonal filters");
    check_filters(final_filters,f,"final seasonal filters");
    // Below the trend length the pinned AsymmetricEndPoints replaces the
    // endpoint region with a mean fallback that is not delivered here.
    if(n<static_cast<std::size_t>(trend_length))
        throw std::invalid_argument("B pass requires at least the trend filter length");
    // Keeps every normalization span a full cycle plus the annual radius
    // (the pinned endpoint restoration indexes one cycle inward).
    if(n<3*t+f) throw std::invalid_argument("B pass requires three annual radii plus one cycle");
    checkpoint(interrupt);
    X11BPass pass;
    pass.b2drop=t;
    pass.first_period=first_period;
    pass.initial_filters=initial_filters; pass.final_filters=final_filters;
    pass.trend_length=trend_length;
    pass.mode=mode;
    pass.b1=to_pinned(series);
    const auto annual=x11_symmetric_filter(series,X11FilterKind::trend,frequency,0,interrupt);
    pass.b2=to_pinned(MaskedSeries{std::vector<double>(annual.values.begin()+t,annual.values.begin()+n-t),
        std::vector<std::uint8_t>(annual.missing.begin()+t,annual.missing.begin()+n-t)});
    MaskedSeries cropped;
    cropped.values.assign(series.values.begin()+t,series.values.begin()+n-t);
    cropped.missing.assign(series.missing.begin()+t,series.missing.begin()+n-t);
    pass.b3=remove_pinned(cropped,pass.b2,mult,interrupt);
    const auto first=static_cast<std::size_t>(first_period-1);
    const auto center_position=(first+t)%f;   // pinned getPosition(b2drop)
    const std::size_t initial_position=first; // pinned getPosition(0)
    checkpoint(interrupt);
    pass.b4a=process_pinned(pass.b3,frequency,center_position,initial_filters,interrupt);
    pass.b4anorm=normalize_pinned(pass.b4a,frequency,initial_position,t,0,final_filters,
        mult,interrupt);
    // X11BStepPseudoAdd.b4d: PseudoAddUtility.removeIrregular(B3, B4anorm).
    pass.b4d=pseudo?pseudo_remove_irregular(pass.b3,pass.b4anorm,interrupt)
        :remove_pinned(pass.b3,pass.b4anorm,mult,interrupt);
    checkpoint(interrupt);
    auto& round1=pass.rounds[0];
    round1.detection=pass.b4d;
    round1.analysis_position=static_cast<std::int64_t>(center_position);
    round1.policy=policy;
    round1.sigma=x11_sigma_weights(repo_view(pass.b4d),frequency,static_cast<int>(center_position)+1,
        policy,0,0,false,lower_sigma,upper_sigma,mult,interrupt);
    round1.replacements=x11_replacements(repo_view(pass.b3),frequency,round1.sigma.weights,
        0,0,false,mult,interrupt);
    pass.b4=MaskedSeries{round1.replacements.replacements,
        missing_of(round1.replacements.replacements)};
    pass.b4g=MaskedSeries{round1.replacements.corrected,
        missing_of(round1.replacements.corrected)};
    checkpoint(interrupt);
    const auto initial5=process_pinned(pass.b4g,frequency,center_position,initial_filters,interrupt);
    pass.b5=normalize_pinned(initial5,frequency,center_position,0,t,final_filters,mult,interrupt);
    // X11BStepPseudoAdd.b6: ratio B1/B5 at the borders, interior
    // B1.drop - B2*(B5.drop - 1) (removeSeasonalAdjustedWithBorders).
    pass.b6=pseudo?pseudo_with_borders(pass.b1,pass.b5,pass.b2,t,interrupt)
        :remove_pinned(pass.b1,pass.b5,mult,interrupt);
    checkpoint(interrupt);
    pass.b7=pinned_trend_filter(pass.b6,frequency,trend_length,interrupt);
    // The pinned b7Step positivity repair of the trend table under a
    // multiplicative context (the one X11BStep makePositivity site).
    if(mult) make_positivity(pass.b7.values);
    pass.b8=remove_pinned(pass.b1,pass.b7,mult,interrupt);
    checkpoint(interrupt);
    const auto final9=process_pinned(pass.b8,frequency,initial_position,final_filters,interrupt);
    const auto seasonal9=normalize_pinned(final9,frequency,initial_position,0,0,final_filters,
        mult,interrupt);
    auto& round2=pass.rounds[1];
    // X11BStepPseudoAdd.b9d: removeIrregular(B8, s9) = B8 - s9 + 1.
    round2.detection=pseudo?pseudo_remove_irregular(pass.b8,seasonal9,interrupt)
        :remove_pinned(pass.b8,seasonal9,mult,interrupt);
    round2.analysis_position=static_cast<std::int64_t>(initial_position);
    round2.policy=policy;
    round2.sigma=x11_sigma_weights(repo_view(round2.detection),frequency,
        static_cast<int>(initial_position)+1,policy,0,0,false,lower_sigma,upper_sigma,mult,
        interrupt);
    round2.replacements=x11_replacements(repo_view(pass.b8),frequency,round2.sigma.weights,
        0,0,false,mult,interrupt);
    pass.b9=MaskedSeries{round2.replacements.replacements,
        missing_of(round2.replacements.replacements)};
    pass.b9g=MaskedSeries{round2.replacements.corrected,
        missing_of(round2.replacements.corrected)};
    checkpoint(interrupt);
    const auto final10=process_pinned(pass.b9g,frequency,initial_position,final_filters,interrupt);
    pass.b10=normalize_pinned(final10,frequency,initial_position,0,0,final_filters,mult,interrupt);
    // X11BStepPseudoAdd.b11: removeSeasonalAdjusted(B1, B10, B7)
    // = B1 - B7*(B10 - 1).
    pass.b11=pseudo?pseudo_remove_seasonal_adjusted(pass.b1,pass.b10,pass.b7,interrupt)
        :remove_pinned(pass.b1,pass.b10,mult,interrupt);
    pass.b13=remove_pinned(pass.b11,pass.b7,mult,interrupt);
    checkpoint(interrupt);
    auto& round3=pass.rounds[2];
    round3.detection=pass.b13;
    round3.analysis_position=static_cast<std::int64_t>(initial_position);
    round3.policy=policy;
    round3.sigma=x11_sigma_weights(repo_view(pass.b13),frequency,
        static_cast<int>(initial_position)+1,policy,0,0,false,lower_sigma,upper_sigma,mult,
        interrupt);
    round3.replacements=x11_replacements(repo_view(pass.b13),frequency,round3.sigma.weights,
        0,0,false,mult,interrupt);
    pass.b17=MaskedSeries{round3.sigma.weights,missing_of(round3.sigma.weights)};
    pass.b20=MaskedSeries{round3.replacements.correction_factors,
        missing_of(round3.replacements.correction_factors)};
    pass.status="available";
    checkpoint(interrupt);
    return pass;
}
} // namespace kalos
