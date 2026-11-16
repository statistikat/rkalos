#include "kalos/x11_selection_pass.hpp"
#include "kalos/detail/x11_pinned.hpp"
#include "kalos/x11_trend_selection.hpp"
#include <algorithm>
#include <stdexcept>
namespace kalos {
namespace {
using namespace detail;
// MsrFilterSelection.doMSR over the D9bis series: drop the explicit
// extensions, trim one trailing partial calendar year at the original
// phase (completeYear), then evaluate the pinned MsrTable over the leading
// window while at least five calendar years remain, removing one year per
// ambiguous band. The qualified item 22b.3c x11_msr_band rule decides and
// the qualified iteration policies stop the loop; an undefined (NaN) table
// selects nothing and marks the stage unsupported, exactly like the
// qualified contract refuses to reproduce the pinned NaN-to-3x9 control.
struct MsrWalk {
    int filter=5;
    std::size_t iterations=0;
    bool using_default=true;
    double global_msr=0.;
    bool unsupported=false;
    X11MsrSelectionStatus status=X11MsrSelectionStatus::available;
    X11MsrSelectionReason reason=X11MsrSelectionReason::unavailable;
    std::vector<X11MsrAttempt> attempts;
};
MsrWalk pinned_msr_walk(const MaskedSeries& d9bis,int frequency,int first_period,
        int backcast,int forecast,X11MsrIterationPolicy policy,
        X11DecompositionMode mode,const Interrupt& interrupt) {
    MsrWalk walk;
    const auto f=static_cast<std::size_t>(frequency);
    const auto back=static_cast<std::size_t>(backcast),fc=static_cast<std::size_t>(forecast);
    std::vector<double> w(d9bis.values.begin()+static_cast<std::ptrdiff_t>(back),
        d9bis.values.end()-static_cast<std::ptrdiff_t>(fc)); // DoubleSeq.drop
    const auto first=(static_cast<std::size_t>(first_period)-1+back)%f; // getPosition(backcast)
    const auto rem=(w.size()+first)%f; // one trailing partial year
    if(rem) w.resize(w.size()-rem);
    std::size_t off=0;
    while(true) {
        checkpoint(interrupt);
        if((w.size()-off)/f<5) { // fewer than five years remain (or never held)
            walk.reason=off==0?X11MsrSelectionReason::insufficient_initial_history
                :X11MsrSelectionReason::insufficient_remaining_history;
            return walk; // method default 3x5, using_default stays true
        }
        std::vector<double> ri,rs;
        std::vector<std::int64_t> counts;
        pinned_msr_table(std::vector<double>(w.begin(),w.begin()+w.size()-off),
            frequency,first,mode,ri,rs,counts); // position unchanged
        double numerator=0,denominator=0;
        for(std::size_t phase=0;phase<f;++phase) {
            numerator=numerator+ri[phase]*static_cast<double>(counts[phase]);
            denominator=denominator+rs[phase]*static_cast<double>(counts[phase]);
        }
        walk.global_msr=numerator/denominator;
        ++walk.iterations;
        X11MsrAttempt row;
        row.count=w.size()-off;
        row.years_removed=off/f;
        row.global_msr=walk.global_msr;
        row.status=X11MsrSelectionStatus::available;
        try {
            row.band=x11_msr_band(walk.global_msr);
        } catch(const std::invalid_argument&) {
            row.band=X11MsrBand::unavailable;
            row.status=X11MsrSelectionStatus::undefined_msr;
            walk.attempts.push_back(row);
            walk.unsupported=true;
            walk.status=X11MsrSelectionStatus::undefined_msr;
            return walk;
        }
        walk.attempts.push_back(row);
        int decided=0;
        switch(row.band) {
        case X11MsrBand::s3x3: decided=3; break;
        case X11MsrBand::s3x5: decided=5; break;
        case X11MsrBand::s3x9: decided=9; break;
        default: break;
        }
        if(decided!=0) {
            walk.filter=decided;
            walk.using_default=false;
            walk.reason=X11MsrSelectionReason::decisive_ratio;
            return walk;
        }
        if(policy==X11MsrIterationPolicy::x11_five_year_cap&&off/f>=5) {
            walk.reason=X11MsrSelectionReason::five_year_cap;
            return walk; // method default 3x5 after six evaluations
        }
        off+=f; // one year removed, leading edge
    }
}
// The pinned b7Step inline branch: ratio >= 1 takes period + 1, anything
// else falls through to AutomaticHenderson.selectFilter (a NaN ratio fails
// the >= 1 test and the qualified rule then refuses it: explicit
// unsupported, the pinned longest-filter control is not reproduced).
int pinned_b7_length(int frequency,double ic_ratio) {
    if(ic_ratio>=1.) return frequency+1;
    return x11_henderson_select(frequency,ic_ratio);
}
} // namespace

X11SelectionPass x11_selection_pass(const MaskedSeries& series,int frequency,int first_period,
        const std::vector<X11SeasonalRequest>& requests,int trend_length,
        X11SigmaPolicy policy,X11MsrIterationPolicy msr_policy,
        double lower_sigma,double upper_sigma,int backcast,int forecast,
        X11DecompositionMode mode,const Interrupt& interrupt) {
    series.validate();
    check_decomposition("selection pass",mode,series);
    const auto mult=x11_mode_is_multiplicative(mode);
    const auto pseudo=mode==X11DecompositionMode::pseudo_additive;
    if(frequency!=4 && frequency!=12)
        throw std::invalid_argument("selection pass supports monthly and quarterly series");
    if(first_period<1 || first_period>frequency)
        throw std::invalid_argument("selection pass first period is within the annual frequency");
    if(trend_length!=0 && (trend_length<5 || trend_length%2==0))
        throw std::invalid_argument("selection pass trend filter is 0 (automatic) or an odd Henderson length >= 5");
    if(msr_policy!=X11MsrIterationPolicy::jdplus && msr_policy!=X11MsrIterationPolicy::x11_five_year_cap)
        throw std::invalid_argument("selection pass MSR iteration policy is unsupported");
    for(const auto value:series.missing)
        if(value) throw std::invalid_argument("selection pass requires complete observations");
    const auto f=static_cast<std::size_t>(frequency),n=series.values.size(),t=f/2;
    const auto first=static_cast<std::size_t>(first_period)-1;
    if(backcast<0 || forecast<0)
        throw std::invalid_argument("selection pass horizons are non-negative");
    if(trend_length!=0 && n<static_cast<std::size_t>(trend_length))
        throw std::invalid_argument("selection pass requires at least the trend filter length");
    if(n<3*t+f) throw std::invalid_argument("selection pass requires three annual radii plus one cycle");
    if(static_cast<std::int64_t>(n)-backcast-forecast<static_cast<std::int64_t>(4*f))
        throw std::invalid_argument("selection pass requires four cycles of original calendar");
    checkpoint(interrupt);
    X11SelectionPass pass;
    // The qualified item 22b.3d plans resolve every pinned stage. Their
    // base arrays are the pinned context arrays: Msr and X11Default take
    // 3x3 at an initial stage and 3x5 at a final stage; a final all-MSR
    // request stays provisional until the D final MSR decision fills it.
    const auto b_initial=x11_seasonal_plan(frequency,requests,X11SeasonalStage::b_initial),
        b_final=x11_seasonal_plan(frequency,requests,X11SeasonalStage::b_final),
        c_initial=x11_seasonal_plan(frequency,requests,X11SeasonalStage::c_initial),
        c_final=x11_seasonal_plan(frequency,requests,X11SeasonalStage::c_final),
        d_initial=x11_seasonal_plan(frequency,requests,X11SeasonalStage::d_initial),
        d_final=x11_seasonal_plan(frequency,requests,X11SeasonalStage::d_final);
    pass.initial_plan=d_initial; pass.final_plan=d_final;
    pass.initial_filters=d_initial.base_filters; pass.final_filters=d_final.base_filters;
    const auto unsupported=[&](X11SelectionPass& out,const char* stage) {
        out.status=std::string("unsupported:")+stage;
        return out;
    };
    // B stage. b1..b6 are selection-independent, so the automatic case
    // probes them with the smallest Henderson length, selects B7 from the
    // pinned I/C ratio over B6 and reruns the finished B pass at the
    // selected length. The pinned inline branch is applied verbatim.
    const auto automatic=trend_length==0;
    const auto probe=automatic?5:trend_length;
    auto b=x11_b_pass(series,frequency,first_period,b_initial.base_filters,
        b_final.base_filters,probe,policy,lower_sigma,upper_sigma,mode,interrupt);
    if(automatic) {
        pass.b7.status=X11TrendStageStatus::selected;
        pass.b7.ic_ratio=pinned_ic_ratio(b.b6.values,frequency,
            static_cast<std::size_t>(backcast),static_cast<std::size_t>(forecast),mult);
        try {
            pass.b7.length=pinned_b7_length(frequency,pass.b7.ic_ratio);
        } catch(const std::invalid_argument&) {
            pass.b7.status=X11TrendStageStatus::unsupported;
            return unsupported(pass,"trend-b7");
        }
        if(pass.b7.length!=probe)
            b=x11_b_pass(series,frequency,first_period,b_initial.base_filters,
                b_final.base_filters,pass.b7.length,policy,lower_sigma,upper_sigma,mode,
                interrupt);
    } else {
        pass.b7.status=X11TrendStageStatus::fixed;
        pass.b7.length=trend_length; // the pinned computes no ratio in fixed mode
    }
    pass.b=b;
    checkpoint(interrupt);
    // C stage. The composition mirrors x11_c_pass.cpp; the automatic C7
    // length may differ from B7 (the pinned kernel never requires them to
    // agree), so the C tables are composed here rather than chained
    // through the fixed C pass, whose validation enforces one length.
    X11CPass c;
    c.mode=mode;
    c.c2drop=t; c.first_period=first_period;
    c.initial_filters=c_initial.base_filters; c.final_filters=c_final.base_filters;
    // The pinned kernel runs the whole pipeline on its (possibly log-
    // transformed) input, so the manual C composition consumes the very
    // sequence the steps see — log(y) under LogAdditive (item 22b.5 M3b-i);
    // the internal fixed-pass calls keep the raw series because those
    // entries apply the same kernel input transform themselves.
    const auto y=to_pinned(mode==X11DecompositionMode::log_additive
        ?log_domain(series):series);
    c.reference_series=y;
    c.c1=pseudo?pseudo_adjust_ref_series(y,b.b20,b.b7,b.b13,interrupt) // X11CStepPseudoAdd
        :remove_pinned(y,b.b20,mult,interrupt);
    const auto annual_c=x11_symmetric_filter(repo_view(c.c1),X11FilterKind::trend,frequency,0,interrupt);
    c.c2=to_pinned(MaskedSeries{std::vector<double>(annual_c.values.begin()+t,annual_c.values.begin()+n-t),
        std::vector<std::uint8_t>(annual_c.missing.begin()+t,annual_c.missing.begin()+n-t)});
    MaskedSeries cropped_c;
    cropped_c.values.assign(c.c1.values.begin()+t,c.c1.values.begin()+n-t);
    cropped_c.missing.assign(c.c1.missing.begin()+t,c.c1.missing.begin()+n-t);
    c.c4=remove_pinned(cropped_c,c.c2,mult,interrupt);
    const auto center=(first+t)%f; // pinned getPosition(c2drop)
    checkpoint(interrupt);
    c.c5a=process_pinned(c.c4,frequency,center,c_initial.base_filters,interrupt);
    c.c5=normalize_pinned(c.c5a,frequency,center,0,t,c_final.base_filters,mult,interrupt);
    c.c6=pseudo?pseudo_with_borders(c.c1,c.c5,c.c2,t,interrupt) // X11CStepPseudoAdd
        :remove_pinned(c.c1,c.c5,mult,interrupt);
    checkpoint(interrupt);
    if(automatic) {
        pass.c7.status=X11TrendStageStatus::selected;
        pass.c7.ic_ratio=pinned_ic_ratio(c.c6.values,frequency,
            static_cast<std::size_t>(backcast),static_cast<std::size_t>(forecast),mult);
        try {
            pass.c7.length=x11_henderson_select(frequency,pass.c7.ic_ratio);
        } catch(const std::invalid_argument&) {
            pass.c7.status=X11TrendStageStatus::unsupported;
            return unsupported(pass,"trend-c7");
        }
    } else {
        pass.c7.status=X11TrendStageStatus::fixed;
        pass.c7.length=trend_length;
    }
    if(n<static_cast<std::size_t>(pass.c7.length))
        throw std::invalid_argument("selection pass requires at least the selected C7 trend length");
    c.trend_length=pass.c7.length;
    c.c7=pinned_trend_filter(c.c6,frequency,pass.c7.length,interrupt);
    if(mult) make_positivity(c.c7.values); // pinned c7Step repair
    c.c9=remove_pinned(c.c1,c.c7,mult,interrupt);
    checkpoint(interrupt);
    c.c10a=process_pinned(c.c9,frequency,first,c_final.base_filters,interrupt);
    c.c10=normalize_pinned(c.c10a,frequency,first,0,0,c_final.base_filters,mult,interrupt);
    c.c11=pseudo?pseudo_remove_seasonal_adjusted(y,c.c10,c.c7,interrupt) // X11CStepPseudoAdd
        :remove_pinned(y,c.c10,mult,interrupt); // remove(refSeries, c10)
    c.c13=remove_pinned(c.c11,c.c7,mult,interrupt);
    checkpoint(interrupt);
    auto& round=c.round;
    round.detection=c.c13;
    round.analysis_position=static_cast<std::int64_t>(first);
    round.policy=policy;
    round.sigma=x11_sigma_weights(repo_view(c.c13),frequency,
        static_cast<int>(first)+1,policy,0,0,false,lower_sigma,upper_sigma,mult,interrupt);
    round.replacements=x11_replacements(repo_view(c.c13),frequency,round.sigma.weights,
        0,0,false,mult,interrupt);
    c.c17=MaskedSeries{round.sigma.weights,missing_of(round.sigma.weights)};
    c.c20=MaskedSeries{round.replacements.correction_factors,
        missing_of(round.replacements.correction_factors)};
    c.status="available";
    pass.c=c;
    checkpoint(interrupt);
    // D stage. The composition mirrors x11_d_pass.cpp with the pinned
    // selection points: the D7 length from D6, the global-MSR D-final
    // seasonal decision filling a D-local copy of the final requests, and
    // the final Henderson length from D11bis (the iCRatio is always
    // computed, also in fixed mode).
    X11DPass d;
    d.mode=mode;
    d.d2drop=t;
    d.backcast=backcast; d.forecast=forecast;
    d.initial_filters=d_initial.base_filters; d.final_filters=d_final.base_filters;
    d.d1=pseudo?pseudo_adjust_ref_series(y,c.c20,c.c7,c.c13,interrupt) // X11DStepPseudoAdd
        :remove_pinned(y,c.c20,mult,interrupt);
    const auto annual_d=x11_symmetric_filter(repo_view(d.d1),X11FilterKind::trend,frequency,0,interrupt);
    d.d2=to_pinned(MaskedSeries{std::vector<double>(annual_d.values.begin()+t,annual_d.values.begin()+n-t),
        std::vector<std::uint8_t>(annual_d.missing.begin()+t,annual_d.missing.begin()+n-t)});
    MaskedSeries cropped_d;
    cropped_d.values.assign(d.d1.values.begin()+t,d.d1.values.begin()+n-t);
    cropped_d.missing.assign(d.d1.missing.begin()+t,d.d1.missing.begin()+n-t);
    d.d4=remove_pinned(cropped_d,d.d2,mult,interrupt);
    checkpoint(interrupt);
    const auto d5raw=process_pinned(d.d4,frequency,center,d_initial.base_filters,interrupt);
    d.d5=normalize_pinned(d5raw,frequency,center,0,t,d_final.base_filters,mult,interrupt);
    d.d6=pseudo?pseudo_with_borders(d.d1,d.d5,d.d2,t,interrupt) // X11DStepPseudoAdd
        :remove_pinned(d.d1,d.d5,mult,interrupt);
    checkpoint(interrupt);
    if(automatic) {
        pass.d7.status=X11TrendStageStatus::selected;
        pass.d7.ic_ratio=pinned_ic_ratio(d.d6.values,frequency,
            static_cast<std::size_t>(backcast),static_cast<std::size_t>(forecast),mult);
        try {
            pass.d7.length=x11_henderson_select(frequency,pass.d7.ic_ratio);
        } catch(const std::invalid_argument&) {
            pass.d7.status=X11TrendStageStatus::unsupported;
            return unsupported(pass,"trend-d7");
        }
    } else {
        pass.d7.status=X11TrendStageStatus::fixed;
        pass.d7.length=trend_length;
    }
    if(n<static_cast<std::size_t>(pass.d7.length))
        throw std::invalid_argument("selection pass requires at least the selected D7 trend length");
    d.d7=pinned_trend_filter(d.d6,frequency,pass.d7.length,interrupt);
    if(mult) make_positivity(d.d7.values); // pinned d7Step repair
    d.d8=pseudo?product_pinned(remove_pinned(d.d1,d.d7,true,interrupt),
        c.c20,interrupt) // X11DStepPseudoAdd: (D1/D7)*C20
        :remove_pinned(y,d.d7,mult,interrupt); // remove(refSeries, d7)
    d.d9bis=remove_pinned(d.d1,d.d7,mult,interrupt);
    d.d9=d9_mask_pinned(d.d9bis,d.d8,interrupt); // pinned dcmpg mask
    checkpoint(interrupt);
    // The pinned dFinalStep: the seasFilter copy of the final requests,
    // filled with the global-MSR decision when the context is all-MSR.
    auto seas=d_final.base_filters;
    pass.msr.engaged=d_final.requires_msr;
    pass.msr.policy=msr_policy;
    if(pass.msr.engaged) {
        const auto walk=pinned_msr_walk(d.d9bis,frequency,first_period,backcast,forecast,
            msr_policy,mode,interrupt);
        pass.msr.iterations=walk.iterations;
        pass.msr.global_msr=walk.global_msr;
        pass.msr.attempts=walk.attempts;
        pass.msr.status=walk.status;
        pass.msr.reason=walk.reason;
        pass.msr.using_default=walk.using_default;
        if(walk.unsupported) return unsupported(pass,"msr-d-final");
        pass.msr.selected_filter=walk.filter;
        std::fill(seas.begin(),seas.end(),walk.filter);
    } else {
        pass.msr.status=X11MsrSelectionStatus::available;
        pass.msr.reason=X11MsrSelectionReason::unavailable;
    }
    d.d10bis=process_pinned(d.d9bis,frequency,first,seas,interrupt);
    // The pinned normalization calls ctx.getFinalSeasonalFilter(), which
    // re-derives the provisional array (Msr and X11Default mapping to 3x5):
    // the seasFilter copy filled by doMSR is a private D-step array (the
    // getter returns a fresh array) and never reaches the normalizer.
    d.d10=normalize_pinned(d.d10bis,frequency,first,0,0,d_final.base_filters,mult,interrupt);
    d.d11bis=pseudo?pseudo_remove_seasonal_adjusted(d.d1,d.d10,d.d7,interrupt) // X11DStepPseudoAdd
        :remove_pinned(d.d1,d.d10,mult,interrupt);
    checkpoint(interrupt);
    pass.final.ic_ratio=pinned_ic_ratio(d.d11bis.values,frequency,
        static_cast<std::size_t>(backcast),static_cast<std::size_t>(forecast),mult); // always stored
    if(automatic) {
        pass.final.status=X11TrendStageStatus::selected;
        try {
            pass.final.length=x11_henderson_select(frequency,pass.final.ic_ratio);
        } catch(const std::invalid_argument&) {
            pass.final.status=X11TrendStageStatus::unsupported;
            return unsupported(pass,"trend-final");
        }
    } else {
        pass.final.status=X11TrendStageStatus::fixed;
        pass.final.length=trend_length;
    }
    if(n<static_cast<std::size_t>(pass.final.length))
        throw std::invalid_argument("selection pass requires at least the selected final trend length");
    d.d12=pinned_trend_filter(d.d11bis,frequency,pass.final.length,interrupt);
    if(mult) make_positivity(d.d12.values); // pinned dFinalStep repair
    d.d11=pseudo?pseudo_remove_seasonal_adjusted(y,d.d10,d.d12,interrupt) // X11DStepPseudoAdd
        :remove_pinned(y,d.d10,mult,interrupt); // remove(refSeries, d10)
    d.d13=remove_pinned(d.d11,d.d12,mult,interrupt);
    checkpoint(interrupt);
    // The pinned MSRE diagnostic runs on the dropped original calendar at
    // the original-calendar phase getPosition(backcast).
    std::vector<double> original(d.d9bis.values.begin()+static_cast<std::ptrdiff_t>(backcast),
        d.d9bis.values.end()-static_cast<std::ptrdiff_t>(forecast));
    std::vector<double> ri,rs;
    std::vector<std::int64_t> counts;
    pinned_msr_table(original,frequency,(first+static_cast<std::size_t>(backcast))%f,
        mode,ri,rs,counts);
    auto& msr=d.msr_table;
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
    d.trend_length=pass.final.length;
    d.d7_henderson_length=pass.d7.length;
    d.final_henderson_length=pass.final.length;
    d.seas_filter=seas;
    d.msr_iterations=pass.msr.engaged?static_cast<int>(pass.msr.iterations):0;
    d.status="available";
    pass.d=d;
    pass.status="available";
    checkpoint(interrupt);
    return pass;
}
} // namespace kalos
