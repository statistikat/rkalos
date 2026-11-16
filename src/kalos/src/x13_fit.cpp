// Item 22c M3-finale: the restricted fixed-model X-13 fit. Pure composition
// of the qualified seams in the pinned X13Kernel order — initialStep, the
// X11Kernel selection run over out with the updateSpec-forced mode, the
// log-additive results-builder stage where the mode is log_additive, and
// finals — with the unsupported-domain guards of docs/contracts/
// x13-fit-contract.md in front.
#include "kalos/x13_fit.hpp"
#include <cmath>
#include <stdexcept>
namespace kalos {
namespace {
// The pinned X13Kernel.updateSpec forcing (item 22c M3c contract): a non-log
// model runs X-11 additively; a log model keeps LogAdditive and
// PseudoAdditive and becomes Multiplicative for anything else.
X11DecompositionMode forceMode(bool log, X11DecompositionMode mode) {
    if(!log) return X11DecompositionMode::additive;
    if(mode==X11DecompositionMode::log_additive
        ||mode==X11DecompositionMode::pseudo_additive) return mode;
    return X11DecompositionMode::multiplicative;
}
// The pinned X11Kernel results-builder stage of the LogAdditive run
// (finalResults -> ratioBiasCorrection, the X11Spec default bias), the
// tables finals consumes. With A=exp(d10step), B=exp(d12step),
// C=exp(d13step), v9=average(A), v13=v9*average(C) over the run's whole
// series (the X11Context period; the whole calendar under the zero-horizon
// contract), the stored tables are D10=A/v9, D12=makePositivity(B*v13) —
// the identity here: exp tables and averages of positive values are
// strictly positive, so that branch is unreachable through this path —
// D11=input/D10 and D13=(input/D10)/D12 with input the X-11 series. The
// other BiasCorrection options are not delivered.
void logAdditiveResults(const std::vector<double>& step10,
    const std::vector<double>& step12, const std::vector<double>& step13,
    const std::vector<double>& input, std::vector<double>& d10,
    std::vector<double>& d11, std::vector<double>& d12, std::vector<double>& d13) {
    const std::size_t n=step10.size();
    std::vector<double> A(n),B(n),C(n);
    double sa=0.,sc=0.;
    for(std::size_t i=0;i<n;++i){
        A[i]=std::exp(step10[i]); B[i]=std::exp(step12[i]); C[i]=std::exp(step13[i]);
        sa+=A[i]; sc+=C[i];
    }
    const double v9=sa/static_cast<double>(n),
        v13=v9*(sc/static_cast<double>(n));
    d10.resize(n); d11.resize(n); d12.resize(n); d13.resize(n);
    for(std::size_t i=0;i<n;++i){
        d10[i]=A[i]/v9;
        d12[i]=B[i]*v13;
        d11[i]=input[i]/d10[i];
        d13[i]=(input[i]/d10[i])/d12[i];
    }
}
} // namespace
X13FitResult x13_fit(const X13FitRequest& request) {
    X13FitResult out;
    // Provenance echoes precede every guard so a refused request stays
    // fully auditable.
    out.frequency=request.frequency; out.first_period=request.first_period;
    out.log=request.log; out.seasonal=request.seasonal;
    out.trend_length=request.trend_length;
    out.sigma_set=request.sigma.has_value();
    if(request.sigma) out.sigma=*request.sigma;
    out.lower_sigma=request.lower_sigma; out.upper_sigma=request.upper_sigma;
    out.msr_policy=request.msr_policy;
    out.mode_requested=request.mode;
    out.mode_effective=forceMode(request.log,request.mode);
    out.backcast=request.backcast; out.forecast=request.forecast;
    // Shape of the request itself: one seasonal request per calendar period
    // (the pinned X11Spec carries one filter per period).
    if(request.seasonal.size()!=static_cast<std::size_t>(request.frequency<0?0:request.frequency))
        throw std::invalid_argument("x13 fit seasonal requests have one entry per calendar period");
    // Unsupported-domain guards, in the pinned order the choices appear in
    // the request: the fixed seasonal contract, the Henderson length, the
    // pinned fixed-path horizon domain {0, 0}, the required sigma policy.
    for(const X11SeasonalRequest r:request.seasonal)
        if(r==X11SeasonalRequest::x11default||r==X11SeasonalRequest::msr){
            out.status="unsupported:x13 fit fixed seasonal requests only";
            return out;
        }
    if(request.trend_length==0){
        out.status="unsupported:x13 fit fixed Henderson trend length";
        return out;
    }
    if(request.backcast!=0||request.forecast!=0){
        // The pinned kernel extends out and a1 through ExactArimaForecasts
        // for any nonzero horizon (negatives become |h| * annual frequency);
        // that route (a1a/a1b, the extended out calendar, the d11a/d11b
        // split finals) is outside the delivered surface.
        out.status="unsupported:x13 fit requires zero horizons";
        return out;
    }
    if(!request.sigma){
        out.status="unsupported:x13 fit requires an explicit sigma policy";
        return out;
    }
    // initialStep: the fixed-coefficient preadjust. Shape violations throw
    // here (and in the seams below), as in the allocation contract.
    X13PreadjustInput pi;
    pi.input=request.input; pi.variables=request.variables; pi.log=request.log;
    out.preadjustment=x13_preadjust(pi);
    // The X11Kernel run over out with the forced mode and the fixed spec.
    MaskedSeries series;
    series.values=out.preadjustment.out;
    series.missing.assign(out.preadjustment.out.size(),0);
    out.selection=x11_selection_pass(series,request.frequency,request.first_period,
        request.seasonal,request.trend_length,*request.sigma,request.msr_policy,
        request.lower_sigma,request.upper_sigma,0,0,out.mode_effective);
    if(out.selection.status!="available"){
        out.status=out.selection.status; // verbatim, stage records kept
        return out;
    }
    // The results tables finals consumes: the raw step tables outside
    // log_additive, the results-builder rescale inside it.
    if(out.mode_effective==X11DecompositionMode::log_additive)
        logAdditiveResults(out.selection.d.d10.values,out.selection.d.d12.values,
            out.selection.d.d13.values,out.preadjustment.out,
            out.d10,out.d11,out.d12,out.d13);
    else {
        out.d10=out.selection.d.d10.values; out.d11=out.selection.d.d11.values;
        out.d12=out.selection.d.d12.values; out.d13=out.selection.d.d13.values;
    }
    out.c17=out.selection.c.c17.values;
    // finals.
    X13RestoreInput ri;
    ri.pre=out.preadjustment;
    ri.d10=out.d10; ri.d11=out.d11; ri.d12=out.d12; ri.d13=out.d13;
    ri.c17=out.c17; ri.mode=out.mode_effective;
    out.finals=x13_restore(ri);
    out.finalSA=out.finals.d11final;
    out.finalTrend=out.finals.d12final;
    out.finalIrregular=out.finals.d13final;
    out.status="available";
    return out;
}
} // namespace kalos
