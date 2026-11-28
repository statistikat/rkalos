// P4 (X-11 diagnostic battery): the pinned Mstatistics battery over a
// finished X13FitResult. Every loop below mirrors a pinned bytecode body,
// operation order included; see the scope header for the input seam.
#include "kalos/x11_mstats.hpp"

#include "kalos/combined_seasonality.hpp"
#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/x13_fit.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace kalos {
namespace {

// X11Utility.calcAbsMeanVariation(seq, lag, mult): |x[i] - x[i-lag]| taken
// BEFORE the multiplicative divide (the pinned operand order, kept exact).
double amv_single(const std::vector<double>& x,int lag,bool mult) {
    const auto n=x.size();
    double sum=0.;
    for(std::size_t i=static_cast<std::size_t>(lag);i<n;++i) {
        double a=std::fabs(x[i]-x[i-lag]);
        if(mult)a/=x[i-lag];
        sum+=a;
    }
    return sum/static_cast<double>(static_cast<long long>(n)-lag);
}

// X11Utility.calcAbsMeanVariations(seq, maxperiod, mult): per-lag table,
// d = x[i] - x[i-lag] divided (when mult) BEFORE the absolute value.
std::vector<double> amv_table(const std::vector<double>& x,int maxperiod,bool mult) {
    const auto n=x.size();
    std::vector<double> g(static_cast<std::size_t>(maxperiod));
    for(int lag=1;lag<=maxperiod;++lag) {
        double sum=0.;
        for(std::size_t i=static_cast<std::size_t>(lag);i<n;++i) {
            double d=x[i]-x[i-lag];
            if(mult)d/=x[i-lag];
            sum+=std::fabs(d);
        }
        g[static_cast<std::size_t>(lag)-1]
            =sum/static_cast<double>(static_cast<long long>(n)-lag);
    }
    return g;
}

// X11Utility.adr(seq, mult): (n-1) / the count of nonzero-sign runs of the
// consecutive differences (relative differences when mult). Zero differences
// keep the previous sign, the pinned transition machine.
double adr(const std::vector<double>& x,bool mult) {
    const auto n=x.size();
    if(n<2)return 0.;
    const auto m=static_cast<std::size_t>(n)-1;
    std::vector<double> d(m);
    for(std::size_t i=0;i<m;++i) {
        d[i]=x[i+1]-x[i];
        if(mult)d[i]/=x[i];
    }
    std::size_t count=0;
    int last=0;
    for(std::size_t i=0;i<m;++i) {
        const int s=d[i]>0.?1:(d[i]<0.?-1:0);
        if(s!=0&&s!=last) {
            ++count;
            last=s;
        }
    }
    return static_cast<double>(m)/static_cast<double>(count);
}

// Mstatistics.variance(seq, mult, zeroMean): the log transform (pinned
// Math.log) when mult, the mean over the finite observations UNLESS
// zeroMean, and the raw sum of squared deviations over the finite
// observations (no division).
double variance(const std::vector<double>& source,bool mult,bool zeroMean) {
    std::vector<double> x=source;
    if(mult)for(double& v:x)v=detail::fdlibm_log(v);
    double mean=0.;
    if(!zeroMean) {
        double sum=0.;
        std::size_t count=0;
        for(double v:x)if(std::isfinite(v)) {
            sum+=v;
            ++count;
        }
        mean=sum/static_cast<double>(count);
    }
    double ssq=0.;
    for(double v:x)if(std::isfinite(v)) {
        const double d=v-mean;
        ssq+=d*d;
    }
    return ssq;
}

std::string grade_of(double value) {
    if(std::isnan(value))return "Undefined";
    if(value>=2.)return "Severe"; // the pinned 2.0 severe threshold
    if(value>=1.)return "Bad";    // the pinned 1.0 bad threshold
    return "Good";
}

} // namespace

X11Mstats x13_mstats(const X13FitResult& fit) {
    if(fit.status!="available")
        throw std::invalid_argument("x13 mstats requires an available fit, got '"+fit.status+"'");
    X11Mstats out;
    out.mode=fit.mode_effective;
    if(fit.mode_effective==X11DecompositionMode::pseudo_additive) {
        out.status="unsupported:pseudo_additive";
        return out;
    }
    const int nf=fit.frequency;
    const bool mult=fit.mode_effective!=X11DecompositionMode::additive;
    const auto& pre=fit.preadjustment;
    const auto n=pre.a1.size();
    if(n==0||nf<2||fit.first_period<1||fit.first_period>nf||fit.backcast<0||fit.forecast<0)
        throw std::invalid_argument("x13 mstats needs a nonempty monthly-or-quarterly fit calendar");
    const auto back=static_cast<std::size_t>(fit.backcast);
    const auto fore=static_cast<std::size_t>(fit.forecast);
    const auto ext=n+back+fore;
    const auto slice=[&](const std::vector<double>& v,const char* what)->std::vector<double> {
        if(v.size()!=ext)
            throw std::invalid_argument(std::string("x13 mstats table ")+what
                +" does not span the fit calendar");
        return std::vector<double>(v.begin()+static_cast<std::ptrdiff_t>(back),
            v.begin()+static_cast<std::ptrdiff_t>(back+n));
    };
    std::vector<double> si_source;
    {
        const auto& d8=fit.selection.d.d8;
        if(d8.values.size()!=ext||d8.missing.size()!=d8.values.size())
            throw std::invalid_argument("x13 mstats D8 does not span the fit calendar");
        si_source.resize(n);
        for(std::size_t j=0;j<n;++j) {
            if(d8.missing[back+j])
                throw std::invalid_argument("x13 mstats D8 carries a missing row");
            // The pinned results builder stores log_additive tables exp'd;
            // calcM7 takes SI.log(), so the battery sees log(exp(value)).
            const double raw=d8.values[back+j];
            si_source[j]=fit.mode_effective==X11DecompositionMode::log_additive
                ?detail::fdlibm_log(detail::fdlibm_exp(raw)):raw;
        }
    }
    const auto Sc=slice(fit.d10,"D10"), Cc=slice(fit.d12,"D12");
    slice(fit.d11,"D11"); // CIc: consumed by the pinned class, not the battery
    const auto Ic=slice(fit.d13,"D13");
    const auto A8t=slice(pre.a8t,"A8t"), A8s=slice(pre.a8s,"A8s");
    (void)A8t;(void)A8s; // Omod/CImod inputs: outside the published battery
    if(fit.finals.e1.size()!=n||fit.finals.e3.size()!=n||fit.finals.d18.size()!=n)
        throw std::invalid_argument("x13 mstats finals tables must sit on the a1 calendar");
    const std::vector<double>& Ome=fit.finals.e1;
    const std::vector<double>& Imod=fit.finals.e3;
    std::vector<double> TD; // the pinned TD field, null when constant
    {
        bool constant=true;
        for(std::size_t i=1;i<n;++i)
            if(!(fit.finals.d18[i]==fit.finals.d18[0])) {constant=false;break;}
        if(!constant)TD=fit.finals.d18;
    }
    const std::vector<double>& P=pre.a8; // the whole extended table
    if(P.size()!=ext)
        throw std::invalid_argument("x13 mstats A8 must span the extended calendar");

    // ---- of(): the D9 record, its presets and bShort ----------------------
    const auto& msr=fit.selection.msr;
    const auto& msrdef=fit.selection.d.msr_table; // the pinned d9msr: X11Utility
    // .defaultMsrTable on the A1-domain D9bis window, built unconditionally in
    // dFinalStep and distinct from the year-removal selection walk (whose last
    // attempt is selection.msr.global_msr). Mstatistics reads getD9Msr(), so
    // the M6 preset is fed by the default table.
    //
    // The M6-in-q gate is s3x5 = !X11Results.isD9default() && d9filter==S3X5.
    // isD9default() reads the X11DStep.d9default field, which is declared and
    // read by the getter but never assigned anywhere in the pinned class (no
    // putfield) — so it is always false; the walk's isUsingDefault() is a
    // separate, unused value. The gate therefore reduces to the fit's D9
    // filter equalling S3X5 (verified: ao reports isD9default=false with the
    // walk's isUsingDefault=true).
    const bool s3x5=msr.engaged&&msr.selected_filter==5;
    if(!msrdef.mean_irregular_evolutions.empty())
        out.raw[5]=0.4*std::fabs(msrdef.global_msr-4.);
    const double icr=fit.selection.final.ic_ratio;
    out.raw[2]=0.5*(icr-1.);
    const bool bShort=n/static_cast<std::size_t>(nf)<6;

    // ---- calcStationaryVariances(): the [1 | 0..n-1] trend partition ------
    std::vector<double> y(n),ones(n),t(n),reg(n);
    for(std::size_t i=0;i<n;++i) {
        y[i]=mult?detail::fdlibm_log(Cc[i]):Cc[i];
        ones[i]=1.;
        t[i]=static_cast<double>(i);
    }
    const auto ols=detail::pinned_ols_fit(y,{ones,t});
    if(!ols) {
        out.status="failed:trend_regression";
        return out;
    }
    for(std::size_t i=0;i<n;++i)
        reg[i]=ols->coefficients[0]+ols->coefficients[1]*t[i]; // LSR.regressionEffect
    std::vector<double> stC(n),stO(n);
    for(std::size_t i=0;i<n;++i) {
        stC[i]=y[i]-reg[i];
        stO[i]=(mult?detail::fdlibm_log(Ome[i]):Ome[i])-reg[i];
    }
    double varO=variance(stO,false,false);
    double varC=variance(stC,false,false);
    double varS=variance(Sc,mult,true);
    double varI=variance(Imod,mult,true);
    double varP=variance(P,mult,false);
    double varTD=TD.empty()?0.:variance(TD,mult,true);
    varP/=varO;varTD/=varO;varS/=varO;varC/=varO;varI/=varO;

    // ---- calcSNorm(): Sc / rms(Sc (minus 1 when mult)) --------------------
    std::vector<double> snorm(n);
    {
        double ssq=0.;
        for(std::size_t i=0;i<n;++i) {
            const double v=mult?Sc[i]-1.:Sc[i];
            ssq+=v*v;
        }
        const double inv=1./std::sqrt(ssq/static_cast<double>(n));
        for(std::size_t i=0;i<n;++i)snorm[i]=Sc[i]*inv;
    }

    // ---- calcEvolutions(): the M5 ratio table ------------------------------
    const auto gIc=amv_table(Ic,nf,mult), gCc=amv_table(Cc,nf,mult);
    const auto smic=[&](int k){return gIc[static_cast<std::size_t>(k)-1]
        /gCc[static_cast<std::size_t>(k)-1];};

    // ---- calcM(lag) ---------------------------------------------------------
    const int lag=nf==12?3:1;
    {
        const double vc=amv_single(Cc,lag,mult),vi0=amv_single(Imod,lag,mult);
        const double vs=amv_single(Sc,lag,mult),vp=amv_single(P,lag,mult);
        const double vt=TD.empty()?0.:amv_single(TD,lag,mult);
        const double vc2=vc*vc,vi=vi0*vi0,vs2=vs*vs,vp2=vp*vp,vt2=vt*vt;
        const double tot=vc2+vi+vs2+vp2+vt2;
        out.raw[0]=10.*(vi/tot)/(1.-vp2/tot);
    }
    out.raw[1]=varP>=1.?3.:10.*varI/(1.-varP);
    {
        const double a=adr(Ic,mult);
        const int ni=static_cast<int>(n); // the pinned O.length() is an int
        const double den=2.577*std::sqrt(1.6*static_cast<double>(ni)-2.9);
        const double num=std::fabs(static_cast<double>(3*(ni-1))/a
            -static_cast<double>(2*ni-1));
        out.raw[3]=num/den;
    }
    {
        const double nfYear=static_cast<double>(12/nf);
        int m=nf;
        while(m>0&&smic(m)<1.)--m;
        if(m<nf)++m;
        double x;
        if(m==1) {
            x=1.+(smic(1)-1.)/(smic(1)-smic(2));
            if(x<0.5)x=0.5;
            if(x>1.)x=1.;
        } else {
            const double diff=smic(m-1)-smic(m);
            if(diff<=0.&&m==nf)x=nfYear*15.5; // the pinned dcmpg keeps NaN on the normal path
            else x=static_cast<double>(m)+(smic(m)-1.)/diff;
        }
        out.raw[4]=(x*nfYear-0.5)/5.;
    }
    {
        // calcM7: the CombinedSeasonality mvalue over the SI slice. The
        // qualified mirror needs three years; shorter calendars (where the
        // pinned class dereferences its null ANOVA) leave the word uncomputed.
        const auto cs=combined_si_seasonality(si_source,nf,fit.first_period,
            fit.mode_effective==X11DecompositionMode::multiplicative);
        if(cs) {
            const double fs=cs->stable_f,fe=cs->moving_f;
            double a=7./fs,b=fs==0.?9.:3.*fe/fs;
            if(a>9.)a=9.;
            if(b>9.)b=9.;
            out.raw[6]=std::sqrt(0.5*(a+b));
        }
    }
    if(!bShort) {
        {
            double sum=0.;
            for(std::size_t i=0;i+static_cast<std::size_t>(nf)<n;++i)
                sum+=std::fabs(snorm[i+static_cast<std::size_t>(nf)]-snorm[i]);
            out.raw[7]=10.*(sum/static_cast<double>(n-static_cast<std::size_t>(nf)));
        }
        {
            double sum=0.;
            for(int p=1;p<=nf;++p) {
                const auto i=static_cast<int>(n)-p;
                const int i12=i%nf;
                const int k=(i-i12)/nf;
                sum+=std::fabs(snorm[static_cast<std::size_t>(i)]
                    -snorm[static_cast<std::size_t>(i12)])/static_cast<double>(k);
            }
            out.raw[8]=sum/static_cast<double>(nf)*10.;
        }
        {
            double sum=0.;
            bool full=true;
            for(int j=1;j<=nf;++j) {
                const int i=static_cast<int>(n)-j-nf;
                const int k=i-3*nf;
                if(k<nf) {full=false;break;}
                for(int tpos=k;tpos<i;tpos+=nf)
                    sum+=std::fabs(snorm[static_cast<std::size_t>(tpos)]
                        -snorm[static_cast<std::size_t>(tpos-nf)]);
            }
            if(full)out.raw[9]=sum/static_cast<double>(3*nf)*10.;
        }
        {
            double sum=0.;
            bool full=true;
            for(int j=1;j<=nf;++j) {
                const int i=static_cast<int>(n)-j-2*nf;
                const int k=i-3*nf;
                if(k<0) {full=false;break;}
                sum+=std::fabs(snorm[static_cast<std::size_t>(i)]
                    -snorm[static_cast<std::size_t>(k)]);
            }
            if(full)out.raw[10]=sum/static_cast<double>(3*nf)*10.;
        }
    }

    // ---- publication: getM clamps, getQ/getQm2, the q/q2 grades ------------
    for(std::size_t k=0;k<out.raw.size();++k) {
        out.used[k]=out.raw[k]!=-1.;
        const double v=out.raw[k];
        out.published[k]=v>3.?3.:(v<0.?0.:v);
    }
    const auto aggregate=[&](bool skip_m2) {
        static const double wt_full[11]={10.,11.,10.,8.,11.,10.,18.,7.,7.,4.,4.};
        static const double wt_short[7]={14.,15.,10.,8.,11.,10.,32.};
        const double* wt=bShort?wt_short:wt_full;
        const std::size_t count=bShort?7:11;
        double num=0.,den=0.;
        for(std::size_t k=0;k<count;++k) {
            if(skip_m2&&k==1)continue;
            const double v=out.raw[k];
            if(v==-1.)continue;
            if(k==5&&!s3x5)continue;
            den+=wt[k];
            if(v>3.)num+=3.*wt[k];
            else if(v>0.)num+=v*wt[k];
        }
        return num/den;
    };
    out.q=aggregate(false);
    out.qm2=aggregate(true);
    out.q_grade=grade_of(out.q);
    out.qm2_grade=grade_of(out.qm2);
    out.status="available";
    return out;
}

} // namespace kalos
