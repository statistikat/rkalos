#ifndef KALOS_DETAIL_X11_PINNED_HPP
#define KALOS_DETAIL_X11_PINNED_HPP
// Internal (detail) pinned-sequence machinery shared by the fixed additive
// B, C and D passes of items 22b.4 S1/S2/S3. These functions reproduce the
// pinned JDemetra+ 3.9.0 X11 arithmetic row for row — the qualified public
// primitives of items 22a/22b.2 document independent endpoint and mean
// conventions and are deliberately not reused for the seasonal phases or the
// Henderson endpoints. Nothing here is public API or stable across releases;
// see docs/contracts/x11-bpass-contract.md, x11-cpass-contract.md and
// x11-dpass-contract.md for the qualification evidence.
#include "../numerics.hpp"
#include "../parallel.hpp"
#include "../x11_filters.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace kalos {
namespace detail {
inline void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }
inline double nan() { return std::numeric_limits<double>::quiet_NaN(); }
// MusgraveFilterFactory.findR: pinned I/C ratio lookup driving the asymmetric
// Henderson endpoints (fixed filters only; quarterly is the four-period case).
inline double pinned_find_r(int length,int period) {
    if(period==4) return length<=5 ? .001 : 4.5;
    if(length<=9) return 1.;
    return length<=13 ? 3.5 : 4.5;
}
inline std::vector<std::uint8_t> missing_of(const std::vector<double>& values) {
    std::vector<std::uint8_t> missing(values.size(),0);
    for(std::size_t row=0;row<values.size();++row) missing[row]=std::isnan(values[row]);
    return missing;
}
// Public MaskedSeries storage keeps finite placeholders with the missing flag
// (numerics.hpp), while the pinned reference carries NaN through exactly the
// arithmetic below. The passes work in the pinned representation and convert
// only where a qualified primitive validates or resolves its own missing
// semantics.
inline MaskedSeries to_pinned(const MaskedSeries& repo) {
    MaskedSeries out=repo;
    for(std::size_t row=0;row<out.values.size();++row)
        if(out.missing[row]) out.values[row]=nan();
    return out;
}
inline MaskedSeries repo_view(const MaskedSeries& pinned) {
    MaskedSeries out=pinned;
    for(std::size_t row=0;row<out.values.size();++row)
        if(!std::isfinite(out.values[row])) {out.values[row]=0.;out.missing[row]=1;}
    return out;
}
// Elementwise difference mirroring the pinned IEEE arithmetic: the missing
// flag is the union of the inputs and masked values carry NaN, exactly the
// NaN the pinned elementwise subtraction produces from NaN-carrying rows.
inline MaskedSeries minus_pinned(const MaskedSeries& left,const MaskedSeries& right,
        const Interrupt& interrupt) {
    MaskedSeries out;
    out.values.assign(left.values.size(),nan());
    out.missing.assign(left.values.size(),1);
    parallel_for(left.values.size(),1,interrupt,[&](std::size_t row) {
        if(!left.missing[row] && !right.missing[row]) {
            const auto value=left.values[row]-right.values[row];
            out.values[row]=value; out.missing[row]=std::isnan(value);
        }
    });
    return out;
}
// DefaultSeasonalNormalizer.normalize(seq, drop, ctx, shift) with the pinned
// coordinates spelled out: `position` is ctx.getPosition(drop), the pinned
// 0-based calendar phase driving the endpoint restoration test, and `shift`
// is the fourth argument offsetting which filter entry governs each phase.
// The restoration set holds the phases whose final request is stable, whose
// span has fewer than five complete cycles, or whose 3x15 span covers fewer
// than twenty calendar years. `extend` embeds the normalized sequence into
// the periodic copy-back of CopyPeriodicEndPoints.
inline MaskedSeries normalize_pinned(const MaskedSeries& seq,int frequency,
        std::size_t position,std::size_t shift,std::size_t extend,
        const std::vector<int>& final_filters,const Interrupt& interrupt) {
    const auto f=static_cast<std::size_t>(frequency),m=seq.values.size();
    const auto t=f/2;
    if(m<f+t) throw std::invalid_argument("normalization span requires a full cycle plus the annual radius");
    const auto annual=x11_symmetric_filter(repo_view(seq),X11FilterKind::trend,frequency,0,interrupt);
    MaskedSeries out;
    out.values.assign(m,nan()); out.missing.assign(m,1);
    for(std::size_t row=t;row<m-t;++row) { out.missing[row]=annual.missing[row];
        out.values[row]=annual.missing[row]?nan():annual.values[row]; }
    // CopyEndPoints: constant extension of the annual band's outer values.
    for(std::size_t row=0;row<t;++row) {
        out.values[row]=out.values[t]; out.missing[row]=out.missing[t];
        out.values[m-1-row]=out.values[m-t-1]; out.missing[m-1-row]=out.missing[m-t-1];
    }
    std::vector<std::uint8_t> restore(f,0);
    const auto cycles=m/f,years=cycles+(m%f?1:0);
    for(std::size_t phase=0;phase<f;++phase) {
        const auto entry=final_filters[(phase+shift)%f]; // pinned (shift+i)%f ids
        restore[phase]=entry==0 || cycles<5 || (entry==15 && years<20);
    }
    for(std::size_t row=0;row<t;++row)
        if(restore[(position+row)%f]) { out.values[row]=out.values[row+f]; out.missing[row]=out.missing[row+f]; }
    for(std::size_t row=m-t;row<m;++row)
        if(restore[(position+row)%f]) { out.values[row]=out.values[row-f]; out.missing[row]=out.missing[row-f]; }
    auto normalized=minus_pinned(seq,out,interrupt);
    if(extend==0) return normalized;
    MaskedSeries padded;
    padded.values.assign(m+2*extend,nan()); padded.missing.assign(m+2*extend,1);
    for(std::size_t row=0;row<m;++row) {
        padded.values[row+extend]=normalized.values[row];
        padded.missing[row+extend]=normalized.missing[row];
    }
    // CopyPeriodicEndPoints in the pinned descending order.
    for(std::int64_t step=static_cast<std::int64_t>(extend)-1;step>=0;--step) {
        padded.values[static_cast<std::size_t>(step)]=
            padded.values[static_cast<std::size_t>(step)+f];
        padded.missing[static_cast<std::size_t>(step)]=
            padded.missing[static_cast<std::size_t>(step)+f];
        const auto right=padded.values.size()-extend+static_cast<std::size_t>(step);
        padded.values[right]=padded.values[right-f];
        padded.missing[right]=padded.missing[right-f];
    }
    return padded;
}
inline void check_filters(const std::vector<int>& filters,std::size_t f,const char* what) {
    if(filters.size()!=f)
        throw std::invalid_argument(std::string(what)+": one fixed seasonal request per calendar period");
    for(const auto choice:filters)
        if(choice!=0 && choice!=1 && choice!=3 && choice!=5 && choice!=9 && choice!=15)
            throw std::invalid_argument(std::string(what)+": unsupported fixed seasonal request");
}
// The passes reproduce the pinned reference's fixed seasonal processing
// arithmetic row for row, which the public x11_seasonal_process primitive
// deliberately does not (its contract documents the independent endpoint and
// mean conventions). Pinned X11SeasonalFiltersFactory static tables: the
// classic asymmetric seasonal endpoint weights, indexed per option and
// available-side points, plus the stable/short-span gates.
inline const std::vector<double>& fc_entry(int option,std::size_t index) {
    static const std::vector<std::vector<double>> tables[]={{
        {{.39,.61}}},                                        // 3x1
        {{.1111111111111111,.25925925925925924,.37037037037037035,.25925925925925924},
         {.18518518518518517,.4074074074074074,.4074074074074074}},                    // 3x3
        {{.06666666666666667,.13333333333333333,.21666666666666667,.21666666666666667,.21666666666666667,.15},
         {.06666666666666667,.18333333333333332,.25,.25,.25},
         {.15,.2833333333333333,.2833333333333333,.2833333333333333}},                 // 3x5
        {{.034,.073,.111,.113,.114,.116,.117,.118,.12,.084},
         {.034,.075,.113,.117,.123,.128,.132,.137,.141},
         {.032,.079,.123,.133,.143,.154,.163,.173},
         {.028,.092,.144,.16,.176,.192,.208},
         {.051,.112,.173,.197,.221,.246}},                                             // 3x9
        {{.0222,.04444,.06667,.06667,.06667,.06667,.06667,.06667,.06667,.06667,.06667,.07111,.07111,.07111,.07111,.04889},
         {.0222,.04444,.06667,.06667,.06667,.06667,.06667,.06667,.06667,.06667,.08,.08,.08,.08,.08},
         {.02222,.04444,.06667,.06667,.06667,.06667,.06667,.06667,.06667,.09333,.09333,.09333,.09333,.09333},
         {.02219,.04444,.06667,.06667,.06667,.06667,.06667,.06667,.10667,.10667,.10667,.10667,.10667},
         {.02221,.04444,.06667,.06667,.06667,.06667,.06667,.12,.12,.12,.12,.12},
         {.02223,.04444,.06667,.06667,.06667,.06667,.13333,.13333,.13333,.13333,.13333},
         {.0222,.04444,.06667,.06667,.06667,.14667,.14667,.14667,.14667,.14667},
         {.02222,.04444,.06667,.06667,.16,.16,.16,.16,.16}}                            // 3x15
    };
    static const int slots[]={0,0,0,1,0,2,0,0,0,3,0,0,0,0,0,4}; // option -> table (1,3,5,9,15)
    const auto&table=tables[slots[option]];
    if(index>=table.size())throw std::logic_error("pinned seasonal endpoint table index");
    return table[index];
}
inline std::size_t fc_count(int option) { return static_cast<std::size_t>((option+2)/2); }
// DataBlock.average(): running sum of the finite values divided by their
// count; NaN only when no value is finite.
inline double block_average(const std::vector<double>& block) {
    double sum=0;std::size_t skipped=0;
    for(const auto value:block)if(std::isfinite(value))sum+=value;else ++skipped;
    if(skipped==block.size())return nan();
    return sum/static_cast<double>(block.size()-skipped);
}
// IFiniteFilter.apply on a window: anchored dot product in ascending tap
// order over the full coefficient array.
inline double anchored_dot(const std::vector<double>& w,const double* window) {
    double value=0;
    for(std::size_t t=0;t<w.size();++t)value+=w[t]*window[t];
    return value;
}
// Arithmetics.sumOfPowers for the orders used by the pinned Henderson kernel.
inline double pinned_sum_powers(int power,long n) {
    const double d=static_cast<double>(n);
    if(power==2) return (d*(d+1.))*(2.*d+1.)/6.;
    if(power==4) return (d*(d+1.))*(2.*d+1.)*((3.*d*d)+3.*d-1.)/30.;
    const double d2=d*d,d3=d2*d,d4=d3*d;
    return (d*(d+1.))*(2.*d+1.)*((((3.*d4)+6.*d3)-3.*d)+1.)/42.;
}
// DiscreteKernel.henderson(radius) plus LocalPolynomialFilters.of2_3: the
// pinned X11Context.trendFilter weights, rebuilt tap for tap in the pinned
// double arithmetic (the classic closed form differs by up to 16 ULP).
// Ascending tap order, center at index radius.
inline std::vector<double> pinned_trend_weights(int length) {
    const auto p=static_cast<std::size_t>(length/2);
    const double a=static_cast<double>(p+1),a2=a*a;
    const double b=static_cast<double>(p+2),b2=b*b;
    const double c=static_cast<double>(p+3),c2=c*c;
    const double s2=pinned_sum_powers(2,static_cast<long>(p));
    const double s4=pinned_sum_powers(4,static_cast<long>(p));
    const double s6=pinned_sum_powers(6,static_cast<long>(p));
    double v=a2*b2; v=v*c2; v=v*static_cast<double>(1+2*static_cast<long>(p));
    double s=a2*b2; s=s+a2*c2; s=s+b2*c2; s=s*s2;
    double q=a2+b2; q=q+c2; q=q*s4;
    const double den=v-2.*((s-q)+s6);
    std::vector<double> kernel(p+1),half(p+1);
    for(std::size_t j=0;j<=p;++j) {
        const double j2=static_cast<double>(j*j);
        double num=a2-j2; num=num*(b2-j2); num=num*(c2-j2);
        kernel[j]=num/den;
    }
    double s1=0,m2=0,m4=0;
    for(std::size_t j=1;j<=p;++j) {
        const double j2=static_cast<double>(j*j);
        s1=s1+kernel[j]; m2=m2+j2*kernel[j]; m4=m4+j2*j2*kernel[j];
    }
    s1=kernel[0]+2.*s1; m2=2.*m2; m4=2.*m4;
    const double curvature=m2/m4,den3=s1-m2*curvature;
    for(std::size_t j=0;j<=p;++j) {
        const double j2=static_cast<double>(j*j);
        half[j]=kernel[j]*(1.-j2*curvature)/den3;
    }
    std::vector<double> w(static_cast<std::size_t>(length));
    for(std::size_t j=0;j<=p;++j) {
        w[p+j]=half[j];
        w[p-j]=half[j];
    }
    return w;
}
// AsymmetricFiltersFactory.mmsreFilter2(tf, side, 0, {dz}, null): the pinned
// endpoint correction of the trend filter. The 1x1 systems reduce to plain
// division (verified against the toolkit's LULinearSystemSolver bit for bit).
// Result has m = p+side+1 taps anchored at lb = -p.
inline std::vector<double> pinned_mmsre_filter(const std::vector<double>& w,int p,int side,
        double dz) {
    const auto m=static_cast<std::size_t>(p+side+1);
    const auto n=static_cast<std::size_t>(2*p+1);
    std::vector<double> out(w.begin(),w.begin()+static_cast<std::ptrdiff_t>(m));
    double xsum=0;
    for(std::size_t t=m;t<n;++t)xsum=xsum+w[t];
    const double x=xsum/static_cast<double>(m);
    // dot(z12, x*1), sign-flipped, plus the whole dot(z13, w[m..n)) added
    // once (DataBlock.chs/addProduct order, not an interleaved accumulation).
    double left=0;
    for(std::size_t r=0;r<m;++r)left=left+static_cast<double>(static_cast<std::ptrdiff_t>(r)-p)*x;
    double moment=-left,right=0;
    for(std::size_t t=m;t<n;++t)right=right+static_cast<double>(static_cast<std::ptrdiff_t>(t)-p)*w[t];
    moment=moment+right;
    std::vector<double> slope(m);
    double s23=0;
    for(std::size_t r=0;r<m;++r) {
        slope[r]=static_cast<double>(static_cast<std::ptrdiff_t>(r)-p)*dz;
        s23=s23+slope[r];
    }
    const double mean=s23/static_cast<double>(m);
    double v27=0;
    for(std::size_t r=0;r<m;++r)
        v27=v27+static_cast<double>(static_cast<std::ptrdiff_t>(r)-p)*(slope[r]-mean);
    const double v30=(moment/(v27*dz+1.))*dz;
    for(std::size_t r=0;r<m;++r) {
        out[r]=out[r]+((slope[r]-mean)*v30);
        out[r]=out[r]+x;
    }
    return out;
}
// X11Context.asymmetricTrendFilters(tf, findR) consumed by AsymmetricEndPoints
// (both sides): interior rows use the IFiniteFilter slow-path ascending dot;
// endpoint rows use the mmsre table stored at array index p-1-side.
inline MaskedSeries pinned_trend_filter(const MaskedSeries& base,int frequency,int length,
        const Interrupt& interrupt) {
    const auto n=base.values.size();
    const auto p=static_cast<std::size_t>(length/2);
    const auto trend=pinned_trend_weights(length);
    const double dz=2./(std::sqrt(3.141592653589793)*pinned_find_r(length,frequency));
    std::vector<std::vector<double>> table(p);
    for(std::size_t s=0;s<p;++s)
        table[s]=pinned_mmsre_filter(trend,static_cast<int>(p),
            static_cast<int>(p-1-s),dz);
    MaskedSeries out;
    out.values.assign(n,nan()); out.missing.assign(n,1);
    parallel_for(n-2*p,static_cast<std::size_t>(length),interrupt,[&](std::size_t offset) {
        checkpoint(interrupt);
        const auto row=p+offset;
        const auto value=anchored_dot(trend,base.values.data()+row-p);
        out.values[row]=value; out.missing[row]=std::isnan(value);
    });
    const auto half=(n+1)/2;
    for(int side=0;side<2;++side) {
        auto avail=std::max<std::int64_t>(0,2*static_cast<std::int64_t>(p)-static_cast<std::int64_t>(n));
        auto count=static_cast<std::int64_t>(p)-avail;
        if(count>static_cast<std::int64_t>(half)) {count=static_cast<std::int64_t>(half);avail=static_cast<std::int64_t>(p)-static_cast<std::int64_t>(half);}
        const auto span=static_cast<std::size_t>(2*static_cast<std::int64_t>(p)-avail);
        if(side==0) {
            std::vector<double> reversed(span);
            for(std::size_t t=0;t<span;++t)reversed[t]=base.values[span-1-t];
            for(auto i=count,shift=static_cast<std::int64_t>(0);i>0;--i,++shift) {
                const auto value=anchored_dot(table[static_cast<std::size_t>(avail)+static_cast<std::size_t>(shift)],
                    reversed.data()+static_cast<std::size_t>(shift));
                out.values[static_cast<std::size_t>(i)-1]=value;
                out.missing[static_cast<std::size_t>(i)-1]=std::isnan(value);
            }
        } else {
            const auto tail=n-span;
            for(auto i=count,shift=static_cast<std::int64_t>(0);i>0;--i,++shift) {
                const auto value=anchored_dot(table[static_cast<std::size_t>(avail)+static_cast<std::size_t>(shift)],
                    base.values.data()+tail+static_cast<std::size_t>(shift));
                out.values[n-static_cast<std::size_t>(i)]=value;
                out.missing[n-static_cast<std::size_t>(i)]=std::isnan(value);
            }
        }
    }
    return out;
}
// X11SeasonalFiltersFactory.DefaultFilter.process plus AsymmetricEndPoints
// with the pinned position 0 (both sides always run).
inline void process_phase_pinned(const std::vector<double>& block,int option,int frequency,
        std::vector<double>& out,std::size_t cycles,std::size_t years) {
    const auto m=block.size(),p=fc_count(option);
    std::fill(out.begin(),out.end(),0.);
    if(option==0 || !(cycles>=5 && !(p>=8 && years<20))) { // StableFilter gate
        const auto mean=block_average(block);
        std::fill(out.begin(),out.end(),mean);
        return;
    }
    const auto definition=x11_filter_definition(X11FilterKind::seasonal,frequency,option);
    // weights are stored in ascending offset order with the center at p.
    const auto&w=definition.weights;
    if(2*p<m)for(std::size_t b=p;b<m-p;++b) {
        auto value=block[b]*w[p];
        for(std::size_t j=1;j<=p;++j)value+=w[p+j]*(block[b-j]+block[b+j]);
        out[b]=value;
    }
    const auto half=(m+1)/2;
    auto fill=2*static_cast<std::int64_t>(p)-static_cast<std::int64_t>(m);
    if(fill<0)fill=0;
    auto k=static_cast<std::int64_t>(p)-fill;
    if(k>static_cast<std::int64_t>(half))k=static_cast<std::int64_t>(half);
    // The StableFilter gate above guarantees a non-negative pinned k.
    if(k<static_cast<std::int64_t>(p)&&k>=0)
        for(auto t=static_cast<std::size_t>(k);t<m-static_cast<std::size_t>(k);++t)
            out[t]=block_average(block);
    // AsymmetricEndPoints.processLeft: reversed prefix window, filters from
    // the available count upward, one oldest point dropped per step.
    auto avail=std::max<std::int64_t>(0,2*static_cast<std::int64_t>(p)-static_cast<std::int64_t>(m));
    auto count=static_cast<std::int64_t>(p)-avail;
    if(count>static_cast<std::int64_t>(half)) {count=static_cast<std::int64_t>(half);avail=static_cast<std::int64_t>(p)-static_cast<std::int64_t>(half);}
    const auto span=static_cast<std::size_t>(2*static_cast<std::int64_t>(p)-avail);
    std::vector<double> reversed(span);
    for(std::size_t t=0;t<span;++t)reversed[t]=block[span-1-t];
    for(auto i=count,shift=static_cast<std::int64_t>(0);i>0;--i,++shift)
        out[static_cast<std::size_t>(i)-1]=anchored_dot(fc_entry(option,static_cast<std::size_t>(avail)+static_cast<std::size_t>(shift)),reversed.data()+static_cast<std::size_t>(shift));
    // processRight: forward suffix window, same filter progression.
    avail=std::max<std::int64_t>(0,2*static_cast<std::int64_t>(p)-static_cast<std::int64_t>(m));
    count=static_cast<std::int64_t>(p)-avail;
    if(count>static_cast<std::int64_t>(half)) {count=static_cast<std::int64_t>(half);avail=static_cast<std::int64_t>(p)-static_cast<std::int64_t>(half);}
    const auto tail=m-static_cast<std::size_t>(2*static_cast<std::int64_t>(p)-avail);
    for(auto i=count,shift=static_cast<std::int64_t>(0);i>0;--i,++shift)
        out[m-static_cast<std::size_t>(i)]=anchored_dot(fc_entry(option,static_cast<std::size_t>(avail)+static_cast<std::size_t>(shift)),block.data()+tail+static_cast<std::size_t>(shift));
}
// X11SeasonalFilterProcessor.process on the fixed per-phase options: phase i
// uses option [(position+i)%f] over the stride block starting at row i.
inline MaskedSeries process_pinned(const MaskedSeries& seq,int frequency,std::size_t position,
        const std::vector<int>& options,const Interrupt& interrupt) {
    const auto f=static_cast<std::size_t>(frequency),n=seq.values.size();
    const auto cycles=n/f,years=cycles+(n%f?1:0);
    MaskedSeries out;out.values.assign(n,nan());out.missing.assign(n,1);
    parallel_for(f,n,interrupt,[&](std::size_t i) {
        checkpoint(interrupt);
        const std::size_t length=1+(n-1-i)/f;
        std::vector<double> block(length),result(length);
        for(std::size_t j=0;j<length;++j)block[j]=seq.values[i+j*f];
        process_phase_pinned(block,options[(position+i)%f],static_cast<int>(f),result,cycles,years);
        for(std::size_t j=0;j<length;++j) {
            out.values[i+j*f]=result[j];
            out.missing[i+j*f]=std::isnan(result[j]);
        }
    });
    return out;
}
// X11FilterFactory.makeSymmetricFilter(n): lag-indexed weights w[0..p],
// lag 0 being the center; all weights are 1/n except the outer lag, which is
// 0.5/n for even n. For odd n the factory's zero outer tap is trimmed by
// SymmetricFilter (makeSymmetricFilter(7) has length 7), so p=(n-1)/2 there;
// the even case has p=n/2. SymmetricFilter.apply(DataBlock,DataBlock)
// accumulates in the pinned pair form over the radius p:
// sum = w[0]x_i + SUM_{j=1..p} w[j](x_{i-j}+x_{i+j}).
inline std::vector<double> pinned_ma_lags(int n) {
    const auto p=static_cast<std::size_t>(n+1)/2;
    std::vector<double> w(p+1);
    const double c=1./static_cast<double>(n);
    for(std::size_t j=0;j<p;++j)w[j]=c;
    w[p]=(n%2==0?0.5:0.)/static_cast<double>(n);
    if(n%2!=0) w.pop_back(); // SymmetricFilter drops the zero outer tap
    return w;
}
// DoubleSeq.average(): plain sum over all values divided by the length.
// Unlike DataBlock.average() this is NOT NaN-skipping: one missing value
// poisons the mean, which is what the D-pass MSRE table pins.
inline double pinned_plain_average(const std::vector<double>& values) {
    double sum=0;
    for(const auto value:values)sum=sum+value;
    return sum/static_cast<double>(values.size());
}
// FilteredMeanEndPoints(ma).process over one stride block: both ends padded
// with the plain averages of the first/last p block values, then the pinned
// pair-form accumulation over the padded window (requires block.size()>=p).
inline void pinned_filtered_mean(const std::vector<double>& block,
        const std::vector<double>& w,std::vector<double>& out) {
    const auto p=w.size()-1,m=block.size();
    std::vector<double> padded(m+2*p,nan());
    for(std::size_t j=0;j<m;++j)padded[p+j]=block[j];
    const std::vector<double> head(block.begin(),block.begin()+p);
    const std::vector<double> tail(block.end()-p,block.end());
    const double left=pinned_plain_average(head),right=pinned_plain_average(tail);
    for(std::size_t j=0;j<p;++j) { padded[j]=left; padded[m+p+j]=right; }
    out.assign(m,nan());
    for(std::size_t k=0;k<m;++k) {
        const auto i=p+k;
        auto value=padded[i]*w[0];
        for(std::size_t j=1;j<=p;++j)
            value=value+w[j]*(padded[i-j]+padded[i+j]);
        out[k]=value;
    }
}
// X11Utility.defaultMsrTable in additive mode over the span `x`: the pinned
// MsrTable difference statistics, indexed per calendar phase (phase numbers
// are absolute; `position` is the phase of the first row). The additive mode
// takes the pinned plain-subtract branch of the decomposition switch — the
// residual is x-mean with no +1 offset (that offset belongs to the other
// decomposition modes and would quantize tiny residuals around 1.0) — and
// the seasonal evolution is the differences of the filtered mean itself;
// both accumulate |difference| in stride order and scale by the pinned
// fis/cs lookup of MsrTable.of.
inline void pinned_msr_table(const std::vector<double>& x,int frequency,
        std::size_t position,std::vector<double>& ri,std::vector<double>& rs,
        std::vector<std::int64_t>& counts) {
    const auto f=static_cast<std::size_t>(frequency),m=x.size();
    static const double ctab[]={1.,1.02584,1.01779,1.01383,1.,3.,1.55291,1.30095};
    const auto fis=[](long n)->double {
        if(n<2)return 1.;
        if(n<6)return ctab[static_cast<std::size_t>(n)-2];
        const double d=static_cast<double>(n);
        return 12.247449*d/(73.239334+12.247449*static_cast<double>(n-6));
    };
    const auto cs=[](long n)->double {
        if(n<2)return 1.;
        if(n<6)return ctab[static_cast<std::size_t>(n)+2];
        const double d=static_cast<double>(n);
        return 1.732051*d/(8.485281+1.732051*static_cast<double>(n-6));
    };
    const auto w=pinned_ma_lags(7); // p=3 (zero outer tap trimmed)
    std::vector<double> mean(m,nan());
    for(std::size_t phase=0;phase<f;++phase) { // stride blocks by row phase
        const std::size_t count=1+(m-1-phase)/f;
        std::vector<double> block(count),filtered;
        for(std::size_t j=0;j<count;++j)block[j]=x[phase+j*f];
        pinned_filtered_mean(block,w,filtered);
        for(std::size_t j=0;j<count;++j)mean[phase+j*f]=filtered[j];
    }
    std::vector<double> residual(m);
    for(std::size_t row=0;row<m;++row)residual[row]=x[row]-mean[row];
    ri.assign(f,nan());rs.assign(f,nan());counts.assign(f,0);
    for(std::size_t phase=0;phase<f;++phase) {
        const auto start=(phase+f-position)%f; // pinned (phase-position+f)%f
        double s_ri=0,s_rs=0;
        long n=0;
        // The pinned MsrTable.of takes the MEAN series (not x) against the
        // residual: seasonal evolutions measure the smoothed series.
        for(std::size_t i1=start,i2=start+f;i2<m;i1=i2,i2+=f,++n) {
            s_ri=s_ri+std::fabs(residual[i2]-residual[i1]);
            s_rs=s_rs+std::fabs(mean[i2]-mean[i1]);
        }
        ri[phase]=(s_ri/static_cast<double>(n))*fis(n);
        rs[phase]=(s_rs/static_cast<double>(n))*cs(n);
        counts[phase]=n;
    }
}
// AutomaticHenderson.calcICR over one complete window: the pilot Henderson
// filter of length period + 1 with the pinned symmetric weights applied to
// the interior span only (the pinned NaN-filled table plus DataBlock.of view
// leaves the outer `side` rows unwritten, so they never contribute), the
// irregular as the aligned interior input minus the pilot, the explicit
// extensions dropped from both windows, and the ratio of the lag-1 mean
// absolute changes of the irregular over the same statistic of the pilot,
// multiplied by 3.0 for quarterly and 6.0 for biannual (pinned constants).
// calcAbsMeanVariation divides the plain sum by length - 1 with no missing
// handling: degenerate dropped windows reproduce the pinned division
// outcomes (NaN from 0/0, signed zero from 0/-1) unchanged.
inline double pinned_ic_ratio(const std::vector<double>& values,int frequency,
        std::size_t backcast,std::size_t forecast) {
    const auto f=static_cast<std::size_t>(frequency),n=values.size();
    const auto length=static_cast<std::size_t>(frequency)+1; // odd: period is even
    const auto side=length/2;
    const auto w=pinned_trend_weights(static_cast<int>(length));
    const auto m=n>2*side?n-2*side:0;
    std::vector<double> pilot(m),irregular(m);
    for(std::size_t row=0;row<m;++row) {
        pilot[row]=anchored_dot(w,values.data()+row);
        irregular[row]=values[side+row]-pilot[row];
    }
    // DoubleSeq.drop(backcast, forecast) = dropFirst then dropLast.
    const auto first=std::min(backcast,m);
    const auto rest=m-first;
    const auto count=rest>forecast?rest-forecast:0;
    double pilot_change=0,irregular_change=0;
    for(std::size_t row=1;row<count;++row) {
        pilot_change=pilot_change+std::fabs(pilot[first+row]-pilot[first+row-1]);
        irregular_change=irregular_change+std::fabs(irregular[first+row]-irregular[first+row-1]);
    }
    const auto divisor=static_cast<double>(count)-1.; // calcAbsMeanVariation's n - 1
    auto ratio=(irregular_change/divisor)/(pilot_change/divisor);
    if(f==4)ratio=ratio*3.;
    else if(f==2)ratio=ratio*6.;
    return ratio;
}
// The D9 mask pair of X11DStep.d9Step with the pinned dcmpg semantics: a NaN
// difference keeps the D9bis value.
inline MaskedSeries d9_mask_pinned(const MaskedSeries& d9bis,const MaskedSeries& d8,
        const Interrupt& interrupt) {
    const auto n=d9bis.values.size();
    MaskedSeries out;out.values.assign(n,nan());out.missing.assign(n,1);
    parallel_for(n,1,interrupt,[&](std::size_t row) {
        const auto difference=std::fabs(d9bis.values[row]-d8.values[row]);
        if(std::isnan(difference)||difference>=1e-9) {
            out.values[row]=d9bis.values[row];
            out.missing[row]=std::isnan(d9bis.values[row]);
        }
    });
    return out;
}
} // namespace detail
} // namespace kalos
#endif
