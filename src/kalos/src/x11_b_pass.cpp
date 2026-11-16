#include "kalos/x11_b_pass.hpp"
#include "kalos/parallel.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace kalos {
namespace {
void checkpoint(const Interrupt& interrupt) { if (interrupt) interrupt(); }
double nan() { return std::numeric_limits<double>::quiet_NaN(); }
// MusgraveFilterFactory.findR: pinned I/C ratio lookup driving the asymmetric
// Henderson endpoints (fixed filters only; quarterly is the four-period case).
double pinned_find_r(int length,int period) {
    if(period==4) return length<=5 ? .001 : 4.5;
    if(length<=9) return 1.;
    return length<=13 ? 3.5 : 4.5;
}
std::vector<std::uint8_t> missing_of(const std::vector<double>& values) {
    std::vector<std::uint8_t> missing(values.size(),0);
    for(std::size_t row=0;row<values.size();++row) missing[row]=std::isnan(values[row]);
    return missing;
}
// Public MaskedSeries storage keeps finite placeholders with the missing flag
// (numerics.hpp), while the pinned reference carries NaN through exactly the
// arithmetic below. The B tables mirror the pinned NaN values, so the pass
// works in the pinned representation and converts only where a qualified
// primitive validates or resolves its own missing semantics.
MaskedSeries to_pinned(const MaskedSeries& repo) {
    MaskedSeries out=repo;
    for(std::size_t row=0;row<out.values.size();++row)
        if(out.missing[row]) out.values[row]=nan();
    return out;
}
MaskedSeries repo_view(const MaskedSeries& pinned) {
    MaskedSeries out=pinned;
    for(std::size_t row=0;row<out.values.size();++row)
        if(!std::isfinite(out.values[row])) {out.values[row]=0.;out.missing[row]=1;}
    return out;
}
// Elementwise difference mirroring the pinned IEEE arithmetic: the missing
// flag is the union of the inputs and masked values carry NaN, exactly the
// NaN the pinned elementwise subtraction produces from NaN-carrying rows.
MaskedSeries minus_pinned(const MaskedSeries& left,const MaskedSeries& right,
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
MaskedSeries normalize_pinned(const MaskedSeries& seq,int frequency,
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
void check_filters(const std::vector<int>& filters,std::size_t f,const char* what) {
    if(filters.size()!=f)
        throw std::invalid_argument(std::string(what)+": one fixed seasonal request per calendar period");
    for(const auto choice:filters)
        if(choice!=0 && choice!=1 && choice!=3 && choice!=5 && choice!=9 && choice!=15)
            throw std::invalid_argument(std::string(what)+": unsupported fixed seasonal request");
}
// The B pass reproduces the pinned reference's fixed seasonal processing
// arithmetic row for row, which the public x11_seasonal_process primitive
// deliberately does not (its contract documents the independent endpoint and
// mean conventions). Pinned X11SeasonalFiltersFactory static tables: the
// classic asymmetric seasonal endpoint weights, indexed per option and
// available-side points, plus the stable/short-span gates.
const std::vector<double>& fc_entry(int option,std::size_t index) {
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
std::size_t fc_count(int option) { return static_cast<std::size_t>((option+2)/2); }
// DataBlock.average(): running sum of the finite values divided by their
// count; NaN only when no value is finite.
double block_average(const std::vector<double>& block) {
    double sum=0;std::size_t skipped=0;
    for(const auto value:block)if(std::isfinite(value))sum+=value;else ++skipped;
    if(skipped==block.size())return nan();
    return sum/static_cast<double>(block.size()-skipped);
}
// IFiniteFilter.apply on a window: anchored dot product in ascending tap
// order over the full coefficient array.
double anchored_dot(const std::vector<double>& w,const double* window) {
    double value=0;
    for(std::size_t t=0;t<w.size();++t)value+=w[t]*window[t];
    return value;
}
// Arithmetics.sumOfPowers for the orders used by the pinned Henderson kernel.
double pinned_sum_powers(int power,long n) {
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
std::vector<double> pinned_trend_weights(int length) {
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
std::vector<double> pinned_mmsre_filter(const std::vector<double>& w,int p,int side,
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
MaskedSeries pinned_trend_filter(const MaskedSeries& b6,int frequency,int length,
        const Interrupt& interrupt) {
    const auto n=b6.values.size();
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
        const auto value=anchored_dot(trend,b6.values.data()+row-p);
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
            for(std::size_t t=0;t<span;++t)reversed[t]=b6.values[span-1-t];
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
                    b6.values.data()+tail+static_cast<std::size_t>(shift));
                out.values[n-static_cast<std::size_t>(i)]=value;
                out.missing[n-static_cast<std::size_t>(i)]=std::isnan(value);
            }
        }
    }
    return out;
}
// X11SeasonalFiltersFactory.DefaultFilter.process plus AsymmetricEndPoints
// with the pinned position 0 (both sides always run).
void process_phase_pinned(const std::vector<double>& block,int option,int frequency,
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
MaskedSeries process_pinned(const MaskedSeries& seq,int frequency,std::size_t position,
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
} // namespace

X11BPass x11_b_pass(const MaskedSeries& series,int frequency,int first_period,
        const std::vector<int>& initial_filters,const std::vector<int>& final_filters,
        int trend_length,X11SigmaPolicy policy,
        double lower_sigma,double upper_sigma,const Interrupt& interrupt) {
    series.validate();
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
    pass.initial_filters=initial_filters; pass.final_filters=final_filters;
    pass.trend_length=trend_length;
    pass.b1=to_pinned(series);
    const auto annual=x11_symmetric_filter(series,X11FilterKind::trend,frequency,0,interrupt);
    pass.b2=to_pinned(MaskedSeries{std::vector<double>(annual.values.begin()+t,annual.values.begin()+n-t),
        std::vector<std::uint8_t>(annual.missing.begin()+t,annual.missing.begin()+n-t)});
    MaskedSeries cropped;
    cropped.values.assign(series.values.begin()+t,series.values.begin()+n-t);
    cropped.missing.assign(series.missing.begin()+t,series.missing.begin()+n-t);
    pass.b3=minus_pinned(cropped,pass.b2,interrupt);
    const auto first=static_cast<std::size_t>(first_period-1);
    const auto center_position=(first+t)%f;   // pinned getPosition(b2drop)
    const std::size_t initial_position=first; // pinned getPosition(0)
    checkpoint(interrupt);
    pass.b4a=process_pinned(pass.b3,frequency,center_position,initial_filters,interrupt);
    pass.b4anorm=normalize_pinned(pass.b4a,frequency,initial_position,t,0,final_filters,interrupt);
    pass.b4d=minus_pinned(pass.b3,pass.b4anorm,interrupt);
    checkpoint(interrupt);
    auto& round1=pass.rounds[0];
    round1.detection=pass.b4d;
    round1.analysis_position=static_cast<std::int64_t>(center_position);
    round1.policy=policy;
    round1.sigma=x11_sigma_weights(repo_view(pass.b4d),frequency,static_cast<int>(center_position)+1,
        policy,0,0,false,lower_sigma,upper_sigma,false,interrupt);
    round1.replacements=x11_replacements(repo_view(pass.b3),frequency,round1.sigma.weights,
        0,0,false,false,interrupt);
    pass.b4=MaskedSeries{round1.replacements.replacements,
        missing_of(round1.replacements.replacements)};
    pass.b4g=MaskedSeries{round1.replacements.corrected,
        missing_of(round1.replacements.corrected)};
    checkpoint(interrupt);
    const auto initial5=process_pinned(pass.b4g,frequency,center_position,initial_filters,interrupt);
    pass.b5=normalize_pinned(initial5,frequency,center_position,0,t,final_filters,interrupt);
    pass.b6=minus_pinned(pass.b1,pass.b5,interrupt);
    checkpoint(interrupt);
    pass.b7=pinned_trend_filter(pass.b6,frequency,trend_length,interrupt);
    pass.b8=minus_pinned(pass.b1,pass.b7,interrupt);
    checkpoint(interrupt);
    const auto final9=process_pinned(pass.b8,frequency,initial_position,final_filters,interrupt);
    const auto seasonal9=normalize_pinned(final9,frequency,initial_position,0,0,final_filters,interrupt);
    auto& round2=pass.rounds[1];
    round2.detection=minus_pinned(pass.b8,seasonal9,interrupt);
    round2.analysis_position=static_cast<std::int64_t>(initial_position);
    round2.policy=policy;
    round2.sigma=x11_sigma_weights(repo_view(round2.detection),frequency,
        static_cast<int>(initial_position)+1,policy,0,0,false,lower_sigma,upper_sigma,false,interrupt);
    round2.replacements=x11_replacements(repo_view(pass.b8),frequency,round2.sigma.weights,
        0,0,false,false,interrupt);
    pass.b9=MaskedSeries{round2.replacements.replacements,
        missing_of(round2.replacements.replacements)};
    pass.b9g=MaskedSeries{round2.replacements.corrected,
        missing_of(round2.replacements.corrected)};
    checkpoint(interrupt);
    const auto final10=process_pinned(pass.b9g,frequency,initial_position,final_filters,interrupt);
    pass.b10=normalize_pinned(final10,frequency,initial_position,0,0,final_filters,interrupt);
    pass.b11=minus_pinned(pass.b1,pass.b10,interrupt);
    pass.b13=minus_pinned(pass.b11,pass.b7,interrupt);
    checkpoint(interrupt);
    auto& round3=pass.rounds[2];
    round3.detection=pass.b13;
    round3.analysis_position=static_cast<std::int64_t>(initial_position);
    round3.policy=policy;
    round3.sigma=x11_sigma_weights(repo_view(pass.b13),frequency,
        static_cast<int>(initial_position)+1,policy,0,0,false,lower_sigma,upper_sigma,false,interrupt);
    round3.replacements=x11_replacements(repo_view(pass.b13),frequency,round3.sigma.weights,
        0,0,false,false,interrupt);
    pass.b17=MaskedSeries{round3.sigma.weights,missing_of(round3.sigma.weights)};
    pass.b20=MaskedSeries{round3.replacements.correction_factors,
        missing_of(round3.replacements.correction_factors)};
    pass.status="available";
    checkpoint(interrupt);
    return pass;
}
} // namespace kalos
