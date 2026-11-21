// Item 22c: the fixed-coefficient X-13 deterministic-effect allocation
// staging. Pure algebra over the given columns and coefficients; every
// accumulation order, associativity and branch follows the pinned
// X13Kernel.initialStep / X13Kernel.finals transcription in
// docs/contracts/x13-allocation-contract.md.
#include "kalos/x13_allocation.hpp"
#include <cmath>
#include <stdexcept>
namespace kalos {
namespace {
// The pinned DataBlock.addAY(alpha, y): alpha == 0 leaves x untouched (the
// column is not read at all), alpha == 1 adds y undamped, otherwise
// x[i] + alpha*y[i] with a plain multiply followed by a plain add (no fma).
void addColumn(std::vector<double>& x, double alpha, const std::vector<double>& y) {
    if (alpha==0.) return;
    if (alpha==1.) {
        for (std::size_t i=0;i<x.size();++i) x[i]+=y[i];
        return;
    }
    for (std::size_t i=0;i<x.size();++i) x[i]=x[i]+alpha*y[i];
}
// Per-row elementwise folds of the pinned TsData.add/subtract/divide/
// multiply; IEEE addition and multiplication are commutative bit-exactly,
// so only the left-associative grouping across more than two tables (the
// v23/v25/v26 chains) needs the pinned order, kept below.
std::vector<double> backTransform(const std::vector<double>& x, bool log) {
    if (!log) return x;
    std::vector<double> v(x.size());
    for (std::size_t i=0;i<v.size();++i) v[i]=std::exp(x[i]);
    return v;
}
bool restoreOpMultiplicative(X11DecompositionMode mode) { // pinned op() test
    return mode==X11DecompositionMode::multiplicative
        ||mode==X11DecompositionMode::pseudo_additive;
}
bool restoreInvOpMultiplicative(X11DecompositionMode mode) {
    // pinned invOp(): isMultiplicative() || mode == PseudoAdditive, and
    // DecompositionMode.isMultiplicative() is (Multiplicative || LogAdditive).
    return mode!=X11DecompositionMode::additive;
}
std::vector<double> op(X11DecompositionMode mode, const std::vector<double>& a,
                       const std::vector<double>& b) {
    const bool mult=restoreOpMultiplicative(mode);
    std::vector<double> v(a.size());
    for (std::size_t i=0;i<v.size();++i) v[i]=mult ? a[i]/b[i] : a[i]-b[i];
    return v;
}
std::vector<double> invOp(X11DecompositionMode mode, const std::vector<double>& a,
                          const std::vector<double>& b) {
    const bool mult=restoreInvOpMultiplicative(mode);
    std::vector<double> v(a.size());
    for (std::size_t i=0;i<v.size();++i) v[i]=mult ? a[i]*b[i] : a[i]+b[i];
    return v;
}
// The pinned X11Utility.correctSeries: a copy of x whose entry i takes y
// exactly where w[i] == 0 (negative zero included, NaN kept).
std::vector<double> correct(const std::vector<double>& x, const std::vector<double>& w,
                            const std::vector<double>& y) {
    std::vector<double> v=x;
    for (std::size_t i=0;i<v.size();++i) if (w[i]==0.) v[i]=y[i];
    return v;
}
std::vector<double> correct(std::vector<double> x, const std::vector<double>& w, double y) {
    for (std::size_t i=0;i<x.size();++i) if (w[i]==0.) x[i]=y;
    return x;
}
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
} // namespace
X13Preadjustment x13_preadjust(const X13PreadjustInput& input) {
    X13Preadjustment out;
    const std::size_t n=input.input.size();
    require(n>0, "x13 preadjust requires a non-empty series");
    require(input.backcast>=0&&input.forecast>=0,
        "x13 preadjust horizons are non-negative row counts");
    const std::size_t nb=static_cast<std::size_t>(input.backcast);
    const std::size_t nf=static_cast<std::size_t>(input.forecast);
    const std::size_t m=n+nb+nf; // the extended calendar d8
    for (std::size_t i=0;i<n;++i) {
        require(std::isfinite(input.input[i]), "x13 preadjust requires complete observations");
        if (input.log) require(input.input[i]>0., "x13 preadjust requires a positive series under a log model");
    }
    for (const X13AllocationVariable& v : input.variables) {
        require(v.coefficient.size()==v.columns.size(),
            "x13 preadjust variable has one coefficient per column");
        for (std::size_t j=0;j<v.columns.size();++j) {
            require(v.columns[j].size()==m, "x13 preadjust variable column spans the extended calendar");
            require(std::isfinite(v.coefficient[j]), "x13 preadjust coefficients are finite");
        }
    }
    std::vector<double> v11(m,0.),v12(m,0.),v13(m,0.),v14(m,0.),v15(m,0.),v16(m,0.),v17(m,0.),
        v18(m,0.),v19(m,0.),v20(m,0.),v21(m,0.),v22(m,0.);
    for (const X13AllocationVariable& v : input.variables) {
        for (std::size_t j=0;j<v.columns.size();++j) {
            const double c=v.coefficient[j];
            const std::vector<double>& col=v.columns[j];
            switch (v.cls) {
            case X13VariableClass::moving_holiday: addColumn(v11, c, col); break;
            case X13VariableClass::days_related: addColumn(v12, c, col); break;
            case X13VariableClass::outlier: // exact regeffect match; an
                // outlier outside the three buckets contributes nowhere.
                switch (v.effect) {
                case X13Component::trend: addColumn(v13, c, col); break;
                case X13Component::seasonal: addColumn(v14, c, col); break;
                case X13Component::irregular: addColumn(v15, c, col); break;
                default: break;
                }
                break;
            case X13VariableClass::user:
                switch (v.effect) {
                case X13Component::trend: addColumn(v16, c, col); break;
                case X13Component::seasonal: addColumn(v17, c, col); break;
                case X13Component::irregular: addColumn(v18, c, col); break;
                case X13Component::seasonally_adjusted: addColumn(v19, c, col); break;
                case X13Component::series: addColumn(v20, c, col); break;
                case X13Component::calendar_effect: addColumn(v22, c, col); break;
                case X13Component::undefined: addColumn(v21, c, col); break;
                }
                break;
            }
        }
    }
    // The pinned left-associative folds of initialStep; v26 excludes a9u and
    // the pinned v27 = v26 + a9u only feeds the ARIMA extension route (the
    // item 21-H.4 branch below, the bytecode v26/v27 split).
    std::vector<double> p13=v13; for (std::size_t i=0;i<m;++i) p13[i]+=v16[i];
    std::vector<double> p14=v14; for (std::size_t i=0;i<m;++i) p14[i]+=v17[i];
    std::vector<double> p15=v15; for (std::size_t i=0;i<m;++i) p15[i]+=v18[i];
    std::vector<double> v23=p13; for (std::size_t i=0;i<m;++i) v23[i]=v23[i]+p14[i];
    for (std::size_t i=0;i<m;++i) v23[i]+=p15[i];
    std::vector<double> v25=v22; for (std::size_t i=0;i<m;++i) v25[i]=v25[i]+v19[i];
    for (std::size_t i=0;i<m;++i) v25[i]+=v20[i];
    std::vector<double> v26=v12; for (std::size_t i=0;i<m;++i) v26[i]=v26[i]+v11[i];
    for (std::size_t i=0;i<m;++i) v26[i]=v26[i]+v23[i];
    for (std::size_t i=0;i<m;++i) v26[i]+=v25[i];
    out.a1=input.input; // the interpolated series of a complete calendar
    out.a6=backTransform(v12, input.log);
    out.a7=backTransform(v11, input.log);
    out.a8=backTransform(v23, input.log);
    out.a8t=backTransform(p13, input.log);
    out.a8s=backTransform(p14, input.log);
    out.a8i=backTransform(p15, input.log);
    out.a9=backTransform(v25, input.log);
    out.a9u=backTransform(v21, input.log);
    out.a9cal=backTransform(v22, input.log);
    out.a9sa=backTransform(v19, input.log);
    out.a9ser=backTransform(v20, input.log);
    const std::vector<double> v26bt=backTransform(v26, input.log);
    if (nb==0&&nf==0) { // the zero-horizon path, bit-for-bit as delivered
        out.out.resize(n);
        for (std::size_t i=0;i<n;++i) out.out[i]=input.log
            ? input.input[i]/v26bt[i] : input.input[i]-v26bt[i];
        out.status="available";
        return out;
    }
    // The initialStep extension branch (item 21-H.4): r = v30 - v27 over the
    // input window, the exact-diffuse ARIMA route, then the a-tables
    // 28 = TsData.of(p9, backcasts) + v27, 29 = TsData.of(p10, forecasts)
    // + v27, backTransform'ed, and out = concat(a1b, a1, a1a) (op) bt(v26).
    require(input.arima.has_value(),
        "x13 preadjust requires the stochastic component for a nonzero horizon");
    std::vector<double> v27=v26;
    for (std::size_t i=0;i<m;++i) v27[i]+=v21[i];
    std::vector<double> r(n);
    for (std::size_t i=0;i<n;++i) {
        const double v30=input.log?std::log(input.input[i]):input.input[i];
        r[i]=v30-v27[nb+i]; // TsData.subtract aligns v27 to the input domain
    }
    const ArimaForecastsResult route=arima_forecasts(*input.arima, r, nb, nf);
    if (route.status!="available") { out.status=route.status; return out; }
    std::vector<double> braw(nb),fraw(nf);
    for (std::size_t i=0;i<nb;++i) braw[i]=route.backcasts[i]+v27[i];
    for (std::size_t i=0;i<nf;++i) fraw[i]=route.forecasts[i]+v27[nb+n+i];
    out.a1b=backTransform(braw, input.log);
    out.a1a=backTransform(fraw, input.log);
    out.out.resize(m);
    for (std::size_t i=0;i<m;++i) {
        const double a1ext=i<nb ? out.a1b[i]
            :(i<nb+n ? input.input[i-nb] : out.a1a[i-nb-n]);
        out.out[i]=input.log?a1ext/v26bt[i]:a1ext-v26bt[i];
    }
    out.status="available";
    return out;
}
X13Finals x13_restore(const X13RestoreInput& input) {
    X13Finals out;
    require(input.pre.status=="available", "x13 restore requires an available preadjustment");
    const std::size_t n=input.pre.a1.size();
    require(n>0, "x13 restore requires a non-empty calendar");
    require(input.backcast>=0&&input.forecast>=0,
        "x13 restore horizons are non-negative row counts");
    const std::size_t nb=static_cast<std::size_t>(input.backcast);
    const std::size_t nf=static_cast<std::size_t>(input.forecast);
    const std::size_t m=n+nb+nf;
    const std::vector<double>* tables[]={&input.pre.a6,&input.pre.a7,&input.pre.a8t,
        &input.pre.a8s,&input.pre.a8i,&input.pre.a9cal,&input.pre.a9sa,
        &input.d10,&input.d11,&input.d12,&input.d13,&input.c17};
    for (const std::vector<double>* t : tables)
        require(t->size()==m, "x13 restore tables span the extended calendar");
    // The pinned TsData.fitToDomain slices: the a1 calendar is the middle
    // window of the extended one, the a1b/a1a domains its two ends.
    auto mid=[&](const std::vector<double>& v){
        return std::vector<double>(v.begin()+nb, v.begin()+nb+n); };
    const X11DecompositionMode mode=input.mode;
    const std::vector<double>& a1=input.pre.a1;
    const std::vector<double> t21=invOp(mode, invOp(mode, input.pre.a6, input.pre.a7),
        input.pre.a9cal);
    const std::vector<double> t16=invOp(mode, invOp(mode, input.d10, input.pre.a8s), t21);
    const std::vector<double> t11=invOp(mode,
        invOp(mode, invOp(mode, input.d11, input.pre.a8t), input.pre.a8i), input.pre.a9sa);
    const std::vector<double> t12=invOp(mode, input.d12, input.pre.a8t);
    const std::vector<double> t13=invOp(mode, input.d13, input.pre.a8i);
    if (nf>0) { // the pinned dom(a1a) != null branch
        out.d11a=std::vector<double>(t11.end()-nf, t11.end());
        out.d12a=std::vector<double>(t12.end()-nf, t12.end());
        out.d16a=std::vector<double>(t16.end()-nf, t16.end());
        out.d18a=std::vector<double>(t21.end()-nf, t21.end());
    }
    if (nb>0) { // the pinned dom(a1b) != null branch
        out.d11b=std::vector<double>(t11.begin(), t11.begin()+nb);
        out.d12b=std::vector<double>(t12.begin(), t12.begin()+nb);
        out.d16b=std::vector<double>(t16.begin(), t16.begin()+nb);
        out.d18b=std::vector<double>(t21.begin(), t21.begin()+nb);
    }
    const std::vector<double> f11=mid(t11),f12=mid(t12),f13=mid(t13),f16=mid(t16),f18=mid(t21);
    out.d18=f18; out.d16=f16; out.d11final=f11; out.d12final=f12; out.d13final=f13;
    const std::vector<double> a8i=mid(input.pre.a8i),c17=mid(input.c17),
        d12=mid(input.d12),d13=mid(input.d13);
    const std::vector<double> v28=op(mode, a1, a8i);
    const std::vector<double> v24b=op(mode, f11, a8i);
    out.e1=correct(v28, c17, op(mode, a1, f13));
    out.e2=correct(v24b, c17, d12);
    out.e3=correct(d13, c17,
        restoreInvOpMultiplicative(mode) ? 1. : 0.); // pinned mean(mode)
    out.e11=correct(v24b, c17, invOp(mode, d12, op(mode, v28, out.e1)));
    out.status="available";
    return out;
}
} // namespace kalos
