#include "kalos/arima_forecasts.hpp"

#include <algorithm>
#include <cmath>

// Faithful port of the JVM ExactArimaForecasts consumer stack as
// X13Kernel.initialStep drives it (item 21-H.4; pinning record
// tmp_opencode/kalos-21h4-pinning.md). Class-to-function map:
//   ExactArimaForecasts.prepare/fcasts -> arima_forecasts + ckms_run
//   SsfArima.ArimaData/ArimaDynamics.TX -> FilterModel + arima_tx
//   SsfArima.stInitialize + lambda$diffuseInitializer$1 +
//   CkmsDiffuseInitializer + CkmsFilter.process -> ckms_run
//   AbstractArimaModel.getAutoCovarianceFunction,
//   AutoCovarianceFunction(+AutoCovarianceComputers.defaultComputer),
//   RationalFunction.coefficients -> the polynomial helpers below.

namespace kalos {
namespace {

// JVM Polynomial/BackFilter convolution (natural order, full result).
std::vector<double> poly_times(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> r(a.size() + b.size() - 1, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        for (std::size_t j = 0; j < b.size(); ++j)
            r[i + j] += a[i] * b[j];
    }
    return r;
}

// (1 - B)^d, iterated convolution (exact integer coefficients).
std::vector<double> differencing_regular(int d)
{
    std::vector<double> r{1.0};
    const std::vector<double> one{1.0, -1.0};
    for (int i = 0; i < d; ++i)
        r = poly_times(r, one);
    return r;
}

// (1 - B^s)^bd, iterated convolution.
std::vector<double> differencing_seasonal(int period, int bd)
{
    std::vector<double> r{1.0};
    std::vector<double> one(static_cast<std::size_t>(period) + 1, 0.0);
    one[0] = 1.0;
    one[static_cast<std::size_t>(period)] = -1.0;
    for (int i = 0; i < bd; ++i)
        r = poly_times(r, one);
    return r;
}

// RationalFunction.coefficients(n): Taylor expansion of num/den at 0
// (long division; den[0] non-zero, monic in all uses here). Trailing
// entries stay 0 when num terminates early (the JVM zero-pads).
std::vector<double> rational_coefficients(const std::vector<double>& num,
    const std::vector<double>& den, std::size_t n)
{
    std::vector<double> c(n, 0.0);
    const double d0 = den[0];
    for (std::size_t i = 0; i < n; ++i)
    {
        double x = i < num.size() ? num[i] : 0.0;
        for (std::size_t k = 1; k <= i && k < den.size(); ++k)
            x -= den[k] * c[i - k];
        c[i] = x / d0;
    }
    return c;
}

// Gaussian elimination with partial pivoting. The JVM solves the Yule-Walker
// system with LinearSystemSolver.robustSolver() (a rank-revealing solver); a
// plain pivoted LU agrees whenever the system is non-singular, and the
// singular case is reported as a failure, matching the "acgf of non stationary
// model" ArimaException path (measured-tier, not bit-locked).
bool solve_dense(std::vector<std::vector<double>>& a, std::vector<double>& b)
{
    const std::size_t n = b.size();
    for (std::size_t c = 0; c < n; ++c)
    {
        std::size_t piv = c;
        double best = std::abs(a[c][c]);
        for (std::size_t r = c + 1; r < n; ++r)
        {
            const double v = std::abs(a[r][c]);
            if (v > best)
            {
                best = v;
                piv = r;
            }
        }
        if (!(best > 0.0) || !std::isfinite(best))
            return false;
        if (piv != c)
        {
            std::swap(a[piv], a[c]);
            std::swap(b[piv], b[c]);
        }
        for (std::size_t r = c + 1; r < n; ++r)
        {
            const double f = a[r][c] / a[c][c];
            if (f == 0.0)
                continue;
            for (std::size_t k = c; k < n; ++k)
                a[r][k] -= f * a[c][k];
            b[r] -= f * b[c];
        }
    }
    for (std::size_t r = n; r-- > 0;)
    {
        for (std::size_t k = r + 1; k < n; ++k)
            b[r] -= a[r][k] * b[k];
        b[r] /= a[r][r];
    }
    return true;
}

// AutoCovarianceFunction(ma, ar, ivar).values(count) through the default
// computer: solve the Whittle-type system of size m = max(p, q) + 1 with
// rhs[i] = sum_{j=i..q} ma[j] psi0[j-i], psi0 = (ma/ar).coefficients(q+1),
// A[i][|i-k|] += ar[k]; scale by ivar; extend k >= m by the AR recursion
// ac[k] = -sum_{i=1..p} ar[i] ac[k-i]. The JVM buffers values rounded up to
// multiples of 36; the recursion tail is identical, so only `count` values
// are materialized here.
std::vector<double> autocovariance(const std::vector<double>& ma,
    const std::vector<double>& ar, double ivar, std::size_t count)
{
    const std::size_t p = ar.size() - 1;
    const std::size_t q = ma.size() - 1;
    const std::size_t m = std::max(p, q) + 1;

    const std::vector<double> psi0 = rational_coefficients(ma, ar, q + 1);
    std::vector<double> w(m, 0.0);
    for (std::size_t i = 0; i <= q && i < m; ++i)
    {
        double x = 0.0;
        for (std::size_t j = i; j <= q; ++j)
            x += ma[j] * psi0[j - i];
        w[i] = x;
    }
    std::vector<std::vector<double>> a(m, std::vector<double>(m, 0.0));
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t k = 0; k <= p; ++k)
        {
            if (ar[k] == 0.0)
                continue;
            a[i][i > k ? i - k : k - i] += ar[k];
        }
    }
    if (!solve_dense(a, w))
        return {};
    if (ivar != 1.0)
    {
        for (std::size_t i = 0; i < m; ++i)
            w[i] *= ivar;
    }
    while (w.size() < count)
    {
        const std::size_t k = w.size();
        double x = 0.0;
        for (std::size_t i = 1; i <= p; ++i)
            x += ar[i] * w[k - i];
        w.push_back(-x);
    }
    w.resize(count);
    return w;
}

// ArimaDynamics.TX(0, x) over a full block of length dim:
// xs = -sum_{j=1..|phi|-1} phi[j] x[dim-j]; x <<= 1; x[dim-1] = xs.
// (deg phi <= dim by construction, so x[dim-j] stays in range.)
void arima_tx(const std::vector<double>& phi, std::size_t dim, double* x)
{
    double xs = 0.0;
    for (std::size_t j = 1; j < phi.size(); ++j)
        xs -= phi[j] * x[dim - j];
    for (std::size_t i = 0; i + 1 < dim; ++i)
        x[i] = x[i + 1];
    if (dim >= 1)
        x[dim - 1] = xs;
}

// stVar(P, d1, d2, var): col0 := d2; P[i+1][i+1] = P[i][i] - d1[i]^2 var;
// P[i+1][j+1] = P[i][j] - d1[i] d1[j] var (both factors from d1, per the
// pinned bytecode); completed by symmetry. m = d1.size().
std::vector<double> st_variance(double var, const std::vector<double>& d1,
    const std::vector<double>& d2)
{
    const std::size_t m = d1.size();
    std::vector<double> p(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        p[i * m] = d2[i];
    for (std::size_t i = 0; i + 1 < m; ++i)
    {
        const double b = d1[i];
        p[(i + 1) * m + (i + 1)] = p[i * m + i] - b * b * var;
        for (std::size_t j = 0; j < i; ++j)
            p[(i + 1) * m + (j + 1)] = p[i * m + j] - b * d1[j] * var;
    }
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = i + 1; j < m; ++j)
            p[i * m + j] = p[j * m + i]; // fromLower: mirror lower into upper
    }
    return p;
}

// sigma(S, dif): lower-triangular S[i][j] = c[i-j], c = 1/dif expansion.
std::vector<double> sigma_matrix(const std::vector<double>& dif, std::size_t m)
{
    const std::vector<double> c = rational_coefficients(std::vector<double>{1.0}, dif, m);
    std::vector<double> s(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
            s[i * m + j] = c[i - j];
    }
    return s;
}

// X S X^T for symmetric S (row-major m x m).
std::vector<double> xsxt(const std::vector<double>& x, const std::vector<double>& s,
    std::size_t m)
{
    std::vector<double> xs(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t k = 0; k < m; ++k)
        {
            double v = 0.0;
            for (std::size_t j = 0; j < m; ++j)
                v += x[i * m + j] * s[j * m + k];
            xs[i * m + k] = v;
        }
    }
    std::vector<double> r(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < m; ++j)
        {
            double v = 0.0;
            for (std::size_t k = 0; k < m; ++k)
                v += xs[i * m + k] * x[j * m + k];
            r[i * m + j] = v;
        }
    }
    return r;
}

struct FilterState
{
    std::vector<double> a;
    std::vector<double> l;
    std::vector<double> m; // UpdateInformation.M
    double f = 0.0;        // innovation variance
    bool invalid = false;  // the JVM "Invalid fast filter" SsfException
};

// The filter model on the (possibly bmean-folded) surface, in ArimaData form.
struct FilterModel
{
    std::vector<double> ar;      // full AR polynomial (ArimaData.phi)
    std::vector<double> ma;      // full MA polynomial (ArimaData surface)
    std::vector<double> sta_ar;  // stationary AR (acgf + stacgf input)
    std::vector<double> sta_ma;  // MA of the stationary transform (acgf input)
    std::vector<double> dif;     // unit-root polynomial (StationaryTransformation)
    double var = 1.0;
    std::size_t dim = 0;
    std::vector<double> psi; // (ma/ar).coefficients(dim), length dim
};

// The full CKMS run: stInitialize or the diffuse initializer chain
// (lambda$diffuseInitializer$1 -> CkmsDiffuseInitializer), then the
// CkmsFilter.process loop. The PredictionErrorDecomposition(false) sink is
// inert for the extraction below, so the save step is omitted.
FilterState ckms_run(const FilterModel& model, const std::vector<double>& data,
    bool diffuse)
{
    const std::size_t dim = model.dim;
    FilterState st;
    st.a.assign(dim, 0.0);
    st.l.assign(dim, 0.0);
    st.m.assign(dim, 0.0);

    std::size_t start = 0;
    if (diffuse)
    {
        const std::size_t r = model.dif.size() - 1;
        // Diffuse state a0: A is the (dim+r) x r power expansion of 1/dif
        // (first r rows = identity); a[i] = row (i+r) of A times data[0..r).
        std::vector<double> big((dim + r) * r, 0.0);
        for (std::size_t j = 0; j < r; ++j)
            big[j * r + j] = 1.0;
        for (std::size_t i = r; i < dim + r; ++i)
        {
            for (std::size_t j = 0; j < r; ++j)
            {
                double x = 0.0;
                for (std::size_t k = 1; k <= r; ++k)
                    x -= model.dif[k] * big[(i - k) * r + j];
                big[i * r + j] = x;
            }
        }
        for (std::size_t i = 0; i < dim; ++i)
        {
            double x = 0.0;
            for (std::size_t j = 0; j < r; ++j)
                x += big[(i + r) * r + j] * data[j];
            st.a[i] = x;
        }
        // P = stVar(stpsi, stacgf, var) sigma(dif) stVar(...)^T with the
        // stationary-model inputs (lambda$diffuseInitializer$1).
        const std::vector<double> stacf =
            autocovariance(model.sta_ma, model.sta_ar, model.var, dim);
        if (stacf.empty())
        {
            st.invalid = true;
            return st;
        }
        const std::vector<double> stpsi = rational_coefficients(
            model.sta_ma, model.sta_ar, dim);
        const std::vector<double> pstat = st_variance(model.var, stpsi, stacf);
        const std::vector<double> s = sigma_matrix(model.dif, dim);
        // SymmetricMatrix.XSXt(A, B) computes B·A·B^T - the sigma matrix is
        // the OUTER factor: P = sigma * stVar * sigma^T (verified against the
        // runtime P[0][0] = stacf[0]).
        std::vector<double> p = xsxt(s, pstat, dim);

        // CkmsDiffuseInitializer: F = P[0][0] (no measurement error),
        // M = P.col(0), l = T M, W = -F (T P T' + V - P - l l'/F), then
        // l = W.col(argmax diag W) / sqrt(W[k][k]) or zero.
        st.f = p[0];
        for (std::size_t i = 0; i < dim; ++i)
        {
            st.m[i] = p[i * dim]; // P.col(0)
            st.l[i] = st.m[i];
        }
        arima_tx(model.ar, dim, st.l.data());
        std::vector<double> tpt(dim * dim, 0.0);
        for (std::size_t j = 0; j < dim; ++j)
        {
            std::vector<double> col(dim, 0.0);
            for (std::size_t i = 0; i < dim; ++i)
                col[i] = p[i * dim + j];
            arima_tx(model.ar, dim, col.data());
            for (std::size_t i = 0; i < dim; ++i)
                tpt[i * dim + j] = col[i];
        }
        // T P T': the transition acts on the second index the same way (T
        // is applied to the state on both sides of P).
        for (std::size_t i = 0; i < dim; ++i)
        {
            std::vector<double> row(dim, 0.0);
            for (std::size_t j = 0; j < dim; ++j)
                row[j] = tpt[i * dim + j];
            arima_tx(model.ar, dim, row.data());
            for (std::size_t j = 0; j < dim; ++j)
                tpt[i * dim + j] = row[j];
        }
        const double invf = st.f != 0.0 ? -1.0 / st.f : 0.0;
        for (std::size_t i = 0; i < dim; ++i)
        {
            for (std::size_t j = 0; j < dim; ++j)
            {
                // V = var psi psi^T (ArmaInitialization.v), not diagonal.
                double v = tpt[i * dim + j] - p[i * dim + j];
                v += model.var * model.psi[i] * model.psi[j];
                v += invf * st.l[i] * st.l[j];
                p[i * dim + j] = -st.f * v;
            }
        }
        std::size_t k = 0;
        double best = p[0];
        for (std::size_t i = 1; i < dim; ++i)
        {
            const double v = p[i * dim + i];
            if (v > best)
            {
                best = v;
                k = i;
            }
        }
        if (best > 0.0)
        {
            const double s = std::sqrt(1.0 / best);
            for (std::size_t i = 0; i < dim; ++i)
                st.l[i] = p[i * dim + k] * s;
        }
        else
        {
            // FastMatrix.isZero(1e-6): max |element| <= tolerance.
            double mx = 0.0;
            for (std::size_t i = 0; i < dim * dim; ++i)
                mx = std::max(mx, std::abs(p[i]));
            if (!(mx <= 1e-6))
            {
                st.invalid = true;
                return st;
            }
            std::fill(st.l.begin(), st.l.end(), 0.0);
        }
        start = r;
    }
    else
    {
        // stInitialize: M = acf(dim), F = acf[0], l = T M, no skip.
        const std::vector<double> acf =
            autocovariance(model.sta_ma, model.sta_ar, model.var, dim);
        if (acf.empty())
        {
            st.invalid = true;
            return st;
        }
        st.m = acf;
        st.f = acf[0];
        st.l = acf;
        arima_tx(model.ar, dim, st.l.data());
    }

    // CkmsFilter.process: eps = 1e-15, neps = eps F fixed after initialize;
    // per step: error -> update (a += v/F M) -> next (deflate L, M while
    // |l[0]| > neps, steady when |l| <= eps, then T l; always T a).
    const double eps = 1e-15;
    const double neps = eps * st.f;
    bool steady = false;
    const std::size_t n = data.size();
    for (std::size_t i = start; i < n; ++i)
    {
        const double v = data[i] - st.a[0]; // error(i); complete data
        const double k = v / st.f;
        for (std::size_t t = 0; t < dim; ++t)
            st.a[t] += k * st.m[t];
        if (!steady)
        {
            const double zx = st.l[0];
            if (std::abs(zx) > neps)
            {
                const double x = zx / st.f;
                const double f2 = st.f - zx * x;
                for (std::size_t t = 0; t < dim; ++t)
                {
                    const double lold = st.l[t];
                    st.l[t] -= st.m[t] * x;
                    st.m[t] -= lold * x;
                }
                st.f = f2;
            }
            else
            {
                double norm = 0.0;
                for (std::size_t t = 0; t < dim; ++t)
                    norm += st.l[t] * st.l[t];
                if (!(std::sqrt(norm) > eps))
                    steady = true;
            }
            arima_tx(model.ar, dim, st.l.data());
        }
        arima_tx(model.ar, dim, st.a.data());
    }
    return st;
}

} // namespace

ArimaForecastsResult arima_forecasts(const ArimaForecastModel& model,
    const std::vector<double>& series, std::size_t backcast_horizon,
    std::size_t forecast_horizon)
{
    ArimaForecastsResult res;
    if (series.empty())
    {
        res.status = "failed:arima_forecasts";
        res.message = "empty series";
        return res;
    }
    for (const double v : series)
    {
        if (!std::isfinite(v))
        {
            res.status = "failed:arima_forecasts";
            res.message = "non-finite series";
            return res;
        }
    }

    // --- polynomials (SarimaModel accessors on the landing point) ---------
    const std::size_t period = static_cast<std::size_t>(model.period);
    std::vector<double> reg_ar(model.p + 1, 0.0);
    reg_ar[0] = 1.0;
    for (int i = 0; i < model.p; ++i)
        reg_ar[static_cast<std::size_t>(i) + 1] = model.phi[static_cast<std::size_t>(i)];
    std::vector<double> sea_ar(period * static_cast<std::size_t>(model.bp) + 1, 0.0);
    sea_ar[0] = 1.0;
    for (int i = 0; i < model.bp; ++i)
        sea_ar[period * (static_cast<std::size_t>(i) + 1)] =
            model.bphi[static_cast<std::size_t>(i)];
    std::vector<double> reg_ma(model.q + 1, 0.0);
    reg_ma[0] = 1.0;
    for (int i = 0; i < model.q; ++i)
        reg_ma[static_cast<std::size_t>(i) + 1] = model.theta[static_cast<std::size_t>(i)];
    std::vector<double> sea_ma(period * static_cast<std::size_t>(model.bq) + 1, 0.0);
    sea_ma[0] = 1.0;
    for (int i = 0; i < model.bq; ++i)
        sea_ma[period * (static_cast<std::size_t>(i) + 1)] =
            model.btheta[static_cast<std::size_t>(i)];

    FilterModel f;
    f.sta_ar = poly_times(reg_ar, sea_ar);
    f.ma = poly_times(reg_ma, sea_ma);
    std::vector<double> ur = differencing_regular(model.d);
    ur = poly_times(ur, differencing_seasonal(model.period, model.bd));
    f.var = 1.0; // SarimaModel innovation variance is 1
    if (model.mean_correction)
    {
        // prepare(m, true): ArimaModel(sta_ar, dif (1-B), ma (1-B), var);
        // ArimaData and the stationary transform (acgf inputs) live on this
        // folded surface - the (1-B) factors do NOT cancel in the JVM's
        // floating-point long division, so keep them folded.
        const std::vector<double> d1{1.0, -1.0};
        f.ma = poly_times(f.ma, d1);
        f.dif = poly_times(ur, d1);
        f.ar = poly_times(f.sta_ar, f.dif);
    }
    else
    {
        f.dif = ur;
        f.ar = poly_times(f.sta_ar, ur);
    }
    f.sta_ma = f.ma; // stationary transform: (ar field, ONE, ma field, var)
    // ArimaData.dim = max(deg ar, deg ma + 1) = max(|ar| - 1, |ma|).
    f.dim = std::max(f.ar.size() - 1, f.ma.size());
    f.psi = rational_coefficients(f.ma, f.ar, f.dim);

    const bool diffuse = f.dif.size() > 1;

    // --- run (ExactArimaForecasts.fcasts + the reverse-based backcasts) ---
    auto run_forecasts = [&](const std::vector<double>& input, std::size_t h,
                             double& mean) -> std::pair<bool, std::vector<double>>
    {
        if (h == 0)
            return {true, {}};
        FilterState st = ckms_run(f, input, diffuse);
        if (st.invalid)
            return {false, {}};
        std::vector<double> out(h, 0.0);
        if (h >= f.dim)
        {
            std::vector<double> a = st.a;
            for (std::size_t i = 0; i < f.dim; ++i)
                out[i] = a[i];
            for (std::size_t i = f.dim; i < h; ++i)
            {
                arima_tx(f.ar, f.dim, a.data());
                out[i] = a[f.dim - 1];
            }
        }
        else
        {
            for (std::size_t i = 0; i < h; ++i)
                out[i] = st.a[i];
        }
        if (model.mean_correction)
        {
            // bar = stationaryAr nonStationaryAr of the ORIGINAL model;
            // mean = last entry of bar applied to the final state.
            const std::vector<double> bar = poly_times(f.sta_ar, ur);
            const std::size_t deg = bar.size() - 1;
            std::vector<double> w(f.dim - deg, 0.0);
            for (std::size_t j = 0; j < w.size(); ++j)
            {
                double x = 0.0;
                for (std::size_t i = 0; i <= deg; ++i)
                    x += bar[i] * st.a[j + deg - i];
                w[j] = x;
            }
            mean = w.back();
        }
        return {true, out};
    };

    double mean = 0.0;
    if (forecast_horizon > 0)
    {
        auto [ok, fc] = run_forecasts(series, forecast_horizon, mean);
        if (!ok)
        {
            res.status = "failed:arima_forecasts";
            res.message = "invalid fast filter";
            return res;
        }
        res.forecasts = std::move(fc);
    }
    if (backcast_horizon > 0)
    {
        // ArimaForecasts default backcasts: reverse, forecast, reverse.
        std::vector<double> rev(series.rbegin(), series.rend());
        double bmean = 0.0;
        auto [ok, bc] = run_forecasts(rev, backcast_horizon, bmean);
        if (!ok)
        {
            res.status = "failed:arima_forecasts";
            res.message = "invalid fast filter";
            return res;
        }
        std::reverse(bc.begin(), bc.end());
        res.backcasts = std::move(bc);
        mean = bmean; // the JVM mean field is whatever fcasts() ran last
    }
    res.mean = mean;
    return res;
}

} // namespace kalos
