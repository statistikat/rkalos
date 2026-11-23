#ifndef KALOS_TRAMO_SCAN_MACHINE_HPP
#define KALOS_TRAMO_SCAN_MACHINE_HPP

// Item 23-F F.2: the shared transcription of the pinned jdplus 3.9.0 TRAMO
// outlier machine, factored out of the round-0 decision unit (item 23-B,
// docs/contracts/tramo-outlier-scan-contract.md) so the multi-round driver
// (tramo_outlier_rounds) and the round-0 scan observe the SAME kernels. The
// 23-B corpus replay (tests/jvm_tramo_outlier_scan.h) guards every move here:
// each helper is the round-0 code unchanged, only the design build and the
// table fill became callable per round instead of a linear tail.
//
// Machine order per FastOutliersDetector.process bytecode pass (FOD.txt
// offsets 0-241, the Scratch23g manualRun port is the 1:1 witness): clear();
// regarima = arg; estimateModel (false returns without a tail); the loop
// head is the estimateModel at 19 - a LOCAL verified flag (true on entry,
// replaced by verifyModel when backwardStep is set) GATES the scan: a
// removal round skips its scan and re-estimates, exit breaks to the tail;
// otherwise SingleOutlierDetector.process (false breaks to the tail),
// round++, and the dcmpg |maxT| < cv gate breaks to the TAIL (it does not
// continue - the first below-cv scan ends the machine), addOutlier(position,
// type, coefficient), the size == maxOutliers break, the round < maxRound
// test; the tail estimateModel runs only when exit or round == maxRound or
// size == maxOutliers, estimationStep := false, then the drain while
// (!verifyModel()) estimateModel, and the catch at 238 maps any
// RuntimeException to a Failed run. The register-D3 note
// applies to the caps: the OutliersDetectionModule builder feeds maxRound and
// maxOutliers through the constructor whose arg3 lands on field maxOutliers
// and arg4 on field maxRound (23-F F.2 decode, javap evidence in
// tmp_opencode/kalos-m5-pinning.md SECTION 27), so the caller hands the
// crossed FIELD values to this machine, never the builder pair.
//
// The coefficient vector follows the ConcentratedLikelihoodWithMissing raw
// b layout: [missing-point column coefficients, constant, user columns] -
// allCoefficients is the raw solve and the tstats view skips the leading
// nmissing entries; verifyModel indexes outliers from the tstats TAIL
// (base = tstats.length - outliers.size()), and removeOutlier drops the
// coefficient at the MEAN-shifted user-column index - which with missings
// present drops a missing-block entry, not the outlier coefficient. That is
// the pinned jar behavior (probe cases mr_missing and the SECTION 27 dump
// agreement), mirrored, not fixed.

#include "kalos/detail/fdlibm.hpp"
#include "kalos/detail/ssq_deviances.hpp"
#include "kalos/estimation_control.hpp"
#include "kalos/hannan_rissanen.hpp"
#include "kalos/outlier_scoring.hpp"
#include "kalos/tramo_outlier_rounds.hpp"
#include "kalos/tramo_outlier_scan.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace kalos::detail {

constexpr int kScanPeriod = 12;

inline double scan_qnan() { return std::numeric_limits<double>::quiet_NaN(); }

// AverageInterpolator.cleanMissings: backward scan; a finite value closing a
// gap (left != -1) fills left + 1..right - 1 with (data[right] +
// data[left]) / 2 (data[left] alone when right >= n), where right is the
// first finite index right of the gap; every non-finite index is recorded
// (descending). The trailing branch fills a leading gap 0..right - 1 with
// data[right] and reports false only when no finite value exists at all.
inline bool clean_missings(std::vector<double>& data, std::vector<int>& missing_at) {
    const int n = static_cast<int>(data.size());
    int left = -1;
    int right = -1;
    for (int i = n - 1; i >= 0; --i) {
        const double v = data[static_cast<std::size_t>(i)];
        if (!std::isfinite(v)) {
            if (right == -1) {
                right = i + 1;
            }
            left = i - 1;
            missing_at.push_back(i);
        } else if (left != -1) {
            const double avg = right >= n
                ? data[static_cast<std::size_t>(left)]
                : (data[static_cast<std::size_t>(right)]
                      + data[static_cast<std::size_t>(left)]) / 2.0;
            for (int j = left + 1; j < right; ++j) {
                data[static_cast<std::size_t>(j)] = avg;
            }
            right = -1;
            left = -1;
        }
    }
    if (right > 0) {
        if (right >= n) {
            return false;
        }
        for (int j = 0; j < right; ++j) {
            data[static_cast<std::size_t>(j)] = data[static_cast<std::size_t>(right)];
        }
    }
    return true;
}

struct ScanOrders {
    int p, d, q, bp, bd, bq;
};

// TramoUtility.autlar(n, orders): the pinned default AR order heuristics
// (note the JVM integer log argument n2 * n2 and the integer divisions).
inline int tramo_autlar(int n, const ScanOrders& o) {
    const int d_total = o.d + kScanPeriod * o.bd;
    const int q_total = o.q + kScanPeriod * o.bq;
    const int p_total = o.p + kScanPeriod * o.bp;
    const int n2 = n - d_total;
    int k = static_cast<int>(std::log(static_cast<double>(n2 * n2)));
    const int kmax = std::max(2 * q_total, k);
    if (kmax > k) {
        k = kmax;
    }
    if (!(k < n2)) {
        k = n2 - n2 / 4;
    }
    if (k > 50) {
        k = 50;
    }
    const int tmp = o.p + o.bp * (1 + o.p) + o.q + o.bq * (1 + o.q);
    return (n2 - k) - std::max(p_total, q_total) - tmp;
}

// OutliersDetectionModule.comatip over the finite count of the raw series:
// -1 rejects the series (make hands back the null driver, module Failed);
// 0 keeps the estimation as requested; 1 forces the maximum-likelihood
// flag of the FastOutliersDetector builder.
inline int tramo_comatip(int n_finite, const ScanOrders& o, bool mean, bool ml) {
    int need = std::max(o.d + o.p + kScanPeriod * (o.bd + o.bp), o.q + kScanPeriod * o.bq);
    need += mean ? 1 : 0;
    need += (15 * n_finite) / 100;
    need += kScanPeriod;
    if (!(n_finite - need > 0)) {
        return -1;
    }
    if (ml) {
        return 1;
    }
    const int v5 = tramo_autlar(n_finite, o);
    int v6 = 0;
    if (o.p + o.bp > 0 && o.q + o.bq > 0) {
        ScanOrders reduced = o;
        reduced.p = 0;
        reduced.bp = 0;
        v6 = tramo_autlar(n_finite - (o.p + kScanPeriod * o.bp), reduced);
    }
    return v5 < 0 || v6 < 0 ? 1 : 0;
}

// RobustStandardDeviationComputer.Mad2(50).compute: sort (JVM Double.compare
// order: NaN last), the plain median, the dcmpg-signed absolute deviations
// (a NaN deviation keeps the "x - med" branch), re-sort and the (n + 1) -
// based centile at 50 divided by Phi^-1(0.75) on the pinned quantile path.
inline double mad2_default(const std::vector<double>& x) {
    const std::size_t n = x.size();
    std::vector<double> a = x;
    std::sort(a.begin(), a.end());
    const std::size_t half = n / 2;
    const double med = 2 * half == n ? (a[half - 1] + a[half]) / 2.0 : a[half];
    for (double& v : a) {
        v = std::isless(v, med) ? med - v : v - med;
    }
    std::sort(a.begin(), a.end());
    const double q = static_cast<double>(n + 1) * 50.0 / 100.0;
    const int idx = static_cast<int>(q);
    const double f = q - static_cast<double>(idx);
    const double v = f < 1e-9 ? a[static_cast<std::size_t>(idx) - 1]
                              : a[static_cast<std::size_t>(idx) - 1] * (1.0 - f)
                                    + a[static_cast<std::size_t>(idx)] * f;
    return v / standard_normal_lower_quantile(0.5 + 0.005 * 50.0);
}

// LinearModel.calcResiduals at lm.mean = false (the scan machine always
// carries an explicit constant column when it carries a constant): a copy of
// the response, then per design column in ascending order the whole column
// with -b(k) (DataBlock.addAY over the shared result block). The coefficient
// sequence is read POSITIONALLY (cursor order), which is how the F.2
// missing-prefix misalignment propagates: b may carry the missing block,
// fewer or shifted entries than the design has columns; the entries the
// walk reaches are applied as they come.
inline std::vector<double> calc_residuals(const std::vector<double>& y,
    const std::vector<std::vector<double>>& xcols, const std::vector<double>& b) {
    std::vector<double> res = y;
    for (std::size_t k = 0; k < xcols.size(); ++k) {
        const double coef = k < b.size() ? -b[k] : -scan_qnan();
        for (std::size_t r = 0; r < res.size(); ++r) {
            res[r] += coef * xcols[k][r];
        }
    }
    return res;
}

// SarimaModel.getRegularAR / getSeasonalAR raw polynomials at JVM-sign
// coordinates (trailing 1e-9 trimming cannot change a power-series quotient,
// so the untrimmed convolution is used directly).
inline std::vector<double> raw_regular(const std::vector<double>& phi) {
    std::vector<double> p(phi.size() + 1, 0.0);
    p[0] = 1.0;
    for (std::size_t i = 0; i < phi.size(); ++i) {
        p[i + 1] = phi[i];
    }
    return p;
}

inline std::vector<double> raw_seasonal(const std::vector<double>& bphi) {
    std::vector<double> p = {1.0};
    for (const double c : bphi) {
        std::vector<double> factor(static_cast<std::size_t>(kScanPeriod) + 1, 0.0);
        factor[0] = 1.0;
        factor[kScanPeriod] = c;
        p = detail::pinned_polynomial_times(p, factor);
    }
    return p;
}

// HannanRissanenStart blocks [ar, sar, ma, sma] (Kalos convention) to the
// JVM coordinates [phi, bphi, theta, btheta] of the pinned machines.
inline std::vector<double> hr_to_jvm(const std::array<std::vector<double>, 4>& c) {
    std::vector<double> jvm;
    for (const double v : c[0]) jvm.push_back(-v);
    for (const double v : c[1]) jvm.push_back(-v);
    for (const double v : c[2]) jvm.push_back(v);
    for (const double v : c[3]) jvm.push_back(v);
    return jvm;
}

// Raw (pre-differencing) column of a factory at a position: AO is the unit
// spike, LS the pinned zero-ended step, TC the repeated-multiply decay
// c^(t - pos) (the probe-pinned floating-point pattern of the convolved
// column requires the repeated multiplication, not pow) and SO the
// PeriodicOutlierFactory(zero-ended) fill: v = -1/(period - 1); walk k = 1 to
// period - 1 while pos > 0 decrementing and writing v, then write 1.0 at the
// decremented position, reset k and repeat until pos reaches 0 (guard: the
// write only lands below the block length). Shared by the pre-registered
// variables and by FastOutliersDetector.addOutlier's factory fill.
inline std::vector<double> factory_raw(
    int kind, double tc_rate, int pos, std::size_t n) {
    std::vector<double> col(n, 0.0);
    if (kind == 0) { // AO
        col[static_cast<std::size_t>(pos)] = 1.0;
    } else if (kind == 1) { // LS
        std::fill(col.begin(), col.begin() + pos, -1.0);
    } else if (kind == 2) { // TC (TransitoryChangeFactory.fill): the ramp is
        // written while the NEXT value stays |v| >= ZERO (1e-15); the dcmpg
        // ifge keeps going on a NaN rate, the strict < cuts the tail to exact
        // zeros below the threshold (bytecode offsets 18-43, TCF.txt).
        double v = 1.0;
        for (std::size_t t = static_cast<std::size_t>(pos); t < n; ++t) {
            col[t] = v;
            v *= tc_rate;
            if (std::abs(v) < 1e-15) {
                break;
            }
        }
    } else { // SO, PeriodicOutlierFactory.fill with zeroEnded = true
        const double v = -1.0 / (kScanPeriod - 1);
        int p = pos;
        int k = 1;
        for (;;) {
            if (k < kScanPeriod && p > 0) {
                --p;
                if (p < static_cast<int>(n)) {
                    col[static_cast<std::size_t>(p)] = v;
                }
            } else if (p > 0) {
                --p;
                if (p < static_cast<int>(n)) {
                    col[static_cast<std::size_t>(p)] = 1.0;
                }
                k = 1;
                continue;
            } else {
                break;
            }
            ++k;
        }
    }
    return col;
}

struct ScanFactory {
    int kind;
    int zone_start;
    int zone_end;
    std::vector<double> filter_den; // numerator is 1 for every rep
    double correction;              // only live at differencing degree 0
};

// The pinned factory registration order (the module's factory mix).
inline std::vector<ScanFactory> register_factories(
    bool ao, bool ls, bool tc, double tc_rate, bool so) {
    std::vector<ScanFactory> factories;
    if (ao) {
        factories.push_back({0, 0, 0, {1.0}, 0.0});
    }
    if (ls) {
        factories.push_back({1, 1, 1, {1.0, -1.0}, -1.0});
    }
    if (tc) {
        factories.push_back({2, 0, 1, {1.0, -tc_rate}, 0.0});
    }
    if (so) {
        std::vector<double> den(static_cast<std::size_t>(kScanPeriod) + 1, 0.0);
        den[0] = 1.0;
        den[kScanPeriod] = -1.0;
        factories.push_back({3, 2 * kScanPeriod, 2 * kScanPeriod, std::move(den), 0.0});
    }
    return factories;
}

// The scan tables of FastOutlierDetector.processOutlier: per factory column
// the pi / (pi * rep) weight column (RationalBackFilter.times with the
// SimplifyingTool exact-cancel step) and the incremental-ssq sweep over the
// live cells (zone, domain and exclusion already folded into the mask). The
// tables use the JVM storage order (row + column * m); dead cells carry
// fill_dead (the round-0 fresh arrays keep NaN where the JVM table stays
// default-zero - decision-invariant for searchMax, pinned per round by the
// probe dumps).
struct ScanTables {
    std::vector<double> t_table;
    std::vector<double> coef_table;
};

inline void build_scan_tables(const std::vector<ScanFactory>& factories,
    const std::vector<double>& jvm, std::size_t p, std::size_t bp, std::size_t q,
    std::size_t /*bq*/, const std::vector<double>& kernel, const std::vector<double>& el,
    double mad, const std::vector<char>& live, bool fill_dead_nan, ScanTables& out) {
    const std::size_t m = el.size() + kernel.size() - 1;
    const std::size_t k = factories.size();
    out.t_table.assign(m * k, fill_dead_nan ? scan_qnan() : 0.0);
    out.coef_table.assign(m * k, fill_dead_nan ? scan_qnan() : 0.0);

    // pi = ar / ma of the full model (the differencing belongs into the
    // numerator, SarimaModel.getAr); processOutlier composes it with the
    // factory filter and reads m weights.
    std::vector<double> phi_jvm(jvm.begin(), jvm.begin() + static_cast<std::ptrdiff_t>(p));
    std::vector<double> bphi_jvm(jvm.begin() + static_cast<std::ptrdiff_t>(p),
        jvm.begin() + static_cast<std::ptrdiff_t>(p + bp));
    std::vector<double> th_jvm(jvm.begin() + static_cast<std::ptrdiff_t>(p + bp),
        jvm.begin() + static_cast<std::ptrdiff_t>(p + bp + q));
    std::vector<double> bth_jvm(jvm.begin() + static_cast<std::ptrdiff_t>(p + bp + q),
        jvm.end());
    std::vector<double> ar_raw = detail::pinned_polynomial_times(
        raw_regular(phi_jvm), raw_seasonal(bphi_jvm));
    ar_raw = detail::pinned_polynomial_times(ar_raw, kernel);
    const std::vector<double> ma_raw = detail::pinned_polynomial_times(
        raw_regular(th_jvm), raw_seasonal(bth_jvm));

    const std::size_t n_el = el.size();
    for (std::size_t o = 0; o < k; ++o) {
        // RationalBackFilter.times(pi, rep.filter): the SimplifyingTool
        // cross-cancels rep.den against the pi numerator through
        // LeastSquaresDivision (simplifyExact) before multiplying. The LS and
        // SO denominators (1 - B, 1 - B^12) divide the unit-root-laden pi
        // numerator exactly, so those columns run on the quotient over the
        // plain MA; the AO constant denominator skips the simplifier and the
        // TC geometric denominator misses its root window (the roots pass
        // tests |pi(1/delta)| < 1e-12, which fails by orders of magnitude),
        // leaving the unsimplified product.
        const std::vector<double>& fden = factories[o].filter_den;
        std::optional<std::vector<double>> quotient;
        if (fden.size() > 1) {
            // The public simplify() runs its private pass on the LONGER
            // polynomial as dividend (swapping the results back afterwards),
            // so a pi numerator shorter than the filter denominator divides
            // INTO it instead, leaving numerator 1 and the quotient in the
            // denominator.
            const bool swapped = fden.size() > ar_raw.size();
            quotient = detail::pinned_ls_divide(swapped ? fden : ar_raw,
                swapped ? ar_raw : fden);
            if (!quotient && !swapped && fden.size() == 2) {
                // Roots pass on a linear denominator: the single root must
                // miss the 1e-12 evaluation window (the TC geometric root
                // 1/c of the pi numerator, off by orders of magnitude).
                const double root = -fden[0] / fden[1];
                double val = 0.0; // real-root Horner (the Complex walk at im = 0)
                for (std::size_t i = ar_raw.size(); i-- > 0;) {
                    val = val * root + ar_raw[i];
                }
                if (std::abs(val) < 1e-12) {
                    throw std::invalid_argument(
                        "tramo scan machine: the simplifier's roots pass "
                        "engaged, which is outside the pinned machine");
                }
            }
            if (!quotient && (swapped || fden.size() > 2)) {
                throw std::invalid_argument("tramo scan machine: a non-exact "
                    "non-linear filter denominator leaves the pinned "
                    "simplifier's roots pass");
            }
        }
        std::vector<double> num, den;
        if (!quotient) {
            num = detail::pinned_polynomial_times(ar_raw, {1.0});
            den = detail::pinned_polynomial_times(ma_raw, fden);
        } else if (fden.size() <= ar_raw.size()) {
            num = detail::pinned_polynomial_times(*quotient, {1.0});
            den = ma_raw; // the MA times ONE
        } else { // swapped exact pass: numerator ONE, quotient in the tail
            num = detail::pinned_polynomial_times(std::vector<double>{1.0}, {1.0});
            den = detail::pinned_polynomial_times(ma_raw, *quotient);
        }
        auto w = detail::pinned_rational_weights(num, den, m);
        // FastOutlierDetector.processOutlier applies the level-shift
        // correction at frequency zero only when the model is stationary.
        double c0 = 0.0;
        if (kernel.size() == 1 && factories[o].correction != 0.0) {
            const auto at_one = [](const std::vector<double>& poly) {
                double x = 0.0;
                for (auto it = poly.rbegin(); it != poly.rend(); ++it) x += *it;
                return x;
            };
            c0 = factories[o].correction * at_one(ar_raw) / at_one(ma_raw);
            for (auto& wi : w) wi += c0;
        }
        double ssq = c0 * c0 * static_cast<double>(n_el);
        for (std::size_t i = 0; i < m; ++i) {
            ssq += w[i] * w[i];
            if (c0 != 0.0) {
                ssq -= c0 * c0;
            }
            std::size_t j = i + 1;
            if (j > n_el) {
                j = n_el;
                ssq -= w[i - n_el] * w[i - n_el];
                if (c0 != 0.0) {
                    ssq += c0 * c0;
                }
            }
            double numv = 0.0;
            for (std::size_t t = 0; t < j; ++t) {
                numv += el[n_el - 1 - t] * w[i - t];
            }
            if (c0 != 0.0) {
                double s = 0.0;
                for (std::size_t t = 0; t < n_el - j; ++t) {
                    s += el[t];
                }
                numv += s * c0;
            }
            const std::size_t row = m - 1 - i;
            if (live[row + o * m]) {
                const double cf = numv / ssq;
                out.coef_table[row + o * m] = cf;
                out.t_table[row + o * m] = cf * std::sqrt(ssq) / mad;
            }
        }
    }
}

// SingleOutlierDetector.searchMax: column-major sweep over the raw storage
// with unit weights, a NaN candidate loses (dcmpl against the running
// maximum) and a column that found nothing at all aborts the search (the
// JVM then trips on the (-1, -1) table read: the driver's process catches
// the RuntimeException and returns false). idx < 0 marks that abort.
inline void search_max(const ScanTables& tables, std::size_t m, int& idx, double& dmax) {
    const std::size_t k = tables.t_table.size() / m;
    idx = -1;
    dmax = 0.0;
    std::size_t storage = 0;
    for (std::size_t o = 0; o < k; ++o) {
        for (std::size_t row = 0; row < m; ++row, ++storage) {
            const double val = std::abs(tables.t_table[storage]) * 1.0;
            if (val > dmax) { // dcmpl ifle: a NaN value never updates
                idx = static_cast<int>(storage);
                dmax = val;
            }
        }
        if (idx == -1) {
            return;
        }
    }
}

// ---------------------------------------------------------------------------
// The multi-round machine state (item 23-F F.2). One struct mirrors the
// FastOutliersDetector fields plus its SingleOutlierDetector, observed by the
// Scratch23g probe dumps (validation/evidence/native/tramo-outlier-rounds/).
// The vector<double> "null" fields carry explicit *_null flags: the JVM
// driver reads them (arraylength, calcResiduals) and the resulting NPE is
// caught by process as a Failed run, which the ports mirror with ok=false.

struct ScanMachine {
    // Pinned request data (post make()/calc(); validation belongs to the
    // caller - both public units share the 23-B static gate).
    ScanOrders orders;
    bool mean = false;
    bool log = false;
    bool mvx = false;
    double precision = 1e-7;
    double tc_rate = 0.7;
    double cv = 0.0;
    // FIELD semantics (crossed relative to the builder, see the header note).
    int max_round = 30;
    int max_outliers = 50;
    int lower_bound = 0;
    int upper_bound = 0;
    std::vector<ScanFactory> factories;

    // Pre-stage results.
    bool failed = false;      // pre-stage dead end (Failed, no scan)
    bool scan_off = false;    // empty mix or collapsed bounds (calc returns false)
    std::vector<double> z;    // interpolated, logged series (raw model length)
    std::vector<double> dly;  // differenced response
    std::vector<double> kernel;
    std::vector<int> missing_at; // ascending raw positions
    std::vector<std::vector<double>> dummy_cols; // differenced point columns
    std::vector<double> jvm;  // ARIMA parameters at JVM coordinates

    // User regressors: the driver's RegArimaModel.getX() (pre-registered
    // columns first, added outlier columns tail it), stored at RAW length.
    std::vector<std::vector<double>> user_raw;

    // FastOutliersDetector fields (constructor defaults, witness `st init`:
    // rflag 1, estimationStep 1, everything else off/empty).
    int round = 0;
    bool rflag = true;
    bool estimation_step = true;
    bool backward_step = false;
    bool exit = false;
    std::vector<double> coeff;
    bool coeff_null = true;
    std::vector<double> tstats;
    bool tstats_null = true;
    std::vector<double> res;
    bool res_null = true;
    std::vector<std::array<int, 2>> outliers; // {position, mix column}
    std::array<int, 2> lastremoved{{-1, -1}};
    bool lastremoved_null = true;

    // Scan-cell masks (storage row + column * m, column = factory ordinal).
    std::size_t m = 0;
    std::vector<char> zone_ok;  // factory zones, immutable
    std::vector<char> excluded; // exclusions, mutated by add/allow
    std::vector<std::pair<std::size_t, std::size_t>> prespec_cells; // (row, column)

    std::size_t par_count() const {
        return static_cast<std::size_t>(
            std::max(0, orders.p) + std::max(0, orders.bp)
            + std::max(0, orders.q) + std::max(0, orders.bq));
    }

    // differencedModel().asLinearModel(): point columns, the constant when
    // the mean is estimated, then the differenced user columns.
    std::vector<std::vector<double>> design() const {
        std::vector<std::vector<double>> xcols = dummy_cols;
        if (mean) {
            xcols.push_back(std::vector<double>(dly.size(), 1.0));
        }
        for (const auto& col : user_raw) {
            xcols.push_back(detail::pinned_backfilter_apply(kernel, col));
        }
        return xcols;
    }
};

// The pre-stage (make + ModelDescription.build + calc): the comatip size
// gate (a negative verdict is Failed and also forces mvx on a positive one),
// the empty-mix / collapsed-bounds scan switch (calc returns false: the
// module reports the empty selection), AverageInterpolator.cleanMissings
// (retaining the positions, ascending) and the log transformation, the
// differenced response, the point columns, the constant handled by
// design(), the pre-registered columns with their own exclusion cells and
// the zone masks. The caller runs the static request gate first (both
// public units share the 23-B validation, including the pre-registered SO
// rejection); the log positivity check is NOT run here (the units throw
// from their own validation, matching the round-0 order). missing_at stays
// empty when interpolation found nothing (or was not requested and none
// exist).
inline void prestage(ScanMachine& sm, const std::vector<double>& y,
    const std::vector<std::pair<int, TramoOutlierKind>>& prespec,
    const char* log_error) {
    if (sm.factories.empty()) {
        sm.scan_off = true;
    }
    if (sm.upper_bound <= sm.lower_bound) {
        sm.scan_off = true;
    }
    sm.z = y;
    std::vector<int> descending;
    {
        bool any_missing = false;
        for (const double v : sm.z) {
            if (!std::isfinite(v)) {
                any_missing = true;
                break;
            }
        }
        if (any_missing) {
            if (!clean_missings(sm.z, descending)) {
                sm.failed = true;
                return;
            }
            std::sort(descending.begin(), descending.end());
            sm.missing_at = descending; // JVM IntList order is descending
        }
    }
    if (sm.log) {
        for (const double v : sm.z) {
            if (!(v > 0.0)) {
                throw std::invalid_argument(log_error);
            }
        }
        for (double& v : sm.z) {
            v = detail::fdlibm_log(v);
        }
    }
    sm.kernel = detail::pinned_differencing_polynomial(
        kScanPeriod, sm.orders.d, sm.orders.bd);
    sm.dly = detail::pinned_backfilter_apply(sm.kernel, sm.z);
    if (sm.dly.empty()) {
        sm.failed = true;
        return;
    }
    for (const int mpos : sm.missing_at) {
        std::vector<double> col(sm.z.size(), 0.0);
        col[static_cast<std::size_t>(mpos)] = 1.0;
        sm.dummy_cols.push_back(detail::pinned_backfilter_apply(sm.kernel, col));
    }
    for (const auto& pre : prespec) {
        // The probe hands the outlier code straight to addOutlier as the
        // type, so the JVM reads it as the MIX COLUMN index (identical to
        // the kind code only in the full mix).
        const int type = static_cast<int>(pre.second);
        if (static_cast<std::size_t>(pre.first) >= y.size()) {
            continue; // excluded by ModelDescription; the module exclude no-ops
        }
        if (type < 0 || static_cast<std::size_t>(type) >= sm.factories.size()) {
            continue; // the factory read is out of bounds: exclude no-ops
        }
        sm.user_raw.push_back(factory_raw(
            sm.factories[static_cast<std::size_t>(type)].kind, sm.tc_rate,
            pre.first, sm.z.size()));
        sm.prespec_cells.emplace_back(
            static_cast<std::size_t>(pre.first),
            static_cast<std::size_t>(type));
    }
    // The SarimaSpec builder fill at JVM coordinates (phi -0.1, theta and
    // btheta -0.2), the model estimateModel keeps when HannanRissanen fails.
    for (int i = 0; i < sm.orders.p + sm.orders.bp; ++i) {
        sm.jvm.push_back(-0.1);
    }
    for (int i = 0; i < sm.orders.q + sm.orders.bq; ++i) {
        sm.jvm.push_back(-0.2);
    }
    sm.m = sm.z.size();
    const std::size_t k = sm.factories.size();
    sm.zone_ok.assign(sm.m * k, 0);
    sm.excluded.assign(sm.m * k, 0);
    for (std::size_t o = 0; o < k; ++o) {
        const long long zs = sm.factories[o].zone_start;
        const long long ze = sm.factories[o].zone_end;
        for (std::size_t row = 0; row < sm.m; ++row) {
            sm.zone_ok[row + o * sm.m] = static_cast<char>(
                static_cast<long long>(row) >= zs
                && static_cast<long long>(row)
                    < static_cast<long long>(sm.m) - ze);
        }
    }
    for (const int mpos : sm.missing_at) {
        for (std::size_t o = 0; o < k; ++o) {
            sm.excluded[static_cast<std::size_t>(mpos) + o * sm.m] = 1;
        }
    }
    for (const auto& cell : sm.prespec_cells) {
        sm.excluded[cell.first + static_cast<std::size_t>(cell.second) * sm.m] = 1;
    }
}

// estimateModel (FOD offsets 0-301): the residual branch (rflag: OLS over
// the current design - a null Ols is "Ols failed", ok=false; else the raw
// coeff vector when present, else the response), the estimationStep
// HannanRissanen block (install at identity, on mvx, or at round 0; the
// round>0 non-identity non-mvx path re-arms rflag without installing), the
// warm re-estimation seam (optimize_sarima_from) on mvx or a non-identity
// landing - which returns through optimizeModel, whose own
// updateLikelihood of the estimation's concentrated likelihood is mirrored
// by the GLS point re-solve over the optimized parameters - and the plain
// updateLikelihood refresh of coeff / tstats / res. A vars==0 design keeps
// the previous res, as the JVM field does.
inline bool estimate_model(ScanMachine& sm) {
    const std::size_t p = static_cast<std::size_t>(sm.orders.p);
    const std::size_t bp = static_cast<std::size_t>(sm.orders.bp);
    const std::size_t q = static_cast<std::size_t>(sm.orders.q);
    const std::size_t bq = static_cast<std::size_t>(sm.orders.bq);
    const std::vector<std::vector<double>> xcols = sm.design();

    if (sm.rflag) {
        if (!xcols.empty()) {
            const auto ols = detail::pinned_ols_coefficients(sm.dly, xcols);
            if (!ols) {
                return false; // "Ols failed"
            }
            sm.res = calc_residuals(sm.dly, xcols, *ols);
        } else {
            sm.res = sm.dly; // LinearModel.getY()
        }
        sm.res_null = false;
    } else if (!sm.coeff_null) {
        sm.res = calc_residuals(sm.dly, xcols, sm.coeff);
        sm.res_null = false;
    } else {
        sm.res = sm.dly;
        sm.res_null = false;
    }
    sm.rflag = false;

    bool identity = true; // SarimaMapping.stabilize returned the HR model
    if (sm.estimation_step && sm.par_count() > 0) {
        const HannanRissanenStart hr = hannan_rissanen_start(sm.res,
            RegressionDesign{0, {}, {}, {}},
            HannanRissanenOrders{sm.orders.p, sm.orders.bp, sm.orders.q,
                sm.orders.bq, kScanPeriod});
        if (hr.source != HannanRissanenStart::Source::Defaulted) {
            const std::vector<double> landing = hr_to_jvm(hr.coefficients);
            identity = !hr.stabilized; // the landing equals the stabilizer input
            if (identity || sm.mvx || sm.round == 0) {
                // offsets 219-245: install at identity, on the
                // maximum-likelihood flag, or at round 0
                sm.jvm = landing;
            } else {
                sm.rflag = true; // offsets 248-250: no install, OLS re-arms
                identity = true; // offsets 253-254: stable := true, which
                                 // gates the optimizeModel below off
            }
        }
    }

    if (sm.mvx || !identity) {
        // optimizeModel (offsets 256-273): regarima :=
        // RegArimaUtility.processor(true, eps).optimize(regarima, mapping)
        // from the CURRENT model parameters - the landing only where the
        // install branch ran - and the estimateModel return happens INSIDE
        // optimizeModel, which feeds updateLikelihood the estimation's own
        // concentrated likelihood. Re-solving the same GLS point over the
        // optimized parameters below reproduces that concentrated
        // likelihood; the warm seam (optimize_sarima_from, transcribed for
        // exactly this re-estimation) covers the regarima surface: the raw
        // level response with its missing mask and the free constant and
        // user columns.
        SarimaEstimationRequest warm;
        warm.series.values = sm.z;
        warm.series.missing.assign(sm.z.size(), 0);
        for (const int mpos : sm.missing_at) {
            warm.series.missing[static_cast<std::size_t>(mpos)] = 1;
        }
        std::vector<std::vector<double>> reg_cols;
        if (sm.mean) {
            reg_cols.push_back(std::vector<double>(sm.z.size(), 1.0));
        }
        for (const auto& col : sm.user_raw) {
            reg_cols.push_back(col);
        }
        if (!reg_cols.empty()) {
            warm.design.columns = reg_cols.size();
            warm.design.values.assign(
                sm.z.size() * reg_cols.size(), 0.0);
            for (std::size_t c = 0; c < reg_cols.size(); ++c) {
                for (std::size_t r = 0; r < sm.z.size(); ++r) {
                    warm.design.values[r + sm.z.size() * c] = reg_cols[c][r];
                }
            }
            warm.design.coefficients.assign(reg_cols.size(), 0.0);
            warm.design.fixed.assign(reg_cols.size(), 0);
        }
        warm.blocks[0].initial.assign(p, 0.0);
        warm.blocks[0].fixed.assign(p, 0);
        warm.blocks[1].initial.assign(q, 0.0);
        warm.blocks[1].fixed.assign(q, 0);
        warm.blocks[2].initial.assign(bp, 0.0);
        warm.blocks[2].fixed.assign(bp, 0);
        warm.blocks[3].initial.assign(bq, 0.0);
        warm.blocks[3].fixed.assign(bq, 0);
        warm.period = kScanPeriod;
        warm.d = sm.orders.d;
        warm.seasonal_d = sm.orders.bd;
        warm.control = EstimationControl{sm.precision, 1e-4, true, false, true,
            true, EstimationStartingPoint::HannanRissanen,
            EstimationSolver::PinnedLevenbergMarquardt};
        std::array<std::vector<double>, 4> start;
        {
            const auto at = [&](std::size_t o) {
                return sm.jvm[static_cast<std::size_t>(o)];
            };
            for (std::size_t i = 0; i < p; ++i) start[0].push_back(-at(i));
            for (std::size_t i = 0; i < bp; ++i) start[2].push_back(-at(p + i));
            for (std::size_t i = 0; i < q; ++i) start[1].push_back(at(p + bp + i));
            for (std::size_t i = 0; i < bq; ++i)
                start[3].push_back(at(p + bp + q + i));
        }
        const SarimaEstimationResult est = optimize_sarima_from(warm, start);
        const auto& c = est.coefficients; // [phi, theta, bphi, btheta], Kalos
        for (const auto& block : c) {
            for (const double value : block) {
                if (!std::isfinite(value)) {
                    return false;
                }
            }
        }
        sm.jvm.clear();
        for (const double v : c[0]) sm.jvm.push_back(-v);
        for (const double v : c[2]) sm.jvm.push_back(-v);
        for (const double v : c[1]) sm.jvm.push_back(v);
        for (const double v : c[3]) sm.jvm.push_back(v);
    }

    if (!xcols.empty()) {
        // updateLikelihood(DEFAULT_COMPUTER.compute(regarima), parCount):
        // allCoefficients is the raw solve b (missing-correction coefficients
        // first), tstats runs over the coefficients() view (skipping the
        // nmissing leading entries) against ser(parCount, true) =
        // sqrt(uc * ssqErr / (ndata - nmissing - nx - parCount)).
        const std::size_t nm = sm.missing_at.size();
        const auto gls = detail::pinned_gls_point(sm.dly, p, bp, q, bq,
            kScanPeriod, xcols, sm.jvm, nm);
        if (!gls) {
            return false; // "Gls failed"
        }
        sm.coeff = gls->b;
        sm.coeff_null = false;
        const std::size_t nx = sm.coeff.size() - nm;
        const double div = static_cast<double>(
            sm.dly.size() - nm - nx - sm.par_count());
        sm.tstats.assign(nx, 0.0);
        for (std::size_t k = 0; k < nx; ++k) {
            const double ser = std::sqrt(gls->uc_diag[k] * gls->ssq_err / div);
            sm.tstats[k] = sm.coeff[nm + k] / ser;
        }
        sm.tstats_null = false;
        sm.res = calc_residuals(sm.dly, xcols, sm.coeff);
        sm.res_null = false;
    }
    return true;
}

// FastOutlierDetector.initmodel + processOutlier + searchMax, the body of
// SingleOutlierDetector.process (its clear(false) keeps the exclusion mask).
// status: CleanFalse mirrors calc rejecting (collapsed bounds / empty mix:
// process returns false, the probe still dumps the sod state); Throw mirrors
// an initmodel failure (the JVM RuntimeException caught by process without
// any dump); Ok carries the winner unless search_max aborted (search_ok
// false: the getMaxTStat table read throws, the driver still dumps).
struct RoundScan {
    enum class Status { CleanFalse, Throw, Ok };
    Status status = Status::CleanFalse;
    bool search_ok = false;
    int position = -1;
    int kind = -1;         // winning column index (the JVM oMax/getMaxType)
    double max_t = 0.0;
    double coefficient = 0.0;
    double mad = 0.0;
    std::vector<double> el;
};

inline RoundScan scan_process(ScanMachine& sm) {
    RoundScan out;
    if (sm.scan_off) {
        return out; // calc found the bounds collapsed or the mix empty
    }
    const std::size_t p = static_cast<std::size_t>(sm.orders.p);
    const std::size_t bp = static_cast<std::size_t>(sm.orders.bp);
    const std::size_t q = static_cast<std::size_t>(sm.orders.q);
    const std::size_t bq = static_cast<std::size_t>(sm.orders.bq);
    const std::vector<std::vector<double>> xcols = sm.design();
    if (!xcols.empty()) {
        const auto gls = detail::pinned_gls_coefficients(
            sm.dly, p, bp, q, bq, kScanPeriod, xcols, sm.jvm);
        if (!gls) {
            out.status = RoundScan::Status::Throw;
            return out; // initmodel throws "Gls failed": process catches
        }
        const std::vector<double> wr = calc_residuals(sm.dly, xcols, *gls);
        const auto innovations = detail::pinned_fast_kalman_residuals(
            wr, p, bp, q, bq, kScanPeriod, sm.jvm);
        if (!innovations) {
            out.status = RoundScan::Status::Throw;
            return out; // a failed whitening throws out of initmodel
        }
        out.el = *innovations;
    } else {
        const auto innovations = detail::pinned_kalman_innovations(
            sm.dly, p, bp, q, bq, kScanPeriod, sm.jvm);
        if (!innovations) {
            out.status = RoundScan::Status::Throw;
            return out; // a failed stationary filter throws out of initmodel
        }
        out.el = *innovations;
    }
    out.status = RoundScan::Status::Ok;
    out.mad = mad2_default(out.el);
    const std::size_t m_el = out.el.size() + sm.kernel.size() - 1;
    if (m_el != sm.m || sm.m == 0) {
        out.status = RoundScan::Status::Throw; // defensive: the pinned
        return out;                            // prepare uses the domain length
    }
    const std::size_t k = sm.factories.size();
    std::vector<char> live(m_el * k, 0);
    for (std::size_t o = 0; o < k; ++o) {
        for (std::size_t row = 0; row < m_el; ++row) {
            live[row + o * m_el] = static_cast<char>(
                sm.zone_ok[row + o * m_el] && !sm.excluded[row + o * m_el]
                && static_cast<long long>(row) >= sm.lower_bound
                && static_cast<long long>(row) < sm.upper_bound);
        }
    }
    ScanTables tables;
    build_scan_tables(sm.factories, sm.jvm, p, bp, q, bq, sm.kernel, out.el,
        out.mad, live, true, tables);
    int idx = -1;
    double dmax = 0.0;
    search_max(tables, m_el, idx, dmax);
    if (idx < 0) {
        return out; // searchMax aborted: the getMaxTStat table read trips
    }
    out.search_ok = true;
    out.position = idx % static_cast<int>(m_el);
    out.kind = idx / static_cast<int>(m_el); // searchMax stores oMax = COLUMN
    out.max_t = tables.t_table[static_cast<std::size_t>(idx)];
    out.coefficient = tables.coef_table[static_cast<std::size_t>(idx)];
    return out;
}

// FastOutliersDetector.addOutlier(position, outlierType, coefficient):
// addOutlier(II) appends the user column (factory fill, RAW length), records
// the outlier, excludes the scan cell, then the coefficient vector EXTENDS
// with the scan coefficient (Doubles.of when null) and backwardStep := true.
// The type is the searchMax oMax, the MIX COLUMN index (SOD.searchMax stores
// idx / rows; getMaxType returns it verbatim) - identical to the kind code
// only in the full AO,LS,TC,SO mix.
inline bool add_outlier(ScanMachine& sm, int position, int type, double coefficient) {
    if (type < 0 || static_cast<std::size_t>(type) >= sm.factories.size()) {
        return false; // the JVM factory read throws: process catches (238)
    }
    sm.user_raw.push_back(factory_raw(
        sm.factories[static_cast<std::size_t>(type)].kind, sm.tc_rate, position,
        sm.z.size()));
    sm.outliers.push_back({position, type});
    if (sm.coeff_null) {
        sm.coeff = {coefficient};
        sm.coeff_null = false;
    } else {
        sm.coeff.push_back(coefficient);
    }
    sm.excluded[static_cast<std::size_t>(position)
        + static_cast<std::size_t>(type) * sm.m] = 1;
    sm.backward_step = true;
    return true;
}

// FastOutliersDetector.allow(position, type): the exclusion mask opens the
// cell (zone membership is re-checked by the table fill, not here).
inline void allow_cell(ScanMachine& sm, int position, int type) {
    if (type < 0 || static_cast<std::size_t>(type) >= sm.factories.size()) {
        return;
    }
    sm.excluded[static_cast<std::size_t>(position)
        + static_cast<std::size_t>(type) * sm.m] = 0;
}

// FastOutliersDetector.removeOutlier(i): j = X.size() - outliers.size() + i
// (the user columns - pre-registered first, outliers tail them), removeX(j),
// outliers.remove(i); the coefficient vector of length 1 (the mean alone)
// goes null, else drops coeff[j + (isMean ? 1 : 0)] - the missing-correction
// block offset is NOT compensated, the pinned off-by-nmissing misalignment
// (SECTION 27), which with missings present drops a block entry instead of
// the outlier coefficient. A null coeff or an out-of-range drop mirrors the
// JVM NPE / index throw caught by process: ok=false (Failed).
inline bool remove_outlier(ScanMachine& sm, std::size_t i) {
    const std::size_t j = sm.user_raw.size() - sm.outliers.size() + i;
    if (j >= sm.user_raw.size()) {
        return false;
    }
    sm.user_raw.erase(sm.user_raw.begin() + static_cast<std::ptrdiff_t>(j));
    sm.outliers.erase(sm.outliers.begin() + static_cast<std::ptrdiff_t>(i));
    if (sm.coeff_null) {
        return false; // the NPE at coeff.length
    }
    if (sm.coeff.size() == 1) {
        sm.coeff.clear();
        sm.coeff_null = true;
        return true;
    }
    std::size_t drop = j + (sm.mean ? 1u : 0u);
    if (drop >= sm.coeff.size()) {
        return false; // the index throw
    }
    std::vector<double> rebuilt;
    rebuilt.reserve(sm.coeff.size() - 1);
    // The JVM walk: copy 0..drop-1 (iload < iload_2 loop), then drop+1..end.
    for (std::size_t k = 0; k < drop; ++k) {
        rebuilt.push_back(sm.coeff[k]);
    }
    for (std::size_t k = drop + 1; k < sm.coeff.size(); ++k) {
        rebuilt.push_back(sm.coeff[k]);
    }
    sm.coeff = std::move(rebuilt);
    return true;
}

// FastOutliersDetector.verifyModel (offsets 0-170): estimationStep := true
// on entry; an empty set is "verified"; the least-significant outlier is the
// argmin of |tstats[base + i]| over the TAIL block (base = length - count,
// the outlier columns always end the view), a NaN candidate never wins the
// argmin (dcmpg) but a NaN incumbent stays best and its dcmpl gate (< cv)
// fires on it, taking the removal path. Removal clears backwardStep and
// estimationStep, opens the cell, removes, and compares against lastremoved
// for the exit flag. Returns 1/0 for the JVM boolean, -1 mirrors the
// exception paths of the null-field reads (the JVM RuntimeException caught
// by process as a Failed run).
inline int verify_model(ScanMachine& sm) {
    sm.estimation_step = true;
    if (sm.outliers.empty()) {
        return 1;
    }
    if (sm.tstats_null) {
        return -1; // the NPE at tstats.length
    }
    if (sm.tstats.size() < sm.outliers.size()) {
        return -1; // JVM arraylength semantics guard (defensive)
    }
    const std::size_t base = sm.tstats.size() - sm.outliers.size();
    std::size_t best = 0;
    for (std::size_t i = 1; i < sm.outliers.size(); ++i) {
        const double a = std::abs(sm.tstats[base + i]);
        const double b = std::abs(sm.tstats[base + best]);
        if (std::isless(a, b)) { // dcmpg ifge: a NaN candidate loses
            best = i;
        }
    }
    const double tmin = std::abs(sm.tstats[base + best]);
    // dcmpl |t|, cv; iflt: |t| < cv OR an unordered (NaN) operand takes the
    // removal path, only an ordered |t| >= cv verifies (offsets 95-100).
    const bool removal_path
        = std::isless(tmin, sm.cv) || std::isnan(tmin) || std::isnan(sm.cv);
    if (!removal_path) {
        return 1;
    }
    sm.backward_step = false;
    sm.estimation_step = false;
    const std::array<int, 2> removed = sm.outliers[best];
    allow_cell(sm, removed[0], removed[1]);
    if (!remove_outlier(sm, best)) {
        return -1; // the NPE at coeff.length / the index throw
    }
    const bool same = !sm.lastremoved_null && sm.lastremoved == removed;
    sm.lastremoved = removed;
    sm.lastremoved_null = false;
    sm.exit = sm.exit || same; // offsets 141-161: exit := equals(lastremoved)
    return 0;
}

// ---------------------------------------------------------------------------
// The probe-shaped trace (public types of kalos/tramo_outlier_rounds.hpp).

inline TramoOutlierRoundState dump_round_state(const ScanMachine& sm) {
    TramoOutlierRoundState d;
    d.round = sm.round;
    d.rflag = sm.rflag;
    d.estimation_step = sm.estimation_step;
    d.backward_step = sm.backward_step;
    d.exit = sm.exit;
    d.outliers = sm.outliers;
    d.lastremoved = !sm.lastremoved_null;
    d.coeff = sm.coeff;
    d.coeff_null = sm.coeff_null;
    d.tstats = sm.tstats;
    d.tstats_null = sm.tstats_null;
    d.res = sm.res;
    d.res_null = sm.res_null;
    d.pars = sm.jvm;
    d.xcount = sm.user_raw.size();
    return d;
}

// The full driver (FastOutliersDetector.process, offsets 0-241, ported
// through the Scratch23g manual-run shape so the trace records line up with
// the witness): see the file header for the loop map. `note` values follow
// the witness grammar. JVM exception paths (the null-field reads, a failing
// initmodel, the searchMax table read) end the run with ret=false and no
// tail, matching the catch at 238.
inline void process_rounds(
    ScanMachine& sm, std::vector<TramoOutlierRoundEvent>& trace) {
    const auto emit = [&trace, &sm](const char* tag) -> TramoOutlierRoundEvent& {
        trace.push_back(TramoOutlierRoundEvent{});
        TramoOutlierRoundEvent& ev = trace.back();
        ev.tag = tag;
        ev.has_state = true;
        ev.state = dump_round_state(sm);
        return ev;
    };
    emit("init");
    bool ret = false;
    std::string note = "run";
    bool reached_tail = false;
    int rounds = 0;
    for (;;) {
        // loop head, offset 19: estimateModel.
        if (!estimate_model(sm)) {
            note = "est"; // offsets 24-28: return false, no tail
            break;
        }
        emit("est");
        bool ans = true;
        if (sm.backward_step) {
            const int v = verify_model(sm);
            if (v < 0) {
                note = "verify-throw"; // the JVM catch at 238
                break;
            }
            ans = v != 0;
            emit("verify");
            if (sm.exit) { // offsets 45-52
                note = "exit";
                reached_tail = true;
                break;
            }
        }
        if (ans) {
            // offsets 60-74: the scan.
            const RoundScan rs = scan_process(sm);
            if (rs.status == RoundScan::Status::Throw) {
                note = "initthrow"; // the JVM catch at 238: no dumps, no tail
                break;
            }
            TramoOutlierRoundEvent& ev = emit("scan");
            ev.has_sod = true;
            TramoOutlierRoundScanDump& sod = ev.sod;
            sod.el_null = rs.status == RoundScan::Status::CleanFalse;
            sod.el_count = rs.el.size();
            sod.mad = rs.mad;
            sod.lower_bound = sm.lower_bound;
            sod.upper_bound = sm.upper_bound;
            sod.t_table_present = true; // allocated by prepare(n) in make
            sod.allowed_present = sm.m > 0;
            sod.allowed_rows = sm.m;
            sod.allowed_cols = sm.factories.size();
            for (std::size_t o = 0; o < sm.factories.size(); ++o) {
                std::size_t cnt = 0;
                for (std::size_t row = 0; row < sm.m; ++row) {
                    if (sm.zone_ok[row + o * sm.m]
                        && !sm.excluded[row + o * sm.m]) {
                        ++cnt;
                    }
                }
                sod.allowed_counts.push_back(cnt);
            }
            sod.maxt_threw = !rs.search_ok;
            sod.max_t = rs.max_t;
            sod.pick = rs.search_ok;
            sod.pick_position = rs.position;
            sod.pick_type = rs.kind;
            sod.coefficient = rs.coefficient;
            if (rs.status == RoundScan::Status::CleanFalse) {
                note = "sodfalse";
                reached_tail = true;
                break;
            }
            if (!rs.search_ok) {
                note = "maxtthrew"; // the (-1, -1) table read, caught at 238
                break;
            }
            // offsets 77-84: round++
            ++sm.round;
            // offsets 87-109: dcmpg(|maxT|, cv) < 0 breaks to the tail; an
            // unordered maxT takes the ADD path (dcmpg returns +1).
            if (std::isless(std::abs(rs.max_t), sm.cv)) {
                note = "gate";
                reached_tail = true;
                break;
            }
            // offsets 112-146: addOutlier(position, type, coefficient).
            if (!add_outlier(sm, rs.position, rs.kind, rs.coefficient)) {
                note = "addthrow"; // the JVM catch at 238: no dumps, no tail
                break;
            }
            TramoOutlierRoundEvent& av = emit("add");
            av.add_position = rs.position;
            av.add_type = rs.kind;
            av.add_coefficient = rs.coefficient;
            // offsets 149-163: outliers.size() == maxOutliers breaks to the
            // tail.
            if (static_cast<int>(sm.outliers.size()) == sm.max_outliers) {
                note = "cap";
                reached_tail = true;
                break;
            }
        }
        // offsets 166-174: round < maxRound back to 19, else the tail.
        if (sm.round >= sm.max_round) {
            note = "roundcap";
            reached_tail = true;
            break;
        }
        if (++rounds > 400) { // defensive pin against a spinning machine
            note = "stepcap";
            break;
        }
    }
    if (reached_tail) {
        // offsets 177-214: the tail estimateModel on exit / round ==
        // maxRound / size == maxOutliers (result popped).
        const bool tail_cond = sm.exit || sm.round == sm.max_round
            || static_cast<int>(sm.outliers.size()) == sm.max_outliers;
        if (tail_cond) {
            estimate_model(sm);
            emit("tailest"); // the witness pairs the bare marker with a dump
        }
        sm.estimation_step = false; // offset 215
        // offsets 220-233: while (!verifyModel()) estimateModel (popped).
        int drain = 0;
        for (;;) {
            const int v = verify_model(sm);
            if (v < 0) {
                note = "verify-throw";
                reached_tail = false;
                break;
            }
            if (v != 0) {
                break;
            }
            emit("drainv");
            estimate_model(sm);
            emit("draine");
            if (++drain > 400) {
                note = "draincap";
                reached_tail = false;
                break;
            }
        }
        if (reached_tail) {
            emit("drainend");
            ret = true;
        }
    }
    trace.push_back(TramoOutlierRoundEvent{});
    TramoOutlierRoundEvent& mv = trace.back();
    mv.tag = "manret";
    mv.manret_returned = ret;
    mv.note = note;
    mv.manret_round = sm.round;
    trace.push_back(TramoOutlierRoundEvent{});
    TramoOutlierRoundEvent& mo = trace.back();
    mo.tag = "manouts";
    mo.state.outliers = sm.outliers; // the final set rides in state.outliers
}

} // namespace kalos::detail
#endif
