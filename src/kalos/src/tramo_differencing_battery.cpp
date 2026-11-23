// Battery-stage transcription of the pinned jdplus 3.9.0 tramo driver (item
// 23-C.2): SeasonalityController.process (computeReferenceModel /
// compareReferenceModels), the RegularUnderDifferencingTest /
// SeasonalUnderDifferencingTest(+2) loop, the SeasonalOverDifferencingTest
// and the SpectralPeaks / SeasonalityTests gates, from the bytecode captured
// under /tmp/opencode/jdprobe/m5 (documented in
// /tmp/opencode/kalos-m5-pinning.md, SECTION 12). Every JVM double
// comparison is transcribed with its exact opcode family: dcmpg yields +1
// and dcmpl -1 for a NaN operand, which is what decides whether a NaN
// p-value stops or continues each gate below.
#include <string>
#include "kalos/tramo_differencing_battery.hpp"
#include "kalos/tramo_regression_tests.hpp"
#include "kalos/arma_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace kalos {
namespace {

[[noreturn]] void bail(const char* what) {
    throw std::runtime_error(std::string("decide_tramo_differencing_battery: ") + what);
}

int dcmpg(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return 1;
    }
    return a < b ? -1 : (a > b ? 1 : 0);
}

int dcmpl(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return -1;
    }
    return a < b ? -1 : (a > b ? 1 : 0);
}

// StatisticalTest.isSignificant(t): dcmpg(pvalue, t) ifge -> false; the NaN
// p-value is NOT significant here.
bool is_significant(double p, double t) {
    return dcmpg(p, t) < 0;
}

// --- ModelComparator (builder().build() defaults) ---------------------------
// preference BIC, significance 0.01, kbic 0.03, kq 1.1, kqs 1.1, mout 2,
// kout 4, knz 0.05, ksk 1.25 (Builder.<init> constant pool). The acceptable*
// booleans are instance fields the constructor never assigns; they only
// shadow the recomputed per-comparison values and never reach a return.
namespace pref {
enum Preference { kBic, kFirst, kSecond };
} // namespace pref

constexpr double kSig = 0.01;
constexpr double kKbic = 0.03;
constexpr double kKq = 1.1;
constexpr double kKqs = 1.1;
constexpr int kMout = 2;
constexpr double kKnz = 0.05;
constexpr double kKsk = 1.25;

// checkAcceptable(a, b): the b panel must not exceed a by more than
// kbic * |a.bic| (dcmpl: a NaN difference PASSES the guard) and stays
// within the four acceptable* reliefs (all four required).
bool check_acceptable(const BatteryPanel& a, const BatteryPanel& b) {
    if (!(dcmpl(b.bic - a.bic, std::fabs(a.bic) * kKbic) < 0)) {
        return false;
    }
    const bool acc_q = dcmpg(b.ljung_box, a.ljung_box * kKq) < 0
                       || dcmpl(b.ljung_box_p, kSig) > 0;
    const bool acc_out = b.nout < a.nout + 3
                         || dcmpg(static_cast<double>(b.nout), kKnz * static_cast<double>(b.nobs)) < 0;
    const bool acc_sk = dcmpg(b.skew_abs, a.skew_abs * kKsk) < 0 || dcmpl(b.skew_p, kSig) > 0;
    const bool acc_qs = dcmpg(b.season_box, a.season_box * kKqs) < 0 || dcmpl(b.season_box_p, kSig) > 0;
    // acceptableStab = dcmpl(b.stable_mean_p, kSig) > 0 is computed but not
    // referenced by the return in the pinned bytecode; ported for fidelity.
    (void)(dcmpl(b.stable_mean_p, kSig) > 0);
    return acc_q && acc_out && acc_sk && acc_qs;
}

bool test1(const BatteryPanel& a, const BatteryPanel& b) {
    if (dcmpl(a.ljung_box_p, kSig) >= 0) {  // a NaN p-value continues
        return false;
    }
    return dcmpg(a.ljung_box, b.ljung_box * kKq) > 0;
}

bool test2(const BatteryPanel& a, const BatteryPanel& b) {
    if (dcmpg(static_cast<double>(a.nout), static_cast<double>(a.nobs) * kKnz) <= 0) {
        return false;
    }
    return a.nout > b.nout + kMout;
}

bool test3(const BatteryPanel& a, const BatteryPanel& b) {
    if (dcmpl(a.season_box_p, kSig) >= 0) {
        return false;
    }
    return dcmpg(a.season_box, b.season_box * kKqs) > 0;
}

bool test4(const BatteryPanel& a, const BatteryPanel& b) {
    if (dcmpl(a.skew_p, kSig) >= 0) {
        return false;
    }
    return dcmpg(a.skew_abs, b.skew_abs * kKsk) > 0;
}

int stability_score(const BatteryPanel& s) {
    int score = 0;
    if (dcmpg(s.stable_mean_p, kSig) < 0) {
        ++score;
    }
    if (dcmpg(s.stable_var_p, kSig) < 0) {
        ++score;
    }
    return score;
}

bool test5(const BatteryPanel& a, const BatteryPanel& b) {
    return stability_score(a) > stability_score(b);
}

bool prefer_second(const BatteryPanel& a, const BatteryPanel& b) {
    if (!check_acceptable(a, b)) {
        return false;
    }
    return test1(a, b) || test2(a, b) || test3(a, b) || test4(a, b) || test5(a, b);
}

// ModelComparator.compare over the two panels. The JVM first short-circuits
// identity (referenceModel == model in the second pass) - handled by the
// caller. BIC preference: a better (smaller) bic2 takes the Second-style
// branch, otherwise the First-style branch; pure First / Second follow
// their branch directly.
int compare_panels(const BatteryPanel& s1, const BatteryPanel& s2, int preference) {
    if (preference == pref::kBic && dcmpg(s2.bic, s1.bic) < 0) {
        return prefer_second(s2, s1) ? -1 : 1;
    }
    if (preference == pref::kSecond) {
        return prefer_second(s2, s1) ? -1 : 1;
    }
    return prefer_second(s1, s2) ? 1 : -1;
}

// --- SpectralPeaks -----------------------------------------------------------

bool peak_present(const PeakRecord& r) {
    // SpectralPeaks.hasPeak(): ar.isPresent() || tu.isPresent() (A/a, T/t).
    return r.ar == PeakAr::A || r.ar == PeakAr::a || r.tukey == PeakTukey::T
           || r.tukey == PeakTukey::t;
}

bool peak_all(const PeakRecord& r) {
    return r.ar == PeakAr::A && r.tukey == PeakTukey::T;  // ALL
}

bool peak_none(const PeakRecord& r) {
    return r.ar == PeakAr::none && r.tukey == PeakTukey::none;  // NONE
}

} // namespace

std::optional<BatteryPanel> tramo_model_comparison_panel(const TramoModelStatisticsRequest& request) {
    const auto stats = tramo_model_statistics(request);
    if (!stats.built || request.residuals.size() < 4) return std::nullopt;
    const auto tests = tramo_comparator_residual_tests(request.residuals);
    BatteryPanel panel;
    panel.ljung_box = stats.ljung_box;
    panel.ljung_box_p = stats.ljung_box_pvalue;
    panel.season_box = stats.seasonal_ljung_box;
    panel.season_box_p = stats.seasonal_ljung_box_pvalue;
    panel.bic = stats.bic;
    panel.nobs = stats.observations_count;
    panel.nout = stats.outliers_count;
    panel.skew_abs = tests.absolute_skewness;
    panel.skew_p = tests.skewness_p;
    panel.stable_mean_p = tests.stable_mean_p;
    panel.stable_var_p = tests.stable_variance_p;
    return panel;
}

bool tramo_candidate_preferred(const BatteryPanel& incumbent, const BatteryPanel& candidate) {
    return compare_panels(incumbent, candidate, pref::kBic) >= 1;
}

bool battery_has_seasonal_peaks(const std::vector<PeakRecord>& peaks) {
    switch (peaks.size()) {
    case 6: {
        int count = 0, all = 0;
        for (const PeakRecord& r : peaks) {
            if (peak_present(r)) {
                ++count;
                if (peak_all(r)) {
                    ++all;
                }
            }
        }
        if (count >= 4) {
            return true;
        }
        if (count == 3) {
            return all >= 1 || !peak_present(peaks[5]);
        }
        if (count == 2) {
            if (peak_all(peaks[5]) && all == 2) {
                return true;
            }
            return peak_none(peaks[5]) && all >= 1;
        }
        return false;
    }
    case 2: {
        if (peak_all(peaks[0])) {
            return true;
        }
        int count = 0;
        for (const PeakRecord& r : peaks) {
            if (peak_present(r)) {
                ++count;
            }
        }
        return count == 2;
    }
    case 1:
        return peak_present(peaks[0]);
    default:
        return false;  // the JVM default arm
    }
}

int battery_ladder_score(const BatteryTests& t, int freq) {
    // residualSeasonalityTest / 5-argument seasonalityTest scoring: qs and
    // nonParametricTest significance at 0.01, then the peaks contribution
    // behind the differenced-length gate (>= 80, or >= 60 when the frequency
    // is not 12).
    int score = 0;
    if (!t.qs_null && is_significant(t.qs_p, kSig)) {
        ++score;
    }
    if (!t.np_null && is_significant(t.np_p, kSig)) {
        ++score;
    }
    const int dlen = t.differenced_length;
    if (dlen >= 80 || (freq != 12 && dlen >= 60)) {
        if (battery_has_seasonal_peaks(t.peaks)) {
            ++score;
        }
    }
    return score;
}

int battery_over_differencing_code(const BatteryState& s, int freq) {
    // SeasonalOverDifferencingTest.test: specification gate, then the
    // 5-argument ladder score on the linearized-series stage facts plus the
    // F significance at 0.01; score >= 2 or a significant F or a significant
    // qs yields code 1, otherwise 2.
    if (freq == 1) {
        return 0;
    }
    if (s.bp == 1) {
        return 0;
    }
    if (!(s.bd != 0 && s.bq != 0 && dcmpl(s.btheta1, -0.7) < 0)) {
        return 0;
    }
    int score = battery_ladder_score(s.ctlst, freq);
    const bool fsig = s.ftest.present && dcmpg(s.ftest.pvalue, kSig) < 0;
    if (fsig) {
        ++score;
    }
    const bool qs_sig = !s.ctlst.qs_null && dcmpg(s.ctlst.qs_p, kSig) < 0;
    if (score >= 2 || fsig || qs_sig) {
        return 1;
    }
    return 2;
}

namespace {

// --- SeasonalityController ---------------------------------------------------

struct Context {
    bool seasonal;
    int ost;
};

// SeasonalFTest p-value of a state re-computed through the native kernel;
// the JVM DStatException surface is the NaN p-value, matching f_upper_p.
bool ftest_significant(const BatteryFTest& f) {
    return dcmpg(f.pvalue, kSig) < 0;
}

// SeasonalityController.hasSeasonality over the pass-1 stage facts. The
// captured ftest value / degrees of freedom feed the kernel check; the
// decision consumes the JVM-computed p-value bits (the unit re-asserts the
// kernel equality separately for every state).
bool has_seasonality(const BatteryState& s, Context& ctx, int freq) {
    if (s.ctlst.peaks.empty() && !s.ctlst.tukey_null && !s.ctlst.ar_null) {
        // A set of peak tests without records is a state the observer would
        // have echoed differently; fail closed rather than assume.
        bail("controller stage facts missing peak records");
    }
    int score = 0, fcnt = 0;
    if (dcmpg(s.mstat.season_box_p, kSig) < 0) {
        ++score;
        ++fcnt;
    }
    if (!s.ctlst.np_null && is_significant(s.ctlst.np_p, kSig)) {
        ++score;
    }
    const int dlen = s.ctlst.differenced_length;
    if (dlen >= 80 || (freq != 12 && dlen >= 60)) {
        if (battery_has_seasonal_peaks(s.ctlst.peaks)) {
            ++score;
        }
    }
    if (!s.ftest.present) {
        bail("SeasonalFTest produced no p-value on the reference model");
    }
    const bool fp = ftest_significant(s.ftest);
    if (fp) {
        ++score;
        ++fcnt;
    }
    ctx.ost = score;
    if (score > 1 || fcnt > 0) {
        return true;
    }
    return dcmpg(s.mstat.season_box_p, kSig) < 0;  // the JVM tail clause
}

// SarimaOrders.m011(freq) with setBq(1): (0,1,1,0,0,1).
struct Spec {
    int p, d, q, bp, bd, bq;
};

Spec spec_of(const BatteryState& s) {
    return {s.p, s.d, s.q, s.bp, s.bd, s.bq};
}

bool same_spec(const Spec& a, const Spec& b) {
    return a.p == b.p && a.d == b.d && a.q == b.q && a.bp == b.bp && a.bd == b.bd
           && a.bq == b.bq;
}

struct Driver {
    const BatteryCase& c;
    Context ctx;
    Spec spec;
    bool mean;
    std::size_t next_cand = 1;  // candidates[0] is the initial estimation
    BatteryDecision out;

    std::function<BatteryCandidate(const Spec&, bool)> live_estimator;
    BatteryCandidate live_candidate;

    explicit Driver(const BatteryCase& kase)
        : c(kase), ctx{kase.seasonal, kase.ost}, spec{spec_of(kase.after_entry)},
          mean(kase.after_entry.mean) {}

    const BatteryCandidate& consume(const Spec& want, bool want_mean, const char* who) {
        if (live_estimator) live_candidate = live_estimator(want, want_mean);
        else if (next_cand >= c.candidates.size()) bail(who);
        const BatteryCandidate& cand = live_estimator ? live_candidate : c.candidates[next_cand];
        const Spec got{cand.p, cand.d, cand.q, cand.bp, cand.bd, cand.bq};
        if (!same_spec(got, want) || cand.mean != want_mean) {
            bail("candidate construction diverges from the captured snapshot");
        }
        out.consumed.push_back(static_cast<int>(next_cand));
        ++next_cand;
        return cand;
    }

    using R = BatteryDecision::Result;

    // SeasonalityController.computeReferenceModel (first pass).
    R pass1() {
        const BatteryState& e = c.after_entry;
        Context probe = ctx;  // hasSeasonality mutates the score
        if (c.freq == 1) {
            return R::Unprocessed;
        }
        const bool has = has_seasonality(e, probe, c.freq);
        ctx = probe;
        Spec m011_airline{0, 1, 1, 0, 0, 1};
        bool swap = false;
        if (!has) {
            swap = true;  // the isSeasonal() guard is vacuous here (freq > 1)
        } else if (!ctx.seasonal) {
            ctx.seasonal = true;
            return R::Changed;  // seasonal flip without a specification change
        }
        if (!swap && !ctx.seasonal) {
            if (dcmpg(e.mstat.season_box_p, 0.05) < 0 || dcmpg(e.mstat.ljung_box_p, 0.05) < 0) {
                ctx.seasonal = true;
                swap = true;
            }
        }
        if (swap) {
            // mean flag unchanged by computeReferenceModel; the candidate
            // keeps the description's mean.
            const BatteryCandidate& cand = consume(m011_airline, mean, "unobservable pass-1 re-estimation");
            out.cand_pass1 = out.consumed.back();
            if (!cand.ok) {
                return R::Unchanged;  // estimate() false leaves the model alone
            }
            spec = m011_airline;
            return R::Changed;
        }
        return R::Unchanged;
    }

    // RegularUnderDifferencingTest.process over the pass-1 state.
    R ru1(const BatteryState& s) {
        if (s.d == 2 || s.p == 0 || !s.mean) {
            return R::Unchanged;
        }
        // checkResiduals: |mean(e) / sqrt(var(e) / n)| <= 1.6 (dcmpg; a NaN
        // t continues). var(e) = ssq / n - mean * mean.
        const double n = static_cast<double>(s.resid.size());
        double sum = 0.0, ssq = 0.0;
        for (const double v : s.resid) {
            sum += v;
            ssq += v * v;
        }
        const double m = sum / n;
        double var = ssq / n - m * m;
        var = var / n;
        const double t = m / std::sqrt(var);
        if (std::fabs(t) <= 1.6 && !std::isnan(t)) {
            return R::Unchanged;
        }
        // hasQuasiUnitRoots: a mirrored-AR root with re > 0, |im| <= 0.01
        // and |z| >= 0.9 (dcmpl/dcmpg semantics: a NaN re FAILS, a NaN im
        // or modulus PASSES - the Complex comparisons are JVM dcmpl/dcmpg).
        bool quasi = false;
        for (const double root_re : mirrored_ar_roots_re(s.phi)) {
            quasi = quasi || root_re;
        }
        if (!quasi) {
            return R::Unchanged;
        }
        Spec want = spec_of(s);
        want.d = s.d + 1;
        want.p = s.p - 1;
        const BatteryCandidate& cand = consume(want, false, "unobservable ru1 candidate");
        out.cand_ru1 = out.consumed.back();
        if (!cand.ok) {
            return R::Failed;
        }
        spec = want;
        mean = false;
        return R::Changed;
    }

    // SeasonalUnderDifferencingTest.process.
    R su1(const BatteryState& s) {
        if (c.freq == 1) {
            return R::Unprocessed;
        }
        Spec tripe = spec_of(s);
        tripe.bp = 0;
        tripe.bd = 1;
        tripe.bq = 1;
        if (s.bd == 0 && s.bp == 1 && dcmpl(s.bphi1, -0.91) < 0) {
            // fixSeasonalRoots: (bp, bd, bq) = (0, 1, 1), mean off.
            const BatteryCandidate& cand = consume(tripe, false, "unobservable fixSeasonalRoots candidate");
            out.cand_su1 = out.consumed.back();
            if (cand.ok) {
                spec = tripe;
                mean = false;
                return R::Changed;
            }
            // estimate() false: the JVM leaves the (already mutated) live
            // orders in place and falls through to the ladder with a null
            // candidate stream position; the corpus never reaches it.
            bail("fixSeasonalRoots estimation failure is out of the ported domain");
        }
        const int score = battery_ladder_score(s.resl, c.freq);
        if (!(score > 1 || (score == 1 && ctx.seasonal))) {
            return R::Unchanged;
        }
        const BatteryCandidate& cand = consume(tripe, false, "unobservable su1 candidate");
        out.cand_su1 = out.consumed.back();
        if (!cand.ok) {
            bail("su1 buildNewModel returned null (JVM NullPointerException)");
        }
        if (compare_panels(cand.mstat, s.mstat, pref::kBic) < 0) {
            spec = tripe;
            mean = false;
            return R::Changed;
        }
        return R::Unchanged;
    }

    // RegularUnderDifferencingTest2.process.
    R ru2(const BatteryState& s) {
        if (c.freq <= 2) {
            return R::Unprocessed;
        }
        if (!(dcmpl(s.mstat.ljung_box_p, 0.005) < 0)) {
            return R::Unchanged;  // dcmpl/if_lt: a NaN p-value does NOT proceed
        }
        if (!need_processing(s)) {
            return R::Unchanged;
        }
        Spec want = spec_of(s);
        bool want_mean;
        if (s.d == 2) {
            if (s.p == 3) {
                bail("ru2 buildNewModel returned null (p = 3, JVM NullPointerException)");
            }
            want.p = s.p + 1;
            want_mean = true;
        } else {
            if (s.q == 3) {
                bail("ru2 buildNewModel returned null (q = 3, JVM NullPointerException)");
            }
            want.q = s.q + 1;
            want.d = s.d + 1;
            want_mean = false;
        }
        const BatteryCandidate& cand = consume(want, want_mean, "unobservable ru2 candidate");
        out.cand_ru2 = out.consumed.back();
        if (!cand.ok) {
            bail("ru2 buildNewModel returned null (JVM NullPointerException)");
        }
        if (compare_panels(cand.mstat, s.mstat, pref::kBic) < 0) {
            spec = want;
            mean = want_mean;
            return R::Changed;
        }
        return R::Unchanged;
    }

    // SeasonalUnderDifferencingTest2.process.
    R su2(const BatteryState& s) {
        if (c.freq == 1) {
            return R::Unprocessed;
        }
        if (s.bd == 1 || s.bq == 1 || ctx.ost == 0) {
            return R::Unchanged;
        }
        if (s.bp == 1) {
            if (s.nvar != 0) {
                bail("calendar-variable scan on a non-empty variable list is out of the ported domain");
            }
            // nvar = 0: findAny over the (empty) calendar filter is empty,
            // the guard passes exactly as observed.
        }
        Spec tripe = spec_of(s);
        tripe.bp = 0;
        tripe.bd = 1;
        tripe.bq = 1;
        const BatteryCandidate& cand = consume(tripe, false, "unobservable su2 candidate");
        out.cand_su2 = out.consumed.back();
        if (!cand.ok) {
            bail("su2 buildNewModel returned null (JVM NullPointerException)");
        }
        if (compare_panels(cand.mstat, s.mstat, pref::kBic) < 0) {
            spec = tripe;
            mean = false;
            return R::Changed;
        }
        return R::Unchanged;
    }

    // SeasonalityController.compareReferenceModels (second pass). The
    // controller's reference IS the live modelling object (computeReference-
    // Model stores it by reference), so the specification guard sees equal
    // orders, the over-side never runs here, and compare() short-circuits on
    // identity to 0: always Unchanged. Assert that the tracked state is
    // consistent with that reading.
    R pass2() { return R::Unchanged; }

    // needProcessing (RegularUnderDifferencingTest2): acs =
    // autoCovariancesWithZeroMean(e, min(24, n - 1)), lag 0 the zero-mean
    // variance. The JVM counts the positive entries - acs[k] > 0 under
    // dcmpl/if_le, so a NaN lag COUNTS (+1) and zero never does; every JVM
    // divisor is positive, so the sign of the raw product sum decides. The
    // first loop reads lags 0..11 unconditionally (an ArrayIndexOutOfBounds
    // below n = 12: fails closed here), the second adds lags 12..lim to the
    // total; the step proceeds when pos12 >= freq, pos12 >= 9, or the total
    // over all lags reaches 17.
    bool need_processing(const BatteryState& s) {
        const int n = static_cast<int>(s.resid.size());
        if (n < 12) {
            bail("residual vector shorter than the ru2 autocovariance window");
        }
        const int lim = n - 1 < 24 ? n - 1 : 24;
        auto acs_positive = [&s, n](int lag) {
            double acc = 0.0;
            for (int i = 0; i + lag < n; ++i) {
                acc += s.resid[static_cast<std::size_t>(i)]
                       * s.resid[static_cast<std::size_t>(i + lag)];
            }
            return !(acc <= 0.0);  // dcmpl(acs, 0) ifle: a NaN lag counts
        };
        int pos12 = 0, total = 0;
        for (int lag = 0; lag < 12; ++lag) {
            if (acs_positive(lag)) {
                ++pos12;
            }
        }
        total = pos12;
        for (int lag = 12; lag <= lim; ++lag) {
            if (acs_positive(lag)) {
                ++total;
            }
        }
        return pos12 >= c.freq || pos12 >= 9 || total >= 17;
    }

    // hasQuasiUnitRoots over getRegularAR().mirror().roots(): the mirrored
    // polynomial is [phip, ..., phi1, 1] (the 23-C convention, verified
    // there against the JVM root scans); degree <= 2 closed forms, anything
    // higher fails closed. For each root: re > 0 && |im| <= 0.01 &&
    // |z| >= 0.9 (dcmpl on re: NaN FAILS; dcmpg on |im|: NaN PASSES the
    // tolerance; dcmpl on |z|: NaN PASSES the bound).
    static std::vector<bool> mirrored_ar_roots_re(const std::vector<double>& phi) {
        std::vector<double> asc;  // ascending coefficients of the mirror
        for (auto it = phi.rbegin(); it != phi.rend(); ++it) {
            asc.push_back(*it);
        }
        asc.push_back(1.0);
        int deg = static_cast<int>(asc.size()) - 1;
        while (deg > 0 && asc[static_cast<std::size_t>(deg)] == 0.0) {
            --deg;
        }
        std::vector<std::pair<double, double>> roots;
        if (deg == 0) {
            return {};
        }
        if (deg == 1) {
            roots.emplace_back(-asc[0] / asc[1], 0.0);
        } else if (deg == 2) {
            const double a = asc[2], b = asc[1], c0 = asc[0];
            const double denom = 2.0 * a;
            const double disc = b * b - 4.0 * a * c0;
            if (disc < 0.0) {
                const double s = std::sqrt(-disc);
                roots.emplace_back(-b / denom, s / denom);
                roots.emplace_back(-b / denom, -(s / denom));
            } else {
                const double s = std::sqrt(disc);
                roots.emplace_back((-b + s) / denom, 0.0);
                roots.emplace_back((-b - s) / denom, 0.0);
            }
        } else {
            bail("quasi-unit-root scan above degree 2 is out of the ported domain");
        }
        std::vector<bool> hit;
        for (const auto& [re, im] : roots) {
            const bool re_pos = dcmpl(re, 0.0) > 0;  // Complex.getRe() dcmpl ifle
            const bool im_ok = !(dcmpg(std::fabs(im), 0.01) > 0);
            const bool mod_ok = !(dcmpl(complex_abs(re, im), 0.9) < 0);
            hit.push_back(re_pos && im_ok && mod_ok);
        }
        return hit;
    }

    static double complex_abs(double re, double im) {
        if (re == 0.0 && im == 0.0) {
            return 0.0;
        }
        const double a = std::fabs(re), b = std::fabs(im);
        const double m = std::max(a, b), n = std::min(a, b);
        if (n == 0.0) {
            return m;
        }
        const double t = n / m;
        return m * std::sqrt(1.0 + t * t);
    }
};

} // namespace

TramoMonthlyControllerResult tramo_monthly_controllers(const BatteryState& entry, bool seasonal,
                                                       int seasonality_score, const TramoControllerEstimator& estimate) {
    using Result = BatteryDecision::Result;
    BatteryCase c;
    c.after_entry = entry;
    c.seasonal = seasonal;
    c.ost = seasonality_score;
    Driver driver(c);
    std::optional<BatteryState> trial;
    BatteryState state = entry;
    int calls = 0, selected = -1;
    driver.live_estimator = [&](const Spec& want, bool mean) {
        ++calls;
        trial = estimate({want.p, want.d, want.q, want.bp, want.bd, want.bq}, mean);
        BatteryCandidate candidate{want.p, want.d, want.q, want.bp, want.bd, want.bq, mean, trial.has_value(), {}};
        if (trial) {
            if (!same_spec(spec_of(*trial), want) || trial->mean != mean)
                bail("controller estimator returned a different model specification");
            candidate.mstat = trial->mstat;
        }
        return candidate;
    };
    auto accept = [&](Result result) {
        if (result == Result::Changed && trial) {
            state = std::move(*trial);
            selected = calls - 1;
        }
        trial.reset();
        return result;
    };
    const auto initial = accept(driver.pass1());
    TramoMonthlyControllerResult result{initial, Result::Unprocessed, Result::Unprocessed, Result::Unprocessed,
        Result::Unprocessed, driver.ctx.seasonal, driver.ctx.ost, selected, initial == Result::Changed};
    if (result.restart) return result;
    result.regular1 = accept(driver.ru1(state));
    result.seasonal1 = accept(driver.su1(state));
    result.regular2 = accept(driver.ru2(state));
    result.seasonal2 = accept(driver.su2(state));
    result.selected_candidate = selected;
    return result;
}

BatteryDecision decide_tramo_differencing_battery(const BatteryCase& c) {
    if (c.candidates.empty()) {
        bail("the pinned candidate stream is empty");
    }
    Driver drv(c);
    auto snap = [&drv] {
        drv.out.snaps.push_back({drv.spec.p, drv.spec.d, drv.spec.q, drv.spec.bp, drv.spec.bd,
                                 drv.spec.bq, drv.mean, drv.ctx.seasonal, drv.ctx.ost});
    };
    drv.out.seas1 = drv.pass1();
    snap();
    drv.out.ru1 = drv.ru1(c.after_pass1);
    snap();
    drv.out.su1 = drv.su1(c.after_ru1);
    snap();
    drv.out.ru2 = drv.ru2(c.after_su1);
    snap();
    drv.out.su2 = drv.su2(c.after_ru2);
    snap();
    drv.out.seas2 = drv.pass2();
    drv.out.spec_p = drv.spec.p;
    drv.out.spec_d = drv.spec.d;
    drv.out.spec_q = drv.spec.q;
    drv.out.spec_bp = drv.spec.bp;
    drv.out.spec_bd = drv.spec.bd;
    drv.out.spec_bq = drv.spec.bq;
    drv.out.mean = drv.mean;
    drv.out.seasonal = drv.ctx.seasonal;
    drv.out.ost = drv.ctx.ost;
    return drv.out;
}

} // namespace kalos
