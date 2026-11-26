import jdplus.toolkit.base.api.arima.SarmaOrders;
import jdplus.toolkit.base.api.arima.SarimaOrders;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.arima.IArimaModel;
import jdplus.toolkit.base.core.arima.estimation.ArmaFilter;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.regarima.RegArimaModel;
import jdplus.toolkit.base.core.regarima.RegArmaModel;
import jdplus.toolkit.base.core.regarima.estimation.ConcentratedLikelihoodComputer;
import jdplus.toolkit.base.core.sarima.SarimaModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

/**
 * Item-29 probe slice: instruments the JVM ARMA-likelihood filter machine that
 * the pinned TRAMO/X-13 estimator chain runs per objective evaluation.
 *
 * Production semantics (bytecode-cited in the dossier): ConcentratedLikelihood
 * Computer.compute(RegArmaModel) constructs a fresh internal KalmanFilter
 * (multiUse = nx > 0) per call; prepare() primes C0 from the analytic
 * autocovariance function and sweeps the gain history (calcC when multiUse),
 * detecting steadyPos; apply() whitens y and every regressor column, reusing
 * the swept rows while t < steadyPos.
 *
 * This probe reproduces that per-evaluation lifecycle (fresh CLC + fresh
 * KalmanFilter per evaluation, wrapped in a counting delegate) at pinned
 * model coefficients, reads the private machine state reflectively (steadyPos,
 * dim, n, s[], C0, the frozen steady row), and separates deterministic output
 * lines (double-capture byte-identical) from TIMING lines.
 *
 * args: caseId yadj.csv x.csv|NULL "p,d,q,P,D,Q,s" "ar;ma;sar;sma" evals
 *       Coefficients follow the kalos model convention
 *       (1-sum ar B^i)(1-sum sar B^{s i}) y = (1+sum ma B^i)(1+sum sma B^{s i}) e;
 *       the sign calibration against this model's internal convention is
 *       settled by the LL line in the capture, not assumed.
 */
public class Item29FilterProbe {

    static final class Probe implements ArmaFilter {
        final ArmaFilter real;
        long prepareNanos, applyNanos, detNanos;
        int applyCalls;

        Probe(boolean multiUse) {
            real = ArmaFilter.kalman(multiUse);
        }

        @Override
        public int prepare(IArimaModel model, int n) {
            long t = System.nanoTime();
            int r = real.prepare(model, n);
            prepareNanos = System.nanoTime() - t;
            return r;
        }

        @Override
        public void apply(DoubleSeq x, jdplus.toolkit.base.core.data.DataBlock u) {
            long t = System.nanoTime();
            real.apply(x, u);
            applyNanos += System.nanoTime() - t;
            ++applyCalls;
        }

        @Override
        public double getLogDeterminant() {
            long t = System.nanoTime();
            double r = real.getLogDeterminant();
            detNanos = System.nanoTime() - t;
            return r;
        }
    }

    static double[] loadCsv(String path, int col, List<Integer> missCol) throws Exception {
        List<Double> v = new ArrayList<>();
        try (BufferedReader in = new BufferedReader(new FileReader(path))) {
            String line = in.readLine();
            while ((line = in.readLine()) != null) {
                String[] f = line.split(",");
                String s = f[col].replaceAll("\"", "").trim();
                if (s.isEmpty() || "NA".equals(s)) {
                    v.add(Double.NaN);
                } else {
                    v.add(Double.parseDouble(s));
                }
                if (missCol != null && f.length > 1 && "1".equals(f[1].trim())) {
                    missCol.add(v.size() - 1);
                }
            }
        }
        double[] a = new double[v.size()];
        for (int i = 0; i < a.length; ++i) {
            a[i] = v.get(i);
        }
        return a;
    }

    static Object rf(Object o, String name) throws Exception {
        Field f = o.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.get(o);
    }

    static String fmt(double[] v, int from, int to) {
        StringBuilder sb = new StringBuilder();
        for (int i = from; i < to; ++i) {
            sb.append(String.format(java.util.Locale.ROOT, "%.17g", v[i])).append(' ');
        }
        return sb.toString();
    }

    static double[] parseList(String s) {
        if (s.isEmpty()) {
            return new double[0];
        }
        String[] ps = s.split(",");
        double[] v = new double[ps.length];
        for (int i = 0; i < ps.length; ++i) {
            v[i] = Double.parseDouble(ps[i]);
        }
        return v;
    }

    public static void main(String[] args) throws Exception {
        String caseId = args[0];
        List<Integer> miss = new ArrayList<>();
        double[] ya = loadCsv(args[1], 0, miss);
        double[] xa = "NULL".equals(args[2]) ? new double[0] : loadCsv(args[2], 0, null);
        String[] o = args[3].split(",");
        int p = Integer.parseInt(o[0]), d = Integer.parseInt(o[1]), q = Integer.parseInt(o[2]);
        int P = Integer.parseInt(o[3]), D = Integer.parseInt(o[4]), Q = Integer.parseInt(o[5]);
        int s = Integer.parseInt(o[6]);
        String[] cparts = args[4].split(";", -1);
        double[] ar = parseList(cparts.length > 0 ? cparts[0] : "");
        double[] ma = parseList(cparts.length > 1 ? cparts[1] : "");
        double[] sar = parseList(cparts.length > 2 ? cparts[2] : "");
        double[] sma = parseList(cparts.length > 3 ? cparts[3] : "");
        int evals = Integer.parseInt(args[5]);

        SarmaOrders so = new SarmaOrders(s);
        so.setP(p);
        so.setQ(q);
        so.setBp(P);
        so.setBq(Q);
        SarimaOrders ord = SarimaOrders.of(so, d, D);
        SarimaModel.Builder b = SarimaModel.builder(ord);
        if (ar.length > 0) {
            b.phi(ar);
        }
        if (ma.length > 0) {
            b.theta(ma);
        }
        if (sar.length > 0) {
            b.bphi(sar);
        }
        if (sma.length > 0) {
            b.btheta(sma);
        }
        SarimaModel mStar = b.build();

        RegArimaModel.Builder<SarimaModel> mb = RegArimaModel.<SarimaModel>builder()
                .y(DoubleSeq.of(ya)).arima(mStar);
        if (xa.length > 0) {
            mb.addX(DoubleSeq.of(xa));
        }
        int[] mi = miss.stream().mapToInt(Integer::intValue).toArray();
        if (mi.length > 0) {
            mb.missing(mi);
        }
        RegArmaModel<SarimaModel> am = mb.build().differencedModel();
        int n = am.getY().length();
        int nx = am.getX().getColumnsCount();
        IArimaModel sta = am.getArma();

        System.out.println("CASE " + caseId + " n=" + n + " nx=" + nx
                + " nmissing=" + am.getMissingCount()
                + " orders=" + ord.toString()
                + " arma=" + String.valueOf(sta)
                + " innovvar=" + String.format(java.util.Locale.ROOT, "%.17g",
                        sta.getInnovationVariance()));

        // Autocovariance initialization cost (analytic prime of the filter).
        long t0 = System.nanoTime();
        Object acf = sta.getAutoCovarianceFunction();
        long acfBuilt = System.nanoTime() - t0;
        t0 = System.nanoTime();
        int dim0 = Math.max(sta.getArOrder(), sta.getMaOrder() + 1);
        double[] acfValues = ((jdplus.toolkit.base.core.arima.AutoCovarianceFunction) acf)
                .values(dim0);
        long acfVals = System.nanoTime() - t0;
        System.out.println("TIMING acf-build-ns=" + acfBuilt + " acf-values-ns=" + acfVals
                + " acf-dim=" + acfValues.length);

        // Per-evaluation lifecycle identical to the production chain: fresh
        // computer + fresh filter per evaluation (bytecode: CLC.compute
        // constructs KalmanFilter(nx>0) per call when its filter field is null).
        for (int e = 0; e < evals; ++e) {
            Probe pr = new Probe(nx > 0);
            ConcentratedLikelihoodComputer cll
                    = new ConcentratedLikelihoodComputer(pr, 1e-12, true, false);
            long tc = System.nanoTime();
            jdplus.toolkit.base.core.stats.likelihood.ConcentratedLikelihoodWithMissing ll
                    = cll.compute(am);
            long evalNanos = System.nanoTime() - tc;
            if (e == evals - 1) {
                Object inner = pr.real;
                int dim = (Integer) rf(inner, "dim");
                int steadyPos = (Integer) rf(inner, "steadyPos");
                double[] s2 = (double[]) rf(inner, "s");
                double[] C0 = (double[]) rf(inner, "C0");
                double[] C = (double[]) rf(inner, "C");
                double h0 = (Double) rf(inner, "h0");
                System.out.println("MACHINE dim=" + dim + " steadyPos=" + steadyPos
                        + " n=" + n + " h0=" + String.format(java.util.Locale.ROOT, "%.17g", h0)
                        + " C0=" + (C0 == null ? "<null>" : fmt(C0, 0, C0.length)));
                if (s2 != null) {
                    System.out.println("STRACE " + fmt(s2, 0, s2.length));
                }
                if (C != null) {
                    System.out.println("CROW0 " + fmt(C, 0, dim));
                    if (steadyPos < n) {
                        System.out.println("CROWFROZEN " + fmt(C, (steadyPos - 1) * dim,
                                (steadyPos - 1) * dim + dim));
                    }
                }
                System.out.println(String.format(java.util.Locale.ROOT,
                        "LL ssq=%.17g ldet=%.17g loglik=%.10f",
                        ll.ssq(), ll.logDeterminant(), ll.logLikelihood()));
            }
            System.out.println("TIMING eval=" + e
                    + " prepare-ns=" + pr.prepareNanos
                    + " apply-ns=" + pr.applyNanos + " apply-calls=" + pr.applyCalls
                    + " det-ns=" + pr.detNanos + " other-ns="
                    + (evalNanos - pr.prepareNanos - pr.applyNanos - pr.detNanos)
                    + " total-ns=" + evalNanos);
        }

        // Golden fidelity path through the production DEFAULT_FULL_COMPUTER.
        jdplus.toolkit.base.core.stats.likelihood.ConcentratedLikelihoodWithMissing lk
                = ConcentratedLikelihoodComputer.DEFAULT_FULL_COMPUTER.compute(am);
        System.out.println(String.format(java.util.Locale.ROOT,
                "FULLCOMPUTER ll=%.10f ssq=%.17g nobs=%d",
                lk.logLikelihood(), lk.ssq(), lk.e().length()));
        // Concentrated GLS coefficients at the pinned model (regression side).
        if (nx > 0) {
            jdplus.toolkit.base.core.math.matrices.FastMatrix xx = am.getX();
            StringBuilder sb = new StringBuilder("GLSCOEF ");
            for (int j = 0; j < nx; ++j) {
                sb.append(String.format(java.util.Locale.ROOT, "%.17g ",
                        lk.coefficient(j)));
            }
            System.out.println(sb.toString());
        }
    }
}
