// 23-F F.1 probe: the live jdplus.tramoseats.base.core.tramo.FinalEstimator
// driven through its own bytecode (split package; reflection ONLY for the
// private ladder methods and the significanceThreshold field injection, the
// DECISION bytecode always runs from the pinned jar).
package jdplus.tramoseats.base.core.tramo;

import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import jdplus.toolkit.base.api.arima.SarimaOrders;
import jdplus.toolkit.base.api.arima.SarimaSpec;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.Complex;
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.math.polynomials.Polynomial;
import jdplus.toolkit.base.core.regsarima.regular.ModelDescription;
import jdplus.toolkit.base.core.regsarima.regular.RegSarimaModelling;

public final class Scratch23f {

  Scratch23f() {
  }

  static String g(double v) {
    return String.format(Locale.ROOT, "%.17g", v);
  }

  static String bits(double v) {
    return String.format(Locale.ROOT, "%016x", Double.doubleToRawLongBits(v));
  }

  static String gb(double v) {
    return g(v) + "#" + bits(v);
  }

  static double tok(String s) {
    if (s.equals("nan")) {
      return Double.NaN;
    }
    if (s.equals("inf")) {
      return Double.POSITIVE_INFINITY;
    }
    if (s.equals("-inf")) {
      return Double.NEGATIVE_INFINITY;
    }
    if (s.startsWith("0x")) {
      return Double.longBitsToDouble(Long.parseUnsignedLong(s.substring(2), 16));
    }
    return Double.parseDouble(s);
  }

  static final class Observer implements ProcessingLog {
    final List<String> evts = new ArrayList<>();
    int depth = 0;
    long seq = 0;

    void emit(String kind, String msg, String obj, Object d) {
      ++seq;
      evts.add("evt " + seq + " " + kind + " " + depth
          + " msg " + (msg == null ? "-" : msg.replace('\n', ' '))
          + " obj " + (obj == null ? "-" : obj)
          + " det " + (d == null ? "-" : String.valueOf(d)));
    }

    public void push(String s) {
      emit("push", s, null, null);
      ++depth;
    }

    public void pop() {
      if (depth > 0) {
        --depth;
      }
      emit("pop", null, null, null);
    }

    public void error(String m, String o, Object d) {
      emit("error", m, o, d);
    }

    public void warning(String m, String o, Object d) {
      emit("warn", m, o, d);
    }

    public void remark(String m, String o, Object d) {
      emit("remark", m, o, d);
    }

    public void info(String m, String o, Object d) {
      emit("info", m, o, d);
    }

    public void step(String m, Object d) {
      emit("step", m, null, d);
    }

    public boolean isVerbose() {
      return true;
    }
  }

  static TsData monthly(double[] v) {
    return TsData.of(TsPeriod.monthly(1990, 1), DoubleSeq.of(v));
  }

  static TsDomain fullDomain(int n) {
    return TsDomain.of(TsPeriod.monthly(1990, 1), n);
  }

  // The SECTION-1 fixed input family (StrictMath only, intrinsic-free).
  static double[] family(String fam, int n) {
    double[] y = new double[n];
    for (int k = 0; k < n; ++k) {
      y[k] = 50 + 0.05 * (k + 1) + 5 * StrictMath.sin(2 * StrictMath.PI * (k + 1) / 12.0)
          + 3 * StrictMath.cos(4 * StrictMath.PI * (k + 1) / 12.0)
          + 0.5 * StrictMath.sin(14 * StrictMath.PI * (k + 1) / 12.0);
    }
    if (fam.equals("ao")) {
      y[30] += 30;
      y[70] -= 25;
    } else if (fam.equals("ls")) {
      for (int k = 48; k < n; ++k) {
        y[k] += 20;
      }
    } else if (fam.equals("out2")) {
      y[10] += 12;
      y[60] -= 18;
    }
    return y;
  }

  static void dumpArimaBits(String tag, ModelDescription desc) {
    jdplus.toolkit.base.core.sarima.SarimaModel ar = desc.arima();
    SarimaOrders o = ar.orders();
    DoubleSeq p = ar.parameters();
    StringBuilder s = new StringBuilder(tag + " orders " + o.getP() + " " + o.getD() + " "
        + o.getQ() + " " + o.getBp() + " " + o.getBd() + " " + o.getBq() + " pars");
    for (int i = 0; i < p.length(); ++i) {
      s.append(' ').append(gb(p.get(i)));
    }
    System.out.println(s);
  }

  // ---- private-ladder accessors (bytecode of the pinned jar) --------------
  static Method searchur, checkUnitRoots, checkCommonRoots, testM;

  static void need() throws Exception {
    if (searchur != null) {
      return;
    }
    searchur = FinalEstimator.class.getDeclaredMethod("searchur", Complex[].class);
    checkUnitRoots = FinalEstimator.class.getDeclaredMethod("checkUnitRoots",
        RegSarimaModelling.class);
    checkCommonRoots = FinalEstimator.class.getDeclaredMethod("checkCommonRoots",
        RegSarimaModelling.class);
    testM = FinalEstimator.class.getDeclaredMethod("test", RegSarimaModelling.class);
    searchur.setAccessible(true);
    checkUnitRoots.setAccessible(true);
    checkCommonRoots.setAccessible(true);
    testM.setAccessible(true);
  }

  static FinalEstimator fe(double eps, double ur, double cancel, int pass, boolean ami,
      boolean outliers) {
    return FinalEstimator.builder().precision(eps).unitRootThreshold(ur).cancel(cancel)
        .pass(pass).ami(ami).outliers(outliers).build();
  }

  // The EXACT live TramoKernel.estimateModel chain: cancel (0.044) and the
  // significanceThreshold (1.0) keep their builder defaults — the kernel
  // never touches them and no setter exists for tsig.
  static FinalEstimator fe_live(double eps, double ur, int pass, boolean ami,
      boolean outliers) {
    return FinalEstimator.builder().precision(eps).unitRootThreshold(ur)
        .pass(pass).ami(ami).outliers(outliers).build();
  }

  // fe_sr: direct private searchur over crafted Complex[] roots.
  static void fe_sr(String[] h, String[] toks) throws Exception {
    need();
    double ur = tok(h[2]);
    System.out.println("case " + h[1] + " searchur ur " + gb(ur));
    FinalEstimator est = fe(1e-7, ur, 0.0, 0, true, false);
    for (int c = 0; c < toks.length; c += 2) {
      double re = tok(toks[c]);
      double im = tok(toks[c + 1]);
      Complex[] roots = {Complex.cart(re, im)};
      int cnt = (Integer) searchur.invoke(est, (Object) roots);
      System.out.println("sr " + c / 2 + " re " + gb(re) + " im " + gb(im) + " count " + cnt);
    }
  }

  // fe_mr: live Polynomial.mirror().roots() bits over crafted phi blocks
  // (degree 1/2 closed forms and the degree >= 3 QR cascade).
  static void fe_mr(String[] h, String[] toks) {
    System.out.println("case " + h[1]);
    double[] c = new double[toks.length];
    for (int i = 0; i < toks.length; ++i) {
      c[i] = tok(toks[i]);
    }
    // phi block in JVM coordinates; regularAR = [1, phi...], mirrored.
    double[] poly = new double[c.length + 1];
    poly[0] = 1.0;
    System.arraycopy(c, 0, poly, 1, c.length);
    Polynomial p = Polynomial.of(poly);
    Complex[] roots = p.mirror().roots();
    StringBuilder s = new StringBuilder("mr roots " + roots.length);
    for (Complex z : roots) {
      s.append(' ').append(gb(z.getRe())).append(',').append(gb(z.getIm()))
          .append(" abs ").append(gb(z.abs()));
    }
    System.out.println(s);
  }

  static SarimaOrders mkOrders(int p, int d, int q, int bp, int bd, int bq) {
    SarimaOrders o = new SarimaOrders(12);
    o.setP(p);
    o.setD(d);
    o.setQ(q);
    o.setBp(bp);
    o.setBd(bd);
    o.setBq(bq);
    return o;
  }

  static jdplus.toolkit.base.api.data.Parameter[] est(double[] v) {
    jdplus.toolkit.base.api.data.Parameter[] a
        = new jdplus.toolkit.base.api.data.Parameter[v.length];
    for (int i = 0; i < v.length; ++i) {
      a[i] = jdplus.toolkit.base.api.data.Parameter.estimated(v[i]);
    }
    return a;
  }

  static ModelDescription crafted(int[] od, double[] phi, double[] bphi, double[] theta,
      double[] bth, double[] y) {
    TsData ts = monthly(y);
    ModelDescription desc = new ModelDescription(ts, fullDomain(y.length));
    desc.setArimaSpec(SarimaSpec.builder().p(od[0]).d(od[1]).q(od[2])
        .bp(od[3]).bd(od[4]).bq(od[5]).period(12)
        .phi(est(phi)).bphi(est(bphi)).theta(est(theta)).btheta(est(bth)).build());
    return desc;
  }

  // fe_ru: private checkUnitRoots over a crafted description.
  static void fe_ru(String[] h, String[] toks) throws Exception {
    need();
    System.out.println("case " + h[1]);
    int[] od = {Integer.parseInt(toks[0]), Integer.parseInt(toks[1]),
        Integer.parseInt(toks[2]), Integer.parseInt(toks[3]), Integer.parseInt(toks[4]),
        Integer.parseInt(toks[5])};
    double ur = tok(h[2]);
    int k = 6;
    int np = Integer.parseInt(toks[k++]);
    double[] phi = new double[np];
    for (int i = 0; i < np; ++i) {
      phi[i] = tok(toks[k++]);
    }
    int nbp = Integer.parseInt(toks[k++]);
    double[] bphi = new double[nbp];
    for (int i = 0; i < nbp; ++i) {
      bphi[i] = tok(toks[k++]);
    }
    int nq = Integer.parseInt(toks[k++]);
    double[] th = new double[nq];
    for (int i = 0; i < nq; ++i) {
      th[i] = tok(toks[k++]);
    }
    int nbq = Integer.parseInt(toks[k++]);
    double[] bth = new double[nbq];
    for (int i = 0; i < nbq; ++i) {
      bth[i] = tok(toks[k++]);
    }
    double[] y = family("none", 96);
    ModelDescription desc = crafted(od, phi, bphi, th, bth, y);
    dumpArimaBits("ru_in", desc);
    Observer log = new Observer();
    RegSarimaModelling m = RegSarimaModelling.of(desc, log);
    FinalEstimator est = fe(1e-7, ur, 0.0, 0, true, false);
    boolean unchanged;
    try {
      unchanged = (Boolean) checkUnitRoots.invoke(est, m);
    } catch (java.lang.reflect.InvocationTargetException e) {
      System.out.println("threw " + e.getCause().getClass().getName());
      for (String e2 : log.evts) {
        System.out.println(e2);
      }
      return;
    }
    System.out.println("ru_unchanged " + (unchanged ? 1 : 0));
    dumpArimaBits("ru_out", m.getDescription());
    for (String e : log.evts) {
      System.out.println(e);
    }
  }

  // fe_cur: private checkCommonRoots over a crafted description (cancel as
  // request data — the live TramoKernel site leaves it 0.0).
  static void fe_cur(String[] h, String[] toks) throws Exception {
    need();
    System.out.println("case " + h[1]);
    int[] od = {Integer.parseInt(toks[0]), Integer.parseInt(toks[1]),
        Integer.parseInt(toks[2]), Integer.parseInt(toks[3]), Integer.parseInt(toks[4]),
        Integer.parseInt(toks[5])};
    double cancel = tok(h[2]);
    int k = 6;
    int np = Integer.parseInt(toks[k++]);
    double[] phi = new double[np];
    for (int i = 0; i < np; ++i) {
      phi[i] = tok(toks[k++]);
    }
    int nbp = Integer.parseInt(toks[k++]);
    double[] bphi = new double[nbp];
    for (int i = 0; i < nbp; ++i) {
      bphi[i] = tok(toks[k++]);
    }
    int nq = Integer.parseInt(toks[k++]);
    double[] th = new double[nq];
    for (int i = 0; i < nq; ++i) {
      th[i] = tok(toks[k++]);
    }
    int nbq = Integer.parseInt(toks[k++]);
    double[] bth = new double[nbq];
    for (int i = 0; i < nbq; ++i) {
      bth[i] = tok(toks[k++]);
    }
    ModelDescription desc = crafted(od, phi, bphi, th, bth, family("none", 96));
    System.out.println("cur cancel " + gb(cancel) + " bphi1 "
        + (nbp > 0 ? gb(bphi[0]) : "-") + " bth1 " + (nbq > 0 ? gb(bth[0]) : "-"));
    Observer log = new Observer();
    RegSarimaModelling m = RegSarimaModelling.of(desc, log);
    FinalEstimator est = fe(1e-7, 0.9, cancel, 0, true, false);
    boolean ok;
    try {
      ok = (Boolean) checkCommonRoots.invoke(est, m);
    } catch (java.lang.reflect.InvocationTargetException e) {
      System.out.println("threw " + e.getCause().getClass().getName());
      return;
    }
    System.out.println("cur_ok " + (ok ? 1 : 0));
    dumpArimaBits("cur_out", m.getDescription());
  }

  // fe_tc: private test(m) after a live estimation, tsig injected through the
  // final field (the Builder exposes NO setter: tsig is the 1.0 builder
  // default in every live build — the injection parameterizes the port).
  static void fe_tc(String[] h, String[] toks) throws Exception {
    need();
    System.out.println("case " + h[1]);
    String fam = h[2];
    double tsig = tok(h[3]);
    int p = Integer.parseInt(h[4]);
    int q = Integer.parseInt(h[5]);
    int bp = Integer.parseInt(h[6]);
    int bq = Integer.parseInt(h[7]);
    double[] y = family(fam, 96);
    ModelDescription desc = new ModelDescription(monthly(y), fullDomain(y.length));
    desc.setMean(true);
    desc.setSpecification(mkOrders(p, 1, q, bp, 1, bq));
    Observer log = new Observer();
    RegSarimaModelling m = RegSarimaModelling.of(desc, log);
    m.estimate(1e-5);
    dumpArimaBits("tc_landing", m.getDescription());
    var cl = m.getEstimation().getConcentratedLikelihood();
    System.out.println("tc_ssq " + gb(cl.ssq()) + " tc_dim " + cl.dim());
    double[] va = m.getEstimation().getMax().asymptoticCovariance().diagonal().toArray();
    System.out.println("tc_vdiag n " + va.length);
    for (int i = 0; i < va.length; ++i) {
      System.out.println("tc_vd " + i + " " + gb(va[i]));
    }
    Field f = FinalEstimator.class.getDeclaredField("tsig");
    f.setAccessible(true);
    FinalEstimator est = fe(1e-7, 0.9, 0.0, 0, true, false);
    f.set(est, tsig);
    int total;
    try {
      total = (Integer) testM.invoke(est, m);
    } catch (java.lang.reflect.InvocationTargetException e) {
      System.out.println("threw " + e.getCause().getClass().getName());
      return;
    }
    System.out.println("tc_total " + total);
    dumpArimaBits("tc_out", m.getDescription());
    for (String e : log.evts) {
      System.out.println(e);
    }
  }

  // fe_est: the LIVE estimate() driver on a real estimation (return, nnsig,
  // orders, landing bits, events).
  static void fe_est(String[] h, String[] toks) {
    System.out.println("case " + h[1]);
    String fam = h[2];
    int n = Integer.parseInt(h[3]);
    boolean mean = h[4].equals("y");
    double eps = tok(h[5]);
    double ur = tok(h[6]);
    int pass = Integer.parseInt(h[7]);
    boolean ami = h[8].equals("y");
    boolean outliers = h[9].equals("y");
    String[] ov = h[10].split(",");
    int p = 0;
    int d = 0;
    int q = 0;
    int bp = 0;
    int bd = 0;
    int bq = 0;
    if (!h[10].equals("air")) {
      p = Integer.parseInt(ov[0]);
      d = Integer.parseInt(ov[1]);
      q = Integer.parseInt(ov[2]);
      bp = Integer.parseInt(ov[3]);
      bd = Integer.parseInt(ov[4]);
      bq = Integer.parseInt(ov[5]);
    }
    System.out.println("in fam " + fam + " n " + n + " mean " + (mean ? 1 : 0)
        + " eps " + gb(eps) + " ur " + gb(ur) + " pass " + pass
        + " ami " + (ami ? 1 : 0) + " outliers " + (outliers ? 1 : 0));
    double[] y = family(fam, n);
    ModelDescription desc = new ModelDescription(monthly(y), fullDomain(n));
    desc.setMean(mean);
    if (h[10].equals("air")) {
      desc.setAirline(true);
    } else {
      desc.setSpecification(mkOrders(p, d, q, bp, bd, bq));
    }
    Observer log = new Observer();
    RegSarimaModelling m = RegSarimaModelling.of(desc, log);
    FinalEstimator est = fe_live(eps, ur, pass, ami, outliers);
    boolean ret;
    try {
      ret = est.estimate(m);
    } catch (RuntimeException err) {
      ret = false;
      System.out.println("threw " + err.getClass().getName());
    }
    System.out.println("est_return " + (ret ? 1 : 0) + " nnsig "
        + est.getChangedParametersCount());
    System.out.println("mean_out " + (m.getDescription().isMean() ? 1 : 0));
    dumpArimaBits("est_out", m.getDescription());
    try {
      var cl = m.getEstimation().getConcentratedLikelihood();
      double[] co = cl.allCoefficients().toArray();
      System.out.println("cl_dim " + cl.dim() + " ssq " + gb(cl.ssq()));
      System.out.println("coeff n " + co.length);
      for (double v : co) {
        System.out.println("coef " + gb(v));
      }
      double[] va = m.getEstimation().getMax().asymptoticCovariance().diagonal().toArray();
      System.out.println("vdiag n " + va.length);
      for (int i = 0; i < va.length; ++i) {
        System.out.println("vd " + i + " " + gb(va[i]));
      }
    } catch (RuntimeException err) {
      System.out.println("cl_unavailable " + err.getClass().getName());
    }
    for (String e : log.evts) {
      System.out.println(e);
    }
  }

  public static void main(String[] args) throws IOException {
    // Case lines arrive on stdin: TYPE id rest... tokens space-separated.
    try (java.io.BufferedReader in = new java.io.BufferedReader(
        new java.io.FileReader(args[0]))) {
        String line;
        while ((line = in.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] toks = line.trim().split("\\s+");
            String[] h = Arrays.copyOf(toks, Math.min(toks.length, 11));
            try {
            if (toks[0].equals("fe_sr")) {
                fe_sr(h, Arrays.copyOfRange(toks, 3, toks.length));
            } else if (toks[0].equals("fe_mr")) {
                fe_mr(h, Arrays.copyOfRange(toks, 3, toks.length));
            } else if (toks[0].equals("fe_ru")) {
                fe_ru(h, Arrays.copyOfRange(toks, 3, toks.length));
            } else if (toks[0].equals("fe_cur")) {
                fe_cur(h, Arrays.copyOfRange(toks, 3, toks.length));
            } else if (toks[0].equals("fe_tc")) {
                fe_tc(h, toks);
            } else if (toks[0].equals("fe_est")) {
                fe_est(h, toks);
            } else {
                throw new IllegalStateException("unknown case " + toks[0]);
            }
            } catch (Throwable err) {
                Throwable c = err instanceof java.lang.reflect.InvocationTargetException
                    && err.getCause() != null ? err.getCause() : err;
                System.out.println("case_threw " + h[1] + " " + c.getClass().getName()
                    + " msg " + String.valueOf(c.getMessage()).replace('\n', ' '));
            }
            System.out.println("end " + h[1]);
        }
    }
  }
}
