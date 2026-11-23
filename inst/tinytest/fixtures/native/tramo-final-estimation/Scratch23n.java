// 23-F F.2 round-step probe: for every fe_est case, (1) a REAL run of
// jdplus.tramoseats.base.core.tramo.FinalEstimator.estimate reproducing the
// fe_est capture grammar (gate: byte-identical to the F.1 capture blocks),
// and (2) a MANUAL replay of the decoded driver (SECTION 29) over the same
// public machinery — per-round computer/process/landing/Point/ladder rows
// that pin the C++ round-step contract. Reflection ONLY for the three
// private ladder methods (access only; the decision bytecode is the jar's).
package jdplus.tramoseats.base.core.tramo;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import jdplus.toolkit.base.api.arima.SarimaOrders;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.arima.estimation.IArimaMapping;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.regsarima.RegSarimaComputer;
import jdplus.toolkit.base.core.regsarima.regular.ModelDescription;
import jdplus.toolkit.base.core.regsarima.regular.RegSarimaModelling;
import jdplus.toolkit.base.core.regarima.RegArimaEstimation;
import jdplus.toolkit.base.core.sarima.SarimaModel;
import jdplus.toolkit.base.core.sarima.estimation.SarimaMapping2;
import jdplus.toolkit.base.core.stats.likelihood.ConcentratedLikelihoodWithMissing;
import jdplus.toolkit.base.core.stats.likelihood.LogLikelihoodFunction;

public final class Scratch23n {

  Scratch23n() {
  }

  static String g(double v) {
    return String.format(Locale.ROOT, "%.17g", v);
  }

  static String gb(double v) {
    return g(v) + "#" + String.format(Locale.ROOT, "%016x", Double.doubleToRawLongBits(v));
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

  static double tok(String s) {
    if (s.equals("nan")) {
      return Double.NaN;
    }
    if (s.startsWith("0x")) {
      return Double.longBitsToDouble(Long.parseUnsignedLong(s.substring(2), 16));
    }
    return Double.parseDouble(s);
  }

  static FinalEstimator fe_live(double eps, double ur, int pass, boolean ami, boolean outliers) {
    return FinalEstimator.builder().precision(eps).unitRootThreshold(ur)
        .pass(pass).ami(ami).outliers(outliers).build();
  }

  static Method ladder(String name) throws Exception {
    Method mm = FinalEstimator.class.getDeclaredMethod(name, RegSarimaModelling.class);
    mm.setAccessible(true);
    return mm;
  }

  static void dumpArimaBits(String tag, ModelDescription desc) {
    SarimaModel ar = desc.arima();
    SarimaOrders o = ar.orders();
    DoubleSeq p = ar.parameters();
    StringBuilder s = new StringBuilder(tag + " orders " + o.getP() + " " + o.getD() + " "
        + o.getQ() + " " + o.getBp() + " " + o.getBd() + " " + o.getBq() + " pars");
    for (int i = 0; i < p.length(); ++i) {
      s.append(' ').append(gb(p.get(i)));
    }
    System.out.println(s);
  }

  static void dumpCl(String tag, RegArimaEstimation<SarimaModel> est) {
    try {
      ConcentratedLikelihoodWithMissing cl = est.getConcentratedLikelihood();
      System.out.println(tag + " cl_dim " + cl.dim() + " ssq " + gb(cl.ssq()) + " ldet "
          + gb(cl.logDeterminant()));
      double[] co = cl.allCoefficients().toArray();
      System.out.println(tag + " coeff n " + co.length);
      for (double v : co) {
        System.out.println(tag + " coef " + gb(v));
      }
      try {
        double[] va = est.getMax().asymptoticCovariance().diagonal().toArray();
        System.out.println(tag + " vdiag n " + va.length);
        for (int i = 0; i < va.length; ++i) {
          System.out.println(tag + " vd " + i + " " + gb(va[i]));
        }
      } catch (NullPointerException err) {
        System.out.println(tag + " vdiag null");
      }
      LogLikelihoodFunction.Point<?, ?> max = est.getMax();
      DoubleSeq mp = max.getParameters();
      System.out.println(tag + " max n " + mp.length());
      for (int i = 0; i < mp.length(); ++i) {
        System.out.println(tag + " maxp " + i + " " + gb(mp.get(i)));
      }
      DoubleSeq sc = max.getScore();
      System.out.println(tag + " score n " + sc.length());
      for (int i = 0; i < sc.length(); ++i) {
        System.out.println(tag + " scor " + i + " " + gb(sc.get(i)));
      }
      FastMatrix inf = max.getInformation();
      System.out.println(tag + " info " + inf.getRowsCount() + "x" + inf.getColumnsCount());
      for (int i = 0; i < inf.getRowsCount(); ++i) {
        StringBuilder s = new StringBuilder(tag + " inf ");
        for (int j = 0; j < inf.getColumnsCount(); ++j) {
          s.append(' ').append(gb(inf.get(i, j)));
        }
        System.out.println(s);
      }
    } catch (NullPointerException err) {
      System.out.println(tag + " cl_unavailable " + err.getClass().getName());
    }
  }

  static final class Cfg {
    String fam;
    int n;
    boolean mean;
    double eps;
    double ur;
    int pass;
    boolean ami;
    boolean outliers;
    String ordertok;
  }

  static ModelDescription buildDesc(Cfg c) {
    ModelDescription desc = new ModelDescription(monthly(family(c.fam, c.n)), fullDomain(c.n));
    desc.setMean(c.mean);
    if (c.ordertok.equals("air")) {
      desc.setAirline(true);
    } else {
      String[] ov = c.ordertok.split(",");
      desc.setSpecification(mkOrders(Integer.parseInt(ov[0]), Integer.parseInt(ov[1]),
          Integer.parseInt(ov[2]), Integer.parseInt(ov[3]), Integer.parseInt(ov[4]),
          Integer.parseInt(ov[5])));
    }
    return desc;
  }

  static void realRun(Cfg c) throws Exception {
    ModelDescription desc = buildDesc(c);
    Observer log = new Observer();
    RegSarimaModelling m = RegSarimaModelling.of(desc, log);
    FinalEstimator est = fe_live(c.eps, c.ur, c.pass, c.ami, c.outliers);
    boolean ret;
    try {
      ret = est.estimate(m);
    } catch (RuntimeException err) {
      ret = false;
      System.out.println("r threw " + err.getClass().getName());
    }
    System.out.println("r est_return " + (ret ? 1 : 0) + " nnsig "
        + est.getChangedParametersCount());
    System.out.println("r mean_out " + (m.getDescription().isMean() ? 1 : 0));
    dumpArimaBits("r est_out", m.getDescription());
    RegArimaEstimation<SarimaModel> e = m.getEstimation();
    if (e == null) {
      System.out.println("r cl_unavailable estimation-null");
    } else {
      dumpCl("r", e);
    }
    for (String ev : log.evts) {
      System.out.println("r " + ev);
    }
  }

  // Manual replay of the decoded driver (SECTION 29) with per-round dumps.
  static void manualRun(Cfg c) throws Exception {
    Method mUr = ladder("checkUnitRoots");
    Method mCr = ladder("checkCommonRoots");
    Method mTst = ladder("test");
    ModelDescription desc = buildDesc(c);
    Observer log = new Observer();
    RegSarimaModelling m = RegSarimaModelling.of(desc, log);
    FinalEstimator fe = fe_live(c.eps, c.ur, c.pass, c.ami, c.outliers);
    System.out.println("m air_spec " + (desc.specification().toString()));
    int nnsig = 0;
    int ret = -1;
    int cont = 0;
    for (int round = 0; round < 12 && ret < 0; ++round) {
      ModelDescription d = m.getDescription();
      dumpArimaBits("m R" + round + " in", d);
      IArimaMapping<SarimaModel> mapping = d.mapping();
      int dim = mapping.getDim();
      System.out.println("m R" + round + " dim " + dim + " pcount "
          + d.specification().getParametersCount() + " mapclass "
          + mapping.getClass().getSimpleName());
      if (dim != 0) {
        // pin the HR start the computer's process() would build
        SarimaModel hr = jdplus.toolkit.base.core.regsarima.internal
            .HannanRissanenInitializer.builder().stabilize(true).useDefaultIfFailed(true)
            .build().initialize(d.regarima().differencedModel());
        StringBuilder hs = new StringBuilder("m R" + round + " hr");
        for (int i = 0; i < hr.parameters().length(); ++i) {
          hs.append(' ').append(gb(hr.parameters().get(i)));
        }
        System.out.println(hs);
        IArimaMapping<SarimaModel> smap = mapping.stationaryMapping();
        DoubleSeq warm = smap.parametersOf(hr);
        StringBuilder ws = new StringBuilder("m R" + round + " warm");
        for (int i = 0; i < warm.length(); ++i) {
          ws.append(' ').append(gb(warm.get(i)));
        }
        System.out.println(ws);
        System.out.println("m R" + round + " validate "
            + smap.validate(jdplus.toolkit.base.core.data.DataBlock.of(warm)).name());
      }
      RegSarimaComputer computer = RegSarimaComputer.builder()
          .minimizer(LevenbergMarquardtMinimizer.builder()).precision(c.eps)
          .computeExactFinalDerivatives(true).build();
      m.estimate(computer);
      RegArimaEstimation<SarimaModel> e = m.getEstimation();
      dumpArimaBits("m R" + round + " land", m.getDescription());
      dumpCl("m R" + round, e);
      if (dim == 0) {
        ret = 1; // dim==0 short-circuit; nnsig untouched
        System.out.println("m R" + round + " shortcircuit");
        break;
      }
      if (!c.ami) {
        ret = 1;
        System.out.println("m R" + round + " notami");
        break;
      }
      boolean cur = (Boolean) mUr.invoke(fe, m);
      System.out.println("m R" + round + " checkur " + (cur ? 1 : 0));
      dumpArimaBits("m R" + round + " afterur", m.getDescription());
      if (!cur) {
        if (++cont >= 6) {
          ret = 0;
          break;
        }
        continue;
      }
      boolean cc = (Boolean) mCr.invoke(fe, m);
      System.out.println("m R" + round + " checkcr " + (cc ? 1 : 0));
      dumpArimaBits("m R" + round + " aftercr", m.getDescription());
      if (cc) {
        int tot = (Integer) mTst.invoke(fe, m);
        nnsig = tot;
        System.out.println("m R" + round + " testtotal " + tot);
        dumpArimaBits("m R" + round + " aftertest", m.getDescription());
      } else {
        nnsig = 2;
      }
      if (nnsig == 0) {
        ret = 1;
        break;
      }
      if (nnsig == 1) {
        if (++cont >= 6) {
          ret = 0;
          break;
        }
        continue;
      }
      if (c.outliers && c.pass <= 1) {
        ret = 0;
        System.out.println("m R" + round + " cancelreturn");
        break;
      }
      if (++cont >= 6) {
        ret = 0;
        break;
      }
    }
    System.out.println("m est_return " + ret + " nnsig " + nnsig);
    System.out.println("m mean_out " + (m.getDescription().isMean() ? 1 : 0));
    dumpArimaBits("m est_out", m.getDescription());
    RegArimaEstimation<SarimaModel> e = m.getEstimation();
    if (e == null) {
      System.out.println("m cl_unavailable estimation-null");
    } else {
      dumpCl("m", e);
    }
    for (String ev : log.evts) {
      System.out.println("m " + ev);
    }
  }

  public static void main(String[] args) throws Exception {
    try (java.io.BufferedReader in = new java.io.BufferedReader(
        new java.io.FileReader(args[0]))) {
      String line;
      while ((line = in.readLine()) != null) {
        if (line.isEmpty() || line.startsWith("#")) {
          continue;
        }
        String[] toks = line.trim().split("\\s+");
        if (!toks[0].equals("fe_est")) {
          continue;
        }
        Cfg c = new Cfg();
        c.fam = toks[2];
        c.n = Integer.parseInt(toks[3]);
        c.mean = toks[4].equals("y");
        c.eps = tok(toks[5]);
        c.ur = tok(toks[6]);
        c.pass = Integer.parseInt(toks[7]);
        c.ami = toks[8].equals("y");
        c.outliers = toks[9].equals("y");
        c.ordertok = toks[10];
        System.out.println("case " + toks[1]);
        System.out.println("in fam " + c.fam + " n " + c.n + " mean " + (c.mean ? 1 : 0)
            + " eps " + gb(c.eps) + " ur " + gb(c.ur) + " pass " + c.pass
            + " ami " + (c.ami ? 1 : 0) + " outliers " + (c.outliers ? 1 : 0));
        realRun(c);
        manualRun(c);
        System.out.println("end " + toks[1]);
      }
    }
  }
}
