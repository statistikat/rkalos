// 23-F.3 probe: the TramoKernel round driver of the pinned jdplus 3.9.0
// stack, observed two ways per case from ONE JVM:
//  (A) the REAL run: TramoKernel.of(spec, ctx).process(data, log) - every
//      ProcessingLog event (seq/kind/depth/message) plus a kernel-field
//      snapshot at every "round N" step event (reflection at the pass
//      boundary: pass, round, pass3, needOutliers, needAutoModelling,
//      curva, pcr, refAuto/refStats set-ness, ctx seasonal/
//      originalSeasonalityTest);
//  (B) the MANUAL run: the DECODED driver line by line (process/ami stage
//      order, the 10-iteration iterate() loop with the bytecode's stage
//      order and gating, estimateModel with the FinalEstimator builder
//      chain and the pass-0 retry ladder, the pass2 refStats gate,
//      lastSolution, reduceVa, testRegression, restore, the round-0 accept
//      rule), writing the kernel fields through reflection EXACTLY where
//      the bytecode's putfields write, and calling the PINNED jar machines
//      as primitives (DifferencingModule/ArmaModule/OutliersDetectionModule
//      via the kernel's own private factories, the pre-stage
//      testSeasonality/testTransformation, tramo FinalEstimator.estimate,
//      ModelStatistics.of, FastRegressionTest.test, SeasonalityController/
//      ModelController.process, m.estimate). Because the manual driver
//      mirrors the putfields, the SAME observer hook snapshots it at the
//      same points. If the decoded loop shape or any state transition were
//      wrong, the manual event trace / pass snapshots / final model bits
//      would diverge from the real run's - each case ends with
//      dump-not-trust identity lines (events, snaps, final bits).
//
// Case grammar: "KER <id> <spec> <freq> <startY> <startM>" header plus
// "values v1 v2 ..." lines (%.17g), mirroring the 23-B/23-F.2 probe style.
// <spec> is a TramoSpec alias (TR0..TR5, TRfull) with optional overrides
// "+oc" (transform outliersCorrection), "+am"/"-am" (usingAutoModel),
// "+accept" (autoModel acceptDefault), "+mean" (regression mean test),
// "-o" (outlier types off = isUsed false), and kv pairs "va=" "pcr="
// "tsig=" "pc=" "fct=".
package jdplus.tramoseats.base.core.tramo;

import jdplus.toolkit.base.api.arima.SarimaOrders;
import jdplus.tramoseats.base.api.tramo.*;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.modelling.TransformationType;
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.regression.ModellingContext;
import jdplus.toolkit.base.api.timeseries.regression.ModellingUtility;
import jdplus.toolkit.base.api.timeseries.regression.Variable;
import jdplus.toolkit.base.core.regsarima.regular.ModelDescription;
import jdplus.toolkit.base.core.regsarima.regular.ProcessingResult;
import jdplus.toolkit.base.core.regsarima.regular.RegSarimaModel;
import jdplus.toolkit.base.core.regsarima.regular.RegSarimaModelling;
import jdplus.toolkit.base.core.stats.likelihood.ConcentratedLikelihood;
import jdplus.toolkit.base.core.stats.likelihood.ConcentratedLikelihoodWithMissing;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

public class Scratch23r {

  // ---------------- JVM comparison primitives ----------------
  static int dcmpl(double x, double y) {
    if (Double.isNaN(x) || Double.isNaN(y)) {
      return 1;
    }
    return Double.compare(x, y);
  }

  static int dcmpg(double x, double y) {
    if (Double.isNaN(x) || Double.isNaN(y)) {
      return -1;
    }
    return Double.compare(x, y);
  }

  static String b(double v) {
    return String.format("%016x", Double.doubleToRawLongBits(v));
  }

  static String bs(double[] v) {
    StringBuilder s = new StringBuilder();
    for (int i = 0; i < v.length; ++i) {
      if (i > 0) {
        s.append(' ');
      }
      s.append(b(v[i]));
    }
    return s.toString();
  }

  static String g17(double v) {
    return String.format("%.17g", v);
  }

  static Object rf(Object target, String name) {
    try {
      Class<?> c = target.getClass();
      while (c != null) {
        try {
          Field f = c.getDeclaredField(name);
          f.setAccessible(true);
          return f.get(target);
        } catch (NoSuchFieldException e) {
          c = c.getSuperclass();
        }
      }
      throw new NoSuchFieldException(name);
    } catch (ReflectiveOperationException e) {
      throw new RuntimeException("field " + name, e);
    }
  }

  static void wf(Object target, String name, Object value) {
    try {
      Field f = TramoKernel.class.getDeclaredField(name);
      f.setAccessible(true);
      f.set(target, value);
    } catch (ReflectiveOperationException e) {
      throw new RuntimeException("field " + name, e);
    }
  }

  static Object call(Object target, String name, Class<?>[] sig, Object... args) {
    try {
      java.util.Deque<Class<?>> q = new java.util.ArrayDeque<>();
      java.util.Set<Class<?>> seen = new java.util.HashSet<>();
      q.add(target.getClass());
      while (!q.isEmpty()) {
        Class<?> c = q.poll();
        if (c == null || !seen.add(c)) {
          continue;
        }
        try {
          Method m = c.getDeclaredMethod(name, sig);
          m.setAccessible(true);
          return m.invoke(target, args);
        } catch (NoSuchMethodException e) {
          q.add(c.getSuperclass());
          for (Class<?> i : c.getInterfaces()) {
            q.add(i);
          }
        }
      }
      throw new NoSuchMethodException(name);
    } catch (ReflectiveOperationException e) {
      throw new RuntimeException("method " + name, e);
    }
  }

  // ---------------- observer with pass-boundary snapshots ----------------
  static final class Observer implements ProcessingLog {
    final List<String> evts = new ArrayList<>();
    final String snapTag;
    TramoKernel kernel; // snapshots taken here at every "round N" step
    int depth = 0;
    long seq = 0;

    Observer(String snapTag) {
      this.snapTag = snapTag;
    }

    void emit(String kind, String msg) {
      ++seq;
      evts.add("evt " + seq + " " + kind + " " + depth + " "
          + (msg == null ? "-" : msg.replace('\n', ' ')));
      if (msg != null && msg.startsWith("round ") && kernel != null) {
        evts.add(snapshotLine(snapTag, kernel));
      }
    }

    public void push(String s) { emit("push", s); ++depth; }
    public void pop() { if (depth > 0) { --depth; } emit("pop", null); }
    public void error(String m, String o, Object d) { emit("error", m); }
    public void warning(String m, String o, Object d) { emit("warn", m); }
    public void remark(String m, String o, Object d) { emit("remark", m); }
    public void info(String m, String o, Object d) { emit("info", m); }
    public void step(String m, Object d) { emit("step", m); }
    public void remark(String m) { emit("remark", m); }
    public void info(String m) { emit("info", m); }
    public void warning(String m) { emit("warn", m); }
    public void error(String m) { emit("error", m); }
    public void step(String m) { emit("step", m); }
    public boolean isVerbose() { return true; }
  }

  static String snapshotLine(String tag, TramoKernel k) {
    StringBuilder s = new StringBuilder(tag);
    try {
      s.append(" pass=").append(rf(k, "pass"))
          .append(" round=").append(rf(k, "round"))
          .append(" pass3=").append((boolean) rf(k, "pass3") ? 1 : 0)
          .append(" needOut=").append((boolean) rf(k, "needOutliers") ? 1 : 0)
          .append(" needAuto=").append((boolean) rf(k, "needAutoModelling") ? 1 : 0)
          .append(" curva=").append(g17((double) rf(k, "curva")))
          .append(" pcr=").append(g17((double) rf(k, "pcr")))
          .append(" refAuto=").append(rf(k, "refAuto") == null ? 0 : 1)
          .append(" refStats=").append(rf(k, "refStats") == null ? 0 : 1);
      Object ctx = rf(k, "context");
      s.append(" seasonal=").append((boolean) rf(ctx, "seasonal") ? 1 : 0)
          .append(" ost=").append(rf(ctx, "originalSeasonalityTest"));
    } catch (Exception e) {
      s.append(" ERR ").append(e);
    }
    return s.toString();
  }

  static String ordersLine(SarimaOrders o) {
    return o.getP() + "," + o.getD() + "," + o.getQ() + ","
        + o.getBp() + "," + o.getBd() + "," + o.getBq() + "@" + o.getPeriod();
  }

  static String outlierNames(Variable[] vars) {
    StringBuilder ov = new StringBuilder();
    for (Variable v : vars) {
      if (ModellingUtility.isOutlier(v, true)) {
        if (ov.length() > 0) {
          ov.append(';');
        }
        ov.append(v.getName());
      }
    }
    return ov.length() == 0 ? "-" : ov.toString();
  }

  static String outlierNames(ModelDescription desc) {
    return outlierNames(desc.variables().toArray(Variable[]::new));
  }

  // arima coefficients of the current description as raw bit words
  static String parsBits(ModelDescription desc) {
    return bs(desc.arima().parameters().toArray());
  }

  // ---------------- model facts ----------------
  static void dumpModel(List<String> out, String tag, RegSarimaModel m) {
    if (m == null) {
      out.add(tag + " null");
      return;
    }
    SarimaOrders o = m.specification();
    jdplus.toolkit.base.core.modelling.GeneralLinearModel.Description<?> desc =
        m.getDescription();
    out.add(tag + " ok spec " + ordersLine(o) + " freq " + m.getAnnualFrequency());
    out.add(tag + " log " + (desc.isLogTransformation() ? 1 : 0)
        + " nvars " + desc.getVariables().length);
    out.add(tag + " outliers " + outlierNames(desc.getVariables()));
    out.add(tag + " pars " + bs(m.arima().parameters().toArray()));
    Object est = m.getEstimation();
    if (est == null) {
      out.add(tag + " est -");
      return;
    }
    jdplus.toolkit.base.core.modelling.GeneralLinearModel.Estimation ge
        = (jdplus.toolkit.base.core.modelling.GeneralLinearModel.Estimation) est;
    out.add(tag + " coefs " + bs(ge.getCoefficients().toArray()));
    jdplus.toolkit.base.api.math.matrices.Matrix cov = ge.getCoefficientsCovariance();
    int n = (int) call(cov, "getRowsCount", new Class<?>[]{});
    StringBuilder cb = new StringBuilder();
    for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
        cb.append(b(cov.get(i, j))).append(' ');
      }
    }
    out.add(tag + " cov " + cb.toString().trim());
    out.add(tag + " params " + bs(ge.getParameters().getValues().toArray()));
    jdplus.toolkit.base.core.stats.likelihood.LikelihoodStatistics ls = ge.getStatistics();
    out.add(tag + " stats obs " + ls.getObservationsCount()
        + " eff " + ls.getEffectiveObservationsCount()
        + " npar " + ls.getEstimatedParametersCount()
        + " ll " + g17(ls.getLogLikelihood())
        + " tadj " + g17(ls.getTransformationAdjustment())
        + " all " + g17(ls.getAdjustedLogLikelihood())
        + " ssq " + g17(ls.getSsqErr())
        + " aic " + g17(ls.getAIC())
        + " aicc " + g17(ls.getAICC())
        + " bic " + g17(ls.getBIC())
        + " bicc " + g17(ls.getBICC())
        + " bic2 " + g17(ls.getBIC2())
        + " hq " + g17(ls.getHannanQuinn()));
    TsData fres = m.fullResiduals();
    out.add(tag + " fres " + fres.getStart().toString() + " "
        + bs(fres.getValues().toArray()));
  }

  static void dumpStats(List<String> out, String site, ModelDescription desc,
      ConcentratedLikelihood cl) {
    ModelStatistics ms = ModelStatistics.of(desc, cl);
    double[] e = new double[cl.e().length()];
    for (int i = 0; i < e.length; ++i) {
      e[i] = cl.e().get(i);
    }
    out.add("stats " + site
        + " freq " + desc.getAnnualFrequency()
        + " k " + desc.getArimaSpec().freeParametersCount()
        + " dim " + cl.dim()
        + " obs " + desc.getEstimationDomain().getLength()
        + " lbp " + g17(ms.getLjungBoxPvalue())
        + " se " + g17(ms.getSe())
        + " lb " + g17(ms.getLjungBox())
        + " bic " + g17(ms.getBic())
        + " oc " + ms.getOutliersCount()
        + " eff " + ms.getEffectiveObservationsCount());
    out.add("stats " + site + " e " + bs(e));
  }

  // ---------------- spec resolution ----------------
  // alias: TR0..TR5/TRfull (TramoSpec.fromString) or "K" (= TR3); flags
  // +oc +accept +mean +am -am -o +noao +nols +notc +noso +ml +td +fs
  // (fs = transform None); kv: va=<double>|calc pcr= tsig= pc= fct= tcrate=
  static TramoSpec resolveSpec(String tok) {
    // split trailing space-separated kv pairs off the alias token
    String[] parts = tok.trim().split("\\s+");
    String base = parts[0];
    Double va = null, pcr = null, tsig = null, pc = null, fct = null, tcrate = null;
    for (int i = 1; i < parts.length; ++i) {
      String kv = parts[i];
      if (kv.startsWith("va=") && !kv.equals("va=calc")) {
        va = Double.parseDouble(kv.substring(3));
      }
      if (kv.startsWith("pcr=")) { pcr = Double.parseDouble(kv.substring(4)); }
      if (kv.startsWith("tsig=")) { tsig = Double.parseDouble(kv.substring(5)); }
      if (kv.startsWith("pc=")) { pc = Double.parseDouble(kv.substring(3)); }
      if (kv.startsWith("fct=")) { fct = Double.parseDouble(kv.substring(4)); }
      if (kv.startsWith("tcrate=")) { tcrate = Double.parseDouble(kv.substring(7)); }
    }
    boolean oc = base.contains("+oc");
    boolean accept = base.contains("+accept");
    boolean mean = base.contains("+mean");
    boolean amSet = false, am = false;
    if (base.contains("+am")) { amSet = true; am = true; }
    if (base.contains("-am")) { amSet = true; am = false; }
    boolean oOff = base.contains("-o");
    boolean noao = base.contains("+noao"), nols = base.contains("+nols");
    boolean notc = base.contains("+notc"), noso = base.contains("+noso");
    boolean ml = base.contains("+ml"), td = base.contains("+td"), fs = base.contains("+fs");
    base = base.replace("+oc", "").replace("+accept", "").replace("+mean", "")
        .replace("+am", "").replace("-am", "").replace("-o", "")
        .replace("+noao", "").replace("+nols", "").replace("+notc", "")
        .replace("+noso", "").replace("+ml", "").replace("+td", "").replace("+fs", "");
    if (base.equals("K")) {
      base = "TR3";
    }
    TramoSpec spec = TramoSpec.fromString(base);
    jdplus.tramoseats.base.api.tramo.TramoSpec.Builder b = spec.toBuilder();
    if (oc || fct != null || fs) {
      jdplus.tramoseats.base.api.tramo.TransformSpec.Builder tb
          = spec.getTransform().toBuilder();
      if (oc) {
        tb.outliersCorrection(true);
      }
      if (fs) {
        tb.function(TransformationType.None);
      }
      if (fct != null) {
        tb.function(TransformationType.Auto).fct(fct);
      }
      b.transform(tb.build());
    }
    if (amSet) {
      b.usingAutoModel(am);
    }
    if (oOff || va != null || noao || nols || notc || noso || ml || tcrate != null) {
      jdplus.tramoseats.base.api.tramo.OutlierSpec.Builder ob
          = spec.getOutliers().toBuilder();
      if (oOff) {
        ob.ao(false).ls(false).tc(false).so(false);
      }
      if (noao) {
        ob.ao(false);
      }
      if (nols) {
        ob.ls(false);
      }
      if (notc) {
        ob.tc(false);
      }
      if (noso) {
        ob.so(false);
      }
      if (ml) {
        ob.maximumLikelihood(true);
      }
      if (tcrate != null) {
        ob.deltaTC(tcrate);
      }
      if (va != null) {
        ob.criticalValue(va);
      }
      b.outliers(ob.build());
    }
    if (td) {
      b.regression(spec.getRegression().toBuilder().calendar(
          spec.getRegression().getCalendar().toBuilder().tradingDays(
              TradingDaysSpec.automatic(
                  jdplus.toolkit.base.api.timeseries.calendars.LengthOfPeriodType.LengthOfPeriod,
                  TradingDaysSpec.AutoMethod.FTEST, 0.7, true)).build()).build());
    }
    if (mean) {
      b.regression(spec.getRegression().toBuilder().mean(
          spec.getRegression().getMean().toBuilder().test(true).build()).build());
    }
    if (accept || pcr != null || tsig != null || pc != null) {
      jdplus.tramoseats.base.api.tramo.AutoModelSpec.Builder ab
          = spec.getAutoModel().toBuilder();
      if (accept) {
        ab.acceptDefault(true);
      }
      if (pcr != null) {
        ab.pcr(pcr);
      }
      if (tsig != null) {
        ab.tsig(tsig);
      }
      if (pc != null) {
        ab.pc(pc);
      }
      b.autoModel(ab.build());
    }
    return b.build();
  }

  // ---------------- case parsing ----------------
  static String[] header;
  static double[] series;

  static void parseCase(String path) throws Exception {
    List<String> lines = java.nio.file.Files.readAllLines(java.nio.file.Paths.get(path));
    List<Double> vals = new ArrayList<>();
    for (String ln : lines) {
      if (ln.isEmpty()) {
        continue;
      }
      String[] t = ln.trim().split("\\s+");
      if (t[0].equals("KER")) {
        header = ln.trim().split("\\s+");
      } else if (t[0].equals("values")) {
        for (int i = 1; i < t.length; ++i) {
          vals.add(t[i].equals("nan") ? Double.NaN : Double.parseDouble(t[i]));
        }
      }
    }
    series = new double[vals.size()];
    for (int i = 0; i < vals.size(); ++i) {
      series[i] = vals.get(i);
    }
  }

  public static void main(String[] args) throws Exception {
    parseCase(args[0]);
    String mode = args.length > 1 ? args[1] : "both";
    System.out.println("mode " + mode);
    String id = header[1];
    // header: KER id spec [extra kv...] freq startY startM
    StringBuilder spectok = new StringBuilder(header[2]);
    for (int i = 3; i < header.length - 3; ++i) {
      spectok.append(' ').append(header[i]);
    }
    TramoSpec spec = resolveSpec(spectok.toString());
    int freq = Integer.parseInt(header[header.length - 3]);
    int y = Integer.parseInt(header[header.length - 2]);
    int mo = Integer.parseInt(header[header.length - 1]);
    TsPeriod start = freq == 12 ? TsPeriod.monthly(y, mo) : TsPeriod.quarterly(y, 1);
    TsData ts = TsData.of(start, DoubleSeq.of(series));

    List<String> out = new ArrayList<>();
    out.add("KER " + id + " specname " + header[2]
        + " freq " + freq + " start " + y + ":" + mo + " n " + series.length);
    out.add("KER " + id + " resolved autoModel " + (spec.isUsingAutoModel() ? 1 : 0)
        + " outliers " + (spec.getOutliers().isUsed() ? 1 : 0)
        + " ao " + (spec.getOutliers().isAo() ? 1 : 0)
        + " ls " + (spec.getOutliers().isLs() ? 1 : 0)
        + " tc " + (spec.getOutliers().isTc() ? 1 : 0)
        + " so " + (spec.getOutliers().isSo() ? 1 : 0)
        + " tcrate " + g17(spec.getOutliers().getDeltaTC())
        + " ml " + (spec.getOutliers().isMaximumLikelihood() ? 1 : 0)
        + " oc " + (spec.getTransform().isOutliersCorrection() ? 1 : 0)
        + " fct " + g17(spec.getTransform().getFct())
        + " accept " + (spec.getAutoModel().isAcceptDefault() ? 1 : 0)
        + " pcr " + g17(spec.getAutoModel().getPcr())
        + " tsig " + g17(spec.getAutoModel().getTsig())
        + " pc " + g17(spec.getAutoModel().getPc())
        + " cancel " + g17(spec.getAutoModel().getCancel())
        + " ub1 " + g17(spec.getAutoModel().getUb1())
        + " ub2 " + g17(spec.getAutoModel().getUb2())
        + " va " + g17(spec.getOutliers().getCriticalValue())
        + " tol " + g17(spec.getEstimate().getTol())
        + " ubp " + g17(spec.getEstimate().getUbp())
        + " meanTest " + (spec.getRegression().getMean().isTest() ? 1 : 0)
        + " tdAuto " + (spec.getRegression().getCalendar().getTradingDays().isAutomatic() ? 1 : 0)
        + " tdTest " + (spec.getRegression().getCalendar().getTradingDays().isTest() ? 1 : 0)
        + " easter " + (spec.getRegression().getCalendar().getEaster().isTest() ? 1 : 0)
        + " tfunc " + spec.getTransform().getFunction().name());

    // ---------- (A) REAL run ----------
    Observer obsR = new Observer("rsnap");
    TramoKernel kR = TramoKernel.of(spec, new ModellingContext());
    obsR.kernel = kR;
    RegSarimaModel realResult = null;
    boolean rthrew = false;
    if (!mode.equals("manual")) {
      try {
        realResult = kR.process(ts, obsR);
      } catch (RuntimeException e) {
        rthrew = true;
        out.add("rthrow " + e.getClass().getName() + " " + e.getMessage());
      }
    }
    if (mode.equals("real2")) { // second REAL run in the same JVM: ordering test
      Observer obsR2 = new Observer("r2snap");
      TramoKernel kR2 = TramoKernel.of(spec, new ModellingContext());
      obsR2.kernel = kR2;
      RegSarimaModel rr2 = null;
      try {
        rr2 = kR2.process(ts, obsR2);
      } catch (RuntimeException e) {
        out.add("r2throw " + e.getClass().getName() + " " + e.getMessage());
      }
      List<String> r2l = new ArrayList<>();
      dumpModel(r2l, "r2final", rr2);
      out.addAll(r2l);
      List<String> rl0 = new ArrayList<>();
      dumpModel(rl0, "r2final", realResult);
      out.add("r2identity " + (rl0.equals(r2l) ? 1 : 0));
    }
    out.addAll(obsR.evts);
    dumpModel(out, "rfinal", realResult);

    // ---------- (B) MANUAL run ----------
    Observer obsM = new Observer("msnap");
    TramoKernel kM = TramoKernel.of(spec, new ModellingContext());
    obsM.kernel = kM;
    Manual man = new Manual(kM, obsM);
    RegSarimaModel manResult = null;
    boolean mthrew = false;
    if (!mode.equals("real")) {
      try {
        manResult = man.run(ts);
      } catch (RuntimeException e) {
        mthrew = true;
        out.add("mthrow " + e.getClass().getName() + " " + e.getMessage());
      }
    }
    out.addAll(obsM.evts);
    out.addAll(man.lines);
    dumpModel(out, "mfinal", manResult);

    // ---------- identity ----------
    boolean evOk = true;
    if (obsR.evts.size() != obsM.evts.size()) {
      evOk = false;
    } else {
      for (int i = 0; i < obsR.evts.size(); ++i) {
        String a = obsR.evts.get(i), bb = obsM.evts.get(i);
        boolean as = a.startsWith("rsnap"), bs2 = bb.startsWith("msnap");
        if (as != bs2 || !(as ? a.substring(6).equals(bb.substring(6)) : a.equals(bb))) {
          evOk = false;
          break;
        }
      }
    }
    StringBuilder rb = new StringBuilder(), mb = new StringBuilder();
    List<String> rl = new ArrayList<>(), ml = new ArrayList<>();
    dumpModel(rl, "f", realResult);
    dumpModel(ml, "f", manResult);
    boolean finOk = rl.equals(ml) && rthrew == mthrew;
    out.add("identity " + id
        + " events " + (evOk ? 1 : 0)
        + " final " + (finOk ? 1 : 0)
        + " all " + ((evOk && finOk) ? 1 : 0));
    for (String s : out) {
      System.out.println(s);
    }
  }

  // ================= the decoded driver =================
  static final class Manual {
    final TramoKernel k;
    final Observer obs;
    final List<String> lines = new ArrayList<>();
    int seam = 0;

    // kernel state (mirrored into k's fields at every bytecode putfield)
    int pass, round;
    boolean pass3, needOutliers, needAutoModelling;
    double curva, pcr;
    RegSarimaModelling refAuto;
    ModelStatistics refStats;

    Object ctx;
    Object options;
    TramoSpec spec;

    Manual(TramoKernel k, Observer obs) {
      this.k = k;
      this.obs = obs;
      this.ctx = rf(k, "context");
      this.options = rf(k, "options");
      this.spec = (TramoSpec) rf(k, "spec");
    }

    void setPass(int v) { pass = v; wf(k, "pass", v); }
    void setRound(int v) { round = v; wf(k, "round", v); }
    void setPass3(boolean v) { pass3 = v; wf(k, "pass3", v); }
    void setNeedOut(boolean v) { needOutliers = v; wf(k, "needOutliers", v); }
    void setNeedAuto(boolean v) { needAutoModelling = v; wf(k, "needAutoModelling", v); }
    void setCurva(double v) { curva = v; wf(k, "curva", v); }
    void setPcr(double v) { pcr = v; wf(k, "pcr", v); }
    void setRefAuto(RegSarimaModelling v) { refAuto = v; wf(k, "refAuto", v); }
    void setRefStats(ModelStatistics v) { refStats = v; wf(k, "refStats", v); }

    boolean isAutoModelling() { return spec.isUsingAutoModel(); }
    boolean isOutliersDetection() { return spec.getOutliers().isUsed(); }

    double opt(String name) { return (double) rf(options, name); }
    boolean optB(String name) { return (boolean) rf(options, name); }
    double intermediatePrecision() { return opt("intermediatePrecision"); }
    double precision() { return opt("precision"); }

    void seam(String kind, String rest) {
      lines.add("seam " + seam + " " + kind + (rest == null ? "" : " " + rest));
      ++seam;
    }

    // process() transcription
    RegSarimaModel run(TsData ts) {
      obs.push("tramo");
      int sfreq = spec.getFrequency();
      if (sfreq != 0 && sfreq != ts.getAnnualFrequency()) {
        obs.error("frequency not allowed for this spec.");
        obs.pop();
        return null;
      }
      ModelDescription desc = (ModelDescription) call(k, "build",
          new Class<?>[]{TsData.class, ProcessingLog.class}, ts, obs);
      if (desc == null) {
        obs.error("initialization failed");
        obs.pop();
        return null;
      }
      RegSarimaModelling m = RegSarimaModelling.of(desc, obs);
      RegSarimaModel r = ami(m);
      obs.pop();
      return r;
    }

    // ami() transcription
    RegSarimaModel ami(RegSarimaModelling m) {
      if (isFullySpecified()) {
        m.getLog().info("fully specified model");
        m.estimate(precision());
        return m.build();
      }
      m.getLog().push("automatic model identification");
      call(k, "testSeasonality", new Class<?>[]{RegSarimaModelling.class}, m);
      boolean transOk = (boolean) call(k, "testTransformation",
          new Class<?>[]{RegSarimaModelling.class}, m);
      if (!transOk) {
        m.getLog().pop();
        return null;
      }
      Object regmod = call(k, "regressionModule", new Class<?>[]{boolean.class},
          m.getDescription().isAdjusted());
      Object regRes = regmod.getClass().getSuperclass() == null ? null : null;
      ProcessingResult rr = (ProcessingResult) invokeIface(regmod, "test",
          new Class<?>[]{RegSarimaModelling.class}, m);
      lines.add("prestate regtest " + rr.name());
      ModelDescription d0 = m.getDescription();
      lines.add("prestate orders " + ordersLine(d0.specification())
          + " mean " + (d0.isMean() ? 1 : 0)
          + " log " + (d0.isLogTransformation() ? 1 : 0)
          + " outliers " + outlierNames(d0));
      int nobs = d0.regarima().getActualObservationsCount();
      initProcessing(nobs);
      for (int i = 1; i <= 10; ++i) {
        if (iterate(m)) {
          break;
        }
      }
      RegSarimaModel r = m.build();
      m.getLog().pop();
      return r;
    }

    static Object invokeIface(Object target, String name, Class<?>[] sig, Object... args) {
      try {
        return target.getClass().getMethod(name, sig).invoke(target, args);
      } catch (ReflectiveOperationException e) {
        throw new RuntimeException(name, e);
      }
    }

    void initProcessing(int nobs) {
      setRound(0);
      setPass(0);
      setPass3(false);
      setNeedOut(isOutliersDetection());
      setNeedAuto(false);
      if (isOutliersDetection()) {
        setCurva(opt("va"));
        if (dcmpl(curva, 0.0) == 0) {
          setCurva(jdplus.tramoseats.base.core.tramo.internal.TramoUtility.calcCv(nobs));
        }
      }
      setPcr(opt("ljungBoxLimit"));
      setRefAuto(null);
      setRefStats(null);
      Object sc = rf(k, "scontroller");
      if (sc != null) {
        ((ModelController) sc).setReferenceModel(null);
      }
      @SuppressWarnings("unchecked")
      List<ModelController> ctrls = (List<ModelController>) rf(k, "controllers");
      for (ModelController c : ctrls) {
        c.setReferenceModel(null);
      }
      lines.add("init curva " + g17(curva) + " pcr " + g17(pcr));
    }

    boolean isFullySpecified() {
      return spec.getTransform().getFunction() != TransformationType.Auto
          && !isAutoModelling() && !isOutliersDetection()
          && !spec.getRegression().getMean().isTest()
          && !spec.getRegression().getCalendar().getTradingDays().isTest()
          && !spec.getRegression().getCalendar().getTradingDays().isAutomatic()
          && !spec.getRegression().getCalendar().getEaster().isTest();
    }

    boolean iterate(RegSarimaModelling m) {
      m.getLog().step("round " + pass);
      if (m.needEstimation()) {
        m.estimate(intermediatePrecision());
        seam("estimate", "taken 1 ip " + b(intermediatePrecision()));
      } else {
        seam("estimate", "taken 0 ip -");
      }
      ModelDescription desc = m.getDescription();
      boolean changed = false;
      SarimaOrders orders0 = desc.specification();
      boolean mean0 = desc.isMean();
      if (needDifferencing(desc)) {
        Object dm = call(k, "differencingModule", new Class<?>[]{});
        ProcessingResult res = (ProcessingResult) invokeIface(dm, "process",
            new Class<?>[]{RegSarimaModelling.class}, m);
        changed = res == ProcessingResult.Changed;
        seam("diff", "res " + res.name() + " orders " + ordersLine(desc.specification())
            + " mean " + (desc.isMean() ? 1 : 0) + " outs " + outlierNames(desc)
            + " pars " + parsBits(desc));
        SarimaOrders orders1 = desc.specification();
        if (pass == 1 && orders1.getDifferenceOrder() != orders0.getDifferenceOrder()) {
          desc.removeVariable(v -> ModellingUtility.isOutlier(v, true));
          m.clearEstimation();
          seam("diffdrop", "outliers cleared");
        }
      }
      if (needAutoModellingGuarded(desc)) {
        Object am = call(k, "armaModule", new Class<?>[]{});
        ProcessingResult res = (ProcessingResult) invokeIface(am, "process",
            new Class<?>[]{RegSarimaModelling.class}, m);
        desc = m.getDescription();
        SarimaOrders orders = desc.specification();
        if (orders.getParametersCount() == 0 && pass >= 3) {
          orders.setQ(1);
          desc.setSpecification(orders);
          seam("armaq1", "orders " + ordersLine(desc.specification()));
        }
        seam("arma", "res " + res.name() + " orders " + ordersLine(desc.specification())
            + " mean " + (desc.isMean() ? 1 : 0));
        changed = !(desc.specification().equals(orders0) && desc.isMean() == mean0);
      }
      if (needOutliersGuarded(desc)) {
        if (desc.removeVariable(v -> ModellingUtility.isOutlier(v, true))) {
          m.clearEstimation();
        }
        Object om = call(k, "outliersModule", new Class<?>[]{});
        ProcessingResult res = (ProcessingResult) invokeIface(om, "process",
            new Class<?>[]{RegSarimaModelling.class, double.class}, m, curva);
        seam("scan", "res " + res.name() + " cv " + g17(curva) + " orders "
            + ordersLine(m.getDescription().specification())
            + " mean " + (m.getDescription().isMean() ? 1 : 0)
            + " outs " + outlierNames(m.getDescription())
            + " pars " + parsBits(m.getDescription()));
        changed = changed || res == ProcessingResult.Changed;
      }
      if (!estimateModel(m)) {
        setNeedOut(isOutliersDetection());
        setNeedAuto(isAutoModelling());
        setRound(round + 1);
        setPass(pass + 1);
        seam("iterate", "estmodelfail");
        return false;
      }
      if (round == 0 && (isAutoModelling() || isOutliersDetection())) {
        setNeedOut(isOutliersDetection());
        setNeedAuto(isAutoModelling());
        setRound(round + 1);
        setPass(pass + 1);
        RegSarimaModelling ra = RegSarimaModelling.copyOf(m);
        setRefAuto(ra);
        ModelStatistics rs = ModelStatistics.of(ra.getDescription(),
            ra.getEstimation().getConcentratedLikelihood());
        setRefStats(rs);
        dumpStats(lines, "refstats", ra.getDescription(),
            ra.getEstimation().getConcentratedLikelihood());
        double lbp = refStats.getLjungBoxPvalue();
        double w = 1.0 - lbp;
        boolean accept = optB("acceptAirline") && dcmpg(w, opt("ljungBoxLimit")) < 0;
        seam("round0", "lbp " + g17(lbp) + " accept " + (accept ? 1 : 0));
        return accept;
      }
      if (pass <= 3 && !pass3 && isAutoModelling()) {
        if (!pass2(!changed, m)) {
          return false;
        }
      }
      if (!testRegression(m, spec.getAutoModel().getTsig(), "iterate")) {
        setPass(4);
        setNeedAuto(false);
        setNeedOut(false);
        seam("iterate", "regfail");
        return false;
      }
      ModelEstimator estimator = new ModelEstimator(precision(), curva,
          (jdplus.tramoseats.base.core.tramo.internal.OutliersDetectionModule)
              call(k, "outliersModule", new Class<?>[]{}));
      if (spec.isUsingAutoModel()) {
        Object scontroller = rf(k, "scontroller");
        ((SeasonalityController) scontroller).setEstimator(estimator);
        ProcessingResult sres = ((SeasonalityController) scontroller)
            .process(m, (TramoContext) ctx);
        seam("scontroller", "res " + sres.name() + " orders "
            + ordersLine(m.getDescription().specification())
            + " mean " + (m.getDescription().isMean() ? 1 : 0)
            + " outs " + outlierNames(m.getDescription())
            + " pars " + parsBits(m.getDescription()));
        if (sres == ProcessingResult.Changed && !pass3) {
          setPass3(true);
          setPass(1);
          setNeedAuto(true);
          setNeedOut(isOutliersDetection());
          seam("iterate", "pass3reset");
          return false;
        }
      }
      @SuppressWarnings("unchecked")
      List<ModelController> ctrls = (List<ModelController>) rf(k, "controllers");
      for (ModelController c : ctrls) {
        c.setEstimator(estimator);
        c.setReferenceModel(refAuto);
        ProcessingResult cres = c.process(m, (TramoContext) ctx);
        seam("controller", "res " + cres.name() + " orders "
            + ordersLine(m.getDescription().specification())
            + " mean " + (m.getDescription().isMean() ? 1 : 0)
            + " outs " + outlierNames(m.getDescription())
            + " pars " + parsBits(m.getDescription()));
      }
      seam("iterate", "done");
      return true;
    }

    boolean needDifferencing(ModelDescription desc) {
      if (!needAutoModelling) {
        return false;
      }
      if (round == 2 && !desc.variables().anyMatch(v -> ModellingUtility.isOutlier(v, true))) {
        return false;
      }
      SarimaOrders o = desc.specification();
      return o.getD() < 2 || o.getBd() >= 1;
    }

    boolean needAutoModellingGuarded(ModelDescription desc) {
      if (!needAutoModelling) {
        return false;
      }
      if (round == 2 && !desc.variables().anyMatch(v -> ModellingUtility.isOutlier(v, true))) {
        return false;
      }
      return true;
    }

    boolean needOutliersGuarded(ModelDescription desc) {
      return isOutliersDetection() && needOutliers;
    }

    boolean reduceVa() {
      if (dcmpl(curva, 2.0) == 0) {
        return false;
      }
      setCurva(Math.max(2.0, curva * (1.0 - opt("reduceVa"))));
      return true;
    }

    boolean estimateModel(RegSarimaModelling m) {
      FinalEstimator est = FinalEstimator.builder()
          .precision(precision())
          .unitRootThreshold(opt("ur"))
          .pass(pass)
          .ami(isAutoModelling())
          .outliers(isOutliersDetection())
          .build();
      int i = 0;
      while (true) {
        boolean ok = est.estimate(m);
        if (!ok) {
          if (pass == 1 && needOutliersGuarded(m.getDescription())
              && needAutoModellingGuarded(m.getDescription())) {
            reduceVa();
            seam("estmodelreduceva", "curva " + g17(curva));
          }
          seam("estmodel", "ret 0 pass " + pass + " nnsig " + est.getChangedParametersCount()
              + " orders " + ordersLine(m.getDescription().specification())
              + " outs " + outlierNames(m.getDescription())
              + " pars " + parsBits(m.getDescription()));
          return false;
        }
        seam("estmodel", "ret 1 pass " + pass + " nnsig " + est.getChangedParametersCount()
            + " orders " + ordersLine(m.getDescription().specification())
            + " outs " + outlierNames(m.getDescription())
            + " pars " + parsBits(m.getDescription()));
        if (pass3 || pass != 0) {
          return true;
        }
        if (i >= 5) {
          return true;
        }
        ++i;
        if (testRegression(m, 1.96, "retry" + i)) {
          return true;
        }
      }
    }

    boolean testRegression(RegSarimaModelling m, double tsig, String site) {
      FastRegressionTest frt = FastRegressionTest.builder()
          .testMean(optB("checkMu"))
          .meanThreshold(tsig)
          .build();
      ProcessingResult res = frt.test(m);
      double[] td = tstatDump(m);
      seam("regtest", "site " + site + " tsig " + g17(tsig) + " res " + res.name()
          + " mean " + (m.getDescription().isMean() ? 1 : 0)
          + " t " + bs(td));
      return res == ProcessingResult.Unchanged;
    }

    double[] tstatDump(RegSarimaModelling m) {
      try {
        jdplus.toolkit.base.core.regarima.RegArimaEstimation<
            jdplus.toolkit.base.core.sarima.SarimaModel> est = m.getEstimation();
        if (est == null) {
          return new double[0];
        }
        ConcentratedLikelihoodWithMissing cl = est.getConcentratedLikelihood();
        int k = m.getDescription().getArimaSpec().freeParametersCount();
        double[] t = cl.tstats(k, true);
        double[] r = new double[t.length];
        for (int i = 0; i < r.length; ++i) {
          r[i] = t[i];
        }
        return r;
      } catch (RuntimeException e) {
        return new double[0];
      }
    }

    boolean pass2(boolean changed, RegSarimaModelling m) {
      ModelDescription desc = m.getDescription();
      SarimaOrders orders = desc.arima().orders();
      ModelStatistics stats = ModelStatistics.of(desc,
          m.getEstimation().getConcentratedLikelihood());
      lines.add("p2dbg orders " + ordersLine(desc.specification())
          + " arima " + ordersLine(orders)
          + " e " + bs(((ConcentratedLikelihoodWithMissing) ((jdplus.toolkit.base.core.regarima
              .RegArimaEstimation<jdplus.toolkit.base.core.sarima.SarimaModel>)
              m.getEstimation()).getConcentratedLikelihood()).e().toArray()).substring(0, 80)
          + " pars " + bs(desc.arima().parameters().toArray()));
      dumpStats(lines, "pass2", desc, m.getEstimation().getConcentratedLikelihood());
      lines.add("statdbg fresh lbp " + g17(stats.getLjungBoxPvalue())
          + " se " + g17(stats.getSe()) + " lb " + g17(stats.getLjungBox()));
      double lbp = 1.0 - stats.getLjungBoxPvalue();
      double se = stats.getSe();
      if (refAuto != null) {
        double rlbp = 1.0 - refStats.getLjungBoxPvalue();
        double rse = refStats.getSe();
        int nb = desc.countRegressors(v -> ModellingUtility.isOutlier(v, true));
        int nr = refAuto.getDescription().countRegressors(v -> ModellingUtility.isOutlier(v, true));
        boolean revert = false;
        String why = "-";
        if (!(nr > nb)) {
          if (dcmpg(lbp, 0.95) < 0 && dcmpg(rlbp, 0.75) < 0 && dcmpg(rse, se) < 0) {
            revert = true;
            why = "A";
          } else if (pass == 1 && dcmpl(lbp, 0.95) >= 0 && dcmpg(rlbp, 0.95) < 0) {
            revert = true;
            why = "B";
          } else if (dcmpg(lbp, 0.95) < 0 && dcmpg(rlbp, 0.75) < 0
              && dcmpg(lbp, rlbp) < 0 && dcmpg(rse, 1.0 * se) >= 0) {
            revert = true;
            why = "C";
          } else if (dcmpl(lbp, 0.95) >= 0 && dcmpg(rlbp, 0.95) < 0
              && dcmpg(rse, 1.0 * se) >= 0) {
            revert = true;
            why = "D";
          } else if (orders.getD() == 0 && orders.getBd() == 1 && orders.getP() == 1
              && dcmpg(desc.arima().phi(1), -0.82) <= 0
              && orders.getQ() <= 1 && orders.getBp() == 0 && orders.getBq() == 1) {
            revert = true;
            why = "E";
          } else if (orders.getD() == 1 && orders.getBd() == 0 && orders.getP() == 0
              && orders.getQ() == 1 && orders.getBp() == 1
              && dcmpg(desc.arima().bphi(1), -0.65) <= 0
              && orders.getBq() <= 1) {
            revert = true;
            why = "F";
          }
        } else {
          why = "skip-nrnb";
        }
        if (!revert) {
          setRefAuto(RegSarimaModelling.copyOf(m));
          setRefStats(stats);
        } else {
          // TramoKernel.restore(refAuto): m.set(copyOf(refAuto.desc), refAuto.est)
          m.set(ModelDescription.copyOf(refAuto.getDescription()),
              refAuto.getEstimation());
          lbp = rlbp;
        }
        seam("pass2gate", "lbp " + g17(lbp) + " se " + g17(se)
            + " rlbp " + g17(rlbp) + " rse " + g17(rse)
            + " nb " + nb + " nr " + nr + " revert " + (revert ? 1 : 0)
            + " why " + why);
      } else {
        setRefAuto(RegSarimaModelling.copyOf(m));
        setRefStats(stats);
        seam("pass2gate", "first lbp " + g17(lbp) + " se " + g17(se));
      }
      if (pass == 1) {
        setPcr(pcr + 0.025);
      } else if (pass >= 2) {
        setPcr(pcr + 0.015);
      }
      if (dcmpg(lbp, pcr) <= 0) {
        seam("pass2", "ret 1 pcr " + g17(pcr));
        return true;
      }
      if (pass == 1 && isOutliersDetection()) {
        reduceVa();
        m.getLog().remark("reduction of the critical value for outliers detection");
      }
      setRound(round + 1);
      setPass(pass + 1);
      if (pass <= 2) {
        setNeedAuto(!changed);
        setNeedOut(isOutliersDetection());
      } else {
        lastSolution(m);
        setNeedAuto(false);
      }
      seam("pass2", "ret 0 pcr " + g17(pcr));
      return false;
    }

    void lastSolution(RegSarimaModelling m) {
      ModelDescription desc = m.getDescription();
      SarimaOrders orders = desc.specification();
      orders.setP(3);
      if (orders.getBd() > 0) {
        orders.setBp(0);
      }
      orders.setQ(1);
      if ((boolean) rf(ctx, "seasonal")) {
        orders.setBq(1);
      }
      desc.removeVariable(v -> ModellingUtility.isOutlier(v, true));
      m.setSpecification(orders);
      m.getLog().warning("last chance model");
      setRound(1);
      setNeedOut(isOutliersDetection());
      setNeedAuto(false);
    }
  }
}
