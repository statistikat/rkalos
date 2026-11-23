// 23-F F.2 probe: the MULTI-ROUND outlier state machine of the pinned
// jdplus 3.9.0 stack, observed two ways per case from ONE JVM:
//  (A) the REAL run: OutliersDetectionModule.process (public entry) - the
//      lines it emits mirror the 23-B probe exactly, so the final state is
//      byte-comparable against the committed tramo-outlier-scan corpus (the
//      jar/inputs identity gate);
//  (B) the MANUAL run: OutliersDetectionModule.make via reflection -> the
//      FastOutliersDetector private machine driven step by step
//      (clear/estimateModel/verifyModel/addOutlier/sod.process) with a full
//      state dump after EVERY call (round, rflag, estimationStep,
//      backwardStep, exit, outliers, lastremoved, coeff, tstats, res, the
//      regarima orders/parameters/X count/mean, and after a scan the sod
//      state el/mad/maxT/posMax/oMax/coef + allowed-cell counts).
// The manual loop is a line-by-line transcription of FastOutliersDetector.
// process's bytecode (tmp_opencode/23f-jd/FOD.txt offsets 0-241, see pinning
// SECTION 26.3 + the 23f2 decode notes): if the loop shape were wrong, the
// manual run's final outliers would diverge from the real run's - the probe
// emits an identity line comparing both (dump-not-trust).
//
// Case grammar: IDENTICAL to the 23-B probe (SCAN header + values lines), so
// validation/evidence/native/tramo-outlier-scan/inputs.csv feeds cases 1:1.
import jdplus.toolkit.base.api.arima.SarimaOrders;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.TimeSelector;
import jdplus.toolkit.base.api.timeseries.regression.Variable;
import jdplus.toolkit.base.core.data.interpolation.AverageInterpolator;
import jdplus.toolkit.base.core.modelling.regression.AdditiveOutlierFactory;
import jdplus.toolkit.base.core.modelling.regression.LevelShiftFactory;
import jdplus.toolkit.base.core.modelling.regression.PeriodicOutlierFactory;
import jdplus.toolkit.base.core.modelling.regression.TransitoryChangeFactory;
import jdplus.toolkit.base.core.regsarima.ami.FastOutliersDetector;
import jdplus.toolkit.base.core.regsarima.regular.IOutliersDetectionModule;
import jdplus.toolkit.base.core.regsarima.regular.ModelDescription;
import jdplus.toolkit.base.core.regsarima.regular.ProcessingResult;
import jdplus.toolkit.base.core.regsarima.regular.RegSarimaModelling;
import jdplus.toolkit.base.core.regarima.RegArimaModel;
import jdplus.toolkit.base.core.sarima.SarimaModel;
import jdplus.tramoseats.base.core.tramo.internal.OutliersDetectionModule;
import jdplus.tramoseats.base.core.tramo.internal.TramoUtility;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

public final class Scratch23g {

  static String g(double v) {
    return String.format(Locale.ROOT, "%.17g", v);
  }

  static long bits(double v) {
    return Double.doubleToRawLongBits(v);
  }

  // %.17g#hexbits - the bit-level token (23-B discipline).
  static String gb(double v) {
    return g(v) + "#" + Long.toHexString(bits(v));
  }

  static double tok(String s) {
    return s.equals("nan") ? Double.NaN : Double.parseDouble(s);
  }

  static final class Observer implements ProcessingLog {
    final List<String> evts = new ArrayList<>();
    int depth = 0;
    long seq = 0;

    static String detail(Object d) {
      if (d == null) {
        return "-";
      }
      if (d instanceof Double) {
        return "num " + g((Double) d);
      }
      if (d instanceof IOutliersDetectionModule.Info) {
        IOutliersDetectionModule.Info i = (IOutliersDetectionModule.Info) d;
        StringBuilder s = new StringBuilder("odmInfo cv=").append(g(i.getVa()))
            .append(" types[");
        String[] t = i.getTypes();
        for (int k = 0; k < t.length; ++k) {
          if (k > 0) {
            s.append(',');
          }
          s.append(t[k]);
        }
        s.append("] outliers[");
        int[][] o = i.getOutliers();
        for (int k = 0; k < o.length; ++k) {
          if (k > 0) {
            s.append(',');
          }
          s.append(Arrays.toString(o[k]));
        }
        return s.append(']').toString();
      }
      if (d instanceof Object[]) {
        StringBuilder s = new StringBuilder("arr[");
        Object[] a = (Object[]) d;
        for (int i = 0; i < a.length; ++i) {
          if (i > 0) {
            s.append(',');
          }
          s.append(String.valueOf(a[i]));
        }
        return s.append(']').toString();
      }
      return "obj " + String.valueOf(d);
    }

    void emit(String kind, String msg, Object details) {
      ++seq;
      evts.add("evt " + seq + " " + kind + " " + depth + " "
          + (msg == null ? "-" : msg.replace('\n', ' ')) + " det " + detail(details));
    }

    public void push(String s) {
      emit("push", s, null);
      ++depth;
    }

    public void pop() {
      if (depth > 0) {
        --depth;
      }
      emit("pop", null, null);
    }

    public void error(String m, String o, Object d) {
      emit("error", m, d);
    }

    public void warning(String m, String o, Object d) {
      emit("warn", m, d);
    }

    public void remark(String m, String o, Object d) {
      emit("remark", m, d);
    }

    public void info(String m, String o, Object d) {
      emit("info", m, d);
    }

    public void step(String m, Object d) {
      emit("step", m, d);
    }

    public boolean isVerbose() {
      return true;
    }
  }

  @SuppressWarnings("unchecked")
  static String varsOf(ModelDescription desc) {
    StringBuilder s = new StringBuilder();
    for (Variable<?> vr : (java.util.List<? extends Variable<?>>) (java.util.List<?>) desc
        .variables().toList()) {
      if (s.length() > 0) {
        s.append(';');
      }
      s.append(vr.getName()).append('{');
      boolean first = true;
      for (java.util.Map.Entry<String, String> e : vr.getAttributes().entrySet()) {
        if (!first) {
          s.append(',');
        }
        first = false;
        s.append(e.getKey()).append('=').append(e.getValue());
      }
      s.append('}');
    }
    return s.length() == 0 ? "-" : s.toString();
  }

  // ------------------------------------------------------------------ build
  static final class Built {
    ModelDescription desc;
    RegSarimaModelling modelling;
    Observer log;
    OutliersDetectionModule mod;
    double cv;
  }

  static Built build(String[] h, double[] y) {
    Built b = new Built();
    int freq = Integer.parseInt(h[2]);
    boolean mean = h[3].equals("y");
    boolean logTr = h[4].equals("y");
    int d = Integer.parseInt(h[5]);
    int D = Integer.parseInt(h[6]);
    int p = Integer.parseInt(h[7]);
    int q = Integer.parseInt(h[8]);
    int bp = Integer.parseInt(h[9]);
    int bq = Integer.parseInt(h[10]);
    String fac = h[11];
    boolean ao = fac.charAt(0) == 'y';
    boolean ls = fac.charAt(1) == 'y';
    boolean tc = fac.charAt(2) == 'y';
    boolean so = fac.charAt(3) == 'y';
    double tcrate = tok(h[12]);
    double prec = tok(h[13]);
    boolean ml = h[14].equals("y");
    int maxround = Integer.parseInt(h[15]);
    int maxout = Integer.parseInt(h[16]);
    String span = h[17];
    b.cv = h[18].equals("calc") ? TramoUtility.calcCv(y.length) : tok(h[18]);
    boolean interp = h[19].equals("y");
    List<String> pres = new ArrayList<>();
    if (h.length > 20 && !h[20].equals("-")) {
      for (String pt : h[20].split("\\+")) {
        pres.add(pt);
      }
    }
    b.log = new Observer();
    TsPeriod start = TsPeriod.monthly(1990, 1);
    TsData ts = TsData.of(start, DoubleSeq.of(y));
    b.desc = new ModelDescription(ts, TsDomain.of(start, y.length));
    b.desc.setMean(mean);
    if (logTr) {
      b.desc.setLogTransformation(true);
    }
    SarimaOrders orders = new SarimaOrders(freq);
    orders.setRegular(p, d, q);
    orders.setSeasonal(bp, D, bq);
    b.desc.setSpecification(orders);
    if (interp) {
      b.desc.interpolate(AverageInterpolator.interpolator());
    }
    TimeSelector sel = TimeSelector.all();
    if (span.startsWith("idx:") || span.startsWith("idxe:")) {
      boolean epoch = span.startsWith("idxe:");
      String[] lr = span.substring(epoch ? 5 : 4).split("-");
      long l = Long.parseLong(lr[0]);
      long h2 = Long.parseLong(lr[1]);
      long base = epoch ? 0 : start.getId();
      sel = TimeSelector.between(start.dateAt(base + l), start.dateAt(base + h2));
    }
    b.mod = OutliersDetectionModule.builder()
        .ao(ao).ls(ls).tc(tc).so(so).maximumLikelihood(ml)
        .tcrate(tcrate).precision(prec).span(sel)
        .maxOutliers(maxout).maxRound(maxround)
        .build();
    for (String pre : pres) {
      String[] pc = pre.split(":");
      java.time.LocalDateTime at = start.dateAt(start.getId() + Long.parseLong(pc[1]));
      jdplus.toolkit.base.api.timeseries.regression.IOutlier core;
      switch (pc[0]) {
        case "AO" -> core = AdditiveOutlierFactory.FACTORY.make(at);
        case "LS" -> core = LevelShiftFactory.FACTORY_ZEROENDED.make(at);
        case "TC" -> core = new TransitoryChangeFactory(tcrate).make(at);
        default -> core = new PeriodicOutlierFactory(freq, true).make(at);
      }
      b.desc.addVariable(Variable.<jdplus.toolkit.base.api.timeseries.regression.IOutlier>builder()
          .name(pc[0] + "." + pc[1]).core(core).build());
    }
    b.modelling = RegSarimaModelling.of(b.desc, b.log);
    return b;
  }

  // ------------------------------------------------------------- reflection
  static Field f(Class<?> c, String name) throws Exception {
    Field fd = c.getDeclaredField(name);
    fd.setAccessible(true);
    return fd;
  }

  static Method m(Class<?> c, String name, Class<?>... args) throws Exception {
    Method mm = c.getDeclaredMethod(name, args);
    mm.setAccessible(true);
    return mm;
  }

  static final Class<?> FOD = FastOutliersDetector.class;
  static Field F_regarima, F_tstats, F_round, F_rflag, F_backwardStep, F_exit,
      F_estimationStep, F_lastremoved, F_coeff, F_res, F_outliers, F_mvx, F_cv,
      F_maxRound, F_maxOutliers, F_sod;
  static Method M_clear, M_estimateModel, M_verifyModel, M_addOutlier, M_removeOutlier;

  // SingleOutlierDetector / FastOutlierDetector internals (state dumps).
  static Field F_el, F_mad, F_stmodel, F_T, F_allowed, F_lbound, F_ubound,
      F_posMax, F_oMax;
  static Method M_getMaxTStat, M_getMaxType, M_getMaxPos, M_coefficient;

  static void initRefl() throws Exception {
    F_regarima = f(FOD, "regarima");
    F_tstats = f(FOD, "tstats");
    F_round = f(FOD, "round");
    F_rflag = f(FOD, "rflag");
    F_backwardStep = f(FOD, "backwardStep");
    F_exit = f(FOD, "exit");
    F_estimationStep = f(FOD, "estimationStep");
    F_lastremoved = f(FOD, "lastremoved");
    F_coeff = f(FOD, "coeff");
    F_res = f(FOD, "res");
    F_outliers = f(FOD, "outliers");
    F_mvx = f(FOD, "mvx");
    F_cv = f(FOD, "cv");
    F_maxRound = f(FOD, "maxRound");
    F_maxOutliers = f(FOD, "maxOutliers");
    F_sod = f(FOD, "sod");
    M_clear = m(FOD, "clear");
    M_estimateModel = m(FOD, "estimateModel",
        jdplus.toolkit.base.core.arima.estimation.IArimaMapping.class);
    M_verifyModel = m(FOD, "verifyModel");
    M_addOutlier = m(FOD, "addOutlier", int.class, int.class, double.class);
    M_removeOutlier = m(FOD, "removeOutlier", int.class);
    Class<?> sodc = Class.forName("jdplus.toolkit.base.core.regarima.outlier."
        + "SingleOutlierDetector");
    Class<?> fodc = Class.forName("jdplus.toolkit.base.core.regarima.outlier."
        + "FastOutlierDetector");
    F_el = f(fodc, "el");
    F_mad = f(fodc, "mad");
    F_stmodel = f(fodc, "stmodel");
    F_T = f(sodc, "T");
    F_allowed = f(sodc, "allowedTable");
    F_lbound = f(sodc, "lbound");
    F_ubound = f(sodc, "ubound");
    F_posMax = f(sodc, "posMax");
    F_oMax = f(sodc, "oMax");
    M_getMaxTStat = m(sodc, "getMaxTStat");
    M_getMaxType = m(sodc, "getMaxOutlierType");
    M_getMaxPos = m(sodc, "getMaxOutlierPosition");
    M_coefficient = m(sodc, "coefficient", int.class, int.class);
  }

  static String vec(double[] v) {
    if (v == null) {
      return "null";
    }
    StringBuilder s = new StringBuilder(Integer.toString(v.length));
    for (double x : v) {
      s.append(' ').append(gb(x));
    }
    return s.toString();
  }

  static String vec(jdplus.toolkit.base.api.data.DoubleSeq v) {
    if (v == null) {
      return "null";
    }
    return vec(v.toArray());
  }

  static String istat(int[] o) throws Exception {
    if (o == null) {
      return "null";
    }
    return Integer.toString(o.length);
  }

  // One full FOD + sod state dump.
  @SuppressWarnings("unchecked")
  static void dumpState(String tag, Object drv, RegArimaModel<?> ra) throws Exception {
    StringBuilder s = new StringBuilder("st ");
    s.append(tag)
        .append(" round ").append(F_round.getInt(drv))
        .append(" rflag ").append(F_rflag.getBoolean(drv) ? 1 : 0)
        .append(" estep ").append(F_estimationStep.getBoolean(drv) ? 1 : 0)
        .append(" bstep ").append(F_backwardStep.getBoolean(drv) ? 1 : 0)
        .append(" exit ").append(F_exit.getBoolean(drv) ? 1 : 0);
    ArrayList<int[]> outs = (ArrayList<int[]>) F_outliers.get(drv);
    s.append(" outs ").append(outs.size());
    for (int[] o : outs) {
      s.append(' ').append(Arrays.toString(o));
    }
    s.append(" lastrem ").append(istat((int[]) F_lastremoved.get(drv)));
    s.append(" coeff ").append(vec((DoubleSeq) F_coeff.get(drv)));
    s.append(" tstats ").append(vec((double[]) F_tstats.get(drv)));
    s.append(" res ").append(vec((DoubleSeq) F_res.get(drv)));
    SarimaModel sm = (SarimaModel) ra.arima();
    s.append(" orders ").append(sm.getP()).append(' ').append(sm.getD()).append(' ')
        .append(sm.getQ()).append(' ').append(sm.getBp()).append(' ')
        .append(sm.getBd()).append(' ').append(sm.getBq());
    s.append(" pars ").append(vec(sm.parameters().toArray()));
    s.append(" xcount ").append(ra.getX().size());
    s.append(" mean ").append(ra.isMean() ? 1 : 0);
    System.out.println(s);
  }

  // The scan-relevant sod state after a scan (or a failed one).
  static void dumpSod(String tag, Object drv) throws Exception {
    Object sod = F_sod.get(drv);
    StringBuilder s = new StringBuilder("sod ");
    s.append(tag);
    double[] el = (double[]) F_el.get(sod);
    s.append(" el ").append(el == null ? "null" : Integer.toString(el.length));
    s.append(" mad ").append(gb(F_mad.getDouble(sod)));
    s.append(" lb ").append(F_lbound.getInt(sod)).append(" ub ")
        .append(F_ubound.getInt(sod));
    Object T = F_T.get(sod);
    s.append(" T ").append(T == null ? "null" : "ok");
    Object at = F_allowed.get(sod);
    if (at == null) {
      s.append(" alw null");
    } else {
      // Allowed-cell counts per factory column (mask persistence pin).
      jdplus.toolkit.base.api.util.TableOfBoolean tb
          = (jdplus.toolkit.base.api.util.TableOfBoolean) at;
      s.append(" alw ").append(tb.getRowsCount()).append('/').append(tb.getColumnsCount());
      for (int c = 0; c < tb.getColumnsCount(); ++c) {
        int cnt = 0;
        for (int r = 0; r < tb.getRowsCount(); ++r) {
          if (tb.get(r, c)) {
            ++cnt;
          }
        }
        s.append(' ').append(cnt);
      }
    }
    String maxt;
    try {
      maxt = gb(((Double) M_getMaxTStat.invoke(sod)).doubleValue());
    } catch (java.lang.reflect.InvocationTargetException e) {
      maxt = "threw:" + e.getCause().getClass().getName();
    }
    s.append(" maxT ").append(maxt);
    try {
      int ty = ((Integer) M_getMaxType.invoke(sod)).intValue();
      int po = ((Integer) M_getMaxPos.invoke(sod)).intValue();
      s.append(" pick ").append(po).append(',').append(ty);
      if (ty >= 0 && po >= 0) {
        s.append(" cf ").append(gb(((Double) M_coefficient.invoke(sod, po, ty)).doubleValue()));
      }
    } catch (java.lang.reflect.InvocationTargetException e) {
      s.append(" pick threw:").append(e.getCause().getClass().getName());
    }
    System.out.println(s);
  }

  // The manual step machine: a line-by-line transcription of
  // FastOutliersDetector.process's bytecode (FOD.txt offsets 0-241, catch at
  // 238 -> return false). Returns the final outliers echo.
  @SuppressWarnings("unchecked")
  static String manualRun(Object drv, RegArimaModel<?> ra0,
      jdplus.toolkit.base.core.arima.estimation.IArimaMapping<?> mapping) throws Exception {
    // offset 0-1: clear()
    M_clear.invoke(drv);
    // offsets 14-17: regarima = arg
    F_regarima.set(drv, ra0);
    dumpState("init", drv, (RegArimaModel<?>) F_regarima.get(drv));
    double cv = F_cv.getDouble(drv);
    int maxRound = F_maxRound.getInt(drv);
    int maxOutliers = F_maxOutliers.getInt(drv);
    Object sod = F_sod.get(drv);
    Method sodProcess = sod.getClass().getMethod("process", RegArimaModel.class);
    boolean ret = false;
    String note = "run";
    boolean reachedTail = false;
    int rounds = 0;
    while (true) {
      // loop head, offset 19: estimateModel
      if (!(Boolean) M_estimateModel.invoke(drv, mapping)) {
        // offsets 24-28: return false (no tail)
        ret = false;
        note = "est";
        break;
      }
      dumpState("est", drv, (RegArimaModel<?>) F_regarima.get(drv));
      // offset 29: ans = true
      // offsets 29-55: verifyModel is a GATE on the scan, not an alternative
      // (bytecode 36 ifeq 55 / fall-through to the if (ans) test at 55).
      boolean ans = true;
      if (F_backwardStep.getBoolean(drv)) {
        ans = (Boolean) M_verifyModel.invoke(drv);
        dumpState("verify", drv, (RegArimaModel<?>) F_regarima.get(drv));
        if (F_exit.getBoolean(drv)) {
          note = "exit";
          reachedTail = true;
          break;
        }
      }
      if (ans) {
        // offsets 60-74: the scan
        boolean ok = (Boolean) sodProcess.invoke(sod, F_regarima.get(drv));
        dumpState("scan", drv, (RegArimaModel<?>) F_regarima.get(drv));
        dumpSod("scan", drv);
        if (!ok) {
          note = "sodfalse";
          reachedTail = true;
          break;
        }
        // offsets 77-84: round++
        F_round.setInt(drv, F_round.getInt(drv) + 1);
        // offsets 87-109: dcmpg(|maxT|, cv) < 0 -> tail
        double maxT;
        try {
          maxT = ((Double) M_getMaxTStat.invoke(sod)).doubleValue();
        } catch (java.lang.reflect.InvocationTargetException e) {
          // searchMax left (-1,-1): the table read throws - the real process
          // catches RuntimeException at 238 and returns false (no tail).
          ret = false;
          note = "maxtthrew:" + e.getCause().getClass().getName();
          break;
        }
        if (dcmpg(Math.abs(maxT), cv) < 0) {
          note = "gate";
          reachedTail = true;
          break;
        }
        // offsets 112-146: addOutlier(position, type, coefficient(pos, type))
        int ty = ((Integer) M_getMaxType.invoke(sod)).intValue();
        int po = ((Integer) M_getMaxPos.invoke(sod)).intValue();
        double cf = ((Double) M_coefficient.invoke(sod, po, ty)).doubleValue();
        M_addOutlier.invoke(drv, po, ty, cf);
        System.out.println("add pos " + po + " type " + ty + " cf " + gb(cf));
        dumpState("add", drv, (RegArimaModel<?>) F_regarima.get(drv));
        // offsets 149-163: outliers.size() == maxOutliers -> tail
        if (((ArrayList<int[]>) F_outliers.get(drv)).size() == maxOutliers) {
          note = "cap";
          reachedTail = true;
          break;
        }
      }
      // offsets 166-174: round < maxRound -> back to 19, else tail
      if (F_round.getInt(drv) >= maxRound) {
        note = "roundcap";
        reachedTail = true;
        break;
      }
      if (++rounds > 400) {
        ret = false;
        note = "stepcap";
        break;
      }
    }
    if (reachedTail) {
      // offsets 177-214: the tail estimateModel on a cap/exit/gate break
      boolean tailCond = F_exit.getBoolean(drv)
          || F_round.getInt(drv) == maxRound
          || ((ArrayList<int[]>) F_outliers.get(drv)).size() == maxOutliers;
      if (tailCond) {
        M_estimateModel.invoke(drv, mapping); // result popped
        System.out.println("tailest 1");
        dumpState("tailest", drv, (RegArimaModel<?>) F_regarima.get(drv));
      }
      // offset 215: estimationStep = false
      F_estimationStep.setBoolean(drv, false);
      // offsets 220-233: while (!verifyModel()) estimateModel (popped)
      int drain = 0;
      while (!(Boolean) M_verifyModel.invoke(drv)) {
        dumpState("drainv", drv, (RegArimaModel<?>) F_regarima.get(drv));
        M_estimateModel.invoke(drv, mapping); // popped
        dumpState("draine", drv, (RegArimaModel<?>) F_regarima.get(drv));
        if (++drain > 400) {
          ret = false;
          note = "draincap";
          reachedTail = false;
          break;
        }
      }
      if (reachedTail) {
        dumpState("drainend", drv, (RegArimaModel<?>) F_regarima.get(drv));
        ret = true;
      }
    }
    System.out.println("manret " + (ret ? 1 : 0) + " note " + note
        + " round " + F_round.getInt(drv));
    int[][] outs = (int[][]) FOD.getMethod("getOutliers").invoke(drv);
    StringBuilder s = new StringBuilder("manouts ");
    for (int k = 0; k < outs.length; ++k) {
      if (k > 0) {
        s.append(',');
      }
      s.append(Arrays.toString(outs[k]));
    }
    System.out.println(s);
    return s.toString();
  }

  static int dcmpg(double a, double b) {
    if (Double.isNaN(a) || Double.isNaN(b)) {
      return 1;
    }
    return Double.compare(a, b);
  }

  public static void main(String[] args) throws Exception {
    initRefl();
    java.nio.file.Path file = java.nio.file.Paths.get(args[0]);
    List<String> all = java.nio.file.Files.readAllLines(file);
    List<String> buf = new ArrayList<>();
    for (String line : all) {
      String t = line.trim();
      if (t.isEmpty()) {
        if (!buf.isEmpty()) {
          runCase(buf);
          buf.clear();
        }
      } else {
        buf.add(t);
      }
    }
    if (!buf.isEmpty()) {
      runCase(buf);
    }
  }

  static void runCase(List<String> lines) {
    try {
      String hdr = null;
      List<String> vals = new ArrayList<>();
      for (String line0 : lines) {
        String line = line0.trim();
        if (line.isEmpty()) {
          continue;
        }
        if (line.startsWith("values ")) {
          for (String t : line.substring(7).trim().split("\\s+")) {
            vals.add(t);
          }
        } else if (hdr == null) {
          hdr = line;
        } else {
          throw new IllegalStateException("unexpected line: " + line);
        }
      }
      if (hdr == null) {
        throw new IllegalStateException("missing header");
      }
      String[] h = hdr.split("\\s+");
      double[] y = new double[vals.size()];
      for (int i = 0; i < vals.size(); ++i) {
        y[i] = vals.get(i).equals("nan") ? Double.NaN : Double.parseDouble(vals.get(i));
      }
      System.out.println("case " + h[1]);
      StringBuilder rq = new StringBuilder("req");
      for (String t : h) {
        rq.append(' ').append(t);
      }
      System.out.println(rq);
      System.out.print("vals");
      for (double v : y) {
        System.out.print(" ");
        System.out.print(gb(v));
      }
      System.out.println();

      // (A) the REAL run - lines mirror the 23-B probe.
      Built b = build(h, y);
      TsDomain ed = b.desc.getEstimationDomain();
      System.out.println("estdom start " + ed.getStartPeriod().toString() + " len "
          + ed.getLength());
      System.out.println("pre_vars " + varsOf(b.desc));
      ProcessingResult r = b.mod.process(b.modelling, b.cv);
      System.out.println("real_status " + r.name());
      System.out.println("real_vars " + varsOf(b.desc));
      System.out.println("real_need_est " + (b.modelling.needEstimation() ? 1 : 0));
      System.out.println("real_spec " + b.desc.specification());
      for (String e : b.log.evts) {
        System.out.println("real " + e);
      }

      // (B) the MANUAL run over a freshly built description.
      Built b2 = build(h, y);
      Method mk = OutliersDetectionModule.class.getDeclaredMethod("make",
          ModelDescription.class, double.class);
      mk.setAccessible(true);
      Object drv = mk.invoke(b2.mod, b2.desc, b2.cv);
      if (drv == null) {
        System.out.println("make null");
      } else {
        StringBuilder s = new StringBuilder("mrd");
        s.append(" cv ").append(gb(F_cv.getDouble(drv)))
            .append(" mvx ").append(F_mvx.getBoolean(drv) ? 1 : 0)
            .append(" maxround ").append(F_maxRound.getInt(drv))
            .append(" maxout ").append(F_maxOutliers.getInt(drv));
        Object sod = F_sod.get(drv);
        s.append(" lb ").append(F_lbound.getInt(sod)).append(" ub ")
            .append(F_ubound.getInt(sod));
        String[] types = (String[]) FOD.getMethod("outlierTypes").invoke(drv);
        s.append(" types[");
        for (int k = 0; k < types.length; ++k) {
          if (k > 0) {
            s.append(',');
          }
          s.append(types[k]);
        }
        s.append(']');
        System.out.println(s);
        RegArimaModel<?> ra0 = b2.desc.regarima();
        manualRun(drv, ra0, b2.desc.mapping());
      }
      System.out.println("end " + h[1]);
    } catch (Throwable t) {
      System.out.println("probe-threw " + t.getClass().getName() + " " + t.getMessage());
      for (StackTraceElement e : t.getStackTrace()) {
        System.out.println("at " + e);
      }
    }
  }
}
