import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.modelling.TransformationType;
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.calendars.LengthOfPeriodType;
import jdplus.toolkit.base.api.timeseries.regression.ModellingContext;
import jdplus.toolkit.base.api.timeseries.regression.Variable;
import jdplus.toolkit.base.core.modelling.GeneralLinearModel;
import jdplus.toolkit.base.core.regsarima.regular.RegSarimaModel;
import jdplus.toolkit.base.core.regsarima.regular.ModelDescription;
import jdplus.toolkit.base.core.sarima.SarimaModel;
import jdplus.toolkit.base.core.stats.likelihood.LikelihoodStatistics;
import jdplus.x13.base.api.regarima.EasterSpec;
import jdplus.x13.base.api.regarima.RegArimaSpec;
import jdplus.x13.base.api.regarima.RegressionSpec;
import jdplus.x13.base.api.regarima.TransformSpec;
import jdplus.x13.base.api.regarima.TradingDaysSpec;
import jdplus.x13.base.api.x13.X13Spec;
import jdplus.x13.base.r.X13;
import java.util.Arrays;

/**
 * m17-td1 stage-1 probe (21-I stage 1): the CEDM ATD step WITHOUT Easter on
 * the fit-matrix corpora (clean/tdeff/air). Templates derive from the pinned
 * rsa4/rsa5c proto bytes with easter cleared:
 *   td2noe   rsa4 minus easter            (TD2, leap year, TEST_REMOVE, autoAdjust)
 *   td7noe   rsa5c minus easter           (TD7, leap year, TEST_REMOVE, autoAdjust)
 *   td2lev   td2noe with transform NONE   (LoP column arm at the landing)
 *   td2noadj td2noe with autoAdjust=false (LoP column arm under log)
 * Emits the full ProcessingInfo stream (ATD infos with the full
 * LikelihoodStatistics[] arrays at %.17g), the landed regression part, the
 * preadjustment state (LoP fork witness) and the final model facts.
 * Mode: run <td2noe|td7noe|td2lev|td2noadj> <clean|tdeff|air>
 */
public class Td1Probe {

  static final String HEX_RSA4
      = "0ab7010a060a0010011801120b080a1900000000000000c01a1f0a040a02414f0a040a024c530a040a025443120029666666666666e63f301e220c18012202100130013a0210012a4c080111666666666666ee3f19000000000000f03f21280af4893c49c23f29cdccccccccccf03f31d151812fae7ef03f392fbae8a28b2ef23f419a9999999999b93f4967202a1dd933f03f580132161a0a08011001300250016200220808011008180152003a0b0a001148afbc9af2d77a3e1224100122010029000000000000f83f31000000000000044038ffffffffffffffffff0160031a16100219000000000000f03f21000000000000f03f2801";
  static final String HEX_RSA5C
      = "0ab7010a060a0010011801120b080a1900000000000000c01a1f0a040a02414f0a040a024c530a040a025443120029666666666666e63f301e220c18012202100130013a0210012a4c080111666666666666ee3f19000000000000f03f21280af4893c49c23f29cdccccccccccf03f31d151812fae7ef03f392fbae8a28b2ef23f419a9999999999b93f4967202a1dd933f03f580132161a0a08061001300250016200220808011008180152003a0b0a001148afbc9af2d77a3e1224100122010029000000000000f83f31000000000000044038ffffffffffffffffff0160031a16100219000000000000f03f21000000000000f03f2801";

  static final double[] AIR = {
    112.0,118.0,132.0,129.0,121.0,135.0,148.0,148.0,136.0,119.0,104.0,118.0,115.0,126.0,141.0,135.0,125.0,149.0,170.0,170.0,158.0,133.0,114.0,140.0,145.0,150.0,178.0,163.0,172.0,178.0,199.0,199.0,184.0,162.0,146.0,166.0,171.0,180.0,193.0,181.0,183.0,218.0,230.0,242.0,209.0,191.0,172.0,194.0,196.0,196.0,236.0,235.0,229.0,243.0,264.0,272.0,237.0,211.0,180.0,201.0,204.0,188.0,235.0,227.0,234.0,264.0,302.0,293.0,259.0,229.0,203.0,229.0,242.0,233.0,267.0,269.0,270.0,315.0,364.0,347.0,312.0,274.0,237.0,278.0,284.0,277.0,317.0,313.0,318.0,374.0,413.0,405.0,355.0,306.0,271.0,306.0,315.0,301.0,356.0,348.0,355.0,422.0,465.0,467.0,404.0,347.0,305.0,336.0,340.0,318.0,362.0,348.0,363.0,435.0,491.0,505.0,404.0,359.0,310.0,337.0,360.0,342.0,406.0,396.0,420.0,472.0,548.0,559.0,463.0,407.0,362.0,405.0,417.0,391.0,419.0,461.0,472.0,535.0,622.0,606.0,508.0,461.0,390.0,432.0
  };

  static String g(double d) {
    return String.format("%.17g", d);
  }

  static byte[] unhex(String h) {
    byte[] b = new byte[h.length() / 2];
    for (int i = 0; i < b.length; ++i) {
      b[i] = (byte) Integer.parseInt(h.substring(2 * i, 2 * i + 2), 16);
    }
    return b;
  }

  static double[] clean(int n) {
    double[] y = new double[n];
    for (int t = 0; t < n; ++t) {
      double k = t + 1;
      y[t] = 50 + 0.05 * k + 5 * StrictMath.sin(2 * StrictMath.PI * k / 12)
          + 3 * StrictMath.cos(4 * StrictMath.PI * k / 12)
          + 0.5 * StrictMath.sin(14 * StrictMath.PI * k / 12);
    }
    return y;
  }

  static double[] tdeff(int n) {
    double[] y = new double[n];
    int year = 1990, month = 1;
    int k = 0;
    for (int t = 0; t < n; ++t) {
      int len = java.time.LocalDate.of(year, month, 1).lengthOfMonth();
      double s = 0;
      for (int d = 1; d <= len; ++d) {
        java.time.LocalDate dt = java.time.LocalDate.of(year, month, d);
        int wd = dt.getDayOfWeek().getValue() - 1;
        double w;
        switch (wd) {
          case 0: w = 0.9; break;
          case 1: w = 1.15; break;
          case 2: w = 0.8; break;
          case 3: w = 1.25; break;
          case 4: w = 0.6; break;
          default: w = 0.0; break;
        }
        double noise = 0.15 * StrictMath.sin(k * 0.7 + 1) + 0.1 * StrictMath.sin(k * 2.3);
        s += (20.0 + 0.02 * k + noise) * w;
        ++k;
      }
      y[t] = s;
      ++month;
      if (month == 13) {
        month = 1;
        ++year;
      }
    }
    return y;
  }

  static TsData ts(String kind) {
    double[] v;
    int y0, m0;
    switch (kind) {
      case "air":
        v = AIR; y0 = 1949; m0 = 1; break;
      case "tdeff":
        v = tdeff(156); y0 = 1990; m0 = 1; break;
      default:
        v = clean(96); y0 = 1990; m0 = 1; break;
    }
    return TsData.of(TsPeriod.monthly(y0, m0), DoubleSeq.of(v));
  }

  static X13Spec spec(String tpl) {
    X13Spec base = X13.specOf(unhex(tpl.equals("td7noe") ? HEX_RSA5C : HEX_RSA4));
    RegArimaSpec reg = base.getRegArima();
    RegressionSpec rg = reg.getRegression();
    TradingDaysSpec td = rg.getTradingDays();
    RegressionSpec.Builder rgp = rg.toBuilder().easter(EasterSpec.none());
    if (tpl.equals("td2noadj")) {
      rgp.tradingDays(TradingDaysSpec.td(td.getTradingDaysType(),
          td.getLengthOfPeriodType(), td.getRegressionTestType(), false));
    }
    RegArimaSpec.Builder newReg = reg.toBuilder().regression(rgp.buildWithoutValidation());
    if (tpl.equals("td2lev")) {
      TransformSpec tf = reg.getTransform().toBuilder()
          .function(TransformationType.None).buildWithoutValidation();
      newReg = newReg.transform(tf);
    }
    return base.toBuilder().regArima(newReg.buildWithoutValidation())
        .buildWithoutValidation();
  }

  static final class Observer implements ProcessingLog {
    int depth = 0;
    long seq = 0;

    void emit(String kind, String msg, Object details) {
      ++seq;
      System.out.print("evt " + seq + " " + kind + " " + depth + " "
          + (msg == null ? "-" : msg.replace('\n', ' ')));
      dumpInfo(details);
      System.out.println();
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

    public boolean isVerbose() {
      return true;
    }
  }

  static void dumpInfo(Object d) {
    if (d == null) {
      System.out.print("det -");
      return;
    }
    if (d instanceof jdplus.sa.base.core.regarima.AutomaticTradingRegressionModule.Info) {
      jdplus.sa.base.core.regarima.AutomaticTradingRegressionModule.Info i
          = (jdplus.sa.base.core.regarima.AutomaticTradingRegressionModule.Info) d;
      System.out.print("det ATD names=" + Arrays.toString(i.getNames()) + " best=" + i.getBest());
      LikelihoodStatistics[] ll = i.getLl();
      for (int k = 0; k < ll.length; ++k) {
        System.out.print(" ll[" + k + "] n=" + ll[k].getObservationsCount()
            + " eff=" + ll[k].getEffectiveObservationsCount()
            + " np=" + ll[k].getEstimatedParametersCount()
            + " ssq=" + g(ll[k].getSsqErr())
            + " bic2=" + g(ll[k].getBIC2())
            + " bicc=" + g(ll[k].getBICC())
            + " aicc=" + g(ll[k].getAICC()));
      }
      return;
    }
    System.out.print("det obj " + String.valueOf(d));
  }

  public static void main(String[] args) throws Exception {
    String tpl = args[0];
    String kind = args[1];
    X13Spec spec = spec(tpl);
    RegArimaSpec reg = spec.getRegArima();
    RegressionSpec rg = reg.getRegression();
    TradingDaysSpec td = rg.getTradingDays();
    System.out.println("specinfo " + tpl + " td=" + td.getTradingDaysType()
        + " lp=" + td.getLengthOfPeriodType() + " test=" + td.getRegressionTestType()
        + " autoAdjust=" + td.isAutoAdjust() + " automatic=" + td.isAutomatic()
        + " easter=" + rg.getEaster().getType() + "/" + rg.getEaster().getTest()
        + " fn=" + reg.getTransform().getFunction()
        + " aicdiff=" + g(rg.getAicDiff()));
    TsData data = ts(kind).cleanExtremities();
    jdplus.x13.base.core.x13.X13Kernel k =
        jdplus.x13.base.core.x13.X13Kernel.of(spec, new ModellingContext());
    Observer obs = new Observer();
    jdplus.x13.base.core.x13.X13Results r;
    long t0 = System.nanoTime();
    try {
      r = k.process(data, obs);
    } catch (Exception ex) {
      System.out.println("run THROW " + ex.getClass().getName() + " " + ex.getMessage());
      return;
    }
    long t1 = System.nanoTime();
    System.out.println("timing " + tpl + " " + kind + " elapsed_ms=" + (t1 - t0) / 1000000);
    if (r == null) {
      System.out.println("run NULL results (kernel failure path)");
      return;
    }
    RegSarimaModel m = r.getPreprocessing();
    GeneralLinearModel.Description<?> dsc = m.getDescription();
    System.out.println("fact " + tpl + "/" + kind + " log " + dsc.isLogTransformation()
        + " lopTransform " + dsc.getLengthOfPeriodTransformation());
    DoubleSeq ser = dsc.getSeries().getValues();
    StringBuilder sh = new StringBuilder("fact " + tpl + "/" + kind + " serieshead");
    for (int i = 0; i < Math.min(6, ser.length()); ++i) {
      sh.append(' ').append(g(ser.get(i)));
    }
    System.out.println(sh);
    Variable<?>[] vars = dsc.getVariables();
    System.out.println("fact " + tpl + "/" + kind + " variables " + vars.length);
    for (int i = 0; i < vars.length; ++i) {
      Variable<?> v = vars[i];
      System.out.println("var " + tpl + "/" + kind + " " + i + " name=" + v.getName()
          + " dim=" + v.dim() + " core=" + v.getCore().getClass().getSimpleName());
    }
    SarimaModel ar = m.arima();
    if (ar != null) {
      DoubleSeq pars = ar.parameters();
      StringBuilder ps = new StringBuilder();
      for (int i = 0; i < pars.length(); ++i) {
        if (i > 0) {
          ps.append(' ');
        }
        ps.append(g(pars.get(i)));
      }
      System.out.println("fact " + tpl + "/" + kind + " armapars " + ps);
      System.out.println("fact " + tpl + "/" + kind + " orders " + ar.orders().getP() + " "
          + ar.orders().getD() + " " + ar.orders().getQ() + " " + ar.orders().getBp() + " "
          + ar.orders().getBd() + " " + ar.orders().getBq());
    }
    GeneralLinearModel.Estimation e = m.getEstimation();
    if (e != null) {
      LikelihoodStatistics ls = e.getStatistics();
      System.out.println("fact " + tpl + "/" + kind + " ll " + g(ls.getLogLikelihood())
          + " adjll " + g(ls.getAdjustedLogLikelihood()) + " aicc " + g(ls.getAICC())
          + " ssq " + g(ls.getSsqErr()) + " npar " + ls.getEstimatedParametersCount());
      DoubleSeq co = e.getCoefficients();
      StringBuilder cs = new StringBuilder("coefs");
      for (int i = 0; i < co.length(); ++i) {
        cs.append(' ').append(g(co.get(i)));
      }
      System.out.println("fact " + tpl + "/" + kind + " " + cs);
    }
  }
}
