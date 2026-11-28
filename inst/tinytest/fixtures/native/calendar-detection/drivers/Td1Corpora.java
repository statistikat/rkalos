import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;

/**
 * m17-td1 replay corpora dump: the exact TsData values the Td1Probe battery
 * consumed (generators copied verbatim), one value per line at %.17g. The
 * native replay probe reads this file so both engines see bit-identical
 * input (StrictMath-generated doubles must not be regenerated in C++).
 */
public class Td1Corpora {

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

  static void dump(String kind, double[] v, int y0, int m0) {
    TsData ts = TsData.of(TsPeriod.monthly(y0, m0), DoubleSeq.of(v)).cleanExtremities();
    System.out.println("series " + kind + " n=" + ts.length() + " start=" + y0 + "-" + m0);
    for (int i = 0; i < ts.length(); ++i) {
      System.out.println(String.format("%.17g", ts.getValues().get(i)));
    }
  }

  public static void main(String[] args) {
    dump("clean", clean(96), 1990, 1);
    dump("tdeff", tdeff(156), 1990, 1);
  }
}
