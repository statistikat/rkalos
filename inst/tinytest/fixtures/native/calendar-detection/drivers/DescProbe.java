// m17-td1 stage-1 introspection: print the initial ModelDescription the
// decoder's IModelBuilder produces for the TD-only templates, so the native
// transformation prepass can replicate the JVM log-test inputs exactly.
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.regression.ModellingContext;
import jdplus.toolkit.base.api.timeseries.regression.Variable;
import jdplus.toolkit.base.core.regsarima.regular.IModelBuilder;
import jdplus.toolkit.base.core.regsarima.regular.ModelDescription;
import jdplus.x13.base.core.x13.regarima.RegArimaKernel;
import java.lang.reflect.Field;
import java.util.List;

public class DescProbe {
  @SuppressWarnings("unchecked")
  public static void main(String[] args) throws Exception {
    String tpl = args[0], kind = args[1];
    TsData ts = Td1Probe.ts(kind);
    RegArimaKernel k =
        RegArimaKernel.of(Td1Probe.spec(tpl).getRegArima(), new ModellingContext());
    Field f = RegArimaKernel.class.getDeclaredField("modelBuilder");
    f.setAccessible(true);
    IModelBuilder mb = (IModelBuilder) f.get(k);
    ModelDescription desc = mb.build(ts, ProcessingLog.dummy());
    System.out.println("tpl=" + tpl + " kind=" + kind + " null=" + (desc == null));
    if (desc == null) return;
    System.out.println("log=" + desc.isLogTransformation() + " has_td=" + desc.contains("td")
        + " has_lp=" + desc.contains("lp") + " has_easter=" + desc.contains("easter"));
    Field fv = ModelDescription.class.getDeclaredField("variables");
    fv.setAccessible(true);
    for (Variable v : (List<Variable>) fv.get(desc))
      System.out.println("var name=" + v.getName() + " core=" + v.getClass().getSimpleName());
    Field fm = ModelDescription.class.getDeclaredField("mean");
    fm.setAccessible(true);
    System.out.println("mean=" + fm.get(desc));
    Field fl = ModelDescription.class.getDeclaredField("lpTransformation");
    fl.setAccessible(true);
    System.out.println("lpTransformation=" + fl.get(desc));
    System.out.println("est_domain=" + desc.getEstimationDomain().getLength());
  }
}
