# Capture provenance for the item 21-F arma-refinement corpus.
# Oracle: the pinned jdplus 3.9.0 x13 AutoModellingModule stage machine
# through the public-API observer ArmaRefinementProbe.java (no decision
# rule there): the current estimation through RegArimaUtility.processor
# (true, 1e-7).process, the linearized data through RegArimaUtility
# .linearizedData, the differencing decision through the decoder-built
# DifferencingModule.process (public; tmean through getTMean), the module
# own BackFilter(UnitRoots.D(d).times(UnitRoots.D(period, bd))) filter
# application, ArmaModule.select plus getPreferedModels (the package-
# private Builder reached by reflection, the RegArmaBic arima field read
# for observation only), and SarimaOrders.of/airline/equals for the final
# comparison. Every stage failure reproduces the module RuntimeException
# catch as a fallback row (RESULT|fallback: SarimaOrders
# .airline(period), mean off). Module scalars are decoder call-site
# values, verified identical across cases through the MODULE line:
#   maxd 2 maxbd 1 ub1 0.96999999999999990 ub2 0.87999999999999990 cancel 0.10000000000000000 precision 1.0000000000000000e-07
#   nmod 5 maxP 2 maxQ 2 maxBp 1 maxBq 1 balanced 0 mixed 1 wn 0 eps 1.0000000000000000e-07
# estimation.csv rows with status=error (Gls failed on annihilated
# regressor designs) are excluded from the sweep: the capture never built
# a case on them; the JVM error text is kept for the record.
# 23 sweep cases (21 normal, 2 airline fallback), 0 estimation-error rows.
# The qualified period domain is {4, 12}; failures are recorded rows.
# rJava JVM 21.0.11 (observer compiled with javac --release 21)
# double capture executed: the installed files are the byte-identical output
# of two independent JVM processes (all csv and this metadata compared with
# cmp: no differing byte), so every captured decision, landing and panel is
# JVM-process reproducible; SHA256SUMS.csv pins the bytes.
