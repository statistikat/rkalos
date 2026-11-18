# Capture provenance for the item 21-E initial-arma-identification corpus.
# Oracle: the pinned jdplus 3.9.0 x13 ArmaModule through the public-API
# observer InitialArmaIdentificationProbe.java (no decision rule there):
# olsResiduals for the cursor, RegArmaBic.sort for one candidate with the
# private fitted SarimaModel read by reflection, the processor(true, eps)
# route for the likelihood statistics, and ArmaModule.select with its
# getPreferedModels panel. The module is assembled through the pinned
# package-private Builder (X13SpecDecoder.readAutoModel builds the same
# object); the decoder call-site state is modelsCount 5, max_p/max_q 2,
# max_bp/max_bq 1, balanced false, mixed true, acceptWhiteNoise false,
# eps 1.0E-7; variants vary those fields to exercise the refinement rule
# groups and the merge eviction path. The stage grids are fed from the
# JVM's own first-BIC-minimum in the module's grid order, so the
# candidates battery witnesses exactly the grids the module builds.
# Candidates report the specification, BIC2, the raw logLikelihood, the
# estimated-parameters count and the effective n (statistics route), and
# the fitted SarimaModel raw parameters [phi, bphi, theta, btheta] under
# the pinned 1 + phi(B) sign convention.
# The qualified period domain is {4, 12}; the pinned period-1 shortcuts are
# out of scope. Failures are recorded as one error line per case.
# double capture: run this script twice into separate directories and
# compare all csv files; identical bytes prove JVM-process reproducibility.
# R 4.5.3
# ground truth: bytecode (FdLibm) semantics; capture JVM runs with
# _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog"
# (guard verified in-process: Math.exp == StrictMath.exp at the probe argument).
# rjd3jars 0.0.6
# rJava JVM 21.0.11 (observer compiled with javac --release 21)
# capture-initial-arma-identification.R sha256 125b16e43eb51538fde83095c936d3b2b11cfc354df07bed31a2a5712900d92b
# InitialArmaIdentificationProbe.java sha256 b2a371b376f36db4329cfbab9c51d1b1d2a69d7b1e26bf7517c75abe1fbf0ac9
captured_at <- "2026-11-18 11:23:06"
