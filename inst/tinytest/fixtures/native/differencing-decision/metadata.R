# Capture provenance for the item 21-D differencing-decision corpus.
# Oracle: the pinned jdplus 3.9.0 x13 DifferencingModule through the
# public-API observer DifferencingDecisionProbe.java (no decision rule there).
# The module is built as X13SpecDecoder.readAutoModel builds it: cancel
# spec.getCancel() = 0.1, ub1 = 1.0 / spec.getUb1() = 1/1.0309278350515465,
# ub2 = 1.0 / spec.getUb2() = 1/1.1363636363636365 (JVM divisions), precision
# 1.0E-7, maxd/maxbd builder defaults 2/1, ml never set; cases.csv stores the
# effective scalars exactly as the module stored them.
# Each case runs the public process(DoubleSeq, int) on a complete finite
# series and records d, bd, mean, tmean plus the reflective end-state
# observation: iter, mlused, useml, rmax, rsmax, c_, final SarimaOrders
# p,q,bp,bq and the final SarimaModel raw parameters [phi, bphi, theta,
# btheta] (pinned 1 + phi(B) convention). Errors are recorded as one line.
# Pure-constant and non-finite responses are out of scope (the pinned raw
# HannanRissanen pipeline keeps running on NaN; the qualified Kalos seam
# rejects them - see the contract).
# double capture: byte-identical csv files against the independent JVM
# capture in /data/home/kowa/seasonaladjustment/kalos/validation/evidence/native/differencing-decision (captured_at differs only in metadata.R).
# R 4.5.3
# ground truth: bytecode (FdLibm) semantics; capture JVM runs with
# _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog"
# (guard verified in-process: Math.exp == StrictMath.exp at the probe argument).
# rjd3jars 0.0.6
# rJava JVM 21.0.11 (observer compiled with javac --release 21)
# capture-differencing-decision.R sha256 d26eaf2d8ba54278ab02fd660fbd2aaee92fa24f83d0a927e6344af265d089a4
# DifferencingDecisionProbe.java sha256 9b4f2f854ec5d038a9dd6ab9b818877a6859880e4f322eec78fe9145efe5aed7
captured_at <- "2026-11-18 08:56:44"
