# Capture provenance for the item 21-H.0 final-information corpus.
# Oracle: pinned jdplus 3.9.0 classes through the observer
# FinalInformationProbe.java (package-private access to RegArmaSsqFunction;
# no stencil, scaling or inverse algorithm is re-implemented there).
# est rows run RegSarimaComputer.optimize with computeExactFinalDerivatives
# (true) and the default LevenbergMarquardtMinimizer, capture the final
# Point, then abort unless a replay of the pinned tail reproduces the
# estimator point/score/information/asymptotic covariance bit for bit.
# inj rows inject the final point (boundary, eps-zero, NaN, negative-nobs
# and singular-information families the optimiser cannot reach on the
# pinned jars, the 21-A precedent).
# Tail: RegArmaSsqFunction.ssqEvaluate -> allCoefficients -> calcResiduals
# -> RegArmaFunction (deviance) -> NumericalDerivatives(fp, false) stencil,
# score grad * ((-0.5 nobs) / fx), information hess * ((0.5 nobs) / fx),
# asymptotic covariance SymmetricMatrix.inverse(information); nobs =
# len(y_d) - xcols - dim unless overridden (nobs column of cases.csv).
# capture.csv: yd/xd are the differenced response/design rows, *_hex rows
# are raw-bit hex of the JVM doubles (row-major flats), status ok/rej:.../
# invfail. nodes.csv holds every stencil node: key = raw-bit hex of the
# node coordinates, value = raw-bit hex or THROW (newval catch -> fx).
# Parameters are KALOS raw convention [phi bphi theta btheta].
# double capture: byte-identical csv files against the independent JVM
# capture in /tmp/opencode/fi-capC (captured_at differs only in metadata.R).
# R 4.5.3
# rjd3jars 0.0.6
# rJava JVM 21.0.11 (observer compiled with javac --release 21)
# ground truth: bytecode (FdLibm) semantics; capture JVM runs with
# _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog"
# (guard verified in-process: Math.exp == StrictMath.exp at the probe argument).
# capture-final-information.R sha256 4df5871fc09fa96af9cc7b56365b6784670a2064efbc7ac4ac5f156df92ee732
# FinalInformationProbe.java sha256 7c58f935e0e59d90253b9166b6c6d0e4ac3749842a2321258aadef0006f3fc36
captured_at <- "2026-11-19 15:40:44"
