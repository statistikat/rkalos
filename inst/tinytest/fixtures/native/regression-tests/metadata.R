# Capture provenance for the item 21-G regression-tests + mean-rules corpus.
# Oracle: the pinned jdplus 3.9.0 stage-R machines through the public-API
# observer RegressionTestsProbe.java (no decision rule there): the design
# built with the pinned factories (GenericTradingDaysVariable over
# GenericTradingDays.contrasts(TD4/TD2), LengthOfPeriod(LeapYear),
# EasterVariable(duration, PreComputed, endPosition -1), UserVariable),
# the estimation through RegSarimaModelling.estimate(1e-7) =
# processor(true, 1e-7), the concentrated-likelihood statistics through
# coefficients/unscaledCovariance/ssq/dim/tstats(0,false), the
# RegressionVariablesTest built by its public builder exactly as
# X13SpecDecoder.readRegression builds it (tdTest(1.96, derived) iff the
# TD spec test is active, movingHolidaysTest(1.96), meanTest 1.96 initial
# / 1.0 final; ftd stays 0: the decoder never calls tdJointTest),
# MeanController(cval).test for cval in {1.96 (round-0 checkMu), 1.6
# (redoEstimation checkMu)} with the mean-added estimation duplicated
# through the same processor(true, 1.0E-5) call for the MEANSTAT rows
# only, and ModelController(0.95).accept/getRTval for the rtval0 the
# final mean rule consumes. OBS rows call the pinned TRegressionTest/
# DerivedRegressionTest accept forms for observation only and never feed
# the process() run (which always executes on a fresh description copy).
# 18 sweep cases, 0 estimation-error rows (recorded, excluded).
# The qualified period domain is {4, 12}; the mean-added MEANSTAT rows
# carry (coef0, unscaledCov00, ssq, dim, tstat(0,0,false)).
# rJava JVM 21.0.11 (observer compiled with javac --release 21)
# double capture: this run was cross-checked against an independent JVM
# capture of the same cases (--repro): all 7 csv files byte-identical: yes
