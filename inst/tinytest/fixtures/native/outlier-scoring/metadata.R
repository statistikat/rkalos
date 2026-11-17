# Capture provenance for the item 21-C outlier-scoring corpus.
# Oracle: pinned jdplus 3.9.0 classes through the public-API observer
# OutlierScoringProbe.java (no scoring, scan or selection algorithm there).
# Grids: Utility.intProbability(x, Lower), Normal.getProbabilityInverse(p,
# Lower) and X13Utility.calcCv(n) next to advancedComputer(0.05).
# Scan: ExactSingleOutlierDetector(Mad(false), Ansley default, MLJ-ML
# residuals), factories AO/LS(zero-ended)/TC(rate(freq,tcrate))/SO(zero-
# ended), prepare(n), setBounds(l, r), exclusions, then process(regarima);
# T/coefficient/allowed tables printed column-major.
# Driver: ExactOutliersDetector.builder with X13Utility.processor(true,
# 1.0E-7), caller cv/maxRound/maxOutliers (module decoder defaults 50/30),
# process(regarima, null) - the fully free map of the qualified 21-A
# convention. changed, kept (position, factory) pairs, the final tstats
# decision vector (reflective read, no public getter), the final design
# column count and final SarimaModel.parameters() are printed.
# Parameters are KALOS convention: stored [ar sar ma sma] with the AR signs
# flipped from the pinned raw parameters; driver cases start from zero
# coefficients at the orders pinned in driver-inputs.csv (p_ar/...).
# double capture: byte-identical csv files against the independent JVM
# capture in /tmp/opencode/recap-a (captured_at differs only in metadata.R).
# R 4.5.3
# rjd3jars 0.0.6
# rJava JVM 21.0.11 (observer compiled with javac --release 21)
# capture-outlier-scoring.R sha256 1395126f2f87cf6c0f06097b432972811847ecb6195eed95fbd64361e89ec356
# OutlierScoringProbe.java sha256 dc2e9b09d347b2f8f9374b36e3e84f7425e1b94b1b4f90002ed3c08d3dc7cb0a
captured_at <- "2026-11-17 13:28:27"
