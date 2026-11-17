# Capture provenance for the item 21-A estimation-control corpus.
# Oracle: pinned jdplus RegSarimaComputer.process through the public-API
# observer EstimationControlProbe.java (no estimation algorithm here).
# R 4.5.3
# rjd3jars 0.0.6
# Pinned control: eps=1e-7 feps=1e-4 ml=true mt=false cdf=true fast=true,
# default Levenberg-Marquardt minimizer builder; startingPoint Zero/Default.
# Regression columns carry the mean/drift (regression-likelihood contract);
# missing responses are zero-filled for the probe and indexed in inputs.csv
# (one-based). Parameters are KALOS convention: stored [ar buar ma buma]
# with AR signs flipped from the pinned SarimaModel.parameters().
# cases.csv orders column is "free_ar/free_ma/d/D".
# capture-estimation-control.R sha256 da35cee9619a6f8772116ffb7d903dee71012262747ca48a1e36b2ab1002ad94
# EstimationControlProbe.java sha256 bcfd7edebab8d1d7628e23e08d13d2c3f1d4fc9e2bc722b1ee6d20a81edf4d3f
captured_at <- "2026-11-16 22:51:35"
