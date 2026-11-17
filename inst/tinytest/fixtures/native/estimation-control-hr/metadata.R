# Capture provenance for the item 21-A.2 estimation-control HR corpus.
# Oracle: pinned jdplus RegSarimaComputer.process through the public-API
# observer EstimationControlProbe.java (no estimation algorithm here).
# R 4.5.3
# rjd3jars 0.0.6
# Pinned control: ml=true mt=false cdf=true fast=true, default
# Levenberg-Marquardt minimizer builder; startingPoint HannanRissanen/
# Multiple; eps/feps per case (1e-7/1e-4 builder defaults, 1e-9 for the
# pinned PROCESSOR case). Complete series and fully free designs only.
# Regression columns carry the mean/drift (regression-likelihood contract).
# Parameters are KALOS convention: stored [ar buar ma buma]
# with AR signs flipped from the pinned SarimaModel.parameters().
# cases.csv orders column is "free_ar/free_ma/d/D".
# capture-estimation-control-hr.R sha256 d495ee6971ed18c40328a1ba1b9d8f15f4fb1ba5a7b650f312668dd5ff666c19
# EstimationControlProbe.java sha256 bcfd7edebab8d1d7628e23e08d13d2c3f1d4fc9e2bc722b1ee6d20a81edf4d3f
captured_at <- "2026-11-16 23:24:23"
