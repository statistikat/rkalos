# Capture provenance for the item 21-B transformation-decision corpus.
# Oracle: pinned jdplus sa-core LogLevelModule.process through the public-
# API observer TransformationDecisionProbe.java (no decision rule or
# estimation algorithm there). ModelDescription: mean off, monthly TsData,
# no estimation domain, user-variable design columns.
# Module call site pinned by X13SpecDecoder: estimationPrecision 1.0E-7,
# preadjust None, outliersCorrection false; AICC correction -2.0
# (TransformSpec DEF_AICDIFF). Decision processor: DemetraUtility.
# processor(true, 1e-7) = ML, Multiple starting point, feps 1e-4 default.
# Pinned rule: log iff aiccLevel > aiccLog + correction (dcmpl; NaN picks
# levels); short-circuit to levels on any value <= 0 (NaN passes the scan).
# Parameters are KALOS convention: stored [ar sar ma sma]
# with AR signs flipped from the pinned SarimaModel.parameters().
# inputs.csv x holds the regarima() design matrix actually consumed
# (post-check), column-major; xcolumns is its column count.
# R 4.5.3
# rjd3jars 0.0.6
# capture-transformation-decision.R sha256 f3e3521cbbcf37f3ee74ac5bd14123e55526aa6ba4e7c61f5eaf4cf1abda6db4
# TransformationDecisionProbe.java sha256 4ac456840bc5c557dbf603046ba995a59d21266db72a3144bd14975b3383fd57
captured_at <- "2026-11-17 00:45:55"
