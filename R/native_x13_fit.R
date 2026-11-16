# Item 22c M3-finale: the restricted fixed-model X-13 fit — one composition
# of the qualified pieces in the pinned X13Kernel.process order: the
# fixed-coefficient preadjust (native_x13_preadjust's seam), the pass-level
# X-11 selection run over `out` with the updateSpec-forced mode, the
# log-additive results-builder stage where the mode is log_additive, and the
# finals restore (native_x13_restore's seam). The delivered surface is
# restricted by the pinned contract, and the boundary mirrors the core's
# fail-closed answers verbatim instead of inventing R-side refusals: only
# fixed seasonal requests (a pinned SeasonalFilterOption per calendar period;
# "x11default"/"msr" answer "unsupported:x13 fit fixed seasonal requests
# only"), a fixed Henderson trend_length (an odd length >= 5; the automatic 0
# answers "unsupported:x13 fit fixed Henderson trend length"), the pinned
# fixed-path horizon domain backcast=forecast=0 (anything else answers
# "unsupported:x13 fit requires zero horizons"), and a REQUIRED sigma policy:
# NULL (or omitted) is passed through untouched — no default is ever applied
# here — and the core answers "unsupported:x13 fit requires an explicit sigma
# policy". Every refusal still echoes the full provenance of the request.
# Estimation, missing values, the length-of-period transform, benchmarking,
# diagnostics and the ARIMA route stay outside; capabilities() reports the
# restricted_fit scope for callers to gate on.
native_x13_fit <- function(series,variables=list(),log=FALSE,seasonal,trend_length,
    sigma=NULL,msr_policy="jdplus",frequency=NULL,start_period=NULL,lower_sigma=1.5,
    upper_sigma=2.5,mode="additive",backcast=0,forecast=0) {
  if(missing(seasonal)||missing(trend_length))
    stop("The per-period seasonal requests and the fixed trend filter length must be stated explicitly.",call.=FALSE)
  if(is.null(frequency)||is.null(start_period)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs explicit frequency and start_period.",call.=FALSE)
    if(is.null(frequency))frequency<-span[3L]
    if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
  }
  f<-as.integer(frequency);periods<-as.character(seq_len(f))
  if(!is.null(names(seasonal))) {
    if(length(seasonal)!=f || anyNA(names(seasonal)) || anyDuplicated(names(seasonal)) || !setequal(names(seasonal),periods))
      stop("Named seasonal requests must identify every calendar period exactly once.",call.=FALSE)
    seasonal<-seasonal[match(periods,names(seasonal))]
  }
  if(!is.list(variables))
    stop("The x13 model variables must be a list of variable specifications.",call.=FALSE)
  tsp_orig<-stats::tsp(series)
  y<-if(!is.null(dim(series)))as.vector(series) else series
  if(!is.numeric(y))y<-as.numeric(y)
  result<-.native_x13_fit(y,variables,log,f,start_period,as.character(seasonal),
    trend_length,sigma,lower_sigma,upper_sigma,msr_policy,mode,backcast,forecast)
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  start<-if(!is.null(tsp_orig))
    c(floor(tsp_orig[1L]),as.integer(round((tsp_orig[1L]-floor(tsp_orig[1L]))*f))+1L)
    else c(1L,as.integer(start_period))
  if(!is.null(result[["preadjustment"]])) { # past the fail-closed guards
    result$preadjustment$tables<-.x13_decorate_tables(result$preadjustment$tables,start,f)
    result$preadjustment$row_basis<-rep(list(seq_along(y)),length(result$preadjustment$tables))
    names(result$preadjustment$row_basis)<-names(result$preadjustment$tables)
    result$preadjustment$available_rows<-lapply(names(result$preadjustment$tables),
      function(name)result$preadjustment$row_basis[[name]][!is.na(result$preadjustment$tables[[name]])])
    names(result$preadjustment$available_rows)<-names(result$preadjustment$tables)
  }
  if(!is.null(result[["selection"]])) { # the same record the selection wrapper delivers
    result$selection<-.x11_decorate_selection_record(result$selection,ordinal,f,length(y))
  }
  if(identical(result$status,"available")) {
    for(name in c("d10","d11","d12","d13","c17"))
      result$results[[name]]<-stats::ts(result$results[[name]],start=start,frequency=f)
    result$finals$tables<-.x13_decorate_tables(result$finals$tables,start,f)
    for(name in c("finalSA","finalTrend","finalIrregular"))
      result[[name]]<-stats::ts(result[[name]],start=start,frequency=f)
  } else {
    result$records_only<-"the provenance echo is complete; the pinned stage refuses any run before that"
  }
  result$fitted<-identical(result$status,"available")
  structure(result,class="persephone_native_x13_fit")
}
