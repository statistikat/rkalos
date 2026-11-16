# Item 22b.4 S4: the pass-level additive X-11 selection flow of the pinned
# X11Kernel — the staged seasonal requests resolved through the qualified
# seasonal plan, the automatic Henderson choices at B7, C7, D7 and the D
# final step through the qualified comparison-ratio rule, and the global-MSR
# D-final seasonal decision through the qualified band rule with the pinned
# attempt loop, composed onto the fixed-pass arithmetic. The sigma policy
# defaults to "legacy": the Kalos maintainer decision of 2026-11-17 promotes
# legacy (bit-exact against the pinned DefaultExtremeValuesCorrector /
# JDemetra behavior) to the product default, applied here at the interface;
# an explicitly stated policy always wins and the raw entry point keeps the
# core's explicitness contract. trend_length is the fixed Henderson length
# or 0 for the pinned automatic selection. Requests are one seasonal request
# per calendar period ("stable", "s3x1", "s3x3", "s3x5", "s3x9", "s3x15",
# "x11default" or "msr", as in the seasonal-plan boundary). Nothing new is
# decided here: a stage whose pinned choice the qualified rules refuse (a
# NaN comparison ratio or a NaN MSR table) yields the explicit status
# "unsupported:<stage>" with records only. Cross-validation, multiplicative
# and non-seasonal modes, the period-specific corrector and ARIMA-regression
# inputs are not computed here.
native_x11_selection_pass <- function(series,requests,policy="legacy",trend_length,
    msr_policy="jdplus",frequency=NULL,start_period=NULL,lower_sigma=1.5,
    upper_sigma=2.5,backcast=0,forecast=0,mode="additive") {
  if(is.null(policy))policy<-"legacy"
  if(missing(requests)||missing(trend_length))
    stop("The per-period seasonal requests and the trend filter length (0 for automatic) must be stated explicitly.",call.=FALSE)
  if(is.null(frequency)||is.null(start_period)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs explicit frequency and start_period.",call.=FALSE)
    if(is.null(frequency))frequency<-span[3L]
    if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
  }
  f<-as.integer(frequency);periods<-as.character(seq_len(f))
  if(!is.null(names(requests))) {
    if(length(requests)!=f || anyNA(names(requests)) || anyDuplicated(names(requests)) || !setequal(names(requests),periods))
      stop("Named requests must identify every calendar period exactly once.",call.=FALSE)
    requests<-requests[match(periods,names(requests))]
  }
  tsp_orig<-stats::tsp(series)
  if(!is.null(dim(series)))series<-as.vector(series)
  if(!is.numeric(series))series<-as.numeric(series)
  result<-.native_x11_selection_pass(series,f,start_period,as.character(requests),
    trend_length,policy,msr_policy,lower_sigma,upper_sigma,backcast,forecast,mode)
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  # The selected run is the composition of the three decorated states;
  # chaining the records through the fixed boundaries reproduces them
  # bit-exactly. An unsupported run keeps only its stage records and plans.
  result<-.x11_decorate_selection_record(result,ordinal,f,length(series),
    if(identical(result$status,"available"))
      list(composition="selected-state B, C and D passes (records chain into native_x11_c_pass/native_x11_d_pass)")
    else
      list(records_only="the stage records and plans are complete; no tables are delivered for an unsupported run"))
  structure(result,class="persephone_native_x11_selection_pass")
}
