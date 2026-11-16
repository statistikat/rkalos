# Item 22b.4 S4: the pass-level additive X-11 selection flow of the pinned
# X11Kernel — the staged seasonal requests resolved through the qualified
# seasonal plan, the automatic Henderson choices at B7, C7, D7 and the D
# final step through the qualified comparison-ratio rule, and the global-MSR
# D-final seasonal decision through the qualified band rule with the pinned
# attempt loop, composed onto the fixed-pass arithmetic. The sigma policy
# must be stated explicitly ("legacy" or "accurate"): which rule becomes the
# product default remains an open Kalos maintainer decision, so this entry
# point takes no default policy. trend_length is the fixed Henderson length
# or 0 for the pinned automatic selection. Requests are one seasonal request
# per calendar period ("stable", "s3x1", "s3x3", "s3x5", "s3x9", "s3x15",
# "x11default" or "msr", as in the seasonal-plan boundary). Nothing new is
# decided here: a stage whose pinned choice the qualified rules refuse (a
# NaN comparison ratio or a NaN MSR table) yields the explicit status
# "unsupported:<stage>" with records only. Cross-validation, multiplicative
# and non-seasonal modes, the period-specific corrector and ARIMA-regression
# inputs are not computed here.
native_x11_selection_pass <- function(series,requests,policy,trend_length,
    msr_policy="jdplus",frequency=NULL,start_period=NULL,lower_sigma=1.5,
    upper_sigma=2.5,backcast=0,forecast=0,mode="additive") {
  if(missing(policy)||is.null(policy))
    stop("The sigma weight policy ('legacy' or 'accurate') must be stated explicitly.",call.=FALSE)
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
  names(result$requests)<-periods
  for(field in c("initial_filters","final_filters"))names(result[[field]])<-periods
  for(field in c("initial_plan","final_plan"))names(result[[field]]$base_filters)<-periods
  if(identical(result$status,"available")) {
    # The selected-state passes keep the fixed-pass calendar decorations; the
    # D9 MSRE table stays the pinned per-phase diagnostic, not a series.
    decorate<-function(pass,center) {
      drop<-if(identical(pass$type,"x11_b_pass"))pass$b2drop else
        if(identical(pass$type,"x11_c_pass"))pass$c2drop else pass$d2drop
      domain<-.native_time_domain(floor(ordinal/f),ordinal%%f+1L,f,length(series))
      domain$start_ordinal<-ordinal
      cropped<-drop+seq_len(length(series)-2L*drop)
      crop_ordinal<-ordinal+drop
      crop<-.native_time_domain(floor(crop_ordinal/f),crop_ordinal%%f+1L,f,length(cropped))
      crop$start_ordinal<-crop_ordinal
      for(name in names(pass$tables))
        pass$tables[[name]]<-stats::ts(pass$tables[[name]],
          start=if(name %in% center) crop$start else domain$start,frequency=f)
      pass$row_basis<-lapply(names(pass$tables),function(name)
        if(name %in% center) cropped else seq_along(series))
      names(pass$row_basis)<-names(pass$tables)
      pass$available_rows<-lapply(names(pass$tables),function(name)
        pass$row_basis[[name]][!is.na(pass$tables[[name]])])
      names(pass$available_rows)<-names(pass$tables)
      pass$input_domain<-domain
      pass$cropped_domain<-crop
      if(!is.null(pass[["rounds"]])) { # the B first round detects on the cropped span
        pass$rounds[[1]]$detection<-stats::ts(pass$rounds[[1]]$detection,start=crop$start,frequency=f)
        for(i in c(2L,3L))
          pass$rounds[[i]]$detection<-stats::ts(pass$rounds[[i]]$detection,start=domain$start,frequency=f)
        for(i in seq_along(pass$rounds))
          pass$rounds[[i]]$analysis_period<-as.integer((ordinal+pass$rounds[[i]]$analysis_position-1L)%%f)+1L
      }
      if(!is.null(pass[["round"]])) { # the C round detects on the full span
        pass$round$detection<-stats::ts(pass$round$detection,start=domain$start,frequency=f)
        pass$round$analysis_period<-as.integer((ordinal+pass$round$analysis_position-1L)%%f)+1L
      }
      if(!is.null(pass[["msr_table"]])) {
        msr<-pass[["msr_table"]]
        pass$msr_table<-data.frame(phase=seq_along(msr$msr)-1L,
          mean_irregular_evolutions=msr$mean_irregular_evolutions,
          mean_seasonal_evolutions=msr$mean_seasonal_evolutions,
          counts=msr$counts,msr=msr$msr)
        pass$global_msr<-msr$global_msr
        original_ordinal<-ordinal+pass$backcast
        original<-.native_time_domain(floor(original_ordinal/f),original_ordinal%%f+1L,f,
          length(series)-pass$backcast-pass$forecast)
        original$start_ordinal<-original_ordinal
        pass$original_domain<-original
        pass$original_length<-length(series)-pass$backcast-pass$forecast
      }
      pass
    }
    result$b<-decorate(result$b,c("B2","B3","B4","B4a","B4anorm","B4d","B4g"))
    result$c<-decorate(result$c,c("C2","C4","C5a"))
    result$d<-decorate(result$d,c("D2","D4"))
    # The selected run is the composition of these three states; chaining the
    # records through the fixed boundaries reproduces them bit-exactly.
    result$composition<-"selected-state B, C and D passes (records chain into native_x11_c_pass/native_x11_d_pass)"
  } else {
    result$records_only<-"the stage records and plans are complete; no tables are delivered for an unsupported run"
  }
  result$msr$attempts<-data.frame(count=result$msr$attempts$count,
    years_removed=result$msr$attempts$years_removed,global_msr=result$msr$attempts$global_msr,
    band=result$msr$attempts$band,status=result$msr$attempts$status)
  result$selected<-identical(result$status,"available")
  structure(result,class="persephone_native_x11_selection_pass")
}
