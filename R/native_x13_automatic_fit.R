# Item 21-H M4: the automatic-spec X-13 fit — the qualified stage-Z pipeline
# (docs/contracts/x13-automatic-pipeline-contract.md) at the R boundary. The
# request arrives as the decoder-equivalent module state in the rsa-template
# shape (transform / automodel / outliers / ami / regression; the X-11 block
# keeps the flat fit conventions). NOTHING is promoted to a default except the
# maintainer sigma-rule default: every module scalar is required request data
# (the core applies no defaults either), a missing or malformed entry is
# refused with an explicit error before anything runs, and the raw entry
# answers the pinned fail-closed statuses verbatim. The X-11 block is the
# decoder automatic surface: one "msr" request per calendar period and the
# automatic Henderson length 0 are the pinned G4 facts, not user options. The
# X-11 horizons (item 21-H.4) are R-qualified over the horizon envelope: month
# counts, a negative count meaning |h| annual frequencies (the core normalizes);
# the initialStep route extends the downstream composition (out, d10..c17) on the
# calendar shifted by the backcast, while the decision chain, the landing model
# and the preadjustment inputs stay horizon-invariant
# (capabilities()$automatic_fit$horizon_envelope). Sigma: NULL (or omitted) receives the
# maintainer default "legacy" — the decision of 2026-11-17, applied here at
# the interface while the core keeps its explicitness contract and the raw
# entry still answers the pinned "unsupported:x13 fit requires an explicit
# sigma policy" for a missing policy; the $sigma_provenance echo reports
# "legacy (default)" against an explicitly stated policy ("legacy (explicit)"
# / "accurate (explicit)"), the schema-v2 restricted_fit semantics extended to
# this surface. Calendar, easter-detection, automatic-trading-days and
# user-variable detection are deferred to 21-I: their flags are exposed and
# the core answers the explicit unsupported status. There is NO silent
# fallback anywhere: outside the scope envelope the result carries the pinned
# unsupported:/failed: status with the complete provenance echo, never a
# restricted fit. capabilities()$automatic_fit reports the envelope.
.auto_required_entry <- function(module, name, kind, value) {
  label<-paste0("automatic spec ",module,"$",name)
  if(is.null(value))
    stop(label," is required request data (the automatic spec promotes no defaults).",call.=FALSE)
  if(length(value)!=1L || is.na(value)) stop(label," must be a single non-NA value.",call.=FALSE)
  if(identical(kind,"num")) {
    if(!is.numeric(value)) stop(label," must be a number.",call.=FALSE)
    return(as.numeric(value))
  }
  if(identical(kind,"int")) {
    if(!is.numeric(value) || as.numeric(value)!=floor(as.numeric(value)))
      stop(label," must be a whole number.",call.=FALSE)
    return(as.integer(value))
  }
  if(identical(kind,"lgl")) {
    if(!is.logical(value)) stop(label," must be TRUE or FALSE.",call.=FALSE)
    return(value)
  }
  if(identical(kind,"str")) {
    if(!is.character(value)) stop(label," must be a character string.",call.=FALSE)
    return(value)
  }
  stop("internal: unknown entry kind ",kind,call.=FALSE)
}
.auto_module <- function(spec, module, shape) {
  if(!is.null(spec) && !is.list(spec))
    stop("automatic spec ",module," must be a named list of module options.",call.=FALSE)
  if(is.null(spec))spec<-list()
  nms<-names(spec)
  if(length(spec)>0L && (is.null(nms) || any(nms=="")))
    stop("automatic spec ",module," entries must be named (",
      paste(names(shape),collapse=", "),").",call.=FALSE)
  unknown<-setdiff(nms,names(shape))
  if(length(unknown))
    stop("automatic spec ",module," has no option(s): ",paste(unknown,collapse=", "),".",call.=FALSE)
  out<-lapply(names(shape),function(name)
    .auto_required_entry(module,name,shape[[name]],spec[[name]]))
  names(out)<-names(shape)
  out
}
native_x13_fit_automatic <- function(series,transform,automodel,outliers,ami,regression,
    sigma=NULL,lower_sigma=1.5,upper_sigma=2.5,msr_policy="jdplus",mode="additive",
    backcast=0,forecast=0,frequency=NULL,start_period=NULL,first_year=NULL,
    calendar=FALSE,easter=FALSE,trading_days=FALSE,user_variables=FALSE) {
  # The rsa-template module blocks; every listed entry is required (the
  # decoder hands the module state as request data — nothing is promoted).
  tr<-.auto_module(transform,"transform",list(aicc="num"))
  am<-.auto_module(automodel,"automodel",list(max_d="int",max_bd="int",ub1="num",ub2="num",
    cancel="num",precision="num",models_count="int",max_p="int",max_q="int",
    max_bp="int",max_bq="int",balanced="lgl",mixed="lgl",accept_white_noise="lgl",eps="num"))
  # types and span are validated below as a set and a window (not scalars);
  # everything else in the module is scalar request data.
  ou_spec<-outliers[setdiff(names(outliers),c('types','span'))]
  ou<-.auto_module(ou_spec,"outliers",list(used="lgl",max_round="int",
    max_outliers="int",tcrate="num"))
  if(!is.null(outliers[["span"]])) {
    sp<-outliers$span
    if(!is.numeric(sp)||length(sp)!=2L||any(is.na(sp))||any(sp!=floor(sp)))
      stop("automatic spec outliers$span must be NULL or a length-2 integer vector.",call.=FALSE)
  } else sp<-NULL
  types<-outliers[["types"]]
  if(!ou$used && is.null(types)) types<-character(0)
  if(!is.character(types))
    stop("automatic spec outliers$types must be a character vector of AO, LS, TC, SO.",call.=FALSE)
  mi<-.auto_module(ami,"ami",list(precision="num",mean_precision="num",va="num",
    reduce_va="num",ljung_box_limit="num",check_mu="lgl",accept_airline="lgl",ur_limit="num"))
  rg<-.auto_module(regression,"regression",list(mean_test_initial="num",
    mean_test_final="num",tsig="num"))
  sigma_defaulted<-is.null(sigma)
  if(sigma_defaulted) sigma<-"legacy" # the maintainer 2026-11-17 product default
  if(!is.character(sigma)||length(sigma)!=1L||is.na(sigma))
    stop('The sigma weight policy must be "legacy" or "accurate" (or NULL for the maintainer default).',call.=FALSE)
  provenance<-if(sigma_defaulted) "legacy (default)" else paste0(sigma," (explicit)")
  if(is.null(frequency)||is.null(start_period)||is.null(first_year)) {
    span<-stats::tsp(series)
    if(is.null(span)&&is.null(frequency))
      stop("A plain numeric series needs explicit frequency, start_period and first_year.",call.=FALSE)
    if(!is.null(span)) {
      if(is.null(frequency))frequency<-span[3L]
      if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
      if(is.null(first_year))first_year<-floor(span[1L])
    } else if(is.null(start_period)||is.null(first_year))
      stop("A plain numeric series needs explicit start_period and first_year.",call.=FALSE)
  }
  f<-as.integer(frequency)
  tsp_orig<-stats::tsp(series)
  y<-if(!is.null(dim(series)))as.vector(series) else series
  if(!is.numeric(y))y<-as.numeric(y)
  # The decoder automatic X-11 surface (contract G4): Msr for every calendar
  # period, automatic Henderson; the refusal for a non-monthly frequency is the
  # core's, answered before any of this state is interpreted.
  seasonal<-rep("msr",f)
  result<-.native_x13_fit_automatic(y,f,as.integer(start_period),as.integer(first_year),
    tr$aicc,am$max_d,am$max_bd,am$ub1,am$ub2,am$cancel,am$precision,am$models_count,
    am$max_p,am$max_q,am$max_bp,am$max_bq,am$balanced,am$mixed,am$accept_white_noise,am$eps,
    ou$used,types,sp,ou$max_round,ou$max_outliers,
    ou$tcrate,mi$precision,mi$mean_precision,mi$va,mi$reduce_va,mi$ljung_box_limit,
    mi$check_mu,mi$accept_airline,mi$ur_limit,rg$mean_test_initial,rg$mean_test_final,
    rg$tsig,calendar,easter,trading_days,user_variables,seasonal,sigma,lower_sigma,
    upper_sigma,msr_policy,mode,as.integer(backcast),as.integer(forecast))
  result$sigma_provenance<-provenance
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  start<-if(!is.null(tsp_orig))
    c(floor(tsp_orig[1L]),as.integer(round((tsp_orig[1L]-floor(tsp_orig[1L]))*f))+1L)
    else c(as.integer(first_year),as.integer(start_period))
  # The horizon calendars (item 21-H.4): the core normalizes a negative count to
  # |h| annual frequencies (kAnnualFrequency at the initialStep route); the
  # extended composition tables (out, d10..c17) carry bc months before the input
  # calendar start, the n-length tables keep the start above unchanged.
  bc_eff<-as.integer(backcast);if(bc_eff<0L)bc_eff<--bc_eff*12L
  year_ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else
    as.integer(first_year)*f+as.integer(start_period)-1L
  ext_start<-if(bc_eff>0L){o<-year_ordinal-bc_eff
    c(floor(o/f),as.integer(o%%f)+1L)} else start
  if(!is.null(result[["preadjustment"]])) { # past the fail-closed guards
    basis<-list();avail<-list()
    for(name in names(result$preadjustment$tables)) {
      tab<-result$preadjustment$tables[[name]]
      result$preadjustment$tables[[name]]<-stats::ts(tab,
        start=if(length(tab)==length(y))start else ext_start,frequency=f)
      basis[[name]]<-seq_along(tab)
      avail[[name]]<-seq_along(tab)[!is.na(tab)]
    }
    result$preadjustment$row_basis<-basis
    result$preadjustment$available_rows<-avail
  }
  if(!is.null(result[["selection"]])) { # the same record the selection wrapper delivers
    result$selection<-.x11_decorate_selection_record(result$selection,ordinal,f,length(y))
  }
  if(identical(result$status,"available")) {
    for(name in c("d10","d11","d12","d13","c17"))
      result$results[[name]]<-stats::ts(result$results[[name]],
        start=if(length(result$results[[name]])==length(y))start else ext_start,frequency=f)
    result$finals$tables<-.x13_decorate_tables(result$finals$tables,start,f)
    for(name in c("finalSA","finalTrend","finalIrregular"))
      result[[name]]<-stats::ts(result[[name]],start=start,frequency=f)
    result$series<-if(!is.null(dim(series)))series else stats::ts(y,start=start,frequency=f)
  } else {
    result$records_only<-if(grepl("^failed:",result$status))
      "the machine threw inside the pinned calc body; the decision chain kept "
    else "the provenance echo is complete; the pinned fail-closed matrix refuses the request"
  }
  result$fitted<-identical(result$status,"available")
  structure(result,class="persephone_native_x13_automatic_fit")
}
