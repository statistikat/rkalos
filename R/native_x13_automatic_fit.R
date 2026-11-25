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
  if(!is.null(dim(value)) || length(value)!=1L || is.na(value)) stop(label," must be a single non-NA value.",call.=FALSE)
  if(identical(kind,"num")) {
    if(!is.numeric(value) || !is.finite(value)) stop(label," must be a number.",call.=FALSE)
    return(as.numeric(value))
  }
  if(identical(kind,"int")) {
    if(!is.numeric(value) || !is.finite(value) || abs(value) > .Machine$integer.max || value!=floor(value))
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
    if(!is.numeric(sp) || !is.null(dim(sp)) || length(sp)!=2L || any(!is.finite(sp)) ||
       any(sp!=floor(sp)) || sp[1L]<1 || sp[2L]<sp[1L] || sp[2L]>length(series))
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
  if (!is.numeric(series) || !is.null(dim(series)) || !length(series))
    stop("series must be a nonempty univariate numeric vector or ts.", call.=FALSE)
  calendar_int <- function(value, name) .auto_required_entry("calendar", name, "int", value)
  if (stats::is.ts(series)) {
    series_frequency <- stats::frequency(series)
    ordinal <- round(stats::tsp(series)[1L] * series_frequency)
    if (abs(stats::tsp(series)[1L] * series_frequency - ordinal) > 1e-6)
      stop("series must start on an exact calendar period.", call.=FALSE)
    calendar_value <- function(value, expected, name) {
      if (!is.null(value) && calendar_int(value, name) != expected)
        stop(name, " disagrees with the ts calendar.", call.=FALSE)
      calendar_int(expected, name)
    }
    frequency <- calendar_value(frequency, series_frequency, "frequency")
    start_period <- calendar_value(start_period, ordinal %% frequency + 1, "start_period")
    first_year <- calendar_value(first_year, floor(ordinal / frequency), "first_year")
  } else {
    frequency <- calendar_int(frequency, "frequency")
    start_period <- calendar_int(start_period, "start_period")
    first_year <- calendar_int(first_year, "first_year")
  }
  f <- frequency
  start <- c(first_year, start_period)
  backcast <- calendar_int(backcast, "backcast")
  forecast <- calendar_int(forecast, "forecast")
  y <- as.numeric(series)
  # The raw selection parser requires a nonempty vector. The monthly decoder
  # always has twelve filters; the core refuses other frequencies first.
  seasonal <- rep("msr", 12L)
  result<-.native_x13_fit_automatic(y,f,start_period,first_year,
    tr$aicc,am$max_d,am$max_bd,am$ub1,am$ub2,am$cancel,am$precision,am$models_count,
    am$max_p,am$max_q,am$max_bp,am$max_bq,am$balanced,am$mixed,am$accept_white_noise,am$eps,
    ou$used,types,sp,ou$max_round,ou$max_outliers,
    ou$tcrate,mi$precision,mi$mean_precision,mi$va,mi$reduce_va,mi$ljung_box_limit,
    mi$check_mu,mi$accept_airline,mi$ur_limit,rg$mean_test_initial,rg$mean_test_final,
    rg$tsig,calendar,easter,trading_days,user_variables,seasonal,sigma,lower_sigma,
    upper_sigma,msr_policy,mode,as.integer(backcast),as.integer(forecast))
  result$sigma_provenance<-provenance
  result <- date_estimated_x13_result(result, series, start, f)
  if (length(result$residuals$tsres)) { # the estimated route's dating convention (native_x13_estimated_fit.R)
    first <- as.double(start[1L]) * f + start[2L] - 1L + result$model$d + f * result$model$bd
    result$residuals$tsres <- stats::ts(result$residuals$tsres,
      start = c(first %/% f, first %% f + 1), frequency = f)
  }
  structure(result, class = "persephone_native_x13_automatic_fit")
}

date_estimated_x13_result <- function(result, series, start, f) {
  y <- as.numeric(series)
  ordinal <- as.double(start[1L]) * f + start[2L] - 1L
  extended_ordinal <- ordinal - result$effective_backcast
  extended_length <- length(y) + result$effective_backcast + result$effective_forecast
  ext_start <- if (result$effective_backcast > 0L)
    c(floor(extended_ordinal / f), extended_ordinal %% f + 1L) else start
  observed_ts <- function(x) {
    out <- stats::ts(x, start = start, frequency = f)
    if (stats::is.ts(series)) attr(out, "tsp") <- stats::tsp(series)
    out
  }
  forecast_ordinal <- ordinal + length(y)
  forecast_start <- c(floor(forecast_ordinal / f), forecast_ordinal %% f + 1)
  if (length(result$series_forecast_stdev))
    result$series_forecast_stdev <- stats::ts(result$series_forecast_stdev, start = forecast_start, frequency = f)
  if(!is.null(result[["preadjustment"]])) { # past the fail-closed guards
    basis<-list();avail<-list()
    for(name in names(result$preadjustment$tables)) {
      tab<-result$preadjustment$tables[[name]]
      table_start <- switch(name, a1a = forecast_start, a1b = ext_start,
                            if (length(tab) == length(y)) start else ext_start)
      result$preadjustment$tables[[name]] <- if (!name %in% c("a1a", "a1b") && length(tab) == length(y))
        observed_ts(tab) else stats::ts(tab, start = table_start, frequency = f)
      basis[[name]]<-seq_along(tab)
      avail[[name]]<-seq_along(tab)[!is.na(tab)]
    }
    result$preadjustment$row_basis<-basis
    result$preadjustment$available_rows<-avail
  }
  if(!is.null(result[["selection"]])) { # the same record the selection wrapper delivers
    result$selection<-.x11_decorate_selection_record(result$selection,extended_ordinal,f,extended_length)
  }
  if(identical(result$status,"available")) {
    for(name in c("d10","d11","d12","d13","c17"))
      result$results[[name]] <- if (length(result$results[[name]]) == length(y))
        observed_ts(result$results[[name]]) else stats::ts(result$results[[name]], start = ext_start, frequency = f)
    for (name in names(result$finals$tables)) {
      table_start <- if (endsWith(name, "a")) forecast_start else if (endsWith(name, "b")) ext_start else start
      tab <- result$finals$tables[[name]]
      result$finals$tables[[name]] <- if (endsWith(name, "a") || endsWith(name, "b"))
        stats::ts(tab, start = table_start, frequency = f) else observed_ts(tab)
    }
    for(name in c("finalSA","finalTrend","finalIrregular"))
      result[[name]] <- observed_ts(result[[name]])
    # Public residual trading-day F-tests (M7 cutover): computed from this
    # fit's own finalSA/finalIrregular on the shared observed domain start
    # (the automatic and estimated routes both land here on their success
    # path; refusals keep the provenance echo and never carry the surface) —
    # on the reference working scale: the fit's preadjustment log flag (the
    # delivered envelope's preadjustment$log) makes the attach log() the
    # components for log fits (the automatic transform decision lands on log
    # for multiplicative fits) and leaves additive fits on the delivered
    # levels (no-op).
    result$diagnostics <- .attach_td_tests(result$diagnostics, result$finalSA, result$finalIrregular, start,
                                           isTRUE(result$preadjustment$log))
    result$series <- observed_ts(y)
  } else {
    result$records_only<-if(grepl("^failed:",result$status))
      "the machine threw inside the pinned calc body; the decision chain kept "
    else "the provenance echo is complete; the pinned fail-closed matrix refuses the request"
  }
  result$fitted<-identical(result$status,"available")
  result
}
