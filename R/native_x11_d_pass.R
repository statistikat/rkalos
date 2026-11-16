# Item 22b.4 S3: the complete fixed additive D pass, composed onto the
# caller's finished fixed additive C pass exactly like the pinned X11Kernel
# (D1 is the input series corrected by the C pass's final factors D1 = Y - C20).
# The C pass is a required argument: this wrapper takes the result object of
# native_x11_c_pass() (or .native_x11_c_pass()) and passes it through for the
# pinned coupling checks (status, setting echoes, round policy, reference-series
# equality, C20 coverage); the C pass is never silently recomputed here. The
# sigma policy must be stated explicitly ("legacy" or "accurate") and must
# match the C pass; the D pass itself performs no extreme-value analysis, so
# no sigma bounds are taken here. backcast/forecast are the pinned context
# horizons of a hand-extended calendar: the tables run over the extended
# series and the D9 MSRE diagnostic drops them onto the original calendar
# (the row span [backcast+1, n-forecast]), exactly like the pinned context.
# The MSRE filter selection, automatic Henderson choice, cross-validation and
# multiplicative modes are not computed here.
native_x11_d_pass <- function(series,policy,initial_filters,final_filters,trend_length,c,
    backcast=0,forecast=0,frequency=NULL,start_period=NULL) {
  if(missing(policy)||is.null(policy))
    stop("The sigma weight policy ('legacy' or 'accurate') must be stated explicitly.",call.=FALSE)
  if(missing(initial_filters)||missing(final_filters)||missing(trend_length))
    stop("Fixed initial and final seasonal filter requests and the trend filter length must be stated explicitly.",call.=FALSE)
  if(missing(c)||is.null(c)||!is.list(c)||is.null(c$tables)||
     is.null(c$round)||!identical(c$type,"x11_c_pass"))
    stop("The D pass requires a C-pass result object from native_x11_c_pass().",call.=FALSE)
  if(is.null(frequency)||is.null(start_period)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs explicit frequency and start_period.",call.=FALSE)
    if(is.null(frequency))frequency<-span[3L]
    if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
  }
  f<-as.integer(frequency);periods<-as.character(seq_len(f))
  align<-function(filters,name) {
    if(!is.null(names(filters))) {
      if(length(filters)!=f || anyNA(names(filters)) || anyDuplicated(names(filters)) || !setequal(names(filters),periods))
        stop("Named ",name," must identify every calendar period exactly once.",call.=FALSE)
      filters<-filters[match(periods,names(filters))]
    }
    filters
  }
  tsp_orig<-stats::tsp(series)
  if(!is.null(dim(series)))series<-as.vector(series)
  if(!is.numeric(series))series<-as.numeric(series)
  result<-.native_x11_d_pass(series,f,start_period,policy,
    align(initial_filters,"initial seasonal filters"),align(final_filters,"final seasonal filters"),
    trend_length,c,backcast,forecast)
  # Extended-calendar time domains: ts inputs keep their (extended) time
  # origin, plain vectors with an explicit start period are placed in year 1.
  # The D2/D4 tables live on the cropped span, every other table on the full
  # extended span; the dropped original calendar carries the MSRE diagnostic.
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  domain<-.native_time_domain(floor(ordinal/f),ordinal%%f+1L,f,length(series))
  domain$start_ordinal<-ordinal
  cropped<-result$d2drop+seq_len(length(series)-2L*result$d2drop)
  crop_ordinal<-ordinal+result$d2drop
  crop<-.native_time_domain(floor(crop_ordinal/f),crop_ordinal%%f+1L,f,length(cropped))
  crop$start_ordinal<-crop_ordinal
  original_ordinal<-ordinal+result$backcast
  original<-.native_time_domain(floor(original_ordinal/f),original_ordinal%%f+1L,f,
    length(series)-result$backcast-result$forecast)
  original$start_ordinal<-original_ordinal
  center<-c("D2","D4")
  for(name in names(result$tables))
    result$tables[[name]]<-stats::ts(result$tables[[name]],
      start=if(name %in% center) crop$start else domain$start,frequency=f)
  result$row_basis<-lapply(names(result$tables),function(name)
    if(name %in% center) cropped else seq_along(series))
  names(result$row_basis)<-names(result$tables)
  result$available_rows<-lapply(names(result$tables),function(name)
    result$row_basis[[name]][!is.na(result$tables[[name]])])
  names(result$available_rows)<-names(result$tables)
  # The pinned MSRE table is per internal phase index (0-based) of the
  # dropped original calendar at getPosition(backcast), exactly like the
  # pinned diagnostic; it is delivered as a table, not a calendar series.
  msr<-result$msr_table
  result$msr_table<-data.frame(phase=seq_along(msr$msr)-1L,
    mean_irregular_evolutions=msr$mean_irregular_evolutions,
    mean_seasonal_evolutions=msr$mean_seasonal_evolutions,
    counts=msr$counts,msr=msr$msr)
  result$global_msr<-msr$global_msr
  result$missing<-is.na(series)
  result$input_length<-length(series)
  result$input_domain<-domain
  result$cropped_domain<-crop
  result$original_domain<-original
  result$original_length<-length(series)-result$backcast-result$forecast
  result$c_coupling<-"checked by the pinned core against the supplied C-pass result"
  result$complete<-all(vapply(result$tables[c("D1","D5","D6","D7","D8","D10","D11","D11bis",
    "D12","D13")],function(x)!anyNA(x),logical(1)))
  structure(result,class="persephone_native_x11_d_pass")
}
