# Item 22b.4 S2: the complete fixed additive C pass, composed onto the
# caller's finished fixed additive B pass exactly like the pinned X11Kernel
# (C1 is the input series corrected by the B pass's final factors C1 = Y - B20).
# The B pass is a required argument: this wrapper takes the result object of
# native_x11_b_pass() (or .native_x11_b_pass()) and passes it through for the
# pinned coupling checks (status, setting echoes, round policies, series
# equality); the B pass is never silently recomputed here. The sigma policy
# must be stated explicitly ("legacy" or "accurate") and must match the B
# pass: which rule becomes the product default remains an open Kalos
# maintainer decision, so neither is promoted and this entry point takes no
# default policy. The D pass, automatic filter choice, cross-validation and
# multiplicative modes are not computed here.
native_x11_c_pass <- function(series,policy,initial_filters,final_filters,trend_length,b,
    frequency=NULL,start_period=NULL,lower_sigma=1.5,upper_sigma=2.5) {
  if(missing(policy)||is.null(policy))
    stop("The sigma weight policy ('legacy' or 'accurate') must be stated explicitly.",call.=FALSE)
  if(missing(initial_filters)||missing(final_filters)||missing(trend_length))
    stop("Fixed initial and final seasonal filter requests and the trend filter length must be stated explicitly.",call.=FALSE)
  if(missing(b)||is.null(b)||!is.list(b)||is.null(b$tables)||
     is.null(b$rounds)||!identical(b$type,"x11_b_pass"))
    stop("The C pass requires a B-pass result object from native_x11_b_pass().",call.=FALSE)
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
  result<-.native_x11_c_pass(series,f,start_period,policy,
    align(initial_filters,"initial seasonal filters"),align(final_filters,"final seasonal filters"),
    trend_length,b,lower_sigma,upper_sigma)
  # Original-calendar time domains, mirroring the B-pass wrapper: the C2/C4/
  # C5a tables live on the cropped span, every other table on the full span.
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  domain<-.native_time_domain(floor(ordinal/f),ordinal%%f+1L,f,length(series))
  domain$start_ordinal<-ordinal
  cropped<-result$c2drop+seq_len(length(series)-2L*result$c2drop)
  crop_ordinal<-ordinal+result$c2drop
  crop<-.native_time_domain(floor(crop_ordinal/f),crop_ordinal%%f+1L,f,length(cropped))
  crop$start_ordinal<-crop_ordinal
  center<-c("C2","C4","C5a")
  for(name in names(result$tables))
    result$tables[[name]]<-stats::ts(result$tables[[name]],
      start=if(name %in% center) crop$start else domain$start,frequency=f)
  result$round$detection<-stats::ts(result$round$detection,start=domain$start,frequency=f)
  result$round$analysis_period<-as.integer((ordinal+result$round$analysis_position-1L)%%f)+1L
  result$row_basis<-lapply(names(result$tables),function(name)
    if(name %in% center) cropped else seq_along(series))
  names(result$row_basis)<-names(result$tables)
  result$available_rows<-lapply(names(result$tables),function(name)
    result$row_basis[[name]][!is.na(result$tables[[name]])])
  names(result$available_rows)<-names(result$tables)
  result$missing<-is.na(series)
  result$input_length<-length(series)
  result$input_domain<-domain
  result$cropped_domain<-crop
  result$b_coupling<-"checked by the pinned core against the supplied B-pass result"
  result$complete<-all(vapply(result$tables,function(x)!anyNA(x),logical(1)))
  structure(result,class="persephone_native_x11_c_pass")
}
