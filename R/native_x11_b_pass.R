# Item 22b.4 S1: the complete fixed additive B pass on a monthly or quarterly
# series with fixed per-phase seasonal requests, a fixed Henderson trend filter
# and three extreme-value correction rounds. The sigma policy must be stated
# explicitly ("legacy" or "accurate"): which rule becomes the product default
# remains an open Kalos maintainer decision, so neither is promoted here and
# this entry point takes no default policy. The C/D passes, automatic filter
# choice, cross-validation and multiplicative modes are not computed here.
native_x11_b_pass <- function(series,policy,initial_filters,final_filters,trend_length,
    frequency=NULL,start_period=NULL,lower_sigma=1.5,upper_sigma=2.5) {
  if(missing(policy)||is.null(policy))
    stop("The sigma weight policy ('legacy' or 'accurate') must be stated explicitly.",call.=FALSE)
  if(missing(initial_filters)||missing(final_filters)||missing(trend_length))
    stop("Fixed initial and final seasonal filter requests and the trend filter length must be stated explicitly.",call.=FALSE)
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
  result<-.native_x11_b_pass(series,f,start_period,policy,
    align(initial_filters,"initial seasonal filters"),align(final_filters,"final seasonal filters"),
    trend_length,lower_sigma,upper_sigma)
  # Original-calendar time domains: ts inputs keep their time origin, plain
  # vectors with an explicit start period are placed in year 1. The ordinal
  # follows the native_regular_domain convention year*f+(period-1).
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  domain<-.native_time_domain(floor(ordinal/f),ordinal%%f+1L,f,length(series))
  domain$start_ordinal<-ordinal
  cropped<-result$b2drop+seq_len(length(series)-2L*result$b2drop)
  crop_ordinal<-ordinal+result$b2drop
  crop<-.native_time_domain(floor(crop_ordinal/f),crop_ordinal%%f+1L,f,length(cropped))
  crop$start_ordinal<-crop_ordinal
  center<-c("B2","B3","B4","B4a","B4anorm","B4d","B4g")
  for(name in names(result$tables))
    result$tables[[name]]<-stats::ts(result$tables[[name]],
      start=if(name %in% center) crop$start else domain$start,frequency=f)
  result$rounds[[1]]$detection<-stats::ts(result$rounds[[1]]$detection,start=crop$start,frequency=f)
  for(i in c(2L,3L))
    result$rounds[[i]]$detection<-stats::ts(result$rounds[[i]]$detection,start=domain$start,frequency=f)
  for(i in seq_along(result$rounds))
    result$rounds[[i]]$analysis_period<-as.integer((ordinal+result$rounds[[i]]$analysis_position-1L)%%f)+1L
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
  result$complete<-all(vapply(result$tables,function(x)!anyNA(x),logical(1)))
  structure(result,class="persephone_native_x11_b_pass")
}
