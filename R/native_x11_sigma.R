# Item 22b.3a: moving five-year dispersion windows and two-pass sigma weights.
# The sigma weight policy defaults to "legacy": the Kalos maintainer decision
# of 2026-11-17 promotes legacy (bit-exact against the pinned
# DefaultExtremeValuesCorrector / JDemetra behavior) to the product default.
# An explicitly stated policy ("legacy" or "accurate") always wins, and the
# core keeps its explicitness contract: the raw entry point still refuses a
# missing policy; applying the default is this interface layer's role.
native_x11_sigma_weights <- function(series,policy="legacy",frequency=NULL,start_period=NULL,
    backcasts=0L,forecasts=0L,exclude_projections=FALSE,
    lower_sigma=1.5,upper_sigma=2.5,multiplicative=FALSE) {
  if(is.null(policy))policy<-"legacy"
  if(!is.null(dim(series)))series<-as.vector(series)
  if(!is.numeric(series))series<-as.numeric(series)
  if(is.null(frequency)||is.null(start_period)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs explicit frequency and start_period.",call.=FALSE)
    if(is.null(frequency))frequency<-span[3L]
    if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
  }
  result<-.native_x11_sigma_weights(series,frequency,start_period,policy,backcasts,forecasts,
    exclude_projections,lower_sigma,upper_sigma,multiplicative)
  result$missing<-is.na(series)
  for(field in c("windows","first_pass_windows")) {
    window<-result[[field]]
    result[[field]]<-data.frame(dispersion=window$dispersions,start=window$starts+1,
      length=window$lengths,usable=window$usable,kind=window$kinds,check.names=FALSE,
      stringsAsFactors=FALSE)
  }
  names(result$weights)<-names(result$first_pass_weights)<-
    if(is.null(names(series)))as.character(seq_along(series)) else names(series)
  result$input_length<-length(series)
  result$replacement_policy<-"not_computed_here"
  structure(result,class="persephone_native_x11_sigma")
}
