# Item 22c M3c: the two deterministic-effect allocation seams of the pinned
# X13Kernel — the fixed-coefficient preadjust (initialStep) and the finals
# restore (finals). The kernel gives no fixed-coefficient entry point, so the
# qualified contract is the exact seam the capture proves: given the complete
# input series and the pinned model-description variable order with each
# variable's pinned regression columns and coefficient values as INPUTS, the
# staging tables are pure algebra. Coefficient and ARIMA estimation, outlier
# detection and the trading-day tests stay outside; the fit of these seams is
# native_x13_fit(). Nothing is decided here: variables arrive in the pinned
# model order and the pinned seam judges their shapes verbatim.
# A variable is a list(class="outlier"|"days_related"|"moving_holiday"|"user",
# effect="undefined"|"series"|"trend"|"seasonal"|"seasonally_adjusted"|
# "irregular"|"calendar_effect", coefficient=<one number per column>,
# columns=<list of numeric vectors on the input calendar>).
native_x13_preadjust <- function(series,variables=list(),log=FALSE,
    frequency=NULL,start_period=NULL) {
  if(!is.list(variables))
    stop("The x13 model variables must be a list of variable specifications.",call.=FALSE)
  if(is.null(frequency)||is.null(start_period)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs explicit frequency and start_period.",call.=FALSE)
    if(is.null(frequency))frequency<-span[3L]
    if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
  }
  f<-as.integer(frequency)
  y<-if(!is.null(dim(series)))as.vector(series) else series
  if(!is.numeric(y))y<-as.numeric(y)
  result<-.native_x13_preadjust(y,log,variables)
  span<-stats::tsp(series)
  start<-if(!is.null(span))
    c(floor(span[1L]),as.integer(round((span[1L]-floor(span[1L]))*f))+1L)
    else c(1L,as.integer(start_period))
  if(identical(result$status,"available")) {
    result$tables<-.x13_decorate_tables(result$tables,start,f)
    result$row_basis<-rep(list(seq_along(y)),length(result$tables))
    names(result$row_basis)<-names(result$tables)
    result$available_rows<-lapply(names(result$tables),function(name)
      result$row_basis[[name]][!is.na(result$tables[[name]])])
    names(result$available_rows)<-names(result$tables)
  }
  structure(result,class="persephone_native_x13_preadjustment")
}
# The finals restore consumes a finished preadjustment result object (the
# pinned seam enforces the coupling on the rebuilt record — the boundary never
# recomputes it) plus the X11Results tables exactly as finals consumes them:
# D10-D13 and C17 on the same calendar. mode is the EFFECTIVE decomposition
# mode of the X-11 run that produced the tables (the pinned X13Kernel.finals
# branches on it); the default "additive" matches the pinned non-log forcing.
native_x13_restore <- function(preadjustment,d10,d11,d12,d13,c17,mode="additive") {
  if(!is.list(preadjustment)||is.null(preadjustment$tables)||
     is.null(preadjustment$status)||!identical(preadjustment$type,"x13_preadjustment"))
    stop("The restore requires a preadjustment result object from native_x13_preadjust().",call.=FALSE)
  span<-stats::tsp(preadjustment$tables$a1)
  if(is.null(span))
    stop("The preadjustment tables carry no calendar; build it from a ts series.",call.=FALSE)
  f<-span[3L]
  start<-c(floor(span[1L]),as.integer(round((span[1L]-floor(span[1L]))*f))+1L)
  tables<-function(x){x<-if(!is.null(dim(x)))as.vector(x) else x;if(!is.numeric(x))x<-as.numeric(x);x}
  result<-.native_x13_restore(preadjustment,tables(d10),tables(d11),tables(d12),
    tables(d13),tables(c17),mode)
  if(identical(result$status,"available")) {
    result$tables<-.x13_decorate_tables(result$tables,start,f)
    result$row_basis<-rep(list(seq_len(length(preadjustment$tables$a1))),length(result$tables))
    names(result$row_basis)<-names(result$tables)
    result$available_rows<-lapply(names(result$tables),function(name)
      result$row_basis[[name]][!is.na(result$tables[[name]])])
    names(result$available_rows)<-names(result$tables)
  }
  structure(result,class="persephone_native_x13_finals")
}
# The ts decoration of the pinned stored tables, all on the input calendar.
.x13_decorate_tables <- function(tables,start,f) {
  for(name in names(tables))
    tables[[name]]<-stats::ts(tables[[name]],start=start,frequency=f)
  tables
}
