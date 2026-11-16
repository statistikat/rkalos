# Item 22b.3b: same-phase replacement values, neighbor provenance, endpoints,
# insufficient-neighbor fallback and correction factors on top of qualified
# item 22b.3a final weights. Only the replacement pass is computed here; the
# weights themselves are supplied by the caller.
native_x11_replacements <- function(series,weights,frequency=NULL,
    backcasts=0L,forecasts=0L,exclude_projections=FALSE,multiplicative=FALSE) {
  if(!is.null(dim(series)))series<-as.vector(series)
  if(!is.numeric(series))series<-as.numeric(series)
  if(is.null(frequency)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs an explicit frequency.",call.=FALSE)
    frequency<-span[3L]
  }
  result<-.native_x11_replacements(series,frequency,as.numeric(weights),backcasts,forecasts,
    exclude_projections,multiplicative)
  neighbor<-result$neighbors
  storage.mode(neighbor)<-"integer"
  neighbor[neighbor==-1]<-NA_integer_
  result$neighbor_positions<-ifelse(is.na(neighbor),NA,neighbor+1L) # 1-based rows
  result$neighbor_positions<-matrix(result$neighbor_positions,nrow=nrow(neighbor),
    dimnames=list(NULL,paste0("slot",seq_len(4))))
  result$missing<-is.na(series)
  result$weights_length<-length(weights)
  result$input_length<-length(series)
  structure(result,class="persephone_native_x11_replacement")
}
