# The selected-state pass calendar decoration, shared verbatim between the
# pass-level selection wrapper and the restricted X-13 fit wrapper: a B, C or
# D record of a selection serialization becomes exactly the decorated pass the
# fixed-pass wrappers deliver, wherever the record came from (ordinal: the
# 1-based ordinal of the series start on the calendar; f: frequency; n: the
# series length the passes were run over). The D9 MSRE table stays the pinned
# per-phase diagnostic, not a series.
.x11_decorate_pass <- function(pass,center,ordinal,f,n) {
  drop<-if(identical(pass$type,"x11_b_pass"))pass$b2drop else
    if(identical(pass$type,"x11_c_pass"))pass$c2drop else pass$d2drop
  domain<-.native_time_domain(floor(ordinal/f),ordinal%%f+1L,f,n)
  domain$start_ordinal<-ordinal
  cropped<-drop+seq_len(n-2L*drop)
  crop_ordinal<-ordinal+drop
  crop<-.native_time_domain(floor(crop_ordinal/f),crop_ordinal%%f+1L,f,length(cropped))
  crop$start_ordinal<-crop_ordinal
  for(name in names(pass$tables))
    pass$tables[[name]]<-stats::ts(pass$tables[[name]],
      start=if(name %in% center) crop$start else domain$start,frequency=f)
  pass$row_basis<-lapply(names(pass$tables),function(name)
    if(name %in% center) cropped else seq_len(n))
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
      n-pass$backcast-pass$forecast)
    original$start_ordinal<-original_ordinal
    pass$original_domain<-original
    pass$original_length<-n-pass$backcast-pass$forecast
  }
  pass
}
# The selection-run record decoration shared the same way: the per-period
# request/filter/plan names, the three pass decorations, optional notes
# (appended after the passes, as the selection wrapper has always ordered
# them), the MSR attempt table and the selected flag, applied to a complete
# selection serialization record.
.x11_decorate_selection_record <- function(selection,ordinal,f,n,notes=NULL) {
  periods<-as.character(seq_len(f))
  names(selection$requests)<-periods
  for(field in c("initial_filters","final_filters"))names(selection[[field]])<-periods
  for(field in c("initial_plan","final_plan"))names(selection[[field]]$base_filters)<-periods
  if(identical(selection$status,"available")) {
    selection$b<-.x11_decorate_pass(selection$b,c("B2","B3","B4","B4a","B4anorm","B4d","B4g"),ordinal,f,n)
    selection$c<-.x11_decorate_pass(selection$c,c("C2","C4","C5a"),ordinal,f,n)
    selection$d<-.x11_decorate_pass(selection$d,c("D2","D4"),ordinal,f,n)
  }
  for(name in names(notes))selection[[name]]<-notes[[name]]
  selection$msr$attempts<-data.frame(count=selection$msr$attempts$count,
    years_removed=selection$msr$attempts$years_removed,global_msr=selection$msr$attempts$global_msr,
    band=selection$msr$attempts$band,status=selection$msr$attempts$status)
  selection$selected<-identical(selection$status,"available")
  selection
}
