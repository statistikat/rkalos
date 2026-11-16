source("helpers/fixtures.R", local = TRUE)
local({
  raws<-getFromNamespace('.native_x11_selection_pass','rkalos')
  rawb<-getFromNamespace('.native_x11_b_pass','rkalos')
  rawc<-getFromNamespace('.native_x11_c_pass','rkalos')
  rawd<-getFromNamespace('.native_x11_d_pass','rkalos')
  pubs<-getFromNamespace('native_x11_selection_pass','rkalos')
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  cases<-read.csv('fixtures/native/x11-selection/cases.csv')
  inputs<-read.csv('fixtures/native/x11-selection/inputs.csv')
  tables<-read.csv('fixtures/native/x11-selection/tables.csv')
  sels<-read.csv('fixtures/native/x11-selection/selections.csv')
  msr<-read.csv('fixtures/native/x11-selection/msr.csv')
  expect_identical(nrow(cases),15L)
  request_names<-function(x) { y<-as.integer(strsplit(x,',')[[1]])
    ifelse(y==0,'stable',ifelse(y==3,'s3x3',ifelse(y==5,'s3x5',
      ifelse(y==9,'s3x9',ifelse(y==-1,'x11default','msr'))))) }
  input_of<-function(id) { rows<-which(inputs$case==id)
    inputs$value[rows][order(as.integer(inputs$row[rows]))] }
  bset<-c('B7','B13','B20'); cset<-c('C7','C13','C20')
  for(i in seq_len(nrow(cases))) {
    ca<-cases[i,];id<-ca$case
    v<-input_of(id);before<-serialize(v,NULL)
    f<-as.integer(ca$frequency)
    req<-request_names(ca$requests)
    r<-raws(v,f,ca$first_period,req,ca$trend,'legacy','jdplus',1.5,2.5,ca$backcast,ca$forecast)
    expect_identical(r$type,'x11_selection_pass')
    expect_identical(as.character(r$requests),req,info=paste(id,'request echo'))
    expect_identical(r$automatic,ca$trend==0,info=paste(id,'automatic echo'))
    # The qualified plans resolve every pinned stage (records exist even for
    # unsupported runs); the delivered arrays echo the D-stage resolutions.
    expect_identical(as.integer(r$initial_filters),r$initial_plan$base_filters,
      info=paste(id,'initial plan echo'))
    expect_identical(as.integer(r$final_filters),r$final_plan$base_filters,
      info=paste(id,'final plan echo'))
    expect_identical(r$final_plan$requires_msr,as.logical(ca$engaged),
      info=paste(id,'requires_msr echoes the pinned isMSR'))
    if(r$status=='available')
      expect_identical(r$msr$engaged,as.logical(ca$engaged),info=paste(id,'MSR engaged'))
    else expect_false(r$msr$engaged,info=paste(id,'unsupported runs stop before the D final step'))
    # Pinned stage records: mode, lengths and the always-computed final ratio.
    srows<-sels[sels$case==id,]
    expect_identical(nrow(srows),4L,info=paste(id,'stage rows'))
    for(j in seq_len(nrow(srows))) {
      row<-srows[j,];st<-row$stage;rec<-r$stages[[st]]
      if(row$mode=='pinned-chooses-longest-not-claimed') {
        # Degenerate pinned controls: the qualified rule refuses the NaN
        # ratio at the first automatic stage and the run stops there; the
        # later stage records stay default-constructed.
        if(st=='b7') {
          expect_identical(rec$status,'unsupported',info=paste(id,st,'refuses NaN'))
          expect_true(is.na(rec$ic_ratio),info=paste(id,st,'NaN refused'))
        }
      } else if(row$mode=='fixed') {
        expect_identical(rec$status,'fixed',info=paste(id,st,'fixed status'))
        expect_identical(rec$length,as.integer(row$pinned_len),info=paste(id,st,'fixed length'))
        # The pinned stores a ratio at fixed stages; the qualified records
        # carry one only where the selection engaged (and at the final step).
        if(st=='final') expect_true(bits(rec$ic_ratio,row$ic_ratio),info=paste(id,'final ratio'))
        else expect_true(is.na(rec$ic_ratio),info=paste(id,st,'fixed mode computes no ratio'))
      } else {
        expect_identical(rec$status,'selected',info=paste(id,st,'selected status'))
        expect_identical(rec$length,as.integer(row$pinned_len),info=paste(id,st,'selected length'))
        expect_true(bits(rec$ic_ratio,row$ic_ratio),info=paste(id,st,'selected ratio'))
        if(st=='b7'&&row$ic_ratio>=1)
          expect_identical(rec$length,f+1L,info=paste(id,'inline ratio>=1 takes period+1'))
      }
    }
    expect_true(bits(r$stages$final$ic_ratio,ca$icratio),info=paste(id,'final iCRatio'))
    if(id %in% c('const-fixed','const-msr')) {
      # Explicit unsupported passthrough: records and plans only, no tables,
      # and the first refused stage names the status.
      expect_identical(r$status,'unsupported:trend-b7',info=paste(id,'unsupported status'))
      expect_true(is.null(r$b)&&is.null(r$c)&&is.null(r$d),info=paste(id,'records only'))
      expect_true(!is.null(r$initial_plan)&&!is.null(r$stages$final),info=paste(id,'records kept'))
      next
    }
    expect_identical(r$status,'available',info=paste(id,'status'))
    expect_identical(r$stages$d7$length,as.integer(ca$d7len),info=paste(id,'D7 length'))
    expect_identical(r$stages$final$length,as.integer(ca$finallen),info=paste(id,'final length'))
    expect_identical(r$b$status,'available')
    expect_identical(r$c$status,'available')
    expect_identical(r$d$status,'available')
    fin<-as.integer(strsplit(gsub('"','',ca$seas),',')[[1]])
    expect_identical(as.integer(r$d$seas_filter),fin,info=paste(id,'seasFilter fill'))
    expect_identical(r$d$msr_iterations,as.integer(ca$d9msriter),info=paste(id,'d9msriter'))
    expect_identical(sum(!is.na(r$b$tables$B20)&r$b$tables$B20!=0),as.integer(ca$b20nz),
      info=paste(id,'B20 corrections'))
    # The global-MSR walk: bit-exact attempt rows plus the final decision.
    arows<-msr[msr$case==id,]
    if(as.logical(ca$engaged)) {
      arows<-arows[order(as.integer(arows$attempt)),]
      evals<-arows[arows$band!='no_history',]
      expect_identical(length(r$msr$attempts$count),nrow(evals),info=paste(id,'attempt count'))
      expect_identical(as.integer(r$msr$attempts$count),as.integer(evals$count),
        info=paste(id,'attempt windows'))
      expect_identical(as.integer(r$msr$attempts$years_removed),as.integer(evals$years_removed),
        info=paste(id,'attempt years removed'))
      expect_true(bits(r$msr$attempts$global_msr,evals$global_msr),info=paste(id,'attempt ratios'))
      expect_identical(as.character(r$msr$attempts$band),as.character(evals$band),
        info=paste(id,'attempt bands'))
      last<-arows[nrow(arows),]
      expect_identical(r$msr$selected_filter,as.integer(last$final_filter),
        info=paste(id,'selected filter'))
      expect_identical(r$msr$using_default,as.logical(last$using_default),
        info=paste(id,'method default flag'))
      expect_identical(r$msr$status,'available',info=paste(id,'MSR status'))
      if(id=='q-exh20')
        expect_identical(r$msr$reason,'insufficient_remaining_history',info=paste(id,'exhausted'))
      if(id=='q-iter0-19') {
        expect_identical(r$msr$reason,'insufficient_initial_history',info=paste(id,'never engaged'))
        expect_identical(r$msr$iterations,0L)
        expect_identical(r$msr$global_msr,0) # the pinned field default
      }
      if(id=='q-longmsr')
        expect_identical(nrow(evals),19L,info=paste(id,'19-evaluation walk'))
    } else {
      expect_identical(r$msr$reason,'unavailable',info=paste(id,'unengaged reason'))
      expect_identical(r$msr$selected_filter,0L,info=paste(id,'undecided'))
      expect_identical(length(r$msr$attempts$count),0L,info=paste(id,'no attempts'))
    }
    # The captured 22-table selected state, bit-exact through the nested passes.
    for(label in unique(tables$table[tables$case==id])) {
      got<-if(label %in% bset) r$b$tables[[label]] else
        if(label %in% cset) r$c$tables[[label]] else r$d$tables[[label]]
      rows<-which(tables$case==id&tables$table==label)
      pv<-tables$value[rows][order(as.integer(tables$row[rows]))]
      expect_true(bits(got,pv),info=paste(id,label))
    }
    again<-raws(v,f,ca$first_period,req,ca$trend,'legacy','jdplus',1.5,2.5,ca$backcast,ca$forecast)
    expect_identical(again,r,info=paste(id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(id,'input untouched'))
  }
  # The five-year-cap policy truncates the 19-evaluation walk against the
  # same pinned trace: the prefix tables are identical, D10bis diverges.
  v<-input_of('q-longmsr')
  jd<-raws(v,4,1,rep('msr',4),13,'legacy','jdplus',1.5,2.5,0,0)
  cap<-raws(v,4,1,rep('msr',4),13,'legacy','x11_five_year_cap',1.5,2.5,0,0)
  expect_identical(cap$status,'available')
  expect_identical(cap$msr$iterations,6L,info='cap stops after five year removals')
  expect_identical(cap$msr$reason,'five_year_cap')
  expect_identical(cap$msr$selected_filter,5L)
  expect_true(cap$msr$using_default)
  expect_identical(cap$msr$attempts$count,jd$msr$attempts$count[1:6])
  expect_identical(cap$msr$attempts$global_msr,jd$msr$attempts$global_msr[1:6],
    info='truncated prefix of the pinned trace')
  prefix<-c('D1','D2','D4','D5','D6','D7','D8','D9','D9bis','D9_g_bis')
  expect_true(all(vapply(prefix,function(label) identical(cap$d$tables[[label]],
    jd$d$tables[[label]]),logical(1))),info='prefix tables identical while the walk continues')
  expect_false(identical(cap$d$tables$D10bis,jd$d$tables$D10bis),info='D10bis diverges')
  # The selected C state chains into the fixed D boundary: without an MSR fill
  # the composition reproduces the selected D exactly (q-mixed13 provisional
  # requests equal the seas fill); with one, only the pre-D-final stages do.
  mixed<-raws(input_of('q-mixed13'),4,1,c('msr','msr','s3x5','s3x5'),13,'legacy','jdplus',
    1.5,2.5,0,0)
  chained<-rawd(input_of('q-mixed13'),4,1,'legacy',mixed$c$initial_filters,
    mixed$c$final_filters,mixed$c$trend_length,mixed$c,0,0)
  expect_identical(lapply(chained$tables,as.numeric),
    lapply(mixed$d$tables,as.numeric),info='fixed chain reproduces the selected D state')
  msr13<-raws(input_of('q-msr13'),4,1,rep('msr',4),13,'legacy','jdplus',1.5,2.5,0,0)
  chained2<-rawd(input_of('q-msr13'),4,1,'legacy',msr13$c$initial_filters,
    msr13$c$final_filters,msr13$c$trend_length,msr13$c,0,0)
  expect_true(all(vapply(prefix,function(label) identical(chained2$tables[[label]],
    msr13$d$tables[[label]]),logical(1))),info='pre-D-final tables chain unchanged')
  expect_false(identical(chained2$tables$D10bis,msr13$d$tables$D10bis),
    info='the all-MSR seas fill only reaches the nested D10bis')
  chained_c<-rawc(input_of('q-msr13'),4,1,'legacy',msr13$c$initial_filters,
    msr13$c$final_filters,msr13$c$trend_length,msr13$b,1.5,2.5)
  expect_identical(chained_c$tables,msr13$c$tables,info='the selected B record chains into the C pass')
  # Core argument validation at the selection boundary.
  v<-as.numeric(1:40)
  expect_error(raws(c(1:39,NA_real_),4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,0,0),
    'complete observations')
  expect_error(raws(v,3,1,rep('s3x5',3),13,'legacy','jdplus',1.5,2.5,0,0),
    'monthly and quarterly')
  expect_error(raws(v,4,5,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,0,0),'first period')
  expect_error(raws(v,4,1,rep('s3x5',4),4,'legacy','jdplus',1.5,2.5,0,0),
    'is 0 .automatic. or an odd Henderson length')
  expect_error(raws(v,4,1,rep('s3x5',3),13,'legacy','jdplus',1.5,2.5,0,0),
    'one request for every calendar period')
  expect_error(raws(v,4,1,c(rep('s3x5',3),'hazard'),13,'legacy','jdplus',1.5,2.5,0,0),
    'unsupported seasonal request')
  expect_error(raws(v,4,1,rep('s3x5',4),13,'modern','jdplus',1.5,2.5,0,0),'legacy or accurate')
  expect_error(raws(v,4,1,rep('s3x5',4),13,'legacy','counting',1.5,2.5,0,0),
    'MSR iteration policy')
  expect_error(raws(v,4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,-1,0),
    'horizons are non-negative')
  expect_error(raws(v,4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,25,0),
    'four cycles of original calendar')
  # Public wrapper: explicit policy and requests, decorated nested passes and
  # the decision records stay intact for unsupported runs.
  ca<-cases[cases$case=='q-msr13',]
  v<-input_of(ca$case)
  series<-stats::ts(v,start=c(2001,ca$first_period),frequency=4)
  out<-pubs(series,rep('msr',4),'legacy',13L)
  expect_true(inherits(out,'persephone_native_x11_selection_pass'),info='wrapper class')
  expect_true(out$selected,info='selected flag')
  expect_identical(stats::tsp(out$d$tables$D11),stats::tsp(series),info='D11 calendar')
  expect_identical(as.integer(stats::cycle(out$c$tables$C2))[1L],3L,info='C2 cropped calendar')
  expect_identical(dim(out$d$msr_table),c(4L,5L),info='MSRE data.frame')
  expect_true(is.data.frame(out$msr$attempts)&&nrow(out$msr$attempts)==out$msr$iterations,
    info='attempt trace data.frame')
  expect_identical(names(out$requests),as.character(1:4),info='named request echo')
  expect_identical(out$stages,raws(v,4,ca$first_period,rep('msr',4),13,'legacy','jdplus',
    1.5,2.5,0,0)$stages,info='stage records identical to raw')
  again<-pubs(series,setNames(rep('msr',4),as.character(1:4)),'legacy',13L)
  expect_identical(again$d$tables,out$d$tables,info='named per-phase requests')
  bad<-pubs(input_of('const-fixed'),rep('s3x5',4),'legacy',0L,frequency=4,start_period=1)
  expect_identical(bad$status,'unsupported:trend-b7',info='wrapper passthrough status')
  expect_false(bad$selected)
  expect_identical(bad$records_only,
    'the stage records and plans are complete; no tables are delivered for an unsupported run')
  dsel<-pubs(input_of('const-fixed'),rep('s3x5',4),trend_length=0L,frequency=4,start_period=1)
  expect_identical(dsel,bad,info='omitted policy is the explicit legacy run')
  expect_error(pubs(series,policy='legacy',trend_length=13L),
    'requests and the trend filter length')
})
