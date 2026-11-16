source("helpers/fixtures.R", local = TRUE)
local({
  rawb<-getFromNamespace('.native_x11_b_pass','rkalos')
  rawc<-getFromNamespace('.native_x11_c_pass','rkalos')
 rawd<-getFromNamespace('.native_x11_d_pass','rkalos')
  pubb<-getFromNamespace('native_x11_b_pass','rkalos')
  pubc<-getFromNamespace('native_x11_c_pass','rkalos')
  pubd<-getFromNamespace('native_x11_d_pass','rkalos')
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  manifest<-read.csv('fixtures/native/x11-dpass/cases.csv')
  inputs<-read.csv('fixtures/native/x11-dpass/inputs.csv')
  tables<-read.csv('fixtures/native/x11-dpass/tables.csv')
  msr<-read.csv('fixtures/native/x11-dpass/msr.csv')
  expect_identical(nrow(manifest),22L)
  expect_identical(sum(manifest$category=="extended"),4L)
  labels<-c('D1','D2','D4','D5','D6','D7','D8','D9','D9bis','D9_g_bis','D10','D10bis',
    'D11','D11bis','D12','D13')
  center<-c('D2','D4')
  num<-function(x)as.integer(strsplit(x,' ')[[1]])
  input_of<-function(id) {
    rows<-which(inputs$case==id);o<-order(as.integer(inputs$row[rows]))
    v<-inputs$value[rows][o];v[inputs$missing[rows][o]]<-NA_real_;v
  }
  for(i in seq_len(nrow(manifest))) {
    ca<-manifest[i,];id<-ca$case
    v<-input_of(id);before<-serialize(v,NULL)
    f<-as.integer(ca$frequency);t<-as.integer(ca$frequency/2)
    initial<-num(ca$initial_filters);final<-num(ca$final_filters)
    # The extended inputs are hand-extended over the original calendar: the
    # dropped horizon rows leave exactly the original case input.
    expect_identical(as.integer(ca$n),length(v),info=paste(id,'input length'))
    if(ca$category=='extended') {
      base<-sub('^ext-','',id)
      expect_true(bits(v[(as.integer(ca$back)+1L):(as.integer(ca$back)+as.integer(ca$orig_n))],
        input_of(base)),info=paste(id,'extends the base input'))
      expect_identical(as.integer(ca$orig_start_period),
        as.integer((ca$start_period-1L+ca$back)%%f)+1L,info=paste(id,'original first period'))
    } else expect_identical(as.integer(ca$orig_n),length(v))
    b<-rawb(v,f,ca$start_period,'legacy',initial,final,ca$trend,ca$lsig,ca$usig)
    c<-rawc(v,f,ca$start_period,'legacy',initial,final,ca$trend,b,ca$lsig,ca$usig)
    r<-rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,c,ca$back,ca$fc)
    expect_identical(r$status,'available',info=paste(id,'status'))
    expect_identical(r$type,'x11_d_pass')
    expect_identical(r$d2drop,t,info=paste(id,'d2drop'))
    expect_identical(r$backcast,as.integer(ca$back),info=paste(id,'backcast echo'))
    expect_identical(r$forecast,as.integer(ca$fc),info=paste(id,'forecast echo'))
    expect_identical(r$trend_length,as.integer(ca$trend),info=paste(id,'trend echo'))
    expect_identical(r$d7_henderson_length,as.integer(ca$trend),info=paste(id,'D7 echo'))
    expect_identical(r$final_henderson_length,as.integer(ca$trend),info=paste(id,'final echo'))
    expect_identical(r$msr_iterations,0L,info=paste(id,'fixed mode never selects'))
    expect_identical(as.integer(r$seas_filter),final,
      info=paste(id,'fixed seas filter equals the final requests'))
    expect_identical(as.integer(r$initial_filters),initial,info=paste(id,'initial echo'))
    expect_identical(as.integer(r$final_filters),final,info=paste(id,'final echo'))
    for(label in labels) {
      rows<-which(tables$case==id & tables$table==label)
      pv<-tables$value[rows][order(as.integer(tables$row[rows]))]
      expect_true(bits(r$tables[[label]],pv),info=paste(id,label))
      expect_identical(length(r$tables[[label]]),
        if(label %in% center) length(v)-2L*t else length(v),info=paste(id,label,'span'))
    }
    expect_true(bits(r$tables$D9_g_bis,r$tables$D9bis),info=paste(id,'D9_g_bis alias'))
    # The pinned D9 MSRE diagnostic on the dropped original calendar.
    mrows<-msr[msr$case==id,]
    expect_identical(nrow(mrows),f,info=paste(id,'MSRE phases'))
    mp<-as.integer(mrows$phase)+1L
    expect_true(bits(r$msr_table$mean_irregular_evolutions[mp],mrows$ri),
      info=paste(id,'MSRE irregular'))
    expect_true(bits(r$msr_table$mean_seasonal_evolutions[mp],mrows$rs),
      info=paste(id,'MSRE seasonal'))
    expect_identical(r$msr_table$counts[mp],as.integer(mrows$count),info=paste(id,'MSRE counts'))
    expect_true(bits(r$msr_table$msr[mp],mrows$msr),info=paste(id,'MSRE ratios'))
    expect_true(bits(r$msr_table$global_msr,ca$global_msr),info=paste(id,'MSRE global'))
    expect_identical(sum(!is.na(r$tables$D9)),as.integer(ca$d9_nonnan),
      info=paste(id,'D9 mask count'))
    expect_identical(as.integer(ca$d9_nonnan),as.integer(ca$corrected_rc),
      info=paste(id,'D9 count equals the C-round corrections'))
    # Additive identities recomputed exactly in R against the pinned input.
    expect_true(bits(r$tables$D1,v-c$tables$C20),info=paste(id,'D1 identity'))
    expect_true(bits(r$tables$D4,r$tables$D1[(t+1L):(length(v)-t)]-r$tables$D2),
      info=paste(id,'D4 identity'))
    expect_true(bits(r$tables$D6,r$tables$D1-r$tables$D5),info=paste(id,'D6 identity'))
    expect_true(bits(r$tables$D8,v-r$tables$D7),info=paste(id,'D8 identity'))
    expect_true(bits(r$tables$D9bis,r$tables$D1-r$tables$D7),info=paste(id,'D9bis identity'))
    expect_true(bits(r$tables$D11,v-r$tables$D10),info=paste(id,'D11 identity'))
    expect_true(bits(r$tables$D11bis,r$tables$D1-r$tables$D10),info=paste(id,'D11bis identity'))
    expect_true(bits(r$tables$D13,r$tables$D11-r$tables$D12),info=paste(id,'D13 identity'))
    # The pinned D9 dcmpg mask, recomputed.
    mask<-ifelse(is.na(r$tables$D9bis-r$tables$D8) | abs(r$tables$D9bis-r$tables$D8)>=1e-9,
      r$tables$D9bis,NA_real_)
    expect_true(bits(r$tables$D9,mask),info=paste(id,'D9 mask rule'))
    again<-rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,c,ca$back,ca$fc)
    expect_identical(again,r,info=paste(id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(id,'input untouched'))
    # The accurate policy runs when the whole composition states it.
    bacc<-rawb(v,f,ca$start_period,'accurate',initial,final,ca$trend,ca$lsig,ca$usig)
    cacc<-rawc(v,f,ca$start_period,'accurate',initial,final,ca$trend,bacc,ca$lsig,ca$usig)
    dacc<-rawd(v,f,ca$start_period,'accurate',initial,final,ca$trend,cacc,ca$back,ca$fc)
    expect_identical(dacc$status,'available',info=paste(id,'accurate status'))
    expect_identical(dacc$policy,'accurate')
    expect_true(bits(dacc$tables$D1,v-cacc$tables$C20),info=paste(id,'accurate D1 identity'))
    # C-pass coupling: the core rejects mismatched companions before computing.
    bad<-c;bad$round$policy<-'accurate'
    expect_error(rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,bad,ca$back,ca$fc),
      'policy must match the C pass',info=paste(id,'policy coupling'))
    bad2<-c;bad2$status<-''
    expect_error(rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,bad2,ca$back,ca$fc),
      'completed C pass',info=paste(id,'status coupling'))
    bad3<-c;bad3$reference_series[1L]<-bad3$reference_series[1L]+1
    expect_error(rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,bad3,ca$back,ca$fc),
      'differs from the C pass input',info=paste(id,'reference-series coupling'))
    bad4<-c;bad4$trend_length<-bad4$trend_length+2
    expect_error(rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,bad4,ca$back,ca$fc),
      'settings must match the C pass',info=paste(id,'settings coupling'))
  }
  # A C companion computed on a shorter series must be refused by coverage.
  v<-as.numeric(1:36);small<-as.numeric(1:24)
  cshort<-rawc(small,4,1,'legacy',rep(5L,4),rep(9L,4),13,
    rawb(small,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5),1.5,2.5)
  expect_error(rawd(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,cshort,0,0),
    'does not cover the series')
  expect_error(rawd(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,list(tables=list()),0,0),
    'C-pass result object')
  # Core argument validation at the D boundary.
  cfull<-rawc(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,
    rawb(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5),1.5,2.5)
  run<-function(...)rawd(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,cfull,...)
  expect_error(run(-1,0),'horizons are non-negative')
  expect_error(run(0,-1),'horizons are non-negative')
  expect_error(run(21,0),'four cycles of original calendar')
  expect_error(rawd(v,3,1,'legacy',rep(5L,3),rep(9L,3),13,cfull,0,0),
    'monthly and quarterly')
  expect_error(rawd(v,4,5,'legacy',rep(5L,4),rep(9L,4),13,cfull,0,0),'first period')
  expect_error(rawd(v,4,1,'legacy',rep(5L,4),rep(9L,4),14,cfull,0,0),'odd Henderson')
  expect_error(rawd(v,4,1,'legacy',rep(5L,4),rep(9L,4),41,cfull,0,0),
    'at least the trend filter length')
  expect_error(rawd(c(1:35,Inf),4,1,'legacy',rep(5L,4),rep(9L,4),13,cfull,0,0),
    'finite or missing')
  expect_error(rawd(v,4,1,'modern',rep(5L,4),rep(9L,4),13,cfull,0,0),'legacy or accurate')
  # Public wrapper on a ts: the extended case keeps the extended calendar for
  # the tables and carries the dropped original calendar as metadata.
  ca<-manifest[manifest$case=='ext-basic-12-108-1-13-1.5-2.5',]
  v<-input_of(ca$case)
  series<-stats::ts(v,start=c(1997,9),frequency=12)
  bb<-pubb(series,'legacy',num(ca$initial_filters),num(ca$final_filters),13L)
  cc<-pubc(series,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,bb)
  dd<-pubd(series,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,cc,
    back=4,forecast=3)
  expect_true(inherits(dd,'persephone_native_x11_d_pass'),info='wrapper class')
  extended_start<-1997+8/12;crop_start<-1998+2/12
  expect_identical(stats::tsp(dd$tables$D1),
    c(extended_start,extended_start+114/12,12),info='D1 extended calendar')
  expect_identical(stats::tsp(dd$tables$D2),
    c(crop_start,crop_start+102/12,12),info='D2 cropped calendar')
  expect_identical(as.integer(stats::cycle(dd$tables$D2))[1L],3L,info='D2 first cycle')
  expect_identical(as.integer(stats::cycle(dd$tables$D4))[1L],3L,info='D4 first cycle')
  expect_identical(dd$original_domain$start,c(1998,1),info='original calendar start')
  expect_identical(dd$original_length,108L,info='original length')
  expect_identical(nrow(dd$msr_table),12L,info='MSRE table rows')
  expect_identical(names(dd$msr_table),c('phase','mean_irregular_evolutions','mean_seasonal_evolutions',
    'counts','msr'))
  raw<-rawd(v,12,9,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,
    rawc(v,12,9,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,
      rawb(v,12,9,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,1.5,2.5),1.5,2.5),
    4,3)
  expect_true(all(vapply(names(dd$tables),function(nm)
    bits(as.numeric(dd$tables[[nm]]),raw$tables[[nm]]),logical(1))),info='wrapper equals raw')
  # A ts-decorated C companion rebuilds identically to the raw one.
  raw_ts<-rawd(v,12,9,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,cc,4,3)
  expect_identical(raw_ts$tables,raw$tables,info='decorated C companion chains')
  again<-pubd(series,'legacy',setNames(as.numeric(num(ca$initial_filters)),as.character(1:12)),
    setNames(as.numeric(num(ca$final_filters)),as.character(1:12)),13L,cc,back=4,forecast=3)
  expect_identical(again$tables,dd$tables,info='named per-phase filters')
  plain<-pubd(v,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,
    rawc(v,12,9,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,
      rawb(v,12,9,'legacy',num(ca$initial_filters),num(ca$final_filters),13L,1.5,2.5),1.5,2.5),
    back=4,forecast=3,frequency=12,start_period=9)
  expect_identical(plain$original_domain$start,c(2,1),info='plain-vector original calendar')
  expect_error(pubd(series,'legacy',num(ca$initial_filters),num(ca$final_filters),13L),
    'requires a C-pass result object')
  dflt<-pubd(series,initial_filters=num(ca$initial_filters),final_filters=num(ca$final_filters),
    trend_length=13L,c=cc,back=4,forecast=3)
  expect_identical(dflt$policy,dd$policy,info='omitted policy inherits the C-pass echo')
  expect_identical(dflt$tables,dd$tables,info='inherited policy is the explicit run')
  expect_error(pubd(series,'accurate',num(ca$initial_filters),num(ca$final_filters),13L,cc,
    back=4,forecast=3),'inherits its sigma weight policy')
})
