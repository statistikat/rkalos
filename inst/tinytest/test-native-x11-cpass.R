source("helpers/fixtures.R", local = TRUE)
local({
  rawb<-getFromNamespace('.native_x11_b_pass','rkalos')
  rawc<-getFromNamespace('.native_x11_c_pass','rkalos')
  pub<-getFromNamespace('native_x11_c_pass','rkalos')
  source_fixture('fixtures/native/x11-bpass-cases.R',local=TRUE)
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  manifest<-read.csv('fixtures/native/x11-cpass/cases.csv')
  inputs<-read.csv('fixtures/native/x11-cpass/inputs.csv')
  tables<-read.csv('fixtures/native/x11-cpass/tables.csv')
  cases<-x11_bpass_cases()  # the C fixture replays the 18 B-pass cases
  expect_identical(sort(names(cases)),sort(unique(manifest$case)))
  expect_identical(nrow(manifest),length(cases))
  labels<-c('C1','C2','C4','C5','C5a','C6','C7','C9','C10','C10a','C11','C13','C17','C20')
  center<-c('C2','C4','C5a')
  for(ca in cases) {
    inp<-x11_bpass_input(ca);v<-inp$values;v[inp$missing]<-NA_real_
    before<-serialize(v,NULL)
    pinned<-manifest[manifest$case==ca$id,]
    expect_identical(as.integer(pinned$frequency),ca$frequency)
    expect_identical(as.integer(pinned$n),length(v))
    rows<-which(inputs$case==ca$id)
    iv<-inputs$value[rows][order(as.integer(inputs$row[rows]))]
    im<-inputs$missing[rows][order(as.integer(inputs$row[rows]))]
    expect_identical(is.na(v),im,info=paste(ca$id,'input mask'))
    expect_true(bits(v[!im],iv[!im]),info=paste(ca$id,'input values'))
    t<-as.integer(ca$frequency/2)
    expect_identical(as.integer(pinned$c2drop),t)
    b<-rawb(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,ca$lsig,ca$usig)
    r<-rawc(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,b,ca$lsig,ca$usig)
    expect_identical(r$status,'available',info=paste(ca$id,'status'))
    expect_identical(r$c2drop,t,info=paste(ca$id,'c2drop'))
    expect_identical(r$trend_length,as.integer(pinned$trend),info=paste(ca$id,'trend length'))
    expect_identical(r$initial_filters,as.integer(ca$initial),info=paste(ca$id,'initial request'))
    expect_identical(r$final_filters,as.integer(ca$final),info=paste(ca$id,'final request'))
    for(label in labels) {
      rows<-which(tables$case==ca$id & tables$table==label)
      pv<-tables$value[rows][order(as.integer(tables$row[rows]))]
      expect_true(bits(r$tables[[label]],pv),info=paste(ca$id,label))
      expect_identical(length(r$tables[[label]]),
        if(label %in% center) length(v)-2L*t else length(v),info=paste(ca$id,label,'span'))
    }
    # Correction history: the pinned original-calendar analysis position, the
    # required legacy policy, the C13/C17/C20 linkage and the flagged and
    # corrected row counts of the final-stage analysis.
    first<-as.integer(ca$start_period-1L)
    expect_identical(r$round$policy,'legacy',info=paste(ca$id,'round policy'))
    expect_identical(r$round$analysis_position,first+1L,info=paste(ca$id,'analysis position'))
    expect_true(bits(r$round$detection,r$tables$C13),info=paste(ca$id,'C13 detection'))
    expect_true(bits(r$round$sigma_weights,r$tables$C17),info=paste(ca$id,'C17 linkage'))
    expect_true(bits(r$round$corrections,r$tables$C20),info=paste(ca$id,'C20 linkage'))
    expect_identical(sum(r$tables$C17!=1),as.integer(pinned$flagged_rc),
      info=paste(ca$id,'flagged rows'))
    expect_identical(sum(r$tables$C20!=0),as.integer(pinned$corrected_rc),
      info=paste(ca$id,'corrected rows'))
    expect_true(all(r$tables$C20[r$tables$C17==1]==0),info=paste(ca$id,'default corrections'))
    # Additive composition identities against the pinned input Y, recomputed
    # exactly in R (C11 and C13 subtract from Y, not from C1).
    expect_true(bits(r$tables$C1,v-b$tables$B20),info=paste(ca$id,'C1 identity'))
    expect_true(bits(r$tables$C4,r$tables$C1[(t+1L):(length(v)-t)]-r$tables$C2),
      info=paste(ca$id,'C4 identity'))
    expect_true(bits(r$tables$C6,r$tables$C1-r$tables$C5),info=paste(ca$id,'C6 identity'))
    expect_true(bits(r$tables$C9,r$tables$C1-r$tables$C7),info=paste(ca$id,'C9 identity'))
    expect_true(bits(r$tables$C11,v-r$tables$C10),info=paste(ca$id,'C11 identity'))
    expect_true(bits(r$tables$C13,r$tables$C11-r$tables$C7),info=paste(ca$id,'C13 identity'))
    again<-rawc(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,b,ca$lsig,ca$usig)
    expect_identical(again,r,info=paste(ca$id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(ca$id,'input untouched'))
    # The accurate policy runs when the whole composition states it.
    bacc<-rawb(v,ca$frequency,ca$start_period,'accurate',ca$initial,ca$final,ca$trend,ca$lsig,ca$usig)
    acc<-rawc(v,ca$frequency,ca$start_period,'accurate',ca$initial,ca$final,ca$trend,bacc,ca$lsig,ca$usig)
    expect_identical(acc$status,'available',info=paste(ca$id,'accurate status'))
    expect_identical(acc$c2drop,r$c2drop)
    expect_identical(acc$policy,'accurate')
    expect_true(bits(acc$tables$C1,v-bacc$tables$B20),info=paste(ca$id,'accurate C1 identity'))
    # B-pass coupling: the core rejects mismatched companions before computing.
    bad<-b;bad$rounds[[2]]$policy<-'accurate'
    expect_error(rawc(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,bad,
      ca$lsig,ca$usig),'policy must match the B pass',info=paste(ca$id,'policy coupling'))
    bad2<-b;bad2$status<-''
    expect_error(rawc(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,bad2,
      ca$lsig,ca$usig),'completed B pass',info=paste(ca$id,'status coupling'))
    bad3<-b;bad3$tables$B1[1L]<-bad3$tables$B1[1L]+1
    expect_error(rawc(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,bad3,
      ca$lsig,ca$usig),'differs from the B pass input',info=paste(ca$id,'series coupling'))
  }
  # One case exercises every remaining core-side validation through a
  # shape-valid B companion (argument checks precede the coupling checks).
  v<-as.numeric(1:40); b<-rawb(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5)
  expect_error(rawc(v,4,1,'modern',rep(5L,4),rep(9L,4),13,b,1.5,2.5),'legacy or accurate')
  expect_error(rawc(v,3,1,'legacy',rep(5L,3),rep(9L,3),13,b,1.5,2.5),'monthly and quarterly')
  expect_error(rawc(v,4,5,'legacy',rep(5L,4),rep(9L,4),13,b,1.5,2.5),'first period')
  expect_error(rawc(v,4,1,'legacy',rep(5L,4),rep(9L,4),14,b,1.5,2.5),'odd Henderson')
  expect_error(rawc(v,4,1,'legacy',rep(5L,4),rep(9L,4),41,b,1.5,2.5),
    'at least the trend filter length')
  expect_error(rawc(as.numeric(1:20),12,1,'legacy',rep(5L,12),rep(9L,12),13,
    rawb(as.numeric(1:20),12,1,'legacy',rep(5L,12),rep(9L,12),13,1.5,2.5),1.5,2.5),
    'three annual radii')
  expect_error(rawc(v,4,1,'legacy',rep(2L,4),rep(9L,4),13,b,1.5,2.5),
    'unsupported fixed seasonal request')
  expect_error(rawc(c(1:19,Inf,21:40),4,1,'legacy',rep(5L,4),rep(9L,4),13,b,1.5,2.5),
    'finite or missing')
  expect_error(rawc(v,4,1,'legacy',rep(5L,4),rep(9L,4),23,b,1.5,2.5),'must match the B pass')
  short<-as.numeric(1:32); short_b<-rawb(short,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5)
  expect_error(rawc(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,short_b,1.5,2.5),
    'does not cover the series')
  expect_error(rawc(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,list(tables=list()),1.5,2.5),
    'B-pass result object')
  # Public wrapper: explicit policy, an X11BPass-shaped companion object and
  # original-calendar tables on the full and cropped spans.
  case<-cases[['basic-4-36-1-13-1.5-2.5']]
  inp<-x11_bpass_input(case);v<-inp$values;v[inp$missing]<-NA_real_
  series<-stats::ts(v,start=c(1998,case$start_period),frequency=4)
  bb<-getFromNamespace('native_x11_b_pass','rkalos')(series,'legacy',case$initial,case$final,13L)
  out<-pub(series,'legacy',case$initial,case$final,13L,bb)
  expect_true(inherits(out,'persephone_native_x11_c_pass'),info='wrapper class')
  expect_identical(stats::tsp(out$tables$C1),c(1998,1998+35/4,4),info='C1 calendar')
  expect_identical(stats::tsp(out$tables$C2),c(1998.5,2006.25,4),info='C2 calendar')
  expect_identical(as.integer(stats::cycle(out$tables$C4))[1L],3L,info='C4 first cycle')
  expect_identical(stats::cycle(out$round$detection),stats::cycle(out$tables$C13),
    info='C13 calendar')
  expect_identical(out$round$analysis_period,1L,info='round analysis period')
  expect_identical(out$row_basis$C4,3:34,info='cropped row basis')
  expect_identical(length(out$available_rows$C20),36L,info='C20 fully available')
  raw<-rawc(v,4,case$start_period,'legacy',case$initial,case$final,13L,
    rawb(v,4,case$start_period,'legacy',case$initial,case$final,13L,1.5,2.5),1.5,2.5)
  expect_true(all(vapply(names(out$tables),function(nm)
    bits(as.numeric(out$tables[[nm]]),raw$tables[[nm]]),logical(1))),info='wrapper equals raw')
  again<-pub(series,'legacy',setNames(as.numeric(case$initial),as.character(1:4)),
    setNames(as.numeric(case$final),as.character(1:4)),13L,bb)
  expect_identical(again$tables,out$tables,info='named per-phase filters')
  plain<-pub(v,'legacy',case$initial,case$final,13L,bb,frequency=4,start_period=case$start_period)
  expect_identical(stats::tsp(plain$tables$C1),c(1,1+35/4,4),info='plain-vector calendar')
  expect_error(pub(series,'legacy',case$initial,case$final,13L),
    'requires a B-pass result object')
  expect_error(pub(series),'stated explicitly')
})
