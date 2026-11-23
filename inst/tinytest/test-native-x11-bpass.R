source("helpers/fixtures.R", local = TRUE)
local({
  raw<-getFromNamespace('.native_x11_b_pass','rkalos')
  pub<-getFromNamespace('native_x11_b_pass','rkalos')
  source_fixture('fixtures/native/x11-bpass-cases.R',local=TRUE)
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  manifest<-read.csv('fixtures/native/x11-bpass/cases.csv')
  inputs<-read_numeric_fixture('fixtures/native/x11-bpass/inputs.csv')
  tables<-read_numeric_fixture('fixtures/native/x11-bpass/tables.csv')
  cases<-x11_bpass_cases()
  expect_identical(sort(names(cases)),sort(unique(manifest$case)))
  expect_identical(nrow(manifest),length(cases))
  labels<-c('B1','B2','B3','B4','B4a','B4anorm','B4d','B4g','B5','B6','B7','B8','B9',
    'B9d','B9g','B10','B11','B13','B17','B20')
  center<-c('B2','B3','B4','B4a','B4anorm','B4d','B4g')
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
    expect_identical(as.integer(pinned$b2drop),t)
    r<-raw(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,ca$lsig,ca$usig)
    expect_identical(r$status,'available',info=paste(ca$id,'status'))
    expect_identical(r$b2drop,t,info=paste(ca$id,'b2drop'))
    expect_identical(r$trend_length,as.integer(pinned$trend),info=paste(ca$id,'trend length'))
    expect_identical(r$initial_filters,as.integer(ca$initial),info=paste(ca$id,'initial request'))
    expect_identical(r$final_filters,as.integer(ca$final),info=paste(ca$id,'final request'))
    # All 20 pinned tables bit-exact on their pinned spans; the B9d detection
    # series is carried by correction round 2.
    for(label in labels) {
      rows<-which(tables$case==ca$id & tables$table==label)
      pv<-tables$value[rows][order(as.integer(tables$row[rows]))]
      got<-if(label=='B9d') r$rounds[[2]]$detection else r$tables[[label]]
      expect_true(bits(got,pv),info=paste(ca$id,label))
      expect_identical(length(got),if(label %in% center) length(v)-2L*t else length(v),
        info=paste(ca$id,label,'span'))
    }
    # Correction history: required legacy policy, pinned 0-based analysis
    # positions ((first+t) %% f, first, first) on the 1-based rows, the
    # detection linkage, the round replacement counts and the final-round
    # B17/B20 linkage.
    first<-as.integer(ca$start_period-1L)
    expect_identical(vapply(r$rounds,function(x)x$policy,character(1)),
      rep('legacy',3L),info=paste(ca$id,'round policies'))
    expect_identical(r$rounds[[1]]$analysis_position,as.integer((first+t)%%ca$frequency)+1L,
      info=paste(ca$id,'round 1 analysis position'))
    expect_identical(r$rounds[[2]]$analysis_position,first+1L,info=paste(ca$id,'round 2 analysis position'))
    expect_identical(r$rounds[[3]]$analysis_position,first+1L,info=paste(ca$id,'round 3 analysis position'))
    expect_identical(r$rounds[[1]]$replacement_count,as.integer(pinned$replacements_r1),
      info=paste(ca$id,'round 1 replacements'))
    expect_identical(r$rounds[[2]]$replacement_count,as.integer(pinned$replacements_r2),
      info=paste(ca$id,'round 2 replacements'))
    expect_identical(sum(!is.na(r$tables$B4)),r$rounds[[1]]$replacement_count,
      info=paste(ca$id,'B4 correction count'))
    expect_identical(sum(!is.na(r$tables$B9)),r$rounds[[2]]$replacement_count,
      info=paste(ca$id,'B9 correction count'))
    expect_true(bits(r$rounds[[1]]$detection,r$tables$B4d),info=paste(ca$id,'round 1 detection'))
    expect_true(bits(r$rounds[[3]]$detection,r$tables$B13),info=paste(ca$id,'round 3 detection'))
    expect_true(bits(r$rounds[[3]]$sigma_weights,r$tables$B17),info=paste(ca$id,'B17 linkage'))
    expect_true(bits(r$rounds[[3]]$corrections,r$tables$B20),info=paste(ca$id,'B20 linkage'))
    expect_identical(length(v)-sum(r$tables$B17==1),as.integer(pinned$flagged_r3),
      info=paste(ca$id,'round 3 flagged rows'))
    expect_true(all(r$tables$B20[r$tables$B17==1]==0),info=paste(ca$id,'default corrections'))
    # Additive composition identities, recomputed exactly in R.
    expect_true(bits(r$tables$B8,r$tables$B1-r$tables$B7),info=paste(ca$id,'B8 identity'))
    expect_true(bits(r$tables$B13,r$tables$B11-r$tables$B7),info=paste(ca$id,'B13 identity'))
    expect_true(bits(r$tables$B6,r$tables$B1-r$tables$B5),info=paste(ca$id,'B6 identity'))
    expect_true(bits(r$tables$B11,r$tables$B1-r$tables$B10),info=paste(ca$id,'B11 identity'))
    expect_true(bits(r$tables$B4d,r$tables$B3-r$tables$B4anorm),info=paste(ca$id,'B4d identity'))
    expect_true(bits(r$tables$B3,r$tables$B1[(t+1L):(length(v)-t)]-r$tables$B2),info=paste(ca$id,'B3 identity'))
    again<-raw(v,ca$frequency,ca$start_period,'legacy',ca$initial,ca$final,ca$trend,ca$lsig,ca$usig)
    expect_identical(again,r,info=paste(ca$id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(ca$id,'input untouched'))
    acc<-raw(v,ca$frequency,ca$start_period,'accurate',ca$initial,ca$final,ca$trend,ca$lsig,ca$usig)
    expect_identical(acc$status,'available',info=paste(ca$id,'accurate status'))
    expect_identical(acc$b2drop,r$b2drop)
    expect_identical(acc$policy,'accurate')
    expect_true(bits(acc$tables$B1,r$tables$B1),info=paste(ca$id,'accurate same input'))
  }
  # Public wrapper: explicit policy, ts domain with original-calendar tables
  # on the full and cropped spans, named per-phase filters, round provenance.
  case<-cases[['basic-4-36-1-13-1.5-2.5']]
  inp<-x11_bpass_input(case);v<-inp$values;v[inp$missing]<-NA_real_
  series<-stats::ts(v,start=c(1998,case$start_period),frequency=4)
  out<-pub(series,'legacy',case$initial,case$final,13L)
  expect_true(inherits(out,'persephone_native_x11_b_pass'),info='wrapper class')
  expect_identical(stats::tsp(out$tables$B1),c(1998,1998+35/4,4),info='B1 calendar')
  expect_identical(stats::tsp(out$tables$B2),c(1998.5,2006.25,4),info='B2 calendar')
  expect_identical(as.integer(stats::cycle(out$tables$B2))[1L],3L,info='B2 first cycle')
  expect_identical(stats::cycle(out$rounds[[2]]$detection),stats::cycle(out$tables$B13),
    info='B9d calendar')
  expect_identical(out$rounds[[1]]$analysis_period,3L,info='round 1 analysis period')
  expect_identical(vapply(out$rounds,`[[`,integer(1),'analysis_period')[2:3],c(1L,1L),
    info='rounds 2/3 analysis period')
  expect_identical(out$row_basis$B2,3:34,info='cropped row basis')
  expect_identical(out$available_rows$B4,out$row_basis$B4[!is.na(out$tables$B4)])
  expect_identical(length(out$available_rows$B4),6L,info='available corrections')
  again<-pub(series,'legacy',setNames(as.numeric(case$initial),as.character(1:4)),
    setNames(as.numeric(case$final),as.character(1:4)),13L)
  expect_identical(again$tables,out$tables,info='named per-phase filters')
  plain<-pub(v,'legacy',case$initial,case$final,13L,frequency=4,start_period=case$start_period)
  expect_identical(stats::tsp(plain$tables$B1),c(1,1+35/4,4),info='plain-vector calendar')
  expect_true(all(vapply(names(out$tables),function(nm)
    bits(plain$tables[[nm]],as.numeric(out$tables[[nm]])),logical(1))),info='plain-vector values')
  dflt<-pub(series,initial_filters=case$initial,final_filters=case$final,trend_length=13L)
  expect_identical(dflt$policy,'legacy',info='legacy default echo')
  expect_identical(dflt$tables,out$tables,info='omitted policy is the explicit legacy run')
  expect_error(pub(series,initial_filters=case$initial,final_filters=case$final),
    'stated explicitly') # the filters/trend guard is untouched
  expect_error(pub(series,'legacy',case$initial,case$final[-1],13L),'per calendar period')
  expect_error(raw(as.numeric(1:40),4,1,'modern',rep(5L,4),rep(9L,4),13,1.5,2.5),
    'legacy or accurate')
  expect_error(raw(as.numeric(1:40),3,1,'legacy',rep(5L,3),rep(9L,3),13,1.5,2.5),
    'monthly and quarterly')
  expect_error(raw(as.numeric(1:40),4,5,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5),'first period')
  expect_error(raw(as.numeric(1:40),4,1,'legacy',rep(5L,4),rep(9L,4),14,1.5,2.5),'odd Henderson')
  expect_error(raw(as.numeric(1:40),4,1,'legacy',rep(5L,4),rep(9L,4),41,1.5,2.5),
    'at least the trend filter length')
  expect_error(raw(as.numeric(1:20),12,1,'legacy',rep(5L,12),rep(9L,12),13,1.5,2.5),
    'three annual radii')
  expect_error(raw(as.numeric(1:40),4,1,'legacy',rep(5L,3),rep(9L,4),13,1.5,2.5),
    'per calendar period')
  expect_error(raw(as.numeric(1:40),4,1,'legacy',rep(2L,4),rep(9L,4),13,1.5,2.5),
    'unsupported fixed seasonal request')
  expect_error(raw(c(1:19,Inf,21:40),4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5),
    'finite or missing')
})
