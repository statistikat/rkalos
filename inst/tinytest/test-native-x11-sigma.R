source("helpers/fixtures.R", local = TRUE)
local({
  raw<-getFromNamespace('.native_x11_sigma_weights','rkalos')
  pub<-getFromNamespace('native_x11_sigma_weights','rkalos')
  source_x11_capture_generator('fixtures/native/x11-sigma-cases.R')
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  manifest<-read.csv('fixtures/native/x11-sigma/cases.csv')
  inputs<-read_numeric_fixture('fixtures/native/x11-sigma/inputs.csv', 'value')
  weights<-read_numeric_fixture('fixtures/native/x11-sigma/weights.csv', 'weight')
  disp1<-read_numeric_fixture('fixtures/native/x11-sigma/dispersions.csv', 'dispersion')
  disp2<-read_numeric_fixture('fixtures/native/x11-sigma/dispersions-pass2.csv', 'dispersion')
  cases<-x11_sigma_cases()
  expect_identical(sort(names(cases)),sort(unique(manifest$case)))
  expect_identical(nrow(manifest),length(cases))
  for(ca in cases) {
    inp<-x11_sigma_input(ca);v<-inp$values;v[inp$missing]<-NA_real_
    pinned<-manifest[manifest$case==ca$id,]
    expect_identical(as.integer(pinned$frequency),ca$frequency)
    expect_identical(as.integer(pinned$n),length(v))
    rows<-which(inputs$case==ca$id)
    iv<-inputs$value[rows][order(as.integer(inputs$row[rows]))]
    im<-inputs$missing[rows][order(as.integer(inputs$row[rows]))]
    expect_identical(is.na(v),im,info=paste(ca$id,'input mask'))
    expect_true(bits(v[!im],iv[!im]),info=paste(ca$id,'regenerated input vs capture'))
    # The native replay consumes the same binary64 inputs as the JVM capture.
    # Keep the independent local generator comparison above.
    v <- iv
    v[im] <- NA_real_
    before <- serialize(v, NULL)
    ro<-order(as.integer(weights$row[weights$case==ca$id]))
    w<-weights$weight[weights$case==ca$id][ro];p<-weights$pass[weights$case==ca$id][ro]
    final<-w[p==max(p)];first<-w[p==1]
    yo<-order(as.integer(disp1$year[disp1$case==ca$id]))
    d1<-disp1$dispersion[disp1$case==ca$id][yo]
    r<-raw(v,ca$frequency,ca$start_period,'legacy',ca$backcasts,ca$forecasts,
      ca$exclude,ca$lsig,ca$usig,ca$multiplicative)
    expect_identical(r$status,'available',info=paste(ca$id,'status'))
    expect_identical(r$passes,as.integer(pinned$passes),info=paste(ca$id,'passes'))
    expect_true(bits(r$weights,final),info=paste(ca$id,'final weights'))
    expect_true(bits(r$first_pass_weights,first),info=paste(ca$id,'first-pass weights'))
    expect_true(bits(r$first_pass_windows$dispersions,d1),info=paste(ca$id,'first-pass dispersions'))
    d2<-if(r$passes==2) {
      disp2$dispersion[disp2$case==ca$id][order(as.integer(disp2$year[disp2$case==ca$id]))]
    } else d1
    expect_true(bits(r$windows$dispersions,d2),info=paste(ca$id,'final dispersions'))
    expect_identical(as.integer(pinned$window_entries),length(d1),info=paste(ca$id,'window count'))
    expect_identical(r$first_pass_zeroed,as.integer(which(r$first_pass_weights==0)),
      info=paste(ca$id,'exclusions'))
    expect_identical(length(r$states),length(v))
    expect_identical(r$states[which(inp$missing)],rep('missing_input',sum(inp$missing)),
      info=paste(ca$id,'masked rows'))
    expect_true(all(r$weights[inp$missing]==1),info=paste(ca$id,'masked weights'))
    expect_true(all(is.finite(r$weights)&r$weights>=0&r$weights<=1),
      info=paste(ca$id,'weight range'))
    expect_identical(names(r$weights),NULL)
    full<-which(r$states=='full');zero<-which(r$states=='zero');ramp<-which(r$states=='ramp')
    expect_true(all(r$weights[full]==1),info=paste(ca$id,'full state'))
    expect_true(all(r$weights[zero]==0),info=paste(ca$id,'zero state'))
    expect_true(all(r$weights[ramp]>0&r$weights[ramp]<1),info=paste(ca$id,'ramp state'))
    again<-raw(v,ca$frequency,ca$start_period,'legacy',ca$backcasts,ca$forecasts,
      ca$exclude,ca$lsig,ca$usig,ca$multiplicative)
    expect_identical(again,r,info=paste(ca$id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(ca$id,'input untouched'))
    acc<-raw(v,ca$frequency,ca$start_period,'accurate',ca$backcasts,ca$forecasts,
      ca$exclude,ca$lsig,ca$usig,ca$multiplicative)
    expect_identical(acc$status,'available',info=paste(ca$id,'accurate status'))
    expect_true(all(is.finite(acc$weights)&acc$weights>=0&acc$weights<=1),
      info=paste(ca$id,'accurate range'))
    expect_true(all(acc$weights[inp$missing]==1),info=paste(ca$id,'accurate masked'))
    expect_false('full_below_dispersion_cutoff' %in% acc$states,
      info=paste(ca$id,'accurate has no legacy cutoff state'))
    expect_true(acc$passes %in% c(1L,2L),info=paste(ca$id,'accurate passes'))
  }
  # Public wrapper: explicit policy, ts domain and tidy provenance frames.
  series<-stats::ts(c(1,2,3,4,1.2,2.4,3.1,4.4,0.9,2.1,3.3,4.2,7.5,2.2,3.4,4.1),
    start=c(2020,1),frequency=4)
  out<-pub(series,'legacy')
  expect_true(inherits(out,'persephone_native_x11_sigma'),info='wrapper class')
  expect_identical(out$frequency,4L)
  expect_identical(out$start_period,1L)
  expect_identical(names(out$weights),as.character(seq_len(16)))
  expect_identical(ncol(out$windows),5L)
  expect_true(all(out$windows$start>=1&out$windows$start+out$windows$length-1<=16))
  dflt<-pub(series) # omitted policy: the maintainer 2026-11-17 legacy default
  expect_identical(dflt$policy,'legacy',info='legacy default echo')
  expect_true(bits(dflt$weights,out$weights),info='legacy default is the explicit legacy run')
  expect_identical(dflt$windows,out$windows)
  expect_error(raw(as.numeric(1:12),4,1,'modern',0L,0L,FALSE,1.5,2.5,FALSE),'legacy or accurate')
  expect_error(raw(c(1,2,3,4,5,6,7,8,9,10,11,Inf),4,1,'legacy',0L,0L,FALSE,1.5,2.5,FALSE),
    'finite or missing')
  expect_error(raw(as.numeric(1:12),4,1,'legacy',-1L,0L,FALSE,1.5,2.5,FALSE),'nonnegative')
  expect_error(raw(as.numeric(1:12),4,1,'legacy',0L,0L,NA,1.5,2.5,FALSE),'TRUE or FALSE')
  expect_error(raw(as.numeric(1:12),4,1,'legacy',0L,0L,FALSE,0.5,2.5,FALSE))
  expect_error(raw(as.numeric(1:12),4,1,'legacy',0L,0L,FALSE,2.5,1.5,FALSE))
})
