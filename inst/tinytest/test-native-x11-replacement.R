source("helpers/fixtures.R", local = TRUE)
local({
  sigma<-getFromNamespace('.native_x11_sigma_weights','rkalos')
  raw<-getFromNamespace('.native_x11_replacements','rkalos')
  pub<-getFromNamespace('native_x11_replacements','rkalos')
  source_fixture('fixtures/native/x11-replacement-cases.R',local=TRUE)
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  kindmap<-c('none','four_neighbors','period_average','outside_span')
  manifest<-read.csv('fixtures/native/x11-replacement/cases.csv')
  inputs<-read.csv('fixtures/native/x11-replacement/inputs.csv')
  weights<-read.csv('fixtures/native/x11-replacement/weights.csv')
  repls<-read.csv('fixtures/native/x11-replacement/replacements.csv')
  kinds<-read.csv('fixtures/native/x11-replacement/kinds.csv')
  neigh<-read.csv('fixtures/native/x11-replacement/neighbors.csv')
  avgs<-read.csv('fixtures/native/x11-replacement/averages.csv')
  factors<-read.csv('fixtures/native/x11-replacement/factors.csv')
  corrected<-read.csv('fixtures/native/x11-replacement/corrected.csv')
  cases<-x11_replacement_cases()
  expect_identical(sort(names(cases)),sort(unique(manifest$case)))
  for(ca in cases) {
    inp<-x11_replacement_input(ca);v<-inp$values;v[inp$missing]<-NA_real_
    before<-serialize(v,NULL)
    pinned<-manifest[manifest$case==ca$id,]
    expect_identical(as.integer(pinned$n),length(v))
    rows<-order(as.integer(inputs$row[inputs$case==ca$id]))
    iv<-inputs$value[inputs$case==ca$id][rows];im<-inputs$missing[inputs$case==ca$id][rows]
    expect_identical(is.na(v),im,info=paste(ca$id,'input mask'))
    expect_true(bits(v[!im],iv[!im]),info=paste(ca$id,'input values'))
    # Chain: the final item 22b.3a legacy weights are exactly the fixture input.
    chained<-sigma(v,ca$frequency,ca$start_period,'legacy',ca$backcasts,ca$forecasts,
      ca$exclude,ca$lsig,ca$usig,ca$multiplicative)
    rows<-order(as.integer(weights$row[weights$case==ca$id]))
    w<-weights$weight[weights$case==ca$id][rows]
    expect_true(bits(chained$weights,w),info=paste(ca$id,'chained sigma weights'))
    r0<-raw(v,ca$frequency,w,ca$backcasts,ca$forecasts,FALSE,ca$multiplicative)
    r1<-raw(v,ca$frequency,w,ca$backcasts,ca$forecasts,TRUE,ca$multiplicative)
    expect_identical(r0$status,'available',info=paste(ca$id,'status'))
    for(mode in 0:1) {
      rr<-if(mode==0) r0 else r1
      sel<-which(repls$case==ca$id & repls$mode==mode)
      e<-repls$replacement[sel][order(as.integer(repls$row[sel]))]
      expect_true(bits(rr$replacements,e),info=paste(ca$id,'replacements mode',mode))
      selv<-which(kinds$case==ca$id & kinds$mode==mode)
      k<-kinds$kind[selv][order(as.integer(kinds$row[selv]))]
      expect_identical(rr$kinds,unname(kindmap[k+1L]),
        info=paste(ca$id,'kinds mode',mode))
      # Cancellation: unit-weight rows keep no replacement and the mode default.
      outside<-mode==1&(seq_along(v)<=ca$backcasts|seq_along(v)>length(v)-ca$forecasts)
      kept<-which(w==1&!outside)
      expect_true(all(is.na(rr$replacements[kept])),info=paste(ca$id,'unit rows unreplaced'))
      expect_identical(rr$kinds[kept],rep('none',length(kept)),info=paste(ca$id,'unit row kinds'))
      # Period-average fallbacks equal the phase average (both NaN allowed).
      fb<-which(rr$kinds=='period_average')
      avg<-rr$period_averages[((fb-1L) %% ca$frequency)+1L]
      expect_true(bits(rr$replacements[fb],avg),info=paste(ca$id,'fallback equals average'))
    }
    # Provenance belongs to the untrimmed (mode 0) pass.
    sel<-which(neigh$case==ca$id & !is.na(neigh$slot))
    mat<-matrix(-1,ncol=4,nrow=length(v))
    if(length(sel)) mat[cbind(as.integer(neigh$row[sel])+1L,as.integer(neigh$slot[sel])+1L)]<-
      neigh$position[sel]
    expect_identical(r0$neighbors,unname(mat),info=paste(ca$id,'neighbor provenance'))
    expect_identical(r0$neighbor_counts,ifelse(rowSums(mat!=-1)>0,4L,0L),
      info=paste(ca$id,'neighbor counts'))
    blends<-which(r0$kinds=='four_neighbors')
    expect_true(all(r0$neighbor_counts[blends]==4L),info=paste(ca$id,'blend counts'))
    for(row in blends) expect_true(all(w[r0$neighbors[row,]+1L]==1),
      info=paste(ca$id,'neighbors fully weighted'))
    sel<-which(avgs$case==ca$id)
    a<-avgs$average[sel][order(as.integer(avgs$phase[sel]))]
    expect_true(bits(r0$period_averages,a),info=paste(ca$id,'period averages'))
    expect_true(bits(r1$period_averages,a),info=paste(ca$id,'mode 1 averages'))
    sel<-which(factors$case==ca$id)
    f<-factors$factor[sel][order(as.integer(factors$row[sel]))]
    expect_true(bits(r0$correction_factors,f),info=paste(ca$id,'correction factors'))
    sel<-which(corrected$case==ca$id)
    c2<-corrected$value[sel][order(as.integer(corrected$row[sel]))]
    expect_true(bits(r0$corrected,c2),info=paste(ca$id,'corrected series'))
    expect_identical(r0$fallbacks,length(which(r0$kinds=='period_average')),
      info=paste(ca$id,'fallback count'))
    # applyCorrections identity: corrected equals the effective value with the
    # non-NaN replacements substituted.
    effective<-ifelse(inp$missing,NA_real_,v)
    manual<-ifelse(is.na(r0$replacements),effective,r0$replacements)
    expect_true(bits(r0$corrected,manual),info=paste(ca$id,'apply identity'))
    again<-raw(v,ca$frequency,w,ca$backcasts,ca$forecasts,TRUE,ca$multiplicative)
    expect_identical(again,r1,info=paste(ca$id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(ca$id,'input untouched'))
  }
  # Public wrapper: 1-based neighbor positions with NA padding and series class.
  w<-rep(1,16);w[c(1,13)]<-0
  series<-stats::ts(c(1,2,3,4,1.2,2.4,3.1,4.4,0.9,2.1,3.3,4.2,7.5,2.2,3.4,4.1),
    start=c(2020,1),frequency=4)
  out<-pub(series,w)
  expect_true(inherits(out,'persephone_native_x11_replacement'),info='wrapper class')
  expect_identical(nrow(out$neighbor_positions),16L)
  expect_identical(out$kinds[1],'period_average')
  expect_true(all(is.na(out$neighbor_positions[1,])))
  expect_error(pub(series,rep(1,15)),'series length')
  expect_error(pub(series,c(rep(1,15),NA_real_)),'finite')
})
