source("helpers/fixtures.R", local = TRUE)
local({
  rawp<-getFromNamespace('.native_x13_preadjust','rkalos')
  rawr<-getFromNamespace('.native_x13_restore','rkalos')
  wpre<-getFromNamespace('native_x13_preadjust','rkalos')
  wres<-getFromNamespace('native_x13_restore','rkalos')
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])
  # The pinned core-qualified protocol (tests/test_x13_allocation.cpp): under a
  # log model the twelve buckets are stored through backTransform's exp, a
  # transcendental boundary between std::exp and Math.exp (1-ulp hits), and are
  # compared at relative tolerance; a1 and out are pure algebra and exact.
  close15<-function(a,b) {a<-as.vector(a);b<-as.vector(b)
    length(a)==length(b) && all(is.na(a)==is.na(b)) &&
      all(abs(a[!is.na(a)]-b[!is.na(b)])<=1e-15*pmax(1,abs(b[!is.na(b)])))}
  cases<-read.csv('fixtures/native/x13-allocation/cases.csv')
  inputs<-read.csv('fixtures/native/x13-allocation/inputs.csv')
  vars<-read.csv('fixtures/native/x13-allocation/vars.csv')
  columns<-read.csv('fixtures/native/x13-allocation/columns.csv')
  tables<-read.csv('fixtures/native/x13-allocation/tables.csv')
  expect_identical(nrow(cases),7L)
  # The boundary vocabulary of the allocation contract, mapped from the
  # capture's pinned Java spellings.
  classes<-c(outlier='outlier',days='days_related',holiday='moving_holiday',user='user')
  effects<-c(Undefined='undefined',Series='series',Trend='trend',Seasonal='seasonal',
    SeasonallyAdjusted='seasonally_adjusted',Irregular='irregular',CalendarEffect='calendar_effect')
  modes<-c(Additive='additive',Multiplicative='multiplicative',LogAdditive='log_additive',
    PseudoAdditive='pseudo_additive')
  anames<-c('a1','a6','a7','a8','a8t','a8s','a8i','a9','a9u','a9cal','a9sa','a9ser','out')
  fnames<-c('d11final','d12final','d13final','d16','d18','e1','e2','e3','e11')
  input_of<-function(id) { rows<-which(inputs$case==id)
    inputs$value[rows][order(as.integer(inputs$row[rows]))] }
  table_of<-function(id,nm) { rows<-which(tables$case==id & tables$table==nm)
    tables$value[rows][order(as.integer(tables$row[rows]))] }
  variables_of<-function(id) {
    vr<-vars[vars$case==id,];out<-list()
    for(k in sort(unique(as.integer(vr$index)))) {
      rows<-vr[as.integer(vr$index)==k,]
      keep<-which(columns$case==id & as.integer(columns$index)==k)
      cols<-list()
      for(cc in sort(unique(as.integer(columns$coef[keep])))) {
        w<-keep[as.integer(columns$coef[keep])==cc]
        cols[[length(cols)+1]]<-columns$value[w][order(as.integer(columns$row[w]))]
      }
      out[[length(out)+1]]<-list(class=classes[[as.character(rows$class[1])]],
        effect=effects[[as.character(rows$effect[1])]],
        coefficient=rows$coefficient,columns=cols)
    }
    out
  }
  for(i in seq_len(nrow(cases))) {
    ca<-cases[i,];id<-ca$case;n<-as.integer(ca$n);f<-as.integer(ca$frequency)
    lg<-as.logical(ca$log)
    v<-input_of(id);before<-serialize(v,NULL)
    variables<-variables_of(id)
    r<-rawp(v,lg,variables)
    expect_identical(r$type,'x13_preadjustment',info=paste(id,'type'))
    expect_identical(r$status,'available',info=paste(id,'status'))
    expect_identical(r$log,lg,info=paste(id,'log echo'))
    expect_identical(r$length,n,info=paste(id,'length echo'))
    expect_identical(names(r$tables),anames,info=paste(id,'table set'))
    expect_identical(as.vector(r$tables$a1),v,info=paste(id,'a1 is the given series'))
    expect_true(bits(r$tables$out,table_of(id,'out')),info=paste(id,'out exact'))
    buckets<-anames[2:12]
    for(nm in buckets)
      if(lg) expect_true(close15(r$tables[[nm]],table_of(id,nm)),
        info=paste(id,nm,'exp boundary at 1e-15'))
      else expect_true(bits(r$tables[[nm]],table_of(id,nm)),info=paste(id,nm,'exact'))
    expect_identical(serialize(v,NULL),before,info=paste(id,'input untouched'))
    # The finals seam over the captured X11Results tables of the same run:
    # pure pinned algebra, exact for every mode.
    m<-modes[[as.character(ca$mode_effective)]]
    s<-rawr(r,table_of(id,'d10'),table_of(id,'d11'),table_of(id,'d12'),
      table_of(id,'d13'),table_of(id,'c17'),m)
    expect_identical(s$type,'x13_finals',info=paste(id,'finals type'))
    expect_identical(s$status,'available',info=paste(id,'finals status'))
    expect_identical(s$mode,m,info=paste(id,'effective mode echo'))
    expect_identical(names(s$tables),fnames,info=paste(id,'finals table set'))
    for(nm in fnames)
      expect_true(bits(s$tables[[nm]],table_of(id,nm)),info=paste(id,nm,'exact'))
    # The seams refuse a broken coupling instead of papering over it.
    if(i==1L) {
      expect_error(rawr(list(),table_of(id,'d10'),table_of(id,'d11'),table_of(id,'d12'),
        table_of(id,'d13'),table_of(id,'c17'),m),'preadjustment result object')
      broken<-r;broken$status<-"tampered"
      expect_error(rawr(broken,table_of(id,'d10'),table_of(id,'d11'),table_of(id,'d12'),
        table_of(id,'d13'),table_of(id,'c17'),m),'available preadjustment')
    }
  }
  # The wrapper path: ts decoration on the input calendar, faithful chaining
  # and the documented refusals (qualified at the first monthly case).
  local({
    id<-'m_alloc';ca<-cases[cases$case==id,]
    v<-ts(input_of(id),start=c(2001,1),frequency=12)
    variables<-variables_of(id)
    w<-wpre(v,variables,as.logical(ca$log))
    expect_true(inherits(w,'persephone_native_x13_preadjustment'))
    expect_identical(stats::start(w$tables$a1),c(2001,1))
    expect_identical(stats::frequency(w$tables$a1),12)
    expect_identical(names(w$row_basis),anames)
    expect_identical(names(w$available_rows),anames)
    w2<-wres(w,as.vector(table_of(id,'d10')),
      as.vector(table_of(id,'d11')),as.vector(table_of(id,'d12')),
      as.vector(table_of(id,'d13')),as.vector(table_of(id,'c17')),'additive')
    expect_true(inherits(w2,'persephone_native_x13_finals'))
    expect_true(bits(w2$tables$d11final,table_of(id,'d11final')))
    expect_identical(stats::start(w2$tables$d11final),c(2001,1))
    expect_error(wpre(as.numeric(v),variables,TRUE),'explicit frequency')
    expect_error(wres(list(1),1,1,1,1,1),'preadjustment result object')
  })
  caps<-getFromNamespace('capabilities','rkalos')
  expect_true('x13_deterministic_allocation' %in% caps()$primitives)
  expect_true(!caps()$fit)
})
