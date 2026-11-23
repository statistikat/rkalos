source("helpers/fixtures.R", local = TRUE)
local({
  fit<-rkalos::native_x13_fit
  raw<-getFromNamespace('.native_x13_fit','rkalos')
  bits<-function(a,b) {a<-as.vector(a);b<-as.vector(b)
    identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])}
  # The pinned tolerance policy for the log-model fit cases (docs/contracts/
  # x13-fit-contract.md, mirroring the core corpus): a-tables buckets and the
  # results/c17 tables cross exp-class boundaries at 1e-15/1e-12 relative;
  # the finals triplet at 1e-9 (std::exp vs Math.exp 1-ulp hits). Everything
  # pure-algebraic stays bit-exact, and every non-log case is bit-exact.
  closen<-function(a,b,tol) {a<-as.vector(a);b<-as.vector(b)
    length(a)==length(b) && all(is.na(a)==is.na(b)) &&
      all(abs(a[!is.na(a)]-b[!is.na(b)])<=tol*pmax(1,abs(b[!is.na(b)])))}
  cases<-read.csv('fixtures/native/x13-fit/cases.csv')
  inputs<-read_numeric_fixture('fixtures/native/x13-fit/inputs.csv')
  vars<-read_numeric_fixture('fixtures/native/x13-fit/vars.csv', "coefficient")
  columns<-read_numeric_fixture('fixtures/native/x13-fit/columns.csv')
  tables<-read_numeric_fixture('fixtures/native/x13-fit/tables.csv')
  horizons<-read.csv('fixtures/native/x13-fit/horizons.csv')
  fast<-read_numeric_fixture('fixtures/native/x13-fit/x13fast.csv')
  expect_identical(nrow(cases),11L)
  classes<-c(outlier='outlier',days='days_related',holiday='moving_holiday',user='user')
  effects<-c(Undefined='undefined',Series='series',Trend='trend',Seasonal='seasonal',
    SeasonallyAdjusted='seasonally_adjusted',Irregular='irregular',CalendarEffect='calendar_effect')
  modes<-c(Additive='additive',Multiplicative='multiplicative',LogAdditive='log_additive',
    PseudoAdditive='pseudo_additive')
  fixed<-c('f_add_full','f_add_s3x1','f_add_s3x15','f_add_stable','f_mult_forced',
    'f_log','f_logla','f_logpa','f_empty')
  expect_identical(sum(cases$route=='fixed'),9L)
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
  requests_of<-function(ca) { y<-as.integer(strsplit(as.character(ca$filters),',')[[1]])
    ifelse(y==0L,'stable',sprintf('s3x%d',y)) }
  for(i in seq_len(nrow(cases))) {
    ca<-cases[i,];id<-ca$case;n<-as.integer(ca$n);f<-as.integer(ca$frequency)
    lg<-as.logical(ca$log);v<-input_of(id);req<-requests_of(ca)
    if(id %in% fixed) {
      r<-fit(v,variables_of(id),lg,req,as.integer(ca$henderson),'legacy','jdplus',
        f,as.integer(ca$first_period),as.numeric(ca$lower),as.numeric(ca$upper),
        modes[[as.character(ca$mode_requested)]],as.integer(ca$backcast),as.integer(ca$forecast))
      expect_identical(r$status,'available',info=paste(id,'fit available'))
      expect_true(inherits(r,'persephone_native_x13_fit'))
      expect_true(r$fitted,info=paste(id,'fitted flag'))
      # The provenance echo with the pinned updateSpec forcing.
      expect_identical(r$type,'x13_fit');expect_identical(r$schema_version,1L)
      expect_identical(r$frequency,f);expect_identical(r$start_period,as.integer(ca$first_period))
      expect_identical(r$log,lg)
      expect_identical(as.character(r$seasonal),req,info=paste(id,'seasonal echo'))
      expect_identical(r$trend_length,as.integer(ca$henderson),info=paste(id,'trend echo'))
      expect_true(r$sigma_set);expect_identical(r$sigma,'legacy')
      expect_identical(r$lower_sigma,as.numeric(ca$lower))
      expect_identical(r$upper_sigma,as.numeric(ca$upper))
      expect_identical(r$msr_policy,'jdplus')
      expect_identical(r$mode_requested,modes[[as.character(ca$mode_requested)]],
        info=paste(id,'mode requested'))
      expect_identical(r$mode_effective,modes[[as.character(ca$mode_effective)]],
        info=paste(id,'mode forcing'))
      expect_identical(r$backcast,0L);expect_identical(r$forecast,0L)
      expect_identical(r$scope,'restricted_fixed_model_x13_fit_only')
      # The composition: preadjustment, the full selection run over out, the
      # results tables as finals consumed them and the pinned finals.
      expect_identical(r$preadjustment$status,'available',info=paste(id,'pre ran'))
      expect_identical(as.vector(r$preadjustment$tables$a1),v,info=paste(id,'a1 given series'))
      # In the fit composition out feeds the X-11 run and the capture carries
      # the exp boundary: the core corpus judges it at 1e-15 under a log model.
      if(lg) expect_true(closen(r$preadjustment$tables$out,table_of(id,'out'),1e-15),
        info=paste(id,'out at 1e-15'))
      else expect_true(bits(r$preadjustment$tables$out,table_of(id,'out')),
        info=paste(id,'out exact'))
      ab<-names(r$preadjustment$tables)[2:12]
      for(nm in ab)
        if(lg) expect_true(closen(r$preadjustment$tables[[nm]],table_of(id,nm),1e-15),
          info=paste(id,nm,'a-bucket at 1e-15'))
        else expect_true(bits(r$preadjustment$tables[[nm]],table_of(id,nm)),
          info=paste(id,nm,'a-bucket exact'))
      expect_identical(r$selection$type,'x11_selection_pass',info=paste(id,'selection record'))
      expect_identical(r$selection$scope,'restricted_fit_selection_run',info=paste(id,'selection scope'))
      expect_identical(r$selection$status,'available',info=paste(id,'selection ran'))
      expect_identical(names(r$selection$requests),as.character(seq_len(f)),
        info=paste(id,'named requests'))
      expect_identical(r$selection$mode,modes[[as.character(ca$mode_effective)]],
        info=paste(id,'selection ran the forced mode'))
      expect_identical(r$selection$b$type,'x11_b_pass')
      expect_identical(r$selection$c$type,'x11_c_pass')
      expect_identical(r$selection$d$type,'x11_d_pass')
      for(nm in c('d10','d11','d12','d13','c17'))
        if(lg) expect_true(closen(r$results[[nm]],table_of(id,nm),1e-12),
          info=paste(id,nm,'results table at 1e-12'))
        else expect_true(bits(r$results[[nm]],table_of(id,nm)),info=paste(id,nm,'results exact'))
      expect_true(grepl(if(modes[[as.character(ca$mode_effective)]]=='log_additive')
          'log-additive' else 'raw step',r$results$note),info=paste(id,'stage note'))
      trio<-c('finalSA','finalTrend','finalIrregular')
      ftab<-c('d11final','d12final','d13final')
      for(j in seq_along(trio)) {
        if(lg) expect_true(closen(r[[trio[j]]],table_of(id,ftab[j]),1e-9),
          info=paste(id,ftab[j],'finals at 1e-9'))
        else expect_true(bits(r[[trio[j]]],table_of(id,ftab[j])),info=paste(id,ftab[j],'exact'))
        expect_identical(r[[trio[j]]],r$finals$tables[[ftab[j]]],
          info=paste(id,'triplet is the finals table'))
      }
      for(nm in names(r$finals$tables))
        if(lg) expect_true(closen(r$finals$tables[[nm]],table_of(id,nm),1e-9),
          info=paste(id,nm,'finals table at 1e-9'))
        else expect_true(bits(r$finals$tables[[nm]],table_of(id,nm)),info=paste(id,nm,'finals exact'))
      if(id=='f_add_full') { # the decorated wrapper path on a ts calendar
        w<-fit(ts(v,start=c(2005,3),frequency=12),variables_of(id),lg,req,13L,'legacy')
        expect_true(inherits(w,'persephone_native_x13_fit'))
        expect_identical(w$status,'available')
        expect_identical(stats::start(w$finalSA),c(2005,3))
        expect_identical(stats::start(w$preadjustment$tables$a1),c(2005,3))
        expect_identical(stats::start(w$results$d10),c(2005,3))
      }
    } else {
      # f_h6 / f_b3h3: the pinned engine routes any nonzero horizon through
      # the ExactArimaForecasts extension (the horizons.csv facts); the fit
      # answers the explicit zero-horizon refusal with the full echo instead.
      hr<-horizons[horizons$case==id,]
      expect_identical(as.integer(hr$backcast),as.integer(ca$backcast))
      expect_identical(as.integer(hr$forecast),as.integer(ca$forecast))
      expect_false(as.integer(hr$out_len)==n,info=paste(id,'pinned engine extended out'))
      expect_identical(as.integer(hr$d11final_len),n)
      r<-fit(v,variables_of(id),lg,req,as.integer(ca$henderson),'legacy','jdplus',
        f,as.integer(ca$first_period),as.numeric(ca$lower),as.numeric(ca$upper),
        modes[[as.character(ca$mode_requested)]],as.integer(ca$backcast),as.integer(ca$forecast))
      expect_identical(r$status,'unsupported:x13 fit requires zero horizons',
        info=paste(id,'refused verbatim'))
      expect_false(r$fitted)
      expect_null(r$preadjustment,info=paste(id,'refused before the seam'))
      expect_null(r$selection)
      expect_identical(r$backcast,as.integer(hr$backcast),info=paste(id,'horizon echoed'))
      expect_identical(r$forecast,as.integer(hr$forecast),info=paste(id,'horizon echoed'))
      expect_identical(r$mode_effective,modes[[as.character(ca$mode_effective)]],
        info=paste(id,'echo complete'))
    }
  }
  # The fail-closed guards in the pinned order, mirrored (never invented):
  # fixed seasonal -> Henderson length -> zero horizons -> required sigma.
  local({
    ca<-cases[cases$case=='f_add_stable',];id<-ca$case
    v<-input_of(id);req<-requests_of(ca);va<-variables_of(id)
    call<-function(seasonal=req,trend=as.integer(ca$henderson),sigma='legacy',
        backcast=0L,forecast=0L)
      fit(v,va,as.logical(ca$log),seasonal,trend,sigma,'jdplus',as.integer(ca$frequency),
        as.integer(ca$first_period),as.numeric(ca$lower),as.numeric(ca$upper),
        'additive',backcast,forecast)
    expect_identical(call(seasonal=rep('x11default',12L))$status,
      'unsupported:x13 fit fixed seasonal requests only')
    expect_identical(call(seasonal=c('msr',req[-1]))$status,
      'unsupported:x13 fit fixed seasonal requests only')
    expect_identical(call(seasonal=rep('msr',12L),trend=0L,sigma=NULL)$status,
      'unsupported:x13 fit fixed seasonal requests only',info='seasonal guard first')
    expect_identical(call(trend=0L)$status,
      'unsupported:x13 fit fixed Henderson trend length')
    expect_identical(call(trend=0L)$status,
      'unsupported:x13 fit fixed Henderson trend length')
    expect_identical(call(backcast=-2L)$status,
      'unsupported:x13 fit requires zero horizons',info='negative horizons too')
    # The maintainer decision of 2026-11-17: the wrapper applies the product
    # default "legacy" to NULL/omitted sigma — bit-identical to explicit legacy.
    no_sig<-call(sigma=NULL)
    expect_identical(no_sig$status,'available',info='NULL sigma gets the legacy default')
    expect_true(no_sig$sigma_set);expect_identical(no_sig$sigma,'legacy')
    explicit<-call(sigma='legacy')
    for(nm in c('finalSA','finalTrend','finalIrregular'))
      expect_true(bits(no_sig[[nm]],explicit[[nm]]),info=paste('legacy default is explicit legacy:',nm))
    for(nm in c('d10','d11','d12','d13','c17'))
      expect_true(bits(no_sig$results[[nm]],explicit$results[[nm]]),
        info=paste('legacy default results equal:',nm))
    expect_identical(no_sig$preadjustment$status,'available')
    # The core keeps its explicitness contract: through the RAW entry a NULL
    # policy still answers the pinned fail-closed status, and the core guard
    # order (seasonal, trend, horizons, sigma) keeps judging it last.
    expect_identical(raw(v,va,FALSE,12L,1L,as.character(req),5L,NULL,1.5,2.5,
      'jdplus','additive',0L,0L)$status,
      'unsupported:x13 fit requires an explicit sigma policy')
    expect_identical(raw(v,va,FALSE,12L,1L,as.character(req),0L,NULL,1.5,2.5,
      'jdplus','additive',0L,0L)$status,
      'unsupported:x13 fit fixed Henderson trend length',info='trend guard before sigma')
    expect_identical(raw(v,va,FALSE,12L,1L,as.character(req),5L,NULL,1.5,2.5,
      'jdplus','additive',0L,6L)$status,
      'unsupported:x13 fit requires zero horizons',info='horizon guard before sigma')
    expect_identical(raw(v,va,FALSE,12L,1L,rep('msr',12L),0L,NULL,1.5,2.5,
      'jdplus','additive',0L,0L)$status,
      'unsupported:x13 fit fixed seasonal requests only',info='seasonal guard first at the core')
  })
  # The X13Fast cross-capture of f_empty (docs: the same pinned run through
  # the X13Kernel export) is byte-identical to the tables the fit reproduces.
  local({
    for(nm in unique(as.character(fast$table)))
      expect_true(bits(fast$value[fast$table==nm][order(as.integer(fast$row[fast$table==nm]))],
        table_of('f_empty',nm)),info=paste('x13fast',nm))
  })
  # Capability honesty: fit stays FALSE; the restricted fit is gated by its
  # own explicit fields (the persephone3 contract, IMPLEMENTATION_STATUS.md).
  local({
    caps<-getFromNamespace('capabilities','rkalos')
    c1<-caps()
    expect_identical(c1$schema_version,2L)
    expect_false(c1$fit)
    expect_true(c1$native_x13_fit_restricted)
    expect_true('x13_restricted_fixed_fit' %in% c1$primitives)
    expect_true(all(c('x11_pass_level_selection','x11_decomposition_modes') %in% c1$primitives))
    rf<-c1$restricted_fit
    expect_true(is.list(rf))
    expect_identical(names(rf),c('schema_version','native_x13_fit_restricted',
      'fixed_coefficients_only','fixed_seasonal_requests_only','fixed_trend_length_only',
      'sigma_policy_default','zero_horizon_only','monthly_only','frequency','modes',
      'refusals','entry_points'))
    expect_true(rf$fixed_coefficients_only && rf$fixed_seasonal_requests_only &&
      rf$fixed_trend_length_only && rf$zero_horizon_only && rf$monthly_only)
    expect_identical(rf$sigma_policy_default,'legacy',info='maintainer 2026-11-17')
    expect_identical(rf$frequency,12L)
    expect_identical(rf$modes,c('additive','multiplicative','log_additive','pseudo_additive'))
    expect_identical(rf$schema_version,2L)
    expect_true(grepl('restricted',c1$reason))
  })
})
