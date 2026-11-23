# Item 21-H M4: the automatic-spec X-13 fit at the R boundary, asserted
# against the pinned 13-case JVM oracle (docs/contracts/
# x13-automatic-pipeline-contract.md, the 21-H.2/H.3 qualification). The
# fixture CSVs are the capture itself (validation/evidence/native/
# x13-automatic, synced through tools/sync-core.py); the projections below
# mirror kalos tools/gen_x13_automatic_golden.py and the C++ oracle
# tests/test_x13_automatic_oracle.cpp one-for-one: decisions BIT-exact
# (status, the final description, variable identities, the differencing
# triples, the arma selections and panel membership, the ODM picks and cv,
# the effective-observation count), with the pinned bands only on the
# quantities the JVM displays (%.6g AICC pair 1e-5, tmean 1e-3, panel BIC2
# 1e-2, landing 2e-3, the ARMA landing coefficient block 1e-2, downstream
# tables 1e-6 relative). No tolerance is ever applied to a decision bit.
# The fail-closed matrix and the maintainer sigma default are asserted too:
# outside the envelope the surface answers the explicit status, NEVER a
# fallback to the restricted fit.
source("helpers/fixtures.R", local = TRUE)
local({
  fit<-getFromNamespace('native_x13_fit_automatic','rkalos')
  raw<-getFromNamespace('.native_x13_fit_automatic','rkalos')
  bits<-function(a,b) {a<-as.vector(a);b<-as.vector(b)
    identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])}
  # The oracle tier policy, transcribed from test_x13_automatic_oracle.cpp:
  # bit-equal (or both NaN) always passes; otherwise the relative deviation
  # against the want value must stay inside the pinned band (absolute when
  # the want is zero) and the got must be finite. The family worsts are
  # recorded so the run also measures the drift it allows.
  worst<-new.env()
  bands<-c('display-aicc'=1e-5,'diff-tmean'=1e-3,'arma-panel'=1e-2,
    'landing-coefficients'=2e-3,'arma-landing'=1e-2,'tables'=1e-6,
    'landing-stats'=2e-3)
  for(f in names(bands))assign(f,0,envir=worst)
  single<-function(got,want,fam) {a<-as.vector(as.numeric(got))
    b<-as.vector(as.numeric(want))
    if(length(a)!=length(b)||!all(is.na(a)==is.na(b)))return(FALSE)
    k<-!is.na(a)
    dev<-ifelse(abs(b[k])>0,abs(a[k]-b[k])/abs(b[k]),abs(a[k]-b[k]))
    if(length(dev)&&max(dev)>get(fam,worst))assign(fam,max(dev),envir=worst)
    all(is.finite(a[k]))&&all(dev<=bands[[fam]])}

  cases<-read.csv('fixtures/native/x13-automatic/cases.csv')
  inputs<-read.csv('fixtures/native/x13-automatic/inputs.csv',colClasses='character')
  events<-read.csv('fixtures/native/x13-automatic/oracle-events.csv',colClasses='character')
  facts<-read.csv('fixtures/native/x13-automatic/oracle-facts.csv',colClasses='character')
  expect_identical(nrow(cases),13L)

  fact<-function(id,key){v<-facts$value[facts$id==id&facts$key==key]
    if(length(v)==0)"" else v[1]}
  facts_of<-function(id,key)facts$value[facts$id==id&facts$key==key]
  reals<-function(field){t<-strsplit(field,';',fixed=TRUE)[[1]]
    fixture_doubles(t[nchar(t)>0])}
  modes<-c(Additive='additive',Multiplicative='multiplicative',
    LogAdditive='log_additive',PseudoAdditive='pseudo_additive',Undefined='additive')
  TYPE<-c(AO=0L,LS=1L,TC=2L,SO=3L)

  # The case request: the req-fact scalars verbatim plus the fixed decoder
  # constants (the C++ mirror_request), ub1/ub2 as the exact IEEE module-space
  # reciprocals the decoder hands the DifferencingModule (the golden
  # generator's repr(1/float(req)) rule).
  case_modules<-function(id) {
    line<-strsplit(fact(id,'req'),' ',fixed=TRUE)[[1]]
    req<-as.list(line[seq(2,length(line),2)])
    names(req)<-line[seq(1,length(line)-1,2)]
    iline<-inputs$line[inputs$id==id]
    inp<-stats::setNames(as.list(sub('^[^ ]* ','',iline)),sub(' .*','',iline))
    y<-reals(inp$values)
    expect_identical(length(y),as.integer(req$n),info=paste(id,'input length'))
    expect_identical(inp$id,id)
    otypes<-if(is.null(inp$otypes))c('AO','LS','TC','SO')
      else strsplit(inp$otypes,',',fixed=TRUE)[[1]]
    sy<-as.integer(strsplit(req$start,':',fixed=TRUE)[[1]])
    list(y=y,first_year=sy[1],start_period=sy[2],req=req,
      used=identical(req$outliers,'on'),
      types=if(identical(req$outliers,'on'))otypes else NULL)
  }
  fit_case<-function(id,sigma='accurate',span=NULL,used=NULL,...) {
    m<-case_modules(id);req<-m$req
    sigma_cut<-as.numeric(strsplit(req$sigma,';',fixed=TRUE)[[1]])
    if(is.null(used))used<-m$used
    types<-if(is.null(span)&&identical(used,m$used))m$types
      else if(used)c('AO','LS','TC','SO') else NULL
    if(!is.null(span))types<-c('AO','LS','TC','SO')
    a<-list(series=m$y,transform=list(aicc=as.numeric(req$aicdiff)),
      automodel=list(max_d=2L,max_bd=1L,ub1=1/as.numeric(req$ub1),
        ub2=1/as.numeric(req$ub2),cancel=as.numeric(req$cancel),precision=1e-7,
        models_count=5L,max_p=2L,max_q=2L,max_bp=1L,max_bq=1L,
        balanced=identical(req$balanced,'true'),mixed=identical(req$mixed,'true'),
        accept_white_noise=FALSE,eps=1e-7),
      outliers=list(used=used,types=types,span=span,
        max_round=as.integer(req$maxiter),max_outliers=30L,
        tcrate=as.numeric(req$tcrate)),
      ami=list(precision=1e-7,mean_precision=1e-5,va=0,reduce_va=as.numeric(req$predcv),
        ljung_box_limit=as.numeric(req$lblimit),check_mu=TRUE,accept_airline=FALSE,
        ur_limit=0.95),
      regression=list(mean_test_initial=1.96,mean_test_final=1,tsig=1),
      sigma=sigma,lower_sigma=sigma_cut[1],upper_sigma=sigma_cut[2],
      msr_policy='jdplus',mode=modes[[as.character(cases$mode_requested[
        match(id,cases$id)])]],frequency=12L,start_period=m$start_period,
      first_year=m$first_year)
    # ... overrides any entry of the assembled call (refusal probes reuse the
    # case module state with one knob turned).
    ov<-list(...);a[names(ov)]<-NULL;c(a,ov)->a;do.call(fit,a)
  }

  # The observer stream projected onto the decision kinds (the golden
  # generator's parse_events, one for one).
  parse_events<-function(id) {
    rows<-events[events$id==id,];out<-list();in_full<-FALSE;full_depth<-NA
    for(k in seq_len(nrow(rows))) {
      d<-rows$detail[k];kind<-rows$kind[k];msg<-rows$msg[k]
      if(kind=='push'&&msg=='full estimation') {in_full<-TRUE;full_depth<-rows$depth[k];next}
      if(kind=='pop'&&in_full&&rows$depth[k]==full_depth)in_full<-FALSE
      grab<-function(pat){m<-regexpr(pat,d,perl=TRUE)
        if(m<0)NA else as.numeric(substring(d,m,m+attr(m,'match.length')-1))}
      if(substr(d,1,17)=='obj AICC on logs=') {
        out[[length(out)+1]]<-list(kind=1L,
          logs=grab('(?<=logs=)[-+0-9.eE]+'),levels=grab('(?<=levels=)[-+0-9.eE]+'),
          pref=grab('(?<=preference\\)=)[-+0-9.eE]+'))
      } else if(substr(d,1,7)=='diff d=') {
        g<-regmatches(d,regexec('diff d=([0-9]+) bd=([0-9]+) mean=(true|false) tmean=([-+0-9.eE]+)',d))[[1]]
        out[[length(out)+1]]<-list(kind=2L,d=as.integer(g[2]),bd=as.integer(g[3]),
          mean=identical(g[4],'true'),tmean=as.numeric(g[5]))
      } else if(substr(d,1,9)=='arma sel=') {
        g<-regmatches(d,regexec('\\((-?[0-9]+),(-?[0-9]+)\\)\\((-?[0-9]+),(-?[0-9]+)\\)',d))[[1]]
        entries<-regmatches(d,gregexpr(
          '\\((-?[0-9]+),(-?[0-9]+)\\)\\((-?[0-9]+),(-?[0-9]+)\\):([-+0-9.eE]+)',d))[[1]]
        panel<-lapply(entries,function(e){
          q<-regmatches(e,regexec('\\((-?[0-9]+),(-?[0-9]+)\\)\\((-?[0-9]+),(-?[0-9]+)\\):([-+0-9.eE]+)',e))[[1]]
          list(p=as.integer(q[2]),q=as.integer(q[3]),bp=as.integer(q[4]),
            bq=as.integer(q[5]),bic=as.numeric(q[6]))})
        out[[length(out)+1]]<-list(kind=3L,p=as.integer(g[2]),q=as.integer(g[3]),
          bp=as.integer(g[4]),bq=as.integer(g[5]),panel=panel)
      } else if(substr(d,1,4)=='odm ') {
        head<-substring(d,5)
        types<-strsplit(strsplit(head,' ',fixed=TRUE)[[1]][1],',',fixed=TRUE)[[1]]
        entries<-regmatches(d,gregexpr('\\[([0-9]+);([0-9]+)\\]',d))[[1]]
        picks<-lapply(entries,function(e){
          q<-regmatches(e,regexec('\\[([0-9]+);([0-9]+)\\]',e))[[1]]
          # The picks carry the FACTORY index into this case's type list;
          # map to the global AO/LS/TC/SO code the driver speaks (and the
          # boundary exports 1-based positions).
          list(position=as.integer(q[2])+1L,type=as.integer(TYPE[types[as.integer(q[3])+1]]))})
        toks<-strsplit(trimws(head),' ')[[1]]
        out[[length(out)+1]]<-list(kind=4L,picks=picks,cv=as.numeric(toks[length(toks)]),
          types=types)
      } else if(in_full&&kind=='info'&&substr(d,1,5)=='obj (') {
        g<-regmatches(d,regexec('obj \\((-?[0-9]+),(-?[0-9]+),(-?[0-9]+)\\)\\((-?[0-9]+),(-?[0-9]+),(-?[0-9]+)\\)',d))[[1]]
        out[[length(out)+1]]<-list(kind=5L,p=as.integer(g[2]),d=as.integer(g[3]),
          q=as.integer(g[4]),bp=as.integer(g[5]),bd=as.integer(g[6]),bq=as.integer(g[7]))
      }
    }
    out
  }

  # The per-case replay against the pinned oracle, decision-bit exact.
  quick_res<-NULL
  for(i in seq_len(nrow(cases))) {
    ca<-cases[i,];id<-as.character(ca$id)
    cat('== case',id,'\n')
    res<-fit_case(id)
    if(id=='quick')quick_res<-res
    expect_true(inherits(res,'persephone_native_x13_automatic_fit'),info=paste(id,'class'))
    expect_true(res$fitted,info=paste(id,'fitted'))
    expect_identical(res$status,'available',info=paste(id,'status'))
    expect_identical(res$type,'x13_automatic_fit')
    expect_identical(res$schema_version,2L)
    expect_identical(res$scope,'automatic_spec_x13_fit_only')
    expect_identical(res$backend,'native')
    expect_identical(res$sigma_provenance,'accurate (explicit)',info=paste(id,'provenance'))
    # The provenance echo: the module state exactly as it was handed in.
    rq<-res$request
    expect_identical(rq$frequency,12L);expect_identical(rq$start_period,1L)
    expect_identical(rq$first_year,1990L)
    expect_identical(rq$x11$sigma,'accurate')
    expect_identical(rq$x11$mode,modes[[as.character(ca$mode_requested)]])
    expect_identical(as.character(rq$x11$seasonal),rep('msr',12L))
    expect_identical(rq$x11$trend_length,0L,info=paste(id,'automatic Henderson echo'))
    expect_false(any(unlist(rq$deferred)))

    # model-facts: the final description, bit-exact against cases.csv.
    md<-res$model
    expect_identical(md$log,identical(as.character(ca$log),'true'),info=paste(id,'log'))
    expect_identical(md$mean,identical(as.character(ca$mean),'true'),info=paste(id,'mean'))
    for(nm in c('p','d','q','bp','bd','bq'))
      expect_identical(md[[nm]],as.integer(ca[[nm]]),info=paste(id,nm))
    expect_identical(res$mode_effective,modes[[as.character(ca$mode_effective)]],
      info=paste(id,'mode effective'))
    m<-case_modules(id)
    sigma_cut<-as.numeric(strsplit(m$req$sigma,';',fixed=TRUE)[[1]])
    expect_identical(res$lower_sigma,sigma_cut[1],info=paste(id,'sigma lower echo'))
    expect_identical(res$upper_sigma,sigma_cut[2],info=paste(id,'sigma upper echo'))
    parse_var<-function(row){t<-strsplit(row,' ',fixed=TRUE)[[1]]
      # "<row> <name...> <ami> <regeffect>" (the generator's toks[-2]/toks[-1])
      list(row=as.integer(t[1]),name=paste(t[-c(1,length(t)-1,length(t))],collapse=' '),
        ami=t[length(t)-1],regeffect=t[length(t)])}
    vs<-lapply(facts_of(id,'var'),parse_var)
    ami<-Filter(function(v)identical(v$ami,'x13'),vs)
    # cases.csv nvars counts every description variable (the const column
    # included); the oracle's exp_nvars is the x13 (ami) subset the fit lands.
    expect_identical(length(ami)+sum(vapply(vs,function(v)identical(v$ami,'-'),TRUE)),
      as.integer(ca$nvars),info=paste(id,'nvars'))
    # The export carries the ARIMA-side variables; the const column (a mean,
    # not a variable) stays in the cases.csv nvars census only.
    expect_identical(length(md$variables$name),length(ami),info=paste(id,'nvars export'))
    regcoef<-reals(fact(id,'regcoef'))
    if(length(ami)) {
      # The golden generator's calendar-offset rule for the JVM variable name
      # ("AO (1994-01-01)"), +1 for the boundary's 1-based position.
      ym<-vapply(ami,function(v)regmatches(v$name,
        regexpr('[0-9]{4}-[0-9]{2}',v$name)),"")
      want_pos<-as.integer(substr(ym,1,4))*12L+as.integer(substr(ym,6,7))-
        (m$first_year*12L+m$start_period)+1L
      want_type<-vapply(ami,function(v)as.integer(TYPE[substr(v$name,1,2)]),integer(1))
      want_row<-vapply(ami,function(v)v$row,integer(1))
      want_name<-vapply(ami,function(v)v$name,"")
      want_eff<-vapply(ami,function(v)if(identical(v$regeffect,'-'))'' else v$regeffect,"")
      ow<-order(want_pos,want_type);og<-order(md$variables$position,md$variables$type_index)
      expect_true(bits(md$variables$position[og],want_pos[ow]),info=paste(id,'var positions'))
      expect_identical(md$variables$type_index[og],want_type[ow],info=paste(id,'var types'))
      expect_identical(md$variables$name[og],want_name[ow],info=paste(id,'var names'))
      expect_identical(md$variables$regeffect[og],want_eff[ow],info=paste(id,'var regeffects'))
      # The pair identity the separate sorts cannot prove on their own.
      expect_identical(sort(paste(md$variables$position[og],md$variables$type_index[og])),
        sort(paste(want_pos[ow],want_type[ow])),info=paste(id,'var pairs'))
      for(k in seq_along(ami))
        expect_true(single(md$variables$coefficient[og][k],regcoef[want_row[ow][k]+1L],
          'landing-coefficients'),info=paste(id,'var coefficient',k))
    }
    # The landed ARMA block (the measured 1e-2 coefficient-basin band).
    flat<-as.vector(unlist(md$coefficients))
    arma<-reals(fact(id,'arma'))
    expect_identical(md$n_coefficients,length(arma),info=paste(id,'narma'))
    expect_true(single(flat,arma,'arma-landing'),info=paste(id,'arma landing'))
    # Landing statistics: neff bit, loglik/AICC in the 2e-3 basin band.
    st<-regmatches(fact(id,'stats'),regexec(
      'nobs=([0-9]+) neff=([0-9]+) npar=([0-9]+) loglik=([-+0-9.eE]+) aicc=([-+0-9.eE]+)',
      fact(id,'stats')))[[1]]
    expect_identical(as.integer(md$effective_observations),as.integer(st[3]),info=paste(id,'neff'))
    expect_true(single(md$log_likelihood,as.numeric(st[5]),'landing-stats'),info=paste(id,'loglik'))
    expect_true(single(md$aicc,as.numeric(st[6]),'landing-stats'),info=paste(id,'aicc'))

    # The decision chain: the golden events against the exported records.
    stages<-vapply(res$chain,function(e)e$stage,"")
    evs<-parse_events(id)
    expect_identical(as.integer(ca$nevents),nrow(events[events$id==id,]),
      info=paste(id,'event stream length'))
    expect_identical(tail(stages,1),'final',info=paste(id,'chain ends at final'))
    idx_trans<-which(stages=='transformation')
    idx_refine<-which(grepl('\\.refine$',stages))
    idx_odm<-which(stages=='outliers0'|grepl('\\.outliers$',stages))
    aicc<-Filter(function(e)e$kind==1L,evs)
    diffs<-Filter(function(e)e$kind==2L,evs)
    armas<-Filter(function(e)e$kind==3L,evs)
    odms<-Filter(function(e)e$kind==4L,evs)
    full<-Filter(function(e)e$kind==5L,evs)
    expect_identical(length(aicc),length(idx_trans),info=paste(id,'aicc count'))
    expect_identical(length(diffs),length(idx_refine),info=paste(id,'refine count'))
    expect_identical(length(armas),length(idx_refine),info=paste(id,'arma count'))
    expect_identical(length(odms),length(idx_odm),info=paste(id,'odm count'))
    # kind 5: the cascade reduced-spec remark equals the final orders (the
    # observer prints it only inside the full-estimation frame: the rw case
    # in this corpus; an empty list is a valid, vacuous comparison).
    for(e in full)for(nm in c('p','d','q','bp','bd','bq'))
      expect_identical(md[[nm]],e[[nm]],info=paste(id,'cascade',nm))
    # kind 1: the %.6g AICC display pair (1e-5 display band).
    for(k in seq_along(aicc)) {
      e<-aicc[[k]];t<-res$chain[[idx_trans[k]]]$transformation
      expect_false(is.null(t),info=paste(id,'transformation payload'))
      expect_true(single(t$aicc_log,e$logs,'display-aicc'),info=paste(id,'aicc logs'))
      expect_true(single(t$aicc_level,e$levels,'display-aicc'),info=paste(id,'aicc levels'))
      expect_true(e$pref==-2,info=paste(id,'aicc correction echo'))
    }
    # kinds 2/3: the differencing triple and the arma selection + panel.
    for(k in seq_along(diffs)) {
      e<-diffs[[k]];rec<-res$chain[[idx_refine[k]]]
      expect_false(is.null(rec$refine),info=paste(id,'refine payload'))
      dd<-rec$refine$differencing
      expect_false(is.null(dd),info=paste(id,'differencing payload'))
      expect_identical(dd$d,e$d,info=paste(id,'diff d'))
      expect_identical(dd$seasonal_d,e$bd,info=paste(id,'diff bd'))
      expect_identical(dd$mean_correction,e$mean,info=paste(id,'diff mean'))
      expect_true(single(dd$t_mean,e$tmean,'diff-tmean'),info=paste(id,'diff tmean'))
      a<-armas[[k]]
      expect_identical(rec$refine$p,a$p,info=paste(id,'arma p'))
      expect_identical(rec$refine$q,a$q,info=paste(id,'arma q'))
      expect_identical(rec$refine$bp,a$bp,info=paste(id,'arma bp'))
      expect_identical(rec$refine$bq,a$bq,info=paste(id,'arma bq'))
      pm<-rec$refine$preferred_models
      expect_identical(length(pm$p),length(a$panel),info=paste(id,'panel size'))
      expect_true(bits(as.numeric(c(pm$p,pm$q,pm$bp,pm$bq)),as.numeric(c(
        vapply(a$panel,function(x)x$p,1),vapply(a$panel,function(x)x$q,1),
        vapply(a$panel,function(x)x$bp,1),vapply(a$panel,function(x)x$bq,1)))),
        info=paste(id,'panel membership'))
      expect_true(single(pm$bic,vapply(a$panel,function(x)x$bic,1),'arma-panel'),
        info=paste(id,'panel bic'))
    }
    # kind 4: the ODM rounds (picks bit-exact, cv bit-identical).
    for(k in seq_along(odms)) {
      e<-odms[[k]];rec<-res$chain[[idx_odm[k]]]
      expect_false(is.null(rec$outliers),info=paste(id,'odm payload'))
      got<-rec$outliers$decisions
      expect_identical(length(got$position),length(e$picks),info=paste(id,'odm count'))
      if(length(e$picks)) {
        # Pair identity per pick (independent sorts of the two components
        # would break the pairing; positions are 1-based on both sides).
        wp<-sort(vapply(e$picks,function(p)paste(p$position,p$type),""))
        expect_identical(sort(paste(got$position,got$type_index)),wp,
          info=paste(id,'odm pairs'))
      }
      expect_true(bits(as.numeric(rec$scalars[['cv']]),e$cv),info=paste(id,'odm cv'))
    }
    # The downstream composition tables (the measured 1e-6 table band).
    expect_identical(res$preadjustment$status,'available',info=paste(id,'preadjustment'))
    expect_identical(res$selection$status,'available',info=paste(id,'selection'))
    tbl<-list(out=res$preadjustment$tables$out,
      d11final=res$finals$tables$d11final,e11=res$finals$tables$e11)
    for(nm in names(tbl)) {
      v<-as.vector(tbl[[nm]]);w<-reals(fact(id,nm))
      expect_identical(length(v),length(w),info=paste(id,nm,'length'))
      expect_true(single(v,w,'tables'),info=paste(id,nm))
    }
    # The boundary decorations the restricted surface also carries.
    expect_true(inherits(res$preadjustment$tables$out,'ts'),info=paste(id,'ts tables'))
    expect_identical(names(res$preadjustment$row_basis),names(res$preadjustment$tables))
    expect_true(inherits(res$finalSA,'ts'),info=paste(id,'ts finalSA'))
    expect_true(inherits(res$series,'ts'),info=paste(id,'ts series'))
  }
  # The measured drift stays inside the pinned bands (recorded, never widened).
  for(f in names(bands))expect_true(get(f,worst)<=bands[[f]],
    info=paste('worst',f,get(f,worst)))
  for(f in names(bands))
    cat(sprintf('tier %-21s worst %.3g (band %.0g)\n',f,get(f,worst),bands[[f]]))

  # The capabilities contract: fit stays FALSE, the automatic fit reports its
  # own block in the rsa-template module shape, the restricted block untouched.
  local({
    caps<-getFromNamespace('capabilities','rkalos')
    c1<-caps()
    expect_identical(c1$schema_version,2L)
    expect_false(c1$fit)
    expect_true(c1$native_x13_fit_automatic)
    expect_true('x13_automatic_spec_fit' %in% c1$primitives)
    af<-c1$automatic_fit
    expect_true(is.list(af))
    expect_identical(names(af),c('schema_version','native_x13_fit_automatic','transform',
      'automodel','outliers','ami','regression','x11','monthly_only','frequency',
      'complete_finite_input_only','deferred','qualification','refusals','entry_points'))
    expect_identical(af$schema_version,3L,info='3L: x11 zero_horizon_only -> horizon_envelope (21-H.4)')
    expect_true(af$native_x13_fit_automatic)
    expect_identical(names(af$transform),c('auto_log_test','aicc_correction_required'))
    expect_true(af$transform$auto_log_test&&af$transform$aicc_correction_required)
    expect_identical(names(af$automodel),c('differencing','arma_grid'))
    expect_identical(af$outliers$types,c('AO','LS','TC','SO'))
    expect_true(grepl('full estimation domain',af$outliers$span))
    expect_identical(names(af$regression),c('mean_tests','final_tsig','user_detection'))
    expect_false(af$regression$user_detection)
    expect_identical(af$x11$sigma_policy_default,'legacy',info='maintainer 2026-11-17')
    expect_identical(af$x11$modes,c('additive','multiplicative','log_additive','pseudo_additive'))
    expect_true(af$x11$msr_requests_pinned&&af$x11$automatic_henderson)
    expect_null(af$x11$zero_horizon_only,info='the {0,0} marker retired with schema 3L')
    expect_true(grepl('21-H.4',af$x11$horizon_envelope)&&
      grepl('annual frequencies',af$x11$horizon_envelope)&&
      grepl('downstream composition only',af$x11$horizon_envelope),
      info='horizon_envelope states the qualified route (21-H.4)')
    expect_true(af$monthly_only);expect_identical(af$frequency,12L)
    expect_true(af$complete_finite_input_only)
    expect_true(grepl('21-I',af$deferred))
    expect_true(grepl('NEVER a fallback',af$refusals))
    expect_identical(af$entry_points,'native_x13_fit_automatic()')
    # The restricted-fit block keeps its exact schema-v2 shape (M3 contract).
    rf<-c1$restricted_fit
    expect_identical(names(rf),c('schema_version','native_x13_fit_restricted',
      'fixed_coefficients_only','fixed_seasonal_requests_only','fixed_trend_length_only',
      'sigma_policy_default','zero_horizon_only','monthly_only','frequency','modes',
      'refusals','entry_points'))
    expect_identical(rf$sigma_policy_default,'legacy')
  })

  # The request schema itself: NO input field carries fixed coefficients or
  # ARMA user variables; the deferred flag is the only user-variable knob and
  # it is refused, never approximated.
  local({
    fm<-names(formals(raw))
    expect_false(any(grepl('fixed|coefficient|arma_variable',fm)),
      info='no fixed-coefficient or ARMA-user-variable input field exists')
    expect_identical(sum(grepl('user_variables',fm)),1L)
    expect_identical(names(formals(fit)),c('series','transform','automodel','outliers',
      'ami','regression','sigma','lower_sigma','upper_sigma','msr_policy','mode',
      'backcast','forecast','frequency','start_period','first_year','calendar',
      'easter','trading_days','user_variables'))
  })

  # The fail-closed matrix, contract order (raw entries answer the pinned
  # statuses verbatim; every refusal keeps the provenance echo).
  local({
    m<-case_modules('quick');req<-m$req
    sigma_cut<-as.numeric(strsplit(req$sigma,';',fixed=TRUE)[[1]])
    raw0<-function(...) {
      a<-list(values=m$y,frequency=12L,start_period=1L,first_year=1990L,
        aicc_log_correction=-2,max_d=2L,max_bd=1L,ub1=1/as.numeric(req$ub1),
        ub2=1/as.numeric(req$ub2),cancel=as.numeric(req$cancel),differencing_precision=1e-7,
        models_count=5L,max_p=2L,max_q=2L,max_bp=1L,max_bq=1L,balanced=FALSE,
        mixed=identical(req$mixed,'true'),accept_white_noise=FALSE,arma_eps=1e-7,
        outliers_used=FALSE,otypes=character(0),span=NULL,outlier_max_round=30L,
        outlier_max_outliers=30L,tcrate=as.numeric(req$tcrate),precision=1e-7,mean_precision=1e-5,
        va=0,reduce_va=as.numeric(req$predcv),ljung_box_limit=as.numeric(req$lblimit),
        check_mu=TRUE,accept_airline=FALSE,ur_limit=0.95,mean_test_initial=1.96,
        mean_test_final=1,tsig=1,calendar=FALSE,easter=FALSE,trading_days=FALSE,
        user_variables=FALSE,seasonal=rep('msr',12L),sigma='accurate',
        lower_sigma=sigma_cut[1],upper_sigma=sigma_cut[2],msr_policy='jdplus',
        mode='additive',backcast=0L,forecast=0L)
      ov<-list(...);a[names(ov)]<-NULL;do.call(raw,c(a,ov))
    }
    r<-raw0(sigma=NULL)
    expect_identical(r$status,'unsupported:x13 fit requires an explicit sigma policy')
    expect_identical(r$schema_version,2L) # the echo is complete on a refusal too
    expect_identical(raw0(frequency=4L,seasonal=rep('msr',4L))$status,
      'unsupported:x13 automatic fit monthly series only')
    # The horizons are no gate (item 21-H.4 delivers the horizon route): a
    # horizon-bearing request flows through to the sigma gate exactly as a
    # zero-horizon one does (the core mirrors this in its t01_status_matrix).
    expect_identical(raw0(backcast=1L,sigma=NULL)$status,
      'unsupported:x13 fit requires an explicit sigma policy',
      info='horizon is no gate; sigma gate next')
    expect_identical(raw0(forecast=3L,sigma=NULL)$status,
      'unsupported:x13 fit requires an explicit sigma policy',
      info='forecast is no gate either')
    for(flag in c('calendar','easter','trading_days','user_variables')) {
      a<-list(FALSE,FALSE,FALSE,FALSE);names(a)<-c('calendar','easter','trading_days','user_variables')
      a[[flag]]<-TRUE
      expect_identical(do.call(raw0,a)$status,
        'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)',
        info=paste('deferred:',flag))
    }
    expect_identical(raw0(outliers_used=TRUE,otypes=c('AO','LS','TC','SO'),
      span=c(2,95))$status,
      'unsupported:x13 automatic fit spans the full series only')
    y<-m$y;y[10]<-NA_real_
    expect_identical(raw0(values=y)$status,
      'unsupported:x13 automatic fit requires a complete finite series')
    # The guard order of the matrix: monthly, finite, sigma, defer, span —
    # the first violation answers (the horizons are no gate since 21-H.4).
    expect_identical(raw0(values=y,frequency=4L,seasonal=rep('msr',4L))$status,
      'unsupported:x13 automatic fit monthly series only',info='monthly first')
    expect_identical(raw0(values=y,backcast=1L,sigma=NULL)$status,
      'unsupported:x13 automatic fit requires a complete finite series',info='finite beats later gates (horizon present)')
    expect_identical(raw0(backcast=1L,calendar=TRUE)$status,
      'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)',
      info='defer gate answers with a horizon present')
    # Shape violations throw instead of answering a status.
    expect_error(raw0(sigma='foolish'),'sigma weight policy',info='sigma policy throws')
    expect_error(raw0(mode='sideways'),'decomposition mode',info='mode throws')
    expect_error(raw0(msr_policy='nope'),'MSR iteration policy',info='msr policy throws')
    expect_error(raw0(otypes=c('AO','XX')),'AO, LS, TC or SO',info='otype throws')
    expect_error(raw0(outliers_used=TRUE,otypes=character(0)),'at least one type',
      info='empty otypes throws')
    expect_error(raw0(outlier_max_round=-1), "outlier_max_round")
    expect_error(raw0(outlier_max_outliers=1.5), "outlier_max_outliers")
    expect_error(raw0(span=c(NA_integer_, 4L)), "span")
    expect_error(raw0(span=c(3L, 2L)), "span")
    expect_error(raw0(forecast=-.Machine$integer.max), "extended series")
    expect_error(raw0(span=c(1,2,3)),'length-2',info='span shape throws')
    expect_error(raw0(start_period=13L),'1-based month index',info='start_period throws')
    expect_error(raw0(values=m$y[1:13]),'too short',info='short series throws')
    expect_error(raw0(values=numeric(0)),'non-empty',info='empty series throws')
  })

  # The wrapper: module validation up front, refusals as explicit statuses
  # (never a restricted fit), and the maintainer sigma default labelled in
  # the provenance while the raw entry keeps the core's explicitness.
  local({
    expect_error(fit_case('quick',transform=list()),'transform\\$aicc is required',
      info='transform validation')
    expect_error(fit_case('quick',ami=list(precision=1e-7,mean_precision=1e-5,va=0,
      reduce_va=0,ljung_box_limit=0.95,check_mu=TRUE,accept_airline=FALSE,
      ur_limit=0.95,extra=1)),'has no option',info='ami option validation')
    expect_error(fit_case('quick',automodel=list(max_d=1.5,max_bd=1L,ub1=1,ub2=1,
      cancel=0.1,precision=1e-7,models_count=5L,max_p=2L,max_q=2L,max_bp=1L,max_bq=1L,
      balanced=FALSE,mixed=TRUE,accept_white_noise=FALSE,eps=1e-7)),'whole number',
      info='automodel integer validation')
    expect_error(fit_case('quick',backcast=1.5),'whole number')
    expect_error(fit_case('quick',first_year=Inf),'whole number')
    expect_error(fit_case('quick',series=matrix(1:96, ncol=2)),'univariate')
    expect_error(fit_case('quick',series=ts(1:96, start=2000, frequency=12)),'calendar')
    quarterly <- fit_case('quick',frequency=.Machine$integer.max)
    expect_identical(quarterly$status, 'unsupported:x13 automatic fit monthly series only')
    expect_identical(quarterly$request$x11$seasonal, character())
    expect_error(fit_case('quick',sigma=42),'legacy',info='sigma label validation')
    # Item 21-H.4 horizon route at this boundary: accepted over the envelope
    # (capabilities()$automatic_fit$horizon_envelope), qualified against the
    # 4-case JVM horizon capture in test-native-x13-horizon.R. This block keeps
    # the shape pins and adds the boundary decorations the horizons drive: the
    # extended composition tables (out, d10..c17) lengthen by bc+fc months and
    # start bc months before the input calendar; the a1/finals calendar stays.
    r<-fit_case('quick',backcast=1,forecast=2)
    expect_identical(r$status,'available',info='horizon route answers (21-H.4)')
    expect_true(r$fitted)
    expect_false(is.null(r$model))
    expect_false(is.null(r$chain))
    expect_identical(as.integer(r$request$x11$backcast),1L)
    expect_identical(as.integer(r$request$x11$forecast),2L,info='horizon echo raw')
    expect_identical(r$sigma_provenance,'accurate (explicit)')
    expect_identical(length(r$preadjustment$tables$out),99L,info='out grows by bc+fc')
    expect_identical(length(r$preadjustment$tables$a1),96L,info='adjustment stays on n')
    expect_identical(stats::tsp(r$preadjustment$tables$out)[1L],1989+11/12,
      info='extended tables start one month before the input calendar')
    expect_identical(stats::tsp(r$preadjustment$tables$a1)[1L],1990,
      info='n-length tables keep the input calendar')
    expect_identical(stats::tsp(r$results$d10)[1L],1989+11/12,info='d10 shares the shift')
    expect_identical(length(r$preadjustment$row_basis$out),99L,info='row basis follows out')
    expect_identical(r$preadjustment$available_rows$out,seq_len(99L),
      info='available rows complete on the extended calendar')
    expect_true(length(r$results$d10)==99L&&length(r$results$c17)==99L)
    expect_identical(length(r$finals$tables$d11final),96L,info='finals stay on n')
    expect_identical(length(r$finalSA),96L)
    expect_false(anyNA(unlist(r$preadjustment$available_rows)),
      info='no index overrun in available_rows')
    # The core normalizes a negative count to |h| annual frequencies; the
    # decoration follows the normalized backcast.
    nh<-fit_case('quick',forecast=-1)
    expect_identical(nh$status,'available',info='negative horizon accepted (21-H.4)')
    expect_identical(length(nh$preadjustment$tables$out),108L,
      info='forecast=-1 means one annual frequency: 12 months')
    expect_identical(as.integer(nh$request$x11$forecast),-1L,info='echo keeps the raw request')
    expect_identical(stats::tsp(nh$preadjustment$tables$out)[1L],1990,
      info='forecast-only horizons do not shift the start')
    expect_identical(fit_case('quick',calendar=TRUE)$status,
      'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)')
    expect_identical(fit_case('quick',used=TRUE,span=c(2,95))$status,
      'unsupported:x13 automatic fit spans the full series only')
    expect_identical(fit_case('quick',frequency=4L)$status,
      'unsupported:x13 automatic fit monthly series only',info='wrapper passes the refusal through')

    # The maintainer default: omitted sigma == explicit legacy, bit for bit,
    # labelled in the provenance; the raw entry still refuses a missing policy.
    d0<-fit_case('quick',sigma=NULL)
    d1<-fit_case('quick',sigma='legacy')
    expect_identical(d0$sigma_provenance,'legacy (default)')
    expect_identical(d1$sigma_provenance,'legacy (explicit)')
    expect_identical(d0$request$x11$sigma,'legacy')
    expect_identical(d0$status,'available')
    for(nm in names(d0$model))expect_identical(d0$model[[nm]],d1$model[[nm]],
      info=paste('default sigma model:',nm))
    for(nm in c('finalSA','finalTrend','finalIrregular'))
      expect_true(bits(d0[[nm]],d1[[nm]]),info=paste('default sigma equals explicit:',nm))
    for(nm in c('d10','d11','d12','d13','c17'))
      expect_true(bits(d0$results[[nm]],d1$results[[nm]]),
        info=paste('default sigma results:',nm))
    expect_true(bits(d0$preadjustment$tables$out,d1$preadjustment$tables$out),
      info='default sigma preadjustment')
    # The wrapper echo and the selection record prove the request carried
    # legacy (the accurate oracle replay is the only facts-pinned policy;
    # the raw refusal of a NULL policy keeps the core's explicitness).
    expect_identical(d0$selection$policy,'legacy')
    expect_identical(d1$selection$policy,'legacy')
    expect_identical(quick_res$selection$policy,'accurate')

    # ts inference: the same run from a proper monthly ts, bit-identical to
    # the explicit-scalar quick run.
    m<-case_modules('quick');req<-m$req
    yts<-stats::ts(m$y,start=c(1990,1),frequency=12)
    t1<-fit(series=yts,transform=list(aicc=as.numeric(req$aicdiff)),
      automodel=list(max_d=2L,max_bd=1L,ub1=1/as.numeric(req$ub1),
        ub2=1/as.numeric(req$ub2),cancel=as.numeric(req$cancel),precision=1e-7,
        models_count=5L,max_p=2L,max_q=2L,max_bp=1L,max_bq=1L,
        balanced=identical(req$balanced,'true'),mixed=identical(req$mixed,'true'),
        accept_white_noise=FALSE,eps=1e-7),
      outliers=list(used=FALSE,types=NULL,span=NULL,max_round=as.integer(req$maxiter),
        max_outliers=30L,tcrate=as.numeric(req$tcrate)),
      ami=list(precision=1e-7,mean_precision=1e-5,va=0,reduce_va=as.numeric(req$predcv),
        ljung_box_limit=as.numeric(req$lblimit),check_mu=TRUE,accept_airline=FALSE,
        ur_limit=0.95),
      regression=list(mean_test_initial=1.96,mean_test_final=1,tsig=1),
      sigma='accurate')
    expect_identical(t1$status,'available')
    for(nm in names(t1$model))expect_identical(t1$model[[nm]],quick_res$model[[nm]],
      info=paste('ts inference model:',nm))
    expect_true(bits(t1$finalSA,quick_res$finalSA),info='ts inference tables bit-identical')
    expect_identical(stats::tsp(t1$series),stats::tsp(yts))
  })
})
