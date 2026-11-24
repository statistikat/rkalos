# Item 21-H.4: the automatic-pipeline X-11 horizons qualified at the R
# boundary against the pinned 4-case JVM horizon capture (capture-
# x13-automatic-horizon.R, vendored under fixtures/native/x13-horizon and
# pinned in core.lock.json; docs/contracts/arima-forecasts-contract.md and
# the "Horizon extension" tier policy of x13-automatic-pipeline-contract.md).
# The tier policy is transcribed from the core acceptance test
# tests/test_x13_automatic_horizon.cpp at f500162: decisions BIT-exact
# (status, final description, variable identities, neff, horizon echo), a1 and
# d18 bit-exact claims, every other boundary-exposed table carries the
# measured band of the kCases matrix (tests/jvm_x13_automatic_horizon.h) —
# copied, NEVER widened (the widest per case, so every case's own band gates
# its rows). The split horizon tables cross the R boundary; d8 remains internal. Landing coefficients and statistics
# reuse the existing zero-horizon contract bands (2e-3 / 1e-2 / 2e-3): the
# landing model is horizon-invariant upstream. The extended composition
# calendar (out, d10..c17 start bc months early, length n+bc+fc) and the
# row-basis decoration are pinned here; zero-horizon behaviour stays with the
# 13-case replay (test-native-x13-automatic.R).
source("helpers/fixtures.R", local = TRUE)
local({
  fit<-rkalos::native_x13_fit_automatic
  bits<-function(a,b) {a<-as.vector(a);b<-as.vector(b)
    identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(b)])}
  # Oracle tier semantics (the main file's single(), one for one): bit-equal or
  # both-NaN rows pass, the relative deviation against the want (absolute when
  # zero) stays inside the GIVEN band, the got is finite. The per-family worst
  # is recorded so the run measures the drift it allows.
  worst<-new.env();CASE_ID<-''
  dev<-function(got,want,famid,band) {
    fam<-paste0(CASE_ID,':',famid)
    if(!exists(fam,worst))assign(fam,0,envir=worst)
    a<-as.vector(as.numeric(got));b<-as.vector(as.numeric(want))
    if(length(a)!=length(b)||!all(is.na(a)==is.na(b)))return(FALSE)
    k<-!is.na(a)
    d<-ifelse(abs(b[k])>0,abs(a[k]-b[k])/abs(b[k]),abs(a[k]-b[k]))
    if(length(d)&&max(d)>get(fam,worst))assign(fam,max(d),envir=worst)
    all(is.finite(a[k]))&&(length(d)==0||all(d<=band))}

  # ---- evidence integrity: the vendored capture is byte-pinned ------------
  ev<-'fixtures/native/x13-horizon'
  sums<-read.csv(file.path(ev,'SHA256SUMS.csv'),colClasses='character')
  have<-sort(list.files(ev)); want<-sort(unique(c(sums$file,'SHA256SUMS.csv')))
  expect_identical(have,want,info='horizon capture file set equals the pinned manifest')
  dig<-tools::sha256sum(file.path(ev,sums$file))
  names(dig)<-sums$file
  for(f in sums$file)
    expect_identical(unname(dig[[f]]),sums$sha256[match(f,sums$file)],
      info=paste('sha256',f))

  # ---- corpus completeness -------------------------------------------------
  cases<-read.csv(file.path(ev,'cases.csv'))
  inputs<-read.csv(file.path(ev,'inputs.csv'),colClasses='character')
  facts<-read.csv(file.path(ev,'oracle-facts.csv'),colClasses='character')
  expect_identical(sort(as.character(cases$id)),
    sort(c('quick_h6','ao_b3','multi_b3h6','so_h6')))
  expect_true(all(cases$n==96L))
  expect_identical(c(cases$bc[match('multi_b3h6',cases$id)],
    cases$fc[match('multi_b3h6',cases$id)]),c(3L,6L))

  fact<-function(id,key){v<-facts$value[facts$id==id&facts$key==key]
    if(length(v)==0)"" else v[1]}
  facts_of<-function(id,key)facts$value[facts$id==id&facts$key==key]
  reals<-function(field){t<-strsplit(field,';',fixed=TRUE)[[1]]
    fixture_doubles(t[nchar(t)>0])}
  modes<-c(Additive='additive',Multiplicative='multiplicative',
    LogAdditive='log_additive',PseudoAdditive='pseudo_additive',Undefined='additive')
  TYPE<-c(AO=0L,LS=1L,TC=2L,SO=3L)

  # The kCases band matrix for the families crossing this surface, transcribed
  # verbatim from tests/jvm_x13_automatic_horizon.h at f500162 (TableIndex
  # order; the zero-band claims, a1 and d18, are the BIT_TABLES).
  BIT_TABLES<-c('a1','d18')
  HZ<-list(
    quick_h6=c(out=2.0e-06,d10=1.0e-05,d11=2.0e-06,d12=2.0e-06,d13=5.0e-03,
      c17=5.0e-04,d11final=5.0e-07,d12final=1.0e-06,d13final=2.0e-03,d16=1.0e-05,
      e1=5.0e-07,e2=5.0e-07,e3=2.0e-03,e11=5.0e-07),
    ao_b3=c(out=5.0e-06,d10=5.0e-07,d11=5.0e-06,d12=2.0e-06,d13=1.0e-06,
      c17=1.0e-03,d11final=5.0e-07,d12final=1.0e-06,d13final=1.0e-06,d16=5.0e-07,
      e1=2.0e-06,e2=1.0e-01,e3=1.0e-06,e11=1.0e-06),
    multi_b3h6=c(out=1.0e-05,d10=2.0e-06,d11=1.0e-05,d12=5.0e-06,d13=2.0e-06,
      c17=5.0e-03,d11final=2.0e-06,d12final=2.0e-06,d13final=2.0e-06,d16=2.0e-06,
      e1=5.0e-06,e2=2.0e-01,e3=2.0e-06,e11=5.0e-06),
    so_h6=c(out=2.0e-05,d10=5.0e-06,d11=1.0e-05,d12=5.0e-06,d13=5.0e-06,
      c17=1.0e-02,d11final=5.0e-06,d12final=5.0e-06,d13final=5.0e-06,d16=5.0e-06,
      e1=1.0e-06,e2=5.0e-06,e3=5.0e-06,e11=5.0e-06))
  SLICE_BANDS <- list(
    quick_h6 = c(a1a=2.0e-06, a1b=0.0e+00, d11a=2.0e-06, d12a=2.0e-06, d16a=1.0e-05, d18a=0.0e+00, d11b=0.0e+00, d12b=0.0e+00, d16b=0.0e+00, d18b=0.0e+00),
    ao_b3 = c(a1a=0.0e+00, a1b=5.0e-06, d11a=0.0e+00, d12a=0.0e+00, d16a=0.0e+00, d18a=0.0e+00, d11b=5.0e-06, d12b=2.0e-06, d16b=5.0e-07, d18b=0.0e+00),
    multi_b3h6 = c(a1a=5.0e-06, a1b=2.0e-06, d11a=2.0e-06, d12a=2.0e-06, d16a=5.0e-07, d18a=0.0e+00, d11b=5.0e-07, d12b=2.0e-07, d16b=2.0e-06, d18b=0.0e+00),
    so_h6 = c(a1a=2.0e-05, a1b=0.0e+00, d11a=1.0e-05, d12a=5.0e-06, d16a=5.0e-06, d18a=0.0e+00, d11b=0.0e+00, d12b=0.0e+00, d16b=0.0e+00, d18b=0.0e+00))
  # Landing bands inherited from the zero-horizon contract (never widened).
  LANDING<-c('landing-coefficients'=2e-3,'arma-landing'=1e-2,'landing-stats'=2e-3)

  # The case request mirror: the req-fact scalars verbatim plus the fixed
  # decoder constants (the projection of the 13-case replay), the horizons of
  # cases.csv arriving as request data.
  case_modules<-function(id) {
    line<-strsplit(fact(id,'req'),' ',fixed=TRUE)[[1]]
    req<-as.list(line[seq(2,length(line),2)])
    names(req)<-line[seq(1,length(line)-1,2)]
    iline<-inputs$line[inputs$id==id]
    inp<-stats::setNames(as.list(sub('^[^ ]* ','',iline)),sub(' .*','',iline))
    y<-reals(inp$values)
    expect_identical(length(y),as.integer(req$n),info=paste(id,'input length'))
    sy<-as.integer(strsplit(req$start,':',fixed=TRUE)[[1]])
    otypes<-strsplit(req$otypes,',',fixed=TRUE)[[1]]
    list(y=y,first_year=sy[1],start_period=sy[2],req=req,
      used=identical(req$outliers,'on'),
      types=if(identical(req$outliers,'on'))otypes else NULL)
  }
  fit_case<-function(id, as_ts=FALSE) {
    m<-case_modules(id);req<-m$req
    ca<-cases[match(id,cases$id),]
    sigma_cut<-as.numeric(strsplit(req$sigma,';',fixed=TRUE)[[1]])
    series <- if (as_ts) stats::ts(m$y, start=c(m$first_year, m$start_period), frequency=12) else m$y
    fit(series=series,transform=list(aicc=as.numeric(req$aicdiff)),
      automodel=list(max_d=2L,max_bd=1L,ub1=1/as.numeric(req$ub1),
        ub2=1/as.numeric(req$ub2),cancel=as.numeric(req$cancel),precision=1e-7,
        models_count=5L,max_p=2L,max_q=2L,max_bp=1L,max_bq=1L,
        balanced=identical(req$balanced,'true'),mixed=identical(req$mixed,'true'),
        accept_white_noise=FALSE,eps=1e-7),
      outliers=list(used=m$used,types=m$types,span=NULL,
        max_round=as.integer(req$maxiter),max_outliers=30L,
        tcrate=as.numeric(req$tcrate)),
      ami=list(precision=1e-7,mean_precision=1e-5,va=0,reduce_va=as.numeric(req$predcv),
        ljung_box_limit=as.numeric(req$lblimit),check_mu=TRUE,accept_airline=FALSE,
        ur_limit=0.95),
      regression=list(mean_test_initial=1.96,mean_test_final=1,tsig=1),
      sigma='accurate',lower_sigma=sigma_cut[1],upper_sigma=sigma_cut[2],
      msr_policy='jdplus',mode=modes[[as.character(req$mode)]],
      frequency=12L,start_period=m$start_period,first_year=m$first_year,
      backcast=as.integer(ca$bc),forecast=as.integer(ca$fc))
  }

  # ---- the envelope the boundary itself states ------------------------------
  local({
    cap<-getFromNamespace('capabilities','rkalos')()
    af<-cap$automatic_fit
    expect_identical(af$schema_version,3L)
    expect_null(af$x11$zero_horizon_only)
    expect_true(grepl('21-H.4',af$x11$horizon_envelope))
    # The restricted surface keeps its {0,0} domain: the horizon route is the
    # automatic initialStep route only.
    expect_true(cap$restricted_fit$zero_horizon_only)
  })

  # ---- per-case JVM horizon replay ------------------------------------------
  for(i in seq_len(nrow(cases))) {
    ca<-cases[i,];id<-as.character(ca$id);CASE_ID<-id
    cat('== horizon case',id,'\n')
    m<-case_modules(id);n<-length(m$y)
    ext<-n+as.integer(ca$bc)+as.integer(ca$fc)
    want_start<-m$first_year+(m$start_period-1)/12
    o<-m$first_year*12L+m$start_period-1L-as.integer(ca$bc)
    ext_start<-floor(o/12)+(o%%12)/12
    res<-fit_case(id)
    expect_identical(res$status,'available',info=paste(id,'status'))
    expect_identical(res$sigma_provenance,'accurate (explicit)',info=paste(id,'provenance'))
    expect_equal(res$selection$start_period, o %% 12 + 1)
    expect_equal(res$selection$backcast, ca$bc)
    expect_equal(res$selection$forecast, ca$fc)
    for (pass_name in c("b", "c", "d")) {
      pass <- res$selection[[pass_name]]
      expect_equal(pass$input_domain$length, ext)
      expect_equal(pass$input_domain$start_ordinal, o)
      first <- paste0(toupper(pass_name), "1")
      expect_equal(stats::tsp(pass$tables[[first]])[1L], ext_start)
      expect_identical(pass$row_basis[[first]], seq_len(ext))
    }
    expect_equal(res$selection$d$original_domain$length, n)
    if (id == "multi_b3h6")
      expect_identical(fit_case(id, as_ts=TRUE)$selection, res$selection)
    # decisions: final description and horizon echo, BIT.
    expect_identical(res$model$log,identical(as.character(ca$log),'true'),info=paste(id,'log'))
    expect_identical(res$model$mean,identical(as.character(ca$mean),'true'),info=paste(id,'mean'))
    for(nm in c('p','d','q','bp','bd','bq'))
      expect_identical(res$model[[nm]],as.integer(ca[[nm]]),info=paste(id,nm))
    expect_identical(res$mode_effective,modes[[as.character(fact(id,'mode_effective'))]],
      info=paste(id,'mode effective'))
    expect_identical(as.integer(res$request$x11$backcast),as.integer(ca$bc),info=paste(id,'bc echo'))
    expect_identical(as.integer(res$request$x11$forecast),as.integer(ca$fc),info=paste(id,'fc echo'))
    # variables: identity bit, landing coefficients banded (2e-3 basin).
    parse_var<-function(row){t<-strsplit(row,' ',fixed=TRUE)[[1]]
      list(row=as.integer(t[1]),name=paste(t[-c(1,length(t)-1,length(t))],collapse=' '),
        ami=t[length(t)-1],regeffect=t[length(t)])}
    vs<-lapply(facts_of(id,'var'),parse_var)
    ami<-Filter(function(v)identical(v$ami,'x13'),vs)
    mv<-res$model$variables
    expect_identical(length(mv$name),length(ami),info=paste(id,'nvars'))
    if(length(ami)) {
      ym<-vapply(ami,function(v)regmatches(v$name,
        regexpr('[0-9]{4}-[0-9]{2}',v$name)),"")
      want_pos<-as.integer(substr(ym,1,4))*12L+as.integer(substr(ym,6,7))-
        (m$first_year*12L+m$start_period)+1L
      want_type<-vapply(ami,function(v)as.integer(TYPE[substr(v$name,1,2)]),integer(1))
      want_name<-vapply(ami,function(v)v$name,"")
      want_eff<-vapply(ami,function(v)if(identical(v$regeffect,'-'))'' else v$regeffect,"")
      ow<-order(want_pos,want_type);og<-order(mv$position,mv$type_index)
      expect_true(bits(mv$position[og],want_pos[ow]),info=paste(id,'var positions'))
      expect_identical(mv$type_index[og],want_type[ow],info=paste(id,'var types'))
      expect_identical(mv$name[og],want_name[ow],info=paste(id,'var names'))
      expect_identical(mv$regeffect[og],want_eff[ow],info=paste(id,'var regeffects'))
      regcoef<-reals(fact(id,'regcoef'))
      for(k in seq_along(ami))
        expect_true(dev(mv$coefficient[og][k],regcoef[ami[[ow[k]]]$row+1L],
          'landing-coefficients',LANDING[['landing-coefficients']]),
          info=paste(id,'var coefficient',k))
    }
    # The landed ARMA block (existing 1e-2 basin band) and landing statistics
    # (neff bit; loglik/aicc in the existing 2e-3 band): horizon-invariant.
    flat<-as.vector(unlist(res$model$coefficients))
    expect_identical(length(flat),length(reals(fact(id,'arma'))),info=paste(id,'narma'))
    expect_true(dev(flat,reals(fact(id,'arma')),'arma-landing',
      LANDING[['arma-landing']]),info=paste(id,'arma landing'))
    st<-regmatches(fact(id,'stats'),regexec(
      'nobs=([0-9]+) neff=([0-9]+) npar=([0-9]+) loglik=([-+0-9.eE]+) aicc=([-+0-9.eE]+)',
      fact(id,'stats')))[[1]]
    expect_identical(as.integer(res$model$effective_observations),as.integer(st[3]),
      info=paste(id,'neff'))
    expect_true(dev(res$model$log_likelihood,as.numeric(st[5]),'landing-stats',
      LANDING[['landing-stats']]),info=paste(id,'loglik'))
    expect_true(dev(res$model$aicc,as.numeric(st[6]),'landing-stats',
      LANDING[['landing-stats']]),info=paste(id,'aicc'))
    # The extended composition: lengths (bit) against the kCases lengths row,
    # then the table values (bit for a1/d18, banded for the rest).
    pa<-res$preadjustment$tables;fi<-res$finals$tables
    tabs<-c(as.list(pa[c('a1','out')]),
      res$results[c('d10','d11','d12','d13','c17')],
      as.list(fi[c('d11final','d12final','d13final','d16','d18','e1','e2','e3','e11')]))
    for(nm in c('a1','d11final','d12final','d13final','d16','d18','e1','e2','e3','e11'))
      expect_identical(length(tabs[[nm]]),n,info=paste(id,nm,'length'))
    expect_identical(length(pa$out),ext,info=paste(id,'out length'))
    for(nm in c('d10','d11','d12','d13','c17'))
      expect_identical(length(tabs[[nm]]),ext,info=paste(id,nm,'length'))
    # the capture carries the same extended lengths (kCases sanity)
    expect_identical(length(reals(fact(id,'out'))),ext,info=paste(id,'capture out length'))
    for(nm in c('d10','d11','d12','d13','c17'))
      expect_identical(length(reals(fact(id,nm))),ext,info=paste(id,'capture',nm,'length'))
    slice_bands <- SLICE_BANDS[[id]]
    expect_identical(rkalos::capabilities()$automatic_fit$x11$horizon_tables, names(slice_bands))
    for (nm in names(slice_bands)) {
      back <- endsWith(nm, "b")
      count <- if (back) res$effective_backcast else res$effective_forecast
      slice <- if (startsWith(nm, "a1")) pa[[nm]] else fi[[nm]]
      if (!count) {
        expect_null(slice, info=paste(id,nm,"absent window"))
        next
      }
      expected_start <- stats::tsp(res$series)[1] + (if (back) -count else n) / 12
      expect_true(stats::is.ts(slice), info=paste(id,nm,"ts"))
      expect_identical(length(slice), count, info=paste(id,nm,"length"))
      expect_equal(stats::tsp(slice), c(expected_start, expected_start + (count-1) / 12, 12),
                   info=paste(id,nm,"calendar"))
      expected <- reals(fact(id,nm))
      if (slice_bands[[nm]] == 0) {
        expect_true(bits(slice,expected),info=paste(id,nm,"bit tier"))
      } else {
        expect_true(dev(slice,expected,nm,slice_bands[[nm]]),info=paste(id,nm,"numerical tier"))
      }
    }
    for(nm in BIT_TABLES)
      expect_true(bits(tabs[[nm]],reals(fact(id,nm))),info=paste(id,nm,'bit tier'))
    for(nm in names(HZ[[id]]))
      expect_true(dev(tabs[[nm]],reals(fact(id,nm)),nm,HZ[[id]][[nm]]),info=paste(id,nm))
    # The boundary decoration over the extended calendar: the extended tables
    # start bc months before the input calendar, the n tables keep it; the row
    # basis follows each table (no index overrun into available_rows).
    expect_identical(stats::tsp(pa$a1)[1L],want_start,info=paste(id,'a1 calendar'))
    expect_identical(stats::tsp(pa$out)[1L],ext_start,
      info=paste(id,'out start shifted by the backcast'))
    for(nm in c('d10','d11','d12','d13','c17'))
      expect_identical(stats::tsp(res$results[[nm]])[1L],ext_start,info=paste(id,nm,'start'))
    for(nm in setdiff(names(fi), names(slice_bands)))
      expect_identical(stats::tsp(fi[[nm]])[1L],want_start,info=paste(id,nm,'finals calendar'))
    rb<-res$preadjustment$row_basis
    for(nm in names(pa))
      expect_identical(rb[[nm]],seq_along(pa[[nm]]),info=paste(id,'row basis',nm))
    expect_identical(names(res$preadjustment$available_rows),names(pa))
    expect_false(anyNA(unlist(res$preadjustment$available_rows)),
      info=paste(id,'available rows complete'))
  }

  # The measured drift stays inside the transcribed bands (recorded, never
  # widened). The per-case rows above already gated every case against its own
  # band; this scoreboard re-checks each per-case worst against that same case
  # band (no cross-case min: a case is never judged by another case's band)
  # and prints the drift the run allows.
  for(key in ls(worst)) {
    spl<-strsplit(key,':',fixed=TRUE)[[1]]
    id<-spl[1];fam<-spl[2]
    band<-if(fam %in% names(LANDING))LANDING[[fam]] else if (fam %in% names(SLICE_BANDS[[id]])) SLICE_BANDS[[id]][[fam]] else HZ[[id]][[fam]]
    expect_true(get(key,worst)<=band,info=paste('worst',key,get(key,worst)))
    cat(sprintf('hz %-27s worst %.3g (band %.0g)\n',key,get(key,worst),band))
  }
})
