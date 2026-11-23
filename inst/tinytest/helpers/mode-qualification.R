# Shared R-boundary qualification for the item 22b.5 decomposition modes.
# One call replays the whole pinned fixture family for a mode through the
# native pass boundaries: the full B- and C-pass sweeps, the full D-pass
# sweep (extended calendars included) and the full selection sweep, plus the
# mode echo/inheritance contract, determinism and the wrapper paths. The
# additive path stays in the per-pass test files; this helper qualifies the
# three non-additive modes on the shared mult-shape input family.
x11_mode_qualification <- function(tag, mode, b20_correct, pinned_label,
                                   unsupported, repair, identities) {
  rawb<-getFromNamespace('.native_x11_b_pass','rkalos')
  rawc<-getFromNamespace('.native_x11_c_pass','rkalos')
  rawd<-getFromNamespace('.native_x11_d_pass','rkalos')
  raws<-getFromNamespace('.native_x11_selection_pass','rkalos')
  pubb<-getFromNamespace('native_x11_b_pass','rkalos')
  pubc<-getFromNamespace('native_x11_c_pass','rkalos')
  pubd<-getFromNamespace('native_x11_d_pass','rkalos')
  pubs<-getFromNamespace('native_x11_selection_pass','rkalos')
  bits<-function(a,b) identical(is.na(a),is.na(b)) && identical(a[!is.na(a)],b[!is.na(a)])
  num<-function(x)as.integer(strsplit(x,' ')[[1]])
  csv <- function(dir, file, columns = NULL) {
    path <- paste0('fixtures/native/x11-', dir, '-', tag, '/', file)
    if (is.null(columns)) return(utils::read.csv(path))
    read_numeric_fixture(path, columns)
  }
  read_input<-function(obj,id) {
    rows<-which(obj$case==id);o<-order(as.integer(obj$row[rows]))
    v<-obj$value[rows][o];if('missing' %in% names(obj))v[obj$missing[rows][o]]<-NA_real_
    v
  }
  blabels<-c('B1','B2','B3','B4','B4a','B4anorm','B4d','B4g','B5','B6','B7','B8','B9','B9g',
    'B10','B11','B13','B17','B20')
  bcenter<-c('B2','B3','B4','B4a','B4anorm','B4d','B4g')
  clabels<-c('C1','C2','C4','C5','C5a','C6','C7','C9','C10','C10a','C11','C13','C17','C20')
  ccenter<-c('C2','C4','C5a')
  dlabels<-c('D1','D2','D4','D5','D6','D7','D8','D9','D9bis','D9_g_bis','D10','D10bis',
    'D11','D11bis','D12','D13')
  dcenter<-c('D2','D4')
  # ---------------- the fixed B-pass sweep (24 cases) ----------------
  bman<-csv('bpass','cases.csv'); binp<-csv('bpass','inputs.csv','value'); btab<-csv('bpass','tables.csv','value')
  expect_identical(nrow(bman),24L)
  expect_identical(unique(bman$mode),pinned_label)
  for(i in seq_len(nrow(bman))) {
    ca<-bman[i,];id<-ca$case
    v<-read_input(binp,id);before<-serialize(v,NULL)
    f<-as.integer(ca$frequency);t<-as.integer(f/2)
    initial<-num(ca$initial_filters);final<-num(ca$final_filters)
    b<-rawb(v,f,ca$start_period,'legacy',initial,final,ca$trend,ca$lsig,ca$usig,mode)
    expect_identical(b$status,'available',info=paste(tag,id,'status'))
    expect_identical(b$mode,mode,info=paste(tag,id,'mode echo'))
    expect_identical(b$scope,paste0('fixed_',mode,'_b_pass_only'),info=paste(tag,id,'scope'))
    expect_identical(b$b2drop,t,info=paste(tag,id,'b2drop'))
    expect_identical(as.integer(b$initial_filters),initial,info=paste(tag,id,'initial echo'))
    expect_identical(as.integer(b$final_filters),final,info=paste(tag,id,'final echo'))
    expect_identical(as.integer(b$final_filters),num(ca$resolved_final),
      info=paste(tag,id,'resolved final'))
    for(label in blabels) {
      rows<-which(btab$case==id & btab$table==label)
      pv<-btab$value[rows][order(as.integer(btab$row[rows]))]
      expect_true(bits(b$tables[[label]],pv),info=paste(tag,id,label))
      expect_identical(length(b$tables[[label]]),
        if(label %in% bcenter) length(v)-2L*t else length(v),info=paste(tag,id,label,'span'))
    }
    expect_identical(b$rounds[[1]]$replacement_count,as.integer(ca$replacements_r1),
      info=paste(tag,id,'round 1 replacements'))
    expect_identical(b$rounds[[2]]$replacement_count,as.integer(ca$replacements_r2),
      info=paste(tag,id,'round 2 replacements'))
    if(repair) # the pinned makePositivity keeps every repaired trend row positive
      expect_true(all(b$tables$B7[!is.na(b$tables$B7)]>0),info=paste(tag,id,'B7 positivity'))
    again<-rawb(v,f,ca$start_period,'legacy',initial,final,ca$trend,ca$lsig,ca$usig,mode)
    expect_identical(again,b,info=paste(tag,id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(tag,id,'input untouched'))
  }
  # ---------------- the fixed C-pass sweep (24 cases) ----------------
  cman<-csv('cpass','cases.csv'); cinp<-csv('cpass','inputs.csv','value'); ctab<-csv('cpass','tables.csv','value')
  expect_identical(nrow(cman),24L)
  expect_identical(unique(cman$mode),pinned_label)
  expect_true(identical(sort(unique(cman$case)),sort(unique(bman$case))),
    info=paste(tag,'cpass replays the bpass cases'))
  for(i in seq_len(nrow(cman))) {
    ca<-cman[i,];id<-ca$case
    v<-read_input(cinp,id);f<-as.integer(ca$frequency);t<-as.integer(f/2)
    initial<-num(ca$initial_filters);final<-num(ca$final_filters)
    if(identical(identities,'log'))L<-log(v)
    b<-rawb(v,f,ca$start_period,'legacy',initial,final,ca$trend,ca$lsig,ca$usig,mode)
    c<-rawc(v,f,ca$start_period,'legacy',initial,final,ca$trend,b,ca$lsig,ca$usig)
    expect_identical(c$status,'available',info=paste(tag,id,'status'))
    expect_identical(c$mode,mode,info=paste(tag,id,'inherited mode echo'))
    expect_identical(c$mode_source,'inherited from the B-pass result',
      info=paste(tag,id,'mode source'))
    expect_identical(c$scope,paste0('fixed_',mode,'_c_pass_only'),info=paste(tag,id,'scope'))
    expect_identical(c$c2drop,t,info=paste(tag,id,'c2drop'))
    for(label in clabels) {
      rows<-which(ctab$case==id & ctab$table==label)
      pv<-ctab$value[rows][order(as.integer(ctab$row[rows]))]
      expect_true(bits(c$tables[[label]],pv),info=paste(tag,id,label))
      expect_identical(length(c$tables[[label]]),
        if(label %in% ccenter) length(v)-2L*t else length(v),info=paste(tag,id,label,'span'))
    }
    expect_identical(c$round$replacement_count,as.integer(ca$corrected_rc),
      info=paste(tag,id,'C round corrections'))
    if(repair)
      expect_true(all(c$tables$C7[!is.na(c$tables$C7)]>0),info=paste(tag,id,'C7 positivity'))
    if(identical(identities,'log')) {
      # the pinned kernel runs the plain additive steps on log(y): the
      # log-domain additive relations survive bit-exactly on those tables.
      expect_true(bits(c$tables$C1,L-b$tables$B20),info=paste(tag,id,'C1 log identity'))
      expect_true(bits(c$tables$C4,c$tables$C1[(t+1L):(length(v)-t)]-c$tables$C2),
        info=paste(tag,id,'C4 identity'))
      expect_true(bits(c$tables$C11,L-c$tables$C10),info=paste(tag,id,'C11 log identity'))
      expect_true(bits(c$tables$C13,c$tables$C11-c$tables$C7),info=paste(tag,id,'C13 identity'))
    }
    again<-rawc(v,f,ca$start_period,'legacy',initial,final,ca$trend,b,ca$lsig,ca$usig)
    expect_identical(again,c,info=paste(tag,id,'determinism'))
  }
  # ---------------- the fixed D-pass sweep (28 cases) ----------------
  dman<-csv('dpass','cases.csv','global_msr'); dinp<-csv('dpass','inputs.csv','value'); dtab<-csv('dpass','tables.csv','value')
  dmsr<-csv('dpass','msr.csv',c('ri', 'rs', 'msr'))
  expect_identical(nrow(dman),28L)
  expect_identical(sum(dman$category=='extended'),4L)
  expect_identical(unique(dman$mode),pinned_label)
  for(i in seq_len(nrow(dman))) {
    ca<-dman[i,];id<-ca$case
    v<-read_input(dinp,id);before<-serialize(v,NULL)
    f<-as.integer(ca$frequency);t<-as.integer(f/2)
    initial<-num(ca$initial_filters);final<-num(ca$final_filters)
    expect_identical(as.integer(ca$n),length(v),info=paste(tag,id,'input length'))
    if(ca$category=='extended') {
      base<-sub('^ext-','',id)
      expect_true(bits(v[(as.integer(ca$back)+1L):(as.integer(ca$back)+as.integer(ca$orig_n))],
        read_input(dinp,base)),info=paste(tag,id,'extends the base input'))
      expect_identical(as.integer(ca$orig_start_period),
        as.integer((ca$start_period-1L+ca$back)%%f)+1L,info=paste(tag,id,'original first period'))
    } else expect_identical(as.integer(ca$orig_n),length(v))
    b<-rawb(v,f,ca$start_period,'legacy',initial,final,ca$trend,ca$lsig,ca$usig,mode)
    c<-rawc(v,f,ca$start_period,'legacy',initial,final,ca$trend,b,ca$lsig,ca$usig)
    d<-rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,c,ca$back,ca$fc)
    expect_identical(d$status,'available',info=paste(tag,id,'status'))
    expect_identical(d$mode,mode,info=paste(tag,id,'inherited mode echo'))
    expect_identical(d$mode_source,'inherited from the C-pass result',
      info=paste(tag,id,'mode source'))
    expect_identical(d$scope,paste0('fixed_',mode,'_d_pass_only'),info=paste(tag,id,'scope'))
    expect_identical(d$d2drop,t,info=paste(tag,id,'d2drop'))
    expect_identical(d$backcast,as.integer(ca$back),info=paste(tag,id,'backcast echo'))
    expect_identical(d$forecast,as.integer(ca$fc),info=paste(tag,id,'forecast echo'))
    expect_identical(d$msr_iterations,0L,info=paste(tag,id,'fixed mode never selects'))
    expect_identical(as.integer(d$seas_filter),final,info=paste(tag,id,'fixed seas filter'))
    for(label in dlabels) {
      rows<-which(dtab$case==id & dtab$table==label)
      pv<-dtab$value[rows][order(as.integer(dtab$row[rows]))]
      expect_true(bits(d$tables[[label]],pv),info=paste(tag,id,label))
      expect_identical(length(d$tables[[label]]),
        if(label %in% dcenter) length(v)-2L*t else length(v),info=paste(tag,id,label,'span'))
    }
    expect_true(bits(d$tables$D9_g_bis,d$tables$D9bis),info=paste(tag,id,'D9_g_bis alias'))
    mrows<-dmsr[dmsr$case==id,]
    expect_identical(nrow(mrows),f,info=paste(tag,id,'MSRE phases'))
    mp<-as.integer(mrows$phase)+1L
    expect_true(bits(d$msr_table$mean_irregular_evolutions[mp],mrows$ri),
      info=paste(tag,id,'MSRE irregular'))
    expect_true(bits(d$msr_table$mean_seasonal_evolutions[mp],mrows$rs),
      info=paste(tag,id,'MSRE seasonal'))
    expect_identical(d$msr_table$counts[mp],as.integer(mrows$count),info=paste(tag,id,'MSRE counts'))
    expect_true(bits(d$msr_table$msr[mp],mrows$msr),info=paste(tag,id,'MSRE ratios'))
    expect_true(bits(d$msr_table$global_msr,ca$global_msr),info=paste(tag,id,'MSRE global'))
    expect_identical(sum(!is.na(d$tables$D9)),as.integer(ca$d9_nonnan),
      info=paste(tag,id,'D9 mask count'))
    expect_identical(as.integer(ca$d9_nonnan),as.integer(ca$corrected_rc),
      info=paste(tag,id,'D9 count equals the C-round corrections'))
    if(repair) { # makePositivity at d7Step and dFinalStep
      expect_true(all(d$tables$D7[!is.na(d$tables$D7)]>0),info=paste(tag,id,'D7 positivity'))
      expect_true(all(d$tables$D12[!is.na(d$tables$D12)]>0),info=paste(tag,id,'D12 positivity'))
    }
    if(identical(identities,'log')) {
      L<-log(v)
      expect_true(bits(d$tables$D1,L-c$tables$C20),info=paste(tag,id,'D1 log identity'))
      expect_true(bits(d$tables$D4,d$tables$D1[(t+1L):(length(v)-t)]-d$tables$D2),
        info=paste(tag,id,'D4 identity'))
      expect_true(bits(d$tables$D6,d$tables$D1-d$tables$D5),info=paste(tag,id,'D6 identity'))
      expect_true(bits(d$tables$D8,L-d$tables$D7),info=paste(tag,id,'D8 log identity'))
      expect_true(bits(d$tables$D11,L-d$tables$D10),info=paste(tag,id,'D11 log identity'))
      expect_true(bits(d$tables$D11bis,d$tables$D1-d$tables$D10),info=paste(tag,id,'D11bis identity'))
      expect_true(bits(d$tables$D13,d$tables$D11-d$tables$D12),info=paste(tag,id,'D13 identity'))
    }
    again<-rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,c,ca$back,ca$fc)
    expect_identical(again,d,info=paste(tag,id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(tag,id,'input untouched'))
    # the companion coupling stays mode-blind: the core rejects mismatches.
    bad<-c;bad$round$policy<-'accurate'
    expect_error(rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,bad,ca$back,ca$fc),
      'policy must match the C pass',info=paste(tag,id,'policy coupling'))
    bad2<-c;bad2$reference_series[1L]<-bad2$reference_series[1L]+1
    expect_error(rawd(v,f,ca$start_period,'legacy',initial,final,ca$trend,bad2,ca$back,ca$fc),
      'differs from the C pass input',info=paste(tag,id,'reference-series coupling'))
  }
  # ---------------- the selection sweep (15 cases) ----------------
  sman<-csv('selection','cases.csv','icratio'); sinp<-csv('selection','inputs.csv','value')
  stab<-csv('selection','tables.csv','value'); ssels<-csv('selection','selections.csv','ic_ratio')
  smsr<-csv('selection','msr.csv','global_msr')
  expect_identical(nrow(sman),15L)
  expect_identical(unique(sman$mode),pinned_label)
  rq<-function(x) {y<-as.integer(strsplit(x,',')[[1]])
    ifelse(y==0,'stable',ifelse(y==3,'s3x3',ifelse(y==5,'s3x5',
      ifelse(y==9,'s3x9',ifelse(y==-1,'x11default','msr')))))}
  bset<-c('B7','B13','B20'); cset<-c('C7','C13','C20')
  for(i in seq_len(nrow(sman))) {
    ca<-sman[i,];id<-ca$case
    v<-read_input(sinp,id);f<-as.integer(ca$frequency)
    req<-rq(ca$requests);before<-serialize(v,NULL)
    r<-raws(v,f,ca$first_period,req,ca$trend,'legacy','jdplus',1.5,2.5,ca$backcast,
      ca$forecast,mode)
    expect_identical(r$mode,mode,info=paste(tag,id,'mode echo'))
    srows<-ssels[ssels$case==id,]
    expect_identical(nrow(srows),4L,info=paste(tag,id,'stage rows'))
    if(id %in% unsupported) {
      expect_identical(r$status,'unsupported:trend-b7',info=paste(tag,id,'unsupported status'))
      expect_true(is.null(r$b)&&is.null(r$c)&&is.null(r$d),info=paste(tag,id,'records only'))
      expect_true(!is.null(r$initial_plan),info=paste(tag,id,'records kept'))
      next
    }
    expect_identical(r$status,'available',info=paste(tag,id,'status'))
    expect_identical(r$b$mode,mode,info=paste(tag,id,'nested b mode'))
    expect_identical(r$c$mode,mode,info=paste(tag,id,'nested c mode'))
    expect_identical(r$d$mode,mode,info=paste(tag,id,'nested d mode'))
    for(j in seq_len(nrow(srows))) {
      row<-srows[j,];st<-row$stage;rec<-r$stages[[st]]
      if(row$mode=='pinned-chooses-longest-not-claimed'&&st=='b7') {
        expect_identical(rec$status,'unsupported',info=paste(tag,id,st,'refuses NaN'))
        expect_true(is.na(rec$ic_ratio),info=paste(tag,id,st,'NaN refused'))
      } else if(row$mode=='fixed') {
        expect_identical(rec$status,'fixed',info=paste(tag,id,st,'fixed status'))
        expect_identical(rec$length,as.integer(row$pinned_len),info=paste(tag,id,st,'fixed length'))
      } else {
        expect_identical(rec$status,'selected',info=paste(tag,id,st,'selected status'))
        expect_identical(rec$length,as.integer(row$pinned_len),info=paste(tag,id,st,'selected length'))
        if(!is.na(row$ic_ratio))
          expect_true(bits(rec$ic_ratio,row$ic_ratio),info=paste(tag,id,st,'selected ratio'))
      }
    }
    expect_true(bits(r$stages$final$ic_ratio,ca$icratio),info=paste(tag,id,'final iCRatio'))
    expect_identical(r$stages$d7$length,as.integer(ca$d7len),info=paste(tag,id,'D7 length'))
    expect_identical(r$stages$final$length,as.integer(ca$finallen),info=paste(tag,id,'final length'))
    fin<-as.integer(strsplit(gsub('"','',ca$seas),',')[[1]])
    expect_identical(as.integer(r$d$seas_filter),fin,info=paste(tag,id,'seasFilter fill'))
    expect_identical(r$d$msr_iterations,as.integer(ca$d9msriter),info=paste(tag,id,'d9msriter'))
    expect_identical(r$final_plan$requires_msr,as.logical(ca$engaged),
      info=paste(tag,id,'requires_msr echoes the pinned isMSR'))
    expect_identical(r$msr$engaged,as.logical(ca$engaged),info=paste(tag,id,'MSR engaged'))
    # B20 counts the pinned corrections on the mode's own neutral value:
    # 0 on the additive/log scale, 1 on the multiplicative factor scale.
    expect_identical(sum(!is.na(r$b$tables$B20)&r$b$tables$B20!=b20_correct),
      as.integer(ca$b20nz),info=paste(tag,id,'B20 corrections'))
    arows<-smsr[smsr$case==id,]
    if(as.logical(ca$engaged)) {
      arows<-arows[order(as.integer(arows$attempt)),]
      evals<-arows[arows$band!='no_history',]
      expect_identical(length(r$msr$attempts$count),nrow(evals),info=paste(tag,id,'attempt count'))
      expect_identical(as.integer(r$msr$attempts$count),as.integer(evals$count),
        info=paste(tag,id,'attempt windows'))
      expect_true(bits(r$msr$attempts$global_msr,evals$global_msr),info=paste(tag,id,'attempt ratios'))
      expect_identical(as.character(r$msr$attempts$band),as.character(evals$band),
        info=paste(tag,id,'attempt bands'))
      last<-arows[nrow(arows),]
      expect_identical(r$msr$selected_filter,as.integer(last$final_filter),
        info=paste(tag,id,'selected filter'))
      expect_identical(r$msr$using_default,as.logical(last$using_default),
        info=paste(tag,id,'method default flag'))
      expect_identical(r$msr$status,'available',info=paste(tag,id,'MSR status'))
    }
    for(label in unique(stab$table[stab$case==id])) {
      got<-if(label %in% bset) r$b$tables[[label]] else
        if(label %in% cset) r$c$tables[[label]] else r$d$tables[[label]]
      rows<-which(stab$case==id&stab$table==label)
      pv<-stab$value[rows][order(as.integer(stab$row[rows]))]
      expect_true(bits(got,pv),info=paste(tag,id,label))
    }
    again<-raws(v,f,ca$first_period,req,ca$trend,'legacy','jdplus',1.5,2.5,ca$backcast,
      ca$forecast,mode)
    expect_identical(again,r,info=paste(tag,id,'determinism'))
    expect_identical(serialize(v,NULL),before,paste(tag,id,'input untouched'))
  }
  # ---------------- mode surface: validation, rejection, legality ----------
  v<-as.numeric(1:40)
  expect_error(rawb(v,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5,'mult'),
    'decomposition mode must be additive, multiplicative, log_additive or pseudo_additive')
  expect_error(raws(v,4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,0,0,'additive2'),
    'selection-pass decomposition mode')
  if(mode %in% c('multiplicative','log_additive')) {
    # the pinned positivity gate fires on zero and negative rows, verbatim.
    zero<-v;zero[7]<-0
    expect_error(rawb(zero,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5,mode),
      'multiplicative decomposition of non positive series')
    neg<-v;neg[3]<--2
    expect_error(rawb(neg,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5,mode),
      'multiplicative decomposition of non positive series')
    expect_error(raws(zero,4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,0,0,mode),
      'multiplicative decomposition of non positive series')
    # flagged-missing rows are skipped by the fixed-pass gate ...
    miss<-v;miss[5]<-NA_real_
    expect_identical(rawb(miss,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5,mode)$status,
      'available',info=paste(tag,'masked rows skip the positivity gate'))
    # ... while the selection pass keeps its separate completeness rejection.
    expect_error(raws(miss,4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,0,0,mode),
      'complete observations')
  } else {
    # pseudo_additive legally runs on zeros and negatives (no positivity gate).
    z<-v;z[c(3,9,20)]<-c(0,-5,-0.25)
    b<-rawb(z,4,1,'legacy',rep(5L,4),rep(9L,4),13,1.5,2.5,mode)
    c<-rawc(z,4,1,'legacy',rep(5L,4),rep(9L,4),13,b,1.5,2.5)
    d<-rawd(z,4,1,'legacy',rep(5L,4),rep(9L,4),13,c,0,0)
    expect_identical(b$status,'available',info='padd accepts zero')
    expect_identical(c$status,'available',info='padd accepts negatives')
    expect_identical(d$status,'available',info='padd chain accepts zero and negatives')
    expect_error(raws(z,4,1,rep('s3x5',4),13,'legacy','jdplus',1.5,2.5,0,0,'multiplicative'),
      'multiplicative decomposition of non positive series',
      info='the same series is refused by the multiplicative gate')
  }
  # mode-blind core argument validation stays verbatim at the mode surfaces.
  expect_error(rawb(v,4,1,'modern',rep(5L,4),rep(9L,4),13,1.5,2.5,mode),'legacy or accurate')
  expect_error(rawb(v,3,1,'legacy',rep(5L,3),rep(9L,3),13,1.5,2.5,mode),
    'monthly and quarterly')
  # ---------------- the public wrappers on this mode ----------------
  ca<-dman[dman$category!='extended',][1,]
  v<-read_input(dinp,ca$case)
  series<-stats::ts(v,start=c(1997,as.integer(ca$start_period)),frequency=as.integer(ca$frequency))
  initial<-num(ca$initial_filters);final<-num(ca$final_filters)
  bb<-pubb(series,'legacy',initial,final,as.integer(ca$trend),lower_sigma=ca$lsig,
    upper_sigma=ca$usig,mode=mode)
  expect_true(inherits(bb,'persephone_native_x11_b_pass'),info=paste(tag,'b wrapper class'))
  expect_identical(bb$mode,mode,info=paste(tag,'b wrapper mode echo'))
  cc<-pubc(series,'legacy',initial,final,as.integer(ca$trend),bb,lower_sigma=ca$lsig,
    upper_sigma=ca$usig)
  expect_identical(cc$mode,mode,info=paste(tag,'c wrapper inherits the mode'))
  expect_error(pubc(series,'legacy',initial,final,as.integer(ca$trend),bb,lower_sigma=ca$lsig,
    upper_sigma=ca$usig,mode='additive'),
    'C pass inherits its decomposition mode from the B-pass result')
  dd<-pubd(series,'legacy',initial,final,as.integer(ca$trend),cc,mode=mode)
  expect_identical(dd$mode,mode,info=paste(tag,'d wrapper mode echo'))
  expect_identical(stats::tsp(dd$tables$D11),stats::tsp(series),info=paste(tag,'D11 calendar'))
  expect_error(pubd(series,'legacy',initial,final,as.integer(ca$trend),cc,mode='additive'),
    'D pass inherits its decomposition mode from the C-pass result')
  raw<-rawd(v,as.integer(ca$frequency),as.integer(ca$start_period),'legacy',initial,final,
    as.integer(ca$trend),
    rawc(v,as.integer(ca$frequency),as.integer(ca$start_period),'legacy',initial,final,
      as.integer(ca$trend),
      rawb(v,as.integer(ca$frequency),as.integer(ca$start_period),'legacy',initial,final,
        as.integer(ca$trend),ca$lsig,ca$usig,mode),ca$lsig,ca$usig),0,0)
  expect_true(all(vapply(dlabels,function(lm)
    bits(as.numeric(dd$tables[[lm]]),raw$tables[[lm]]),logical(1))),
    info=paste(tag,'wrapper tables equal raw'))
  sca<-sman[sman$case=='q-s3x9-auto',]
  sv<-read_input(sinp,sca$case)
  sout<-pubs(stats::ts(sv,start=c(2001,as.integer(sca$first_period)),
    frequency=as.integer(sca$frequency)),
    rq(sca$requests),'legacy',as.integer(sca$trend),mode=mode)
  expect_true(inherits(sout,'persephone_native_x11_selection_pass'),
    info=paste(tag,'selection wrapper class'))
  expect_identical(sout$mode,mode,info=paste(tag,'selection wrapper mode echo'))
  expect_identical(sout$d$mode,mode,info=paste(tag,'nested d mode through the wrapper'))
  # the additive default is unchanged: the same call without a mode replays
  # the additive family, not this mode.
  bdef<-pubb(series,'legacy',initial,final,as.integer(ca$trend),lower_sigma=ca$lsig,
    upper_sigma=ca$usig)
  expect_identical(bdef$mode,'additive',info=paste(tag,'wrapper defaults to additive'))
  expect_false(identical(as.numeric(bdef$tables$B7),as.numeric(bb$tables$B7)),
    info=paste(tag,'the additive default differs from this mode'))
  invisible(TRUE)
}
