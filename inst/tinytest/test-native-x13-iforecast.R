# M15 (1e2e3c8 re-pin): the published irregular-component forecast i_f — the
# constant placeholder the JVM X13Extractor manufactures on the finals
# trend-forecast (D12a) domain — at the R boundary, mirrored on BOTH envelopes
# the core entry lands on: the automatic lane (native_x13_fit_automatic) and
# the estimated-template lane (native_x13_estimated_fit), established from the
# synced core header (include/kalos/x13_automatic.hpp, the free entry
# kalos::x13_automatic_irregular_forecast over the X13AutomaticFitResult both
# routes return) by a live compile check of the vendored headers, not by prose.
# The ENGINE's claim against the JVM lives in the core: tests/
# test_x13_automatic_iforecast.cpp replays the thirteen-case double-JVM panel
# (twelve stock fits plus the zero-horizon RSA3 arm) with BIT claims for the
# mode echo, the log flag, the domain shape and all 144 published i_f words
# (12/12 word-exact), pinned against the X13Extractor bytecode dossier (the
# lambda$new$28 chain getFinals/getD12a/size/getMode/isMultiplicative): i_f is
# NOT model-based — 1.0 under the JVM DecompositionMode predicate
# (Multiplicative OR LogAdditive), +0.0 otherwise (Additive, PseudoAdditive),
# no mean correction (rsa3_ao answers 1.0), no back-transform (the log fit
# answers 1.0 verbatim), and the D12a words are DOMAIN-only, never part of i_f.
# What this pin asserts is PLUMBING FIDELITY to that core entry: the two gold
# envelopes deliver the entry's output BIT-exactly as recorded live from the
# installed 9035 lanes (the serial and OMP recordings bit-identical), under
# the names the entry gives its two members (status, values), plus the live
# constant semantics by mode (mode_effective echo and the recorded words for
# multiplicative, log_additive — answering 1.0 — and additive and
# pseudo_additive — answering +0.0, the JVM predicate, unlike the X11 context's
# isMultiplicative, excluding PseudoAdditive), both reachable status arms
# ("available", "no_requested_predictions" with the synced null-row echo; the
# finals-null "unavailable" arm is the entry's defensive answer, unreachable
# here and asserted honestly as never fabricated), the refusal absence, and
# the cross-check against the synced iforecast evidence rows: hex-BIT for the
# i_f words, since the reference publishes constants (not a tolerance claim).
source("helpers/fixtures.R", local = TRUE)
local({
  fit_automatic <- getFromNamespace('native_x13_fit_automatic', 'rkalos')
  fit_estimated <- getFromNamespace('native_x13_estimated_fit', 'rkalos')
  tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]

  # The synced double-JVM battery: twelve stock cases with the i_f row (start
  # year, month, count, words) and the zero-horizon RSA3 arm with null words.
  panel <- local({
    lines <- readLines('fixtures/native/x13-automatic-iforecast/iforecast-x86_64.txt')
    out <- list(); id <- NULL
    for (line in lines) {
      tk <- tokens(line)
      if (tk[1] == 'case') { id <- tk[2]; out[[id]] <- list() }
      else if (tk[1] == 'info') out[[id]]$mode <- tk[2]
      else if (tk[1] == 'i_f') out[[id]]$i_f <- tk[-1]
    }
    out
  })
  expect_equal(length(panel), 13L, info = 'synced iforecast panel complete')
  expect_identical(panel[['stock_rsa3_air_h0']]$i_f, 'null',
    info = 'the zero-horizon JVM row is null (the no_requested_predictions gate)')

  # The mirrored surface: the core entry's two members, name-for-name.
  shape <- c('status', 'values')
  multiplicative_modes <- c('Multiplicative', 'LogAdditive')

  # The golden "ao" module state (the stock_rsa3_* requests' base), transcribed
  # as in test-native-x11-mstats.R, over the three battery calendars.
  facts <- read.csv('fixtures/native/x13-automatic/oracle-facts.csv', colClasses = 'character')
  inputs <- read.csv('fixtures/native/x13-automatic/inputs.csv', colClasses = 'character')
  fact <- function(id, key) { v <- facts$value[facts$id == id & facts$key == key]; if (length(v) == 0) "" else v[1] }
  line <- strsplit(fact('ao', 'req'), ' ', fixed = TRUE)[[1]]
  req <- as.list(line[seq(2, length(line), 2)]); names(req) <- line[seq(1, length(line) - 1, 2)]
  iline <- inputs$line[inputs$id == 'ao']
  inp <- stats::setNames(as.list(sub('^[^ ]* ', '', iline)), sub(' .*', '', iline))
  otypes <- if (is.null(inp$otypes)) c('AO', 'LS', 'TC', 'SO') else strsplit(inp$otypes, ',', fixed = TRUE)[[1]]
  sigma_cut <- as.numeric(strsplit(req$sigma, ';', fixed = TRUE)[[1]])
  modes <- c(Additive = 'additive', Multiplicative = 'multiplicative',
    LogAdditive = 'log_additive', PseudoAdditive = 'pseudo_additive', Undefined = 'additive')
  read_calendar <- function(stem) {
    input <- readLines(paste0('fixtures/native/x13-refusal-jvm-profiles/', stem, '.csv'))
    start <- as.integer(strsplit(input[1L], ',', fixed = TRUE)[[1L]])
    stats::ts(fixture_doubles(input[-1L]), start = start, frequency = 12)
  }
  cal <- stats::setNames(lapply(c('air', 'ao', 'seas'), read_calendar), c('air', 'ao', 'seas'))
  rsa3 <- function(series, ...) {
    a <- list(series = series, transform = list(aicc = as.numeric(req$aicdiff)),
      automodel = list(max_d = 2L, max_bd = 1L, ub1 = 1 / as.numeric(req$ub1),
        ub2 = 1 / as.numeric(req$ub2), cancel = as.numeric(req$cancel), precision = 1e-7,
        models_count = 5L, max_p = 2L, max_q = 2L, max_bp = 1L, max_bq = 1L,
        balanced = identical(req$balanced, 'true'), mixed = identical(req$mixed, 'true'),
        accept_white_noise = FALSE, eps = 1e-7),
      outliers = list(used = identical(req$outliers, 'on'),
        types = if (identical(req$outliers, 'on')) otypes else NULL, span = NULL,
        max_round = as.integer(req$maxiter), max_outliers = 30L, tcrate = as.numeric(req$tcrate)),
      ami = list(precision = 1e-7, mean_precision = 1e-5, va = 0, reduce_va = as.numeric(req$predcv),
        ljung_box_limit = as.numeric(req$lblimit), check_mu = TRUE, accept_airline = FALSE,
        ur_limit = 0.95),
      regression = list(mean_test_initial = 1.96, mean_test_final = 1, tsig = 1),
      sigma = 'accurate', lower_sigma = sigma_cut[1], upper_sigma = sigma_cut[2],
      msr_policy = 'jdplus', mode = modes[['Undefined']], forecast = 12L)
    ov <- list(...); a[names(ov)] <- NULL
    do.call(fit_automatic, c(a, ov))
  }
  landed <- list()
  battery_of <- function(id) {
    res <- if (startsWith(id, 'stock_rsa3_')) rsa3(cal[[sub('^.*_', '', id)]]) else
      fit_estimated(cal[[sub('^.*_', '', id)]],
        template = sub('^stock_(rsa[0-9]c?)_.*', '\\1', id))
    expect_identical(res$status, 'available', info = paste(id, 'fit status'))
    expect_false(is.null(res$irregular_forecast), info = paste(id, 'forecast present'))
    expect_identical(names(res$irregular_forecast), shape, info = paste(id, 'forecast names'))
    expect_identical(res$irregular_forecast$status, 'available', info = paste(id, 'forecast status'))
    landed[[id]] <<- res
    res
  }

  # Gold envelope 1 (automatic lane, stock_rsa3_air): the twelve i_f words are
  # the core entry's output recorded live from the installed 9035 lanes
  # (serial and OMP recordings bit-identical; the recording reads, verbatim,
  # twelve times "3ff0000000000000"), delivered BIT-exact under the entry's
  # own member names — plumbing fidelity; the constant 1.0 verbatim, not
  # exponentiated, on this log fit (mode_effective multiplicative).
  m <- battery_of('stock_rsa3_air')
  expect_identical(m$mode_effective, 'multiplicative', info = 'rsa3_air mode echo')
  expect_identical(m$model$log, TRUE, info = 'rsa3_air is a log fit')
  recorded_rsa3_air <- c(
    "3ff0000000000000", "3ff0000000000000", "3ff0000000000000", "3ff0000000000000",
    "3ff0000000000000", "3ff0000000000000", "3ff0000000000000", "3ff0000000000000",
    "3ff0000000000000", "3ff0000000000000", "3ff0000000000000", "3ff0000000000000")
  expect_identical(writeBin(as.numeric(m$irregular_forecast$values), raw()),
    writeBin(fixture_hex_doubles(recorded_rsa3_air), raw()),
    info = 'automatic i_f words are the core entry verbatim')
  expect_identical(as.numeric(m$irregular_forecast$values), rep(1, 12L),
    info = 'the log fit answers the constant 1.0 (no back-transform)')

  # Gold envelope 2 (estimated lane, stock_rsa0_air): the same pin on the
  # estimated-template route — the twelve recorded words are twelve times
  # "0000000000000000", the byte pattern of +0.0 (a negative zero would read
  # 8000000000000000): the additive-level landing answers the reference's
  # +0.0 sign as published.
  me <- battery_of('stock_rsa0_air')
  expect_identical(me$mode_effective, 'additive', info = 'rsa0_air mode echo')
  recorded_rsa0_air <- c(
    "0000000000000000", "0000000000000000", "0000000000000000", "0000000000000000",
    "0000000000000000", "0000000000000000", "0000000000000000", "0000000000000000",
    "0000000000000000", "0000000000000000", "0000000000000000", "0000000000000000")
  expect_identical(writeBin(as.numeric(me$irregular_forecast$values), raw()),
    writeBin(fixture_hex_doubles(recorded_rsa0_air), raw()),
    info = 'estimated i_f words are the core entry verbatim (+0.0 sign included)')

  # The twelve-case cross-check against the synced evidence the core replay
  # pins: the i_f words are BIT-identical for every case (constants — the
  # reference publishes 3ff0000000000000 twelve times on the multiplicative
  # domain and 0000000000000000 twelve times otherwise; measured live, the
  # twelve-case words came out identical in both lanes), the count equals the
  # evidence D12a-domain count, and the answered constant matches the JVM
  # DecompositionMode predicate read off the evidence info line.
  for (id in setdiff(names(panel), 'stock_rsa3_air_h0')) {
    b <- if (identical(id, 'stock_rsa3_air')) m$irregular_forecast else
      if (identical(id, 'stock_rsa0_air')) me$irregular_forecast else
      battery_of(id)$irregular_forecast
    row <- panel[[id]]$i_f
    expect_identical(length(b$values), as.integer(row[3]), info = paste(id, 'D12a-domain count'))
    expect_identical(writeBin(as.numeric(b$values), raw()),
      writeBin(fixture_hex_doubles(row[-seq_len(3L)]), raw()),
      info = paste(id, 'i_f words bit-identical to the evidence row'))
    want <- if (panel[[id]]$mode %in% multiplicative_modes) 1 else 0
    expect_identical(as.numeric(b$values), rep(want, length(b$values)),
      info = paste(id, 'constant follows the JVM DecompositionMode predicate'))
  }

  # Constant semantics by mode, established LIVE on the automatic lane (both
  # lane recordings identical): the explicit log_additive request lands with
  # the log-additive results-builder note and answers 1.0 — the JVM predicate
  # counts LogAdditive — while the explicit pseudo_additive request lands with
  # mode_effective pseudo_additive and answers +0.0: the JVM
  # DecompositionMode.isMultiplicative() predicate, NOT the X11 context's
  # isMultiplicative() (which counts PseudoAdditive), measured to separate at
  # exactly this boundary.
  loga <- rsa3(cal$air, mode = modes[['LogAdditive']])
  expect_identical(loga$mode_effective, 'log_additive', info = 'log_additive mode echo')
  expect_identical(as.numeric(loga$irregular_forecast$values), rep(1, 12L),
    info = 'log_additive answers 1.0 (JVM predicate counts LogAdditive)')
  padd <- rsa3(cal$air, mode = modes[['PseudoAdditive']])
  expect_identical(padd$mode_effective, 'pseudo_additive', info = 'pseudo_additive mode echo')
  expect_identical(as.numeric(padd$irregular_forecast$values), rep(0, 12L),
    info = 'pseudo_additive answers +0.0 (the fit lands; only i_f refuses the mode)')
  # The level (mixed of logs out) arm: suppressing the log decision on the same
  # calendar lands additive and answers +0.0 — the 0.0 half of the battery on
  # the automatic lane itself.
  level <- rsa3(cal$air, transform = list(aicc = 1e9))
  expect_identical(level$model$log, FALSE, info = 'level arm is not a log fit')
  expect_identical(level$mode_effective, 'additive', info = 'level arm mode echo')
  expect_identical(as.numeric(level$irregular_forecast$values), rep(0, 12L),
    info = 'the level arm answers +0.0 on the automatic lane')

  # The zero-horizon arm passes through verbatim: d12a empty answers
  # "no_requested_predictions" with an empty vector (the JVM zero-horizon row
  # in the synced battery answers i_f null, pinned above).
  zero <- rsa3(cal$air, forecast = 0L)
  expect_identical(zero$status, 'available')
  expect_identical(zero$irregular_forecast$status, 'no_requested_predictions')
  expect_identical(zero$irregular_forecast$values, numeric())

  # The "unavailable" arm (the JVM finals-null answer) is the entry's
  # defensive half: an available landing always carries its finals, so NO
  # boundary-reachable fit produces it. Asserted honestly: every landed
  # envelope here carries its finals tables, and the arm is never fabricated
  # — the only arms the boundary shows are the reachable two.
  arms <- vapply(landed, function(res) res$irregular_forecast$status, character(1))
  expect_true(all(vapply(landed, function(res) !is.null(res$finals), logical(1))),
    info = 'every available landing carries its finals')
  expect_identical(unique(arms), 'available',
    info = 'the boundary shows only reachable arms; "unavailable" is never fabricated')

  # The refusal half of the absence contract: a documented non-available status
  # never carries the surface (the core entry throws outside the available
  # domain; the boundary keeps the field missing, never estimated on the side).
  refusal <- rsa3(cal$air, calendar = TRUE)
  expect_identical(refusal$status,
    'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)')
  expect_null(refusal$irregular_forecast, info = 'absent on refusal')
})
