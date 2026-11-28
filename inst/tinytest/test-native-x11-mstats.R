# M14 (c3706ea re-pin): the X-11 diagnostic battery — m1..m11 through the
# pinned getM clamp, q, qm2 and the two MDiagnostics grades — at the R boundary,
# mirrored on BOTH routes the pinned X13Document publishes it on: the automatic
# lane (native_x13_fit_automatic) and the estimated-template lane
# (native_x13_estimated_fit), established from the synced core header
# (include/kalos/x11_mstats.hpp, the free entry kalos::x13_mstats over the
# finished X13FitResult both routes return) by a live compile check of the
# vendored headers, not by prose. The ENGINE's claim against the JVM lives in
# the core: tests/test_x11_mstats.cpp replays the twelve-case double-JVM panel
# (the three stock RSA3 template fits and the nine stock RSA0/RSA1/RSA2c
# estimated fits on the x13-refusal-jvm-profiles calendars) with bit claims for
# m4 (12/12), both grades and the M6 preset invariant against the fit's own
# D9 default-table global, and a max(1e-11, 5e-12*|want|) tolerance tier for
# the twelve regression-seam-coupled words (measured worst 2.4e-13, m6 9.4e-13,
# q/qm2 5.1e-14 — residual ulp noise inherited from the upstream A8/B1 seam,
# not the Mstatistics layer). What this pin asserts is PLUMBING FIDELITY to
# that core entry: the two gold envelopes deliver the core entry's output
# BIT-exactly as recorded live from the installed 9034 lanes (serial and OMP
# recordings bit-identical), under the names the core uses (status, m1..m11,
# q, qm2, q_grade, qm2_grade). The twelve-case cross-check against the synced
# evidence battery (the same file the core replay reads) re-derives the core's
# tiers at this boundary — bit for m4 and the grades, tolerance for the rest —
# and is not an independent bit claim against the JVM. The refusal arms pin the
# core entry's explicit refuse-answers passing through verbatim (the
# pseudo-additive status with the untouched default battery, never a partial
# one) and the absence contract: a non-available fit carries no battery at all.
source("helpers/fixtures.R", local = TRUE)
local({
  fit_automatic <- getFromNamespace('native_x13_fit_automatic', 'rkalos')
  fit_estimated <- getFromNamespace('native_x13_estimated_fit', 'rkalos')
  tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  hexof <- function(x) vapply(as.numeric(x), function(d)
    paste(rev(as.character(writeBin(d, raw()))), collapse = ""), character(1))

  # The synced double-JVM battery: twelve cases, each with the thirteen
  # published words and the two grades.
  panel <- local({
    lines <- readLines('fixtures/native/x11-mstats/battery-x86_64.txt')
    out <- list(); id <- NULL
    for (line in lines) {
      tk <- tokens(line)
      if (tk[1] == 'case') { id <- tk[2]; out[[id]] <- list() }
      else if (tk[1] == 'info') out[[id]]$mode <- tk[2]
      else if (tk[1] == 'mstats') out[[id]]$words <- tk[-1]
      else if (tk[1] == 'grade') out[[id]]$grades <- tk[-1]
    }
    out
  })
  expect_equal(length(panel), 12L, info = 'synced battery panel complete')

  # The mirrored surface: the core entry's publication, name-for-name.
  fields <- c(paste0('m', 1:11), 'q', 'qm2')
  shape <- c('status', fields, 'q_grade', 'qm2_grade')

  # The golden "ao" module state (the stock_rsa3_* requests' base), transcribed
  # as in test-native-x13-forecast-errors.R, over the three battery calendars.
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
  battery_of <- function(id) {
    res <- if (startsWith(id, 'stock_rsa3_')) rsa3(cal[[sub('^.*_', '', id)]]) else
      fit_estimated(cal[[sub('^.*_', '', id)]],
        template = sub('^stock_(rsa[0-9]c?)_.*', '\\1', id))
    expect_identical(res$status, 'available', info = paste(id, 'fit status'))
    expect_false(is.null(res$x11_mstats), info = paste(id, 'battery present'))
    expect_identical(names(res$x11_mstats), shape, info = paste(id, 'battery names'))
    expect_identical(res$x11_mstats$status, 'available', info = paste(id, 'battery status'))
    res$x11_mstats
  }

  # Gold envelope 1 (automatic lane, stock_rsa3_air): the fourteen words are
  # the core entry's output recorded live from the installed 9034 lanes,
  # delivered BIT-exact (plumbing fidelity; the core replay certifies the
  # engine behind them — bit for m4 and the grades, tier 2.4e-13 for the rest).
  m <- battery_of('stock_rsa3_air')
  recorded_rsa3_air <- c(
    "3fa40f77cf735d97", "3fa3c7911fde7d48", "0000000000000000", "3febfdc8d332eebf",
    "3fd106a55fbbe161", "3fe528e3c9a30fcd", "3fc9e73ec5e0b286", "3fdbfbfb4b0e4c0b",
    "3fd691e71be655f5", "3fdd0eeec1efe956", "3fda07b032b9ccf1", "3fd09bd347baae56",
    "3fd293b84cb4da4e")
  expect_identical(writeBin(as.numeric(m[fields]), raw()),
    writeBin(fixture_hex_doubles(recorded_rsa3_air), raw()),
    info = 'automatic battery words are the core entry verbatim')
  expect_identical(m$q_grade, 'Good')
  expect_identical(m$qm2_grade, 'Good')

  # Gold envelope 2 (estimated lane, stock_rsa0_air): the same pin on the
  # estimated-template route.
  me <- battery_of('stock_rsa0_air')
  recorded_rsa0_air <- c(
    "3fb0300a7c5c82da", "3fa937adce853070", "3fada5bb222383e0", "3fd567b7b0905c38",
    "3fd4560fa4c660b3", "3ff2cf7f57e301c7", "3fe925f965812fb6", "3fec80f44c98af71",
    "3feaae7c10d613c9", "3ff316fd746ccbbb", "3ff2b976f6015174", "3fdf004c07d196e9",
    "3fe1707ea81ff2e8")
  expect_identical(writeBin(as.numeric(me[fields]), raw()),
    writeBin(fixture_hex_doubles(recorded_rsa0_air), raw()),
    info = 'estimated battery words are the core entry verbatim')

  # The twelve-case tier cross-check against the synced evidence the core
  # replay pins: m4 and the grades bit-exact, the other twelve words inside
  # the core's max(1e-11, 5e-12*|want|) tier (worst measured at this boundary:
  # 3.226e-13, stock_rsa0_ao m6 — inside the core's own m6 reading of 9.4e-13).
  for (id in names(panel)) {
    b <- if (identical(id, 'stock_rsa3_air')) m else if (identical(id, 'stock_rsa0_air')) me else
      battery_of(id)
    want <- fixture_hex_doubles(panel[[id]]$words)
    got <- as.numeric(b[fields])
    expect_identical(hexof(b$m4), panel[[id]]$words[4], info = paste(id, 'm4 bit anchor'))
    expect_identical(c(b$q_grade, b$qm2_grade), panel[[id]]$grades, info = paste(id, 'grades'))
    tol <- pmax(1e-11, 5e-12 * abs(want))
    ok <- abs(got - want) <= tol
    expect_true(all(ok[-4]), info = paste(id, 'tier words within the core tier'))
  }

  # The pseudo-additive refuse-answer passes through verbatim: the fit lands,
  # and the battery answers its explicit status with every word at the
  # untouched default — never a partial battery.
  padd <- rsa3(cal$air, mode = modes[['PseudoAdditive']])
  expect_identical(padd$status, 'available')
  expect_identical(padd$x11_mstats$status, 'unsupported:pseudo_additive')
  expect_identical(as.numeric(padd$x11_mstats[fields]),
    c(rep(-1, 11), NaN, NaN), info = 'refused battery stays untouched')
  expect_identical(padd$x11_mstats$q_grade, '')
  expect_identical(padd$x11_mstats$qm2_grade, '')

  # The absence half of the contract: a documented non-available status never
  # carries the battery (the core entry throws outside the available domain;
  # the boundary keeps the field missing, never estimated on the side).
  refusal <- rsa3(cal$air, calendar = TRUE)
  expect_identical(refusal$status,
    'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)')
  expect_null(refusal$x11_mstats, info = 'absent on refusal')
})
