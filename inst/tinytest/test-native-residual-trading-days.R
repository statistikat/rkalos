# M7 (item 27 cutover gap): the residual trading-day F-tests at the R
# boundary. The BIT-claim layering is the honest line, in two halves.
# (1) The primitive native_residual_trading_days_tests() is replayed BIT-EXACT
# against the 24-case double-JVM dashboard oracle
# (fixtures/native/residual-trading-days/public-x86_64.txt, synced through
# tools/sync-core.py; the kalos C++ tier replays the same panel in
# tests/test_residual_trading_days.cpp, 192 doubles) on each case's own
# bit-frozen component series: every case feeds the fixture's embedded sa and
# i domains and observations — the inputs are bit-defined, so all four public
# names are compared as raw binary64 words, the convention of the native pin
# tests. The oracle was captured under the bridge flag profile (capture-domain
# note). (2) The PRODUCT fit routes surface diagnostics[["td.ftest.sa"]] /
# diagnostics[["td.ftest.i"]] computed from the fit's OWN finalSA /
# finalIrregular series: their component series are products of qualified but
# tolerance-qualified compositions (estimation decisions, filters and
# extraction chains pinned at bands, not bit-frozen like this fixture's
# inputs), so NO fixture-bit claim is made for product-fit scalars — the route
# wiring is pinned identical (bit-identical doubles) to applying the primitive
# to that same fit's own component series, plus finiteness and consistency.
# The pinned capture option is flast = 0: -all and -last share one fit.
source("helpers/fixtures.R", local = TRUE)
local({
  # (1) The bit-exact primitive replay: all 24 cases, all four names.
  primitive <- getFromNamespace('native_residual_trading_days_tests', 'rkalos')
  tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  panel <- fixture_named_rsa0_cases('fixtures/native/residual-trading-days/public-x86_64.txt')
  expect_identical(length(panel), 24L, info = '24-case oracle panel complete')
  words <- 0L
  for (id in names(panel)) {
    block <- panel[[id]]
    line <- function(prefix) tokens(block[startsWith(block, paste0(prefix, " "))])
    # Each component domain: year, one-based month and the monthly frequency,
    # followed by the bit-frozen observations of that component. The
    # parameters/likelihood/coefficients witness lines are never fed back.
    component <- function(tag) {
      at <- which(startsWith(block, paste0("domain ", tag, " ")))
      stopifnot(length(at) == 1L)
      domain <- tokens(block[at])
      observations <- tokens(block[at + 1L])
      stopifnot(observations[1L] == "observations")
      list(year = as.integer(domain[3L]), month = as.integer(domain[4L]),
           frequency = as.integer(domain[5L]), values = fixture_hex_doubles(observations[-1L]))
    }
    sa <- component("sa"); ir <- component("i")
    expect_identical(c(sa$year, sa$month, sa$frequency, ir$year, ir$month, ir$frequency),
      c(ir$year, ir$month, 12L, ir$year, ir$month, 12L), info = paste(id, 'shared monthly domain'))
    res <- primitive(sa$values, ir$values, sa$year, sa$month)
    expect_identical(names(res), c('td-sa-all', 'td-sa-last', 'td-i-all', 'td-i-last'),
      info = paste(id, 'four public names'))
    for (name in names(res)) {
      reference <- line(name)
      expect_identical(names(res[[name]]), c('value', 'pvalue'), info = paste(id, name, 'shape'))
      expect_identical(writeBin(c(res[[name]]$value, res[[name]]$pvalue), raw()),
        writeBin(fixture_hex_doubles(reference[2:3]), raw()), info = paste(id, name, 'bit'))
      words <- words + 2L
    }
    # flast = 0: the JVM publishes the one fit under both names of a component.
    expect_identical(res[['td-sa-last']], res[['td-sa-all']], info = paste(id, 'sa -last = -all'))
    expect_identical(res[['td-i-last']], res[['td-i-all']], info = paste(id, 'i -last = -all'))
  }
  expect_identical(words, 192L, info = 'all 192 oracle doubles replayed bit-exact')
  # The primitive rides the capability surface.
  caps <- rkalos::capabilities()
  expect_true('residual_trading_days_tests' %in% caps$primitives)
  expect_identical(caps$td_tests$fields, c('td.ftest.sa', 'td.ftest.i'))
  expect_identical(caps$td_tests$names, c('td-sa-all', 'td-sa-last', 'td-i-all', 'td-i-last'))
  # Fail-closed guards (never a fallback): at or below the eight-observation
  # degrees-of-freedom bound, non-finite data, and an invalid month.
  expect_error(primitive(rep(1, 8), rep(1, 8), 1949, 1), 'degrees of freedom')
  expect_error(primitive(c(1, NA_real_), c(1, 2), 1949, 1), 'finite')
  expect_error(primitive(as.numeric(AirPassengers), as.numeric(AirPassengers),
                         1949, 13), 'monthly domain')
})
local({
  # (2) Product-route wiring: the attached pair equals the primitive applied
  # to the SAME fit's own component series, bit for bit. Bit-exactness of the
  # scalars against the fixture is NOT claimed here: the product composition
  # inputs (the fitted decomposition) are tolerance-qualified, not the
  # bit-frozen series this fixture embeds — only the wiring is pinned.
  primitive <- getFromNamespace('native_residual_trading_days_tests', 'rkalos')
  air <- AirPassengers
  check_route <- function(fit, label) {
    expect_identical(names(fit$diagnostics), c('combined_seasonality', 'td.ftest.sa', 'td.ftest.i'),
      info = paste(label, 'fields attached'))
    start <- as.integer(stats::start(fit$finalSA))
    expect_identical(start, as.integer(stats::start(air)), info = paste(label, 'component domain'))
    expected <- primitive(as.numeric(fit$finalSA), as.numeric(fit$finalIrregular),
                          start[1L], start[2L])
    for (pair in list(c('td.ftest.sa', 'td-sa-all'), c('td.ftest.i', 'td-i-all'))) {
      attached <- fit$diagnostics[[pair[1L]]]
      reference <- expected[[pair[2L]]]
      expect_identical(attached[c('value', 'pvalue')], reference,
        info = paste(label, pair[1L], 'wiring identical'))
      expect_true(is.finite(attached$value) && attached$value >= 0 &&
                    is.finite(attached$pvalue) && attached$pvalue >= 0 && attached$pvalue <= 1,
        info = paste(label, pair[1L], 'finite and consistent'))
    }
  }
  x13 <- rkalos::native_x13_estimated_fit(air, 'RSA0')
  expect_identical(x13$status, 'available')
  check_route(x13, 'x13_estimated_rsa0')
  full <- rkalos::native_tramoseats_fit(air, 'RSAfull')
  expect_identical(full$status, 'ok')
  check_route(full, 'tramoseats_rsafull')
  restricted <- rkalos::native_x13_fit(air, seasonal = rep('s3x5', 12L), trend_length = 13L)
  expect_identical(restricted$status, 'available')
  check_route(restricted, 'x13_restricted')
  # The refusal half: a non-available status never carries the surface. The
  # documented 21-I calendar deferral stays provenance-only (the same golden
  # request shape test-native-x13-automatic-diagnostics.R drives).
  facts <- read.csv('fixtures/native/x13-automatic/oracle-facts.csv', colClasses = 'character')
  line <- strsplit(facts$value[facts$id == 'ao' & facts$key == 'req'][1L], ' ', fixed = TRUE)[[1L]]
  req <- as.list(line[seq(2, length(line), 2)]); names(req) <- line[seq(1, length(line) - 1, 2)]
  refusal <- rkalos::native_x13_fit_automatic(series = air,
    transform = list(aicc = as.numeric(req$aicdiff)),
    automodel = list(max_d = 2L, max_bd = 1L, ub1 = 1 / as.numeric(req$ub1),
      ub2 = 1 / as.numeric(req$ub2), cancel = as.numeric(req$cancel), precision = 1e-7,
      models_count = 5L, max_p = 2L, max_q = 2L, max_bp = 1L, max_bq = 1L,
      balanced = identical(req$balanced, 'true'), mixed = identical(req$mixed, 'true'),
      accept_white_noise = FALSE, eps = 1e-7),
    outliers = list(used = identical(req$outliers, 'on'),
      types = if (identical(req$outliers, 'on')) c('AO', 'LS', 'TC', 'SO') else NULL,
      span = NULL, max_round = as.integer(req$maxiter), max_outliers = 30L,
      tcrate = as.numeric(req$tcrate)),
    ami = list(precision = 1e-7, mean_precision = 1e-5, va = 0,
      reduce_va = as.numeric(req$predcv), ljung_box_limit = as.numeric(req$lblimit),
      check_mu = TRUE, accept_airline = FALSE, ur_limit = 0.95),
    regression = list(mean_test_initial = 1.96, mean_test_final = 1, tsig = 1),
    calendar = TRUE)
  expect_identical(refusal$status,
    'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)')
  expect_false(refusal$fitted)
  expect_null(refusal$diagnostics)
})
local({
  # The omission half (not the refusal half): the automatic route delivers an
  # available fit whose published composition collapses — the 'quick' golden
  # request under mode 'pseudo_additive' exports identically-zero finalSA and
  # finalIrregular (the route's own published tables, witnessed by the D8
  # panel of test-native-x13-automatic.R). The pinned machine declines to fit
  # such a component (no usable residual variance), and the attach records
  # the decline: the fields stay absent, the fit is delivered untouched —
  # never a fabricated diagnostic and never a failed fit.
  attach <- getFromNamespace('.attach_td_tests', 'rkalos')
  zero <- stats::ts(rep(0, 96), start = c(1990, 1), frequency = 12)
  expect_identical(attach(list(a = 1), zero, zero, c(1990L, 1L)), list(a = 1),
    info = 'declined components leave the diagnostics untouched')
  short <- stats::ts(as.numeric(AirPassengers)[1:8], start = c(1949, 1), frequency = 12)
  expect_identical(attach(list(a = 1), short, short, c(1949L, 1L)), list(a = 1),
    info = 'components inside the degrees-of-freedom guard omit the fields')
  facts <- read.csv('fixtures/native/x13-automatic/oracle-facts.csv', colClasses = 'character')
  inputs <- read.csv('fixtures/native/x13-automatic/inputs.csv', colClasses = 'character')
  line <- strsplit(facts$value[facts$id == 'quick' & facts$key == 'req'][1L], ' ', fixed = TRUE)[[1L]]
  req <- as.list(line[seq(2, length(line), 2)]); names(req) <- line[seq(1, length(line) - 1, 2)]
  iline <- inputs$line[inputs$id == 'quick']
  inp <- stats::setNames(as.list(sub('^[^ ]* ', '', iline)), sub(' .*', '', iline))
  start <- as.integer(strsplit(req$start, ':', fixed = TRUE)[[1]])
  fit <- rkalos::native_x13_fit_automatic(series = fixture_doubles(strsplit(inp$values, ';',
    fixed = TRUE)[[1]][nchar(strsplit(inp$values, ';', fixed = TRUE)[[1]]) > 0]),
    transform = list(aicc = as.numeric(req$aicdiff)),
    automodel = list(max_d = 2L, max_bd = 1L, ub1 = 1 / as.numeric(req$ub1),
      ub2 = 1 / as.numeric(req$ub2), cancel = as.numeric(req$cancel), precision = 1e-7,
      models_count = 5L, max_p = 2L, max_q = 2L, max_bp = 1L, max_bq = 1L,
      balanced = identical(req$balanced, 'true'), mixed = identical(req$mixed, 'true'),
      accept_white_noise = FALSE, eps = 1e-7),
    outliers = list(used = identical(req$outliers, 'on'),
      types = if (identical(req$outliers, 'on')) c('AO', 'LS', 'TC', 'SO') else NULL,
      span = NULL, max_round = as.integer(req$maxiter), max_outliers = 30L,
      tcrate = as.numeric(req$tcrate)),
    ami = list(precision = 1e-7, mean_precision = 1e-5, va = 0,
      reduce_va = as.numeric(req$predcv), ljung_box_limit = as.numeric(req$lblimit),
      check_mu = TRUE, accept_airline = FALSE, ur_limit = 0.95),
    regression = list(mean_test_initial = 1.96, mean_test_final = 1, tsig = 1),
    mode = 'pseudo_additive', frequency = 12L, start_period = start[2L], first_year = start[1L])
  expect_identical(fit$status, 'available', info = 'the zero-composition fit is still delivered')
  expect_true(all(as.numeric(fit$finalSA) == 0) && all(as.numeric(fit$finalIrregular) == 0),
    info = 'the route publishes the collapsed composition the attach declines')
  expect_null(fit$diagnostics$td.ftest.sa)
  expect_null(fit$diagnostics$td.ftest.i)
})
