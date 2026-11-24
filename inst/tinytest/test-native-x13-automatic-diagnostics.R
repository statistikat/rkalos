# M7 (item 27 cutover): the automatic route's public residual diagnostics at
# the R boundary, asserted BIT-exact against the nine-case double-JVM oracle
# (validation/evidence/native/x13-automatic-diagnostics/public-x86_64.txt,
# synced through tools/sync-core.py). The kalos C++ tier replays all nine
# blocks (tests/test_x13_automatic_diagnostics.cpp); this pin replays the
# stock_rsa3_air block through the decoder-equivalent RSA3 automatic request
# — the golden "ao" spec of the x13-automatic fixtures (the request shape
# test-native-x13-automatic.R builds) over the frozen airline capture — and
# claims the residual series (dates, count, every binary64 word), the
# standard error, the Ljung-Box triple (df 22 mirroring the (0,1,1)(0,1,1)
# landing against the 24-lag statistic) and the Doornik-Hansen triple (df 2).
# The Doornik-Hansen arm of this route is the BRIDGE capture domain
# (tramo_residual_normality, the domain the file's words match per the kalos
# ledger); the estimated route keeps its portable-arm surface untouched.
# The refusal half of the pin: a non-available status never carries a
# residual surface — the documented 21-I deferral answers verbatim.
source("helpers/fixtures.R", local = TRUE)
local({
  fit <- getFromNamespace('native_x13_fit_automatic', 'rkalos')
  tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  facts <- read.csv('fixtures/native/x13-automatic/oracle-facts.csv', colClasses = 'character')
  inputs <- read.csv('fixtures/native/x13-automatic/inputs.csv', colClasses = 'character')
  cases <- read.csv('fixtures/native/x13-automatic/cases.csv')
  fact <- function(id, key) { v <- facts$value[facts$id == id & facts$key == key]; if (length(v) == 0) "" else v[1] }
  # The golden "ao" module state, transcribed as in test-native-x13-automatic.R
  # (the C++ mirror_request of the oracle test, one for one).
  line <- strsplit(fact('ao', 'req'), ' ', fixed = TRUE)[[1]]
  req <- as.list(line[seq(2, length(line), 2)]); names(req) <- line[seq(1, length(line) - 1, 2)]
  iline <- inputs$line[inputs$id == 'ao']
  inp <- stats::setNames(as.list(sub('^[^ ]* ', '', iline)), sub(' .*', '', iline))
  otypes <- if (is.null(inp$otypes)) c('AO', 'LS', 'TC', 'SO') else strsplit(inp$otypes, ',', fixed = TRUE)[[1]]
  sigma_cut <- as.numeric(strsplit(req$sigma, ';', fixed = TRUE)[[1]])
  modes <- c(Additive = 'additive', Multiplicative = 'multiplicative',
    LogAdditive = 'log_additive', PseudoAdditive = 'pseudo_additive', Undefined = 'additive')
  # The airline input: the frozen capture the kalos oracle reads for the
  # stock_rsa3_air block (the same CSV test-native-x13-estimated.R consumes).
  input <- readLines('fixtures/native/x13-refusal-jvm-profiles/air.csv')
  start <- as.integer(strsplit(input[1L], ',', fixed = TRUE)[[1L]])
  series <- stats::ts(fixture_doubles(input[-1L]), start = start, frequency = 12)
  rsa3 <- function(...) {
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
      msr_policy = 'jdplus',
      mode = modes[[as.character(cases$mode_requested[match('ao', cases$id)])]])
    ov <- list(...); a[names(ov)] <- NULL; do.call(fit, c(a, ov))
  }

  panel <- fixture_named_rsa0_cases('fixtures/native/x13-automatic-diagnostics/public-x86_64.txt')
  block <- panel[['stock_rsa3_air']]
  expect_identical(length(panel), 9L, info = 'nine-case oracle panel complete')
  item <- function(name) tokens(block[startsWith(block, paste0(name, ' '))])

  res <- rsa3()
  expect_identical(res$status, 'available')
  expect_true(res$fitted)
  # The final description (the landing the diagnostics ride on), bit-exact.
  model <- item('model')
  expect_identical(c(res$model$log, res$model$mean), model[2:3] == 'true', info = 'description flags')
  expect_identical(unname(unlist(res$model[c('p', 'd', 'q', 'bp', 'bd', 'bq')])),
    as.integer(model[4:9]), info = 'description orders')
  expect_identical(as.integer(res$model$effective_observations),
    as.integer(item('likelihood')[2L]), info = 'effective observations')

  # The residual series: dates, count and every binary64 word bit-exact.
  residual <- item('residuals')
  tsres <- res$residuals$tsres
  expect_identical(as.integer(stats::start(tsres)), as.integer(residual[2:3]), info = 'residual dates')
  expect_identical(stats::frequency(tsres), 12, info = 'residual frequency')
  expect_identical(length(tsres), as.integer(residual[4L]), info = 'residual count')
  expect_identical(writeBin(as.numeric(tsres), raw()),
    writeBin(fixture_hex_doubles(residual[-seq_len(4L)]), raw()), info = 'residual words')
  # The test record: ser, the Ljung-Box triple and the Doornik-Hansen triple.
  expect_identical(res$residuals$ser, fixture_hex_doubles(item('residuals.ser')[2L]), info = 'ser')
  for (name in c('lb', 'doornikhansen')) {
    ref <- item(paste0('residuals.', name))
    actual <- res$residuals[[name]]
    expect_identical(c(actual$value, actual$pvalue), fixture_hex_doubles(ref[2:3]),
      info = paste('bit', name))
    expect_identical(as.numeric(actual$df), as.numeric(ref[6L]), info = paste('df', name))
  }
  expect_identical(res$residuals$lb$df, 22, info = 'lb df mirrors p+q+bp+bq of the landing')
  expect_identical(res$residuals$doornikhansen$df, 2, info = 'dh df')
  # The capability the pin justifies.
  expect_identical(rkalos::capabilities()$automatic_fit$residuals,
    c('tsres', 'ser', 'lb', 'doornikhansen'))

  # The refusal half: a documented non-available status keeps the provenance
  # echo and never carries a residual surface (no fabricated domain).
  refusal <- rsa3(calendar = TRUE)
  expect_identical(refusal$status,
    'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)')
  expect_identical(refusal$schema_version, 2L, info = 'the provenance echo stays complete')
  expect_false(refusal$fitted)
  expect_null(refusal$residuals)
})
