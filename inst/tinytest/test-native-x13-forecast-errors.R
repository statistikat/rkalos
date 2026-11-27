# M13 (69f83be re-pin): the automatic route's public forecast standard errors
# (y_ef) at the R boundary, mirroring the estimated route's series_forecast_stdev
# pair name-for-name. The ENGINE's claim against the JVM lives in the core:
# kalos tests/test_x13_automatic_diagnostics.cpp replays the nine-case panel's
# y_ef rows (144 predictions, h in {1,7,12,24}) through
# kalos::x13_automatic_forecast_errors at a tolerance tier of 1e-12 relative
# against the dual-JVM capture words (measured worst 1.9e-14; the coefficients
# rows join the bit claims there, the y_ef rows explicitly do not). What this
# pin asserts is PLUMBING FIDELITY to that core entry: the RSA3 AirPassengers
# envelope carries the core entry's certified words bit-exactly, exactly as
# recorded below from the core-verified output, with the pair's dating and the
# entry's own status strings passed through verbatim. The tolerance-tier
# cross-check against the synced evidence y_ef words (same file the core
# replay reads) re-derives the core's tier at this boundary; it is not an
# independent bit claim against the JVM. The refusal half pins the absence
# contract: outside the complete-landing domain the fields never appear.
source("helpers/fixtures.R", local = TRUE)
local({
  fit <- getFromNamespace('native_x13_fit_automatic', 'rkalos')
  tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  facts <- read.csv('fixtures/native/x13-automatic/oracle-facts.csv', colClasses = 'character')
  inputs <- read.csv('fixtures/native/x13-automatic/inputs.csv', colClasses = 'character')
  cases <- read.csv('fixtures/native/x13-automatic/cases.csv')
  fact <- function(id, key) { v <- facts$value[facts$id == id & facts$key == key]; if (length(v) == 0) "" else v[1] }
  # The golden "ao" module state (the stock_rsa3_air block's request), transcribed
  # as in test-native-x13-automatic-diagnostics.R.
  line <- strsplit(fact('ao', 'req'), ' ', fixed = TRUE)[[1]]
  req <- as.list(line[seq(2, length(line), 2)]); names(req) <- line[seq(1, length(line) - 1, 2)]
  iline <- inputs$line[inputs$id == 'ao']
  inp <- stats::setNames(as.list(sub('^[^ ]* ', '', iline)), sub(' .*', '', iline))
  otypes <- if (is.null(inp$otypes)) c('AO', 'LS', 'TC', 'SO') else strsplit(inp$otypes, ',', fixed = TRUE)[[1]]
  sigma_cut <- as.numeric(strsplit(req$sigma, ';', fixed = TRUE)[[1]])
  modes <- c(Additive = 'additive', Multiplicative = 'multiplicative',
    LogAdditive = 'log_additive', PseudoAdditive = 'pseudo_additive', Undefined = 'additive')
  # The airline input: the frozen capture of the stock_rsa3_air block.
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

  res <- rsa3(forecast = 24L)
  expect_identical(res$status, 'available')
  errors <- res$series_forecast_stdev
  expect_identical(res$series_forecast_stdev_status, 'available')
  # The core entry's certified words (the recording this boundary ships): the
  # envelope delivers them BIT-exact, unmodified (plumbing fidelity; the core
  # replay certifies the engine behind them against the JVM panel at 1e-12).
  recorded <- c(
    "4030bf7d86640372", "403273be89c3c1ed", "403758a6aa0a5a1b", "403a64bfd12485cd",
    "403d8eca41471ba8", "4042255613ae830d", "4046258677542ef0", "4047472ed083cc80",
    "404474be85975272", "40430dbb86331352", "404129e7d5fa61cb", "4043cb0430b5755a",
    "4046bd9cfc57c9c1", "4046cd666eeb7a19", "404b0dafcf2e05d6", "404d2e0e37f1a96c",
    "404f85e3c749a8b9", "4052cf17ca5806d7", "40566e708bd57856", "4057216ca3084a08",
    "4054007533d54dd2", "405260640ba1b276", "40505be317aa0f1e", "4052abeb3fe6165c")
  expect_identical(writeBin(as.numeric(errors), raw()),
    writeBin(fixture_hex_doubles(recorded), raw()), info = 'y_ef words are the core entry verbatim')
  expect_identical(stats::frequency(errors), 12)
  expect_identical(length(errors), 24L)
  # Tolerance-tier cross-check against the synced evidence the core replay pins:
  # every word agrees with the JVM capture at 1e-12 relative (the core tier;
  # its measured worst over the panel is 1.9e-14), never claimed bit-exact.
  panel <- fixture_named_rsa0_cases('fixtures/native/x13-automatic-diagnostics/public-x86_64.txt')
  block <- panel[['stock_rsa3_air']]
  row <- tokens(block[startsWith(block, 'y_ef_24 ')])
  expect_identical(as.integer(stats::start(errors)), as.integer(row[2:3]), info = 'evidence dates')
  evidence <- fixture_hex_doubles(row[-seq_len(4L)])
  expect_identical(length(errors), length(evidence), info = 'evidence count')
  expect_true(all(abs(as.numeric(errors) - evidence) <= 1e-12 * abs(evidence)),
    info = 'y_ef within the core replay tolerance tier of the JVM words')

  # The entry's zero-horizon status passes through verbatim (the estimated
  # route's convention, mirrored): no predictions requested, empty vector.
  zero <- rsa3(forecast = 0L)
  expect_identical(zero$series_forecast_stdev_status, 'no_requested_predictions')
  expect_identical(zero$series_forecast_stdev, numeric())

  # The refusal half of the absence contract: a documented non-available status
  # never carries the pair (the fields are missing, never estimated on the side).
  refusal <- rsa3(calendar = TRUE)
  expect_identical(refusal$status,
    'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)')
  expect_null(refusal$series_forecast_stdev, info = 'absent on refusal')
  expect_null(refusal$series_forecast_stdev_status, info = 'absent on refusal')
})
