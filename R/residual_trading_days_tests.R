# M7 (item 27 cutover gap): the public residual trading-day F-tests. The pure
# primitive native_residual_trading_days_tests() serves the four JVM dashboard
# scalars diagnostics.td-sa-all / td-sa-last / td-i-all / td-i-last from the
# seasonally-adjusted and irregular component series and their shared monthly
# domain start: ARIMA(0,1,1)[12] (nonseasonal) over the six GenericTradingDays
# TD7 contrast columns (no mean), the pinned JointTest F(6, dim-7). The pinned
# capture option is flast = 0, so the "-all" and "-last" names of a component
# share one fit and agree bit for bit. Bit-exactness claim layering: the
# primitive is BIT-exact against the 24-case double-JVM oracle
# (fixtures/native/residual-trading-days/public-x86_64.txt) whenever the
# component inputs are the oracle's own bit-frozen series; product fits supply
# tolerance-qualified component series, so their td scalars are wiring-exact,
# not fixture-bit-exact (see test-native-residual-trading-days.R).
# .attach_td_tests() is the product-route surfacing: the fit routes attach
# diagnostics[["td.ftest.sa"]] / diagnostics[["td.ftest.i"]] (rjd3-side shape
# list(value, pvalue), the resolution target of Persephone's
# metric("diagnostics.td-sa-all", c("diagnostics","td.ftest.sa"))) computed
# from the fit's OWN decorated finalSA / finalIrregular series on their
# published observed domain, on the success path only. The surface is
# diagnostic, never a fit precondition: a route without one of the components
# — or with a component the pinned machine cannot test (at or below the core's
# eight-observation degrees-of-freedom guard, or with no usable residual
# variance, as for the zero-composition pseudo-additive automatic fits) —
# omits the fields; never a fabricated component, domain or number.
.residual_td_whole_number <- function(value, name) {
  if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) || length(value) != 1L ||
      !is.finite(value) || value != floor(value) || abs(value) > .Machine$integer.max)
    stop(name, " must be a whole number in the native range.", call. = FALSE)
  as.integer(value)
}
native_residual_trading_days_tests <- function(sa, irregular, start_year, start_month) {
  finite_vector <- function(value, name) {
    if (!is.numeric(value) || is.complex(value) || !length(value) || any(!is.finite(as.numeric(value))))
      stop(name, " must be a nonempty finite numeric vector.", call. = FALSE)
    as.numeric(value)
  }
  year <- .residual_td_whole_number(start_year, "start_year")
  month <- .residual_td_whole_number(start_month, "start_month")
  if (year < 1L || year > 9999L || month < 1L || month > 12L)
    stop("the residual trading-day test requires a valid monthly domain (year 1-9999, month 1-12).",
         call. = FALSE)
  .native_residual_trading_days_tests(finite_vector(sa, "sa"), finite_vector(irregular, "irregular"),
                                      year, month)
}
# The product-route attach: both components must be present on their shared
# monthly domain start and longer than the core's eight-observation
# positive-degrees-of-freedom guard; otherwise the attach is a no-op and the
# fields stay absent. When the pinned machine itself declines a present,
# long-enough component (the residual-variance and covariance-factorisation
# guards of the core fit — e.g. the automatic route publishes identically-zero
# finalSA/finalIrregular for pseudo-additive requests whose additive
# decomposition collapses), the test is simply not computable and the attach
# likewise omits the fields rather than failing the delivered fit: the JVM
# dashboard collects these scalars as optional diagnostics, and absence is
# reported as absence, never as a fabricated number.
.attach_td_tests <- function(diagnostics, sa, irregular, start) {
  if (is.null(sa) || is.null(irregular) || length(sa) <= 8L || length(irregular) <= 8L)
    return(diagnostics)
  tests <- tryCatch(native_residual_trading_days_tests(sa, irregular, start[1L], start[2L]),
                    error = function(e) NULL)
  if (is.null(tests))
    return(diagnostics)
  diagnostics[["td.ftest.sa"]] <- tests[["td-sa-all"]]
  diagnostics[["td.ftest.i"]] <- tests[["td-i-all"]]
  diagnostics
}
