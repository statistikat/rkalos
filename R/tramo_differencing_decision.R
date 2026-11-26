#' Choose regular and seasonal differencing for TRAMO
#'
#' Assess how much regular and seasonal differencing is needed and whether
#' to include a mean term in the differenced model.
#'
#' @param series Complete finite monthly or quarterly \code{ts}. Fill gaps
#'   and remove any outlier effects before calling. Supply the original
#'   scale when requesting a log transformation.
#' @param orders Named vector of initial SARIMA orders, as returned by
#'   \code{\link{tramo_prepare_model}}. Regular orders \code{p}, \code{d},
#'   and \code{q} range from zero through two. Seasonal orders \code{bp},
#'   \code{bd}, and \code{bq} range from zero through one.
#' @param mean Whether the initial model includes a mean.
#' @param log Whether to take the natural logarithm. Requires positive data.
#' @param seasonal Whether to consider seasonal differencing.
#' @param initial_round Start the search at zero differencing orders
#'   (\code{TRUE}), or at the supplied orders (\code{FALSE}).
#' @return A list with the selected orders and mean setting, the mean test
#'   and statistic, flags indicating a default or airline model or failed
#'   selection, and the sequence of decisions. Input or estimation errors
#'   can prevent a decision.
tramo_differencing_decision <- function(series, orders = c(p = 0L, d = 0L, q = 0L,
                                        bp = 0L, bd = 0L, bq = 0L),
                                        mean = TRUE, log = FALSE, seasonal = TRUE,
                                        initial_round = TRUE) {
  domain <- native_regular_domain(series)
  if (length(series) > .Machine$integer.max || any(!is.finite(series)))
    stop("series must contain complete finite values within the native index range.", call. = FALSE)
  keys <- c("p", "d", "q", "bp", "bd", "bq")
  if (!is.numeric(orders) || is.complex(orders) || !is.null(dim(orders)) ||
      length(orders) != length(keys) || !setequal(names(orders), keys) ||
      any(!is.finite(orders)) || any(orders != floor(orders)) ||
      any(orders < 0) || any(orders[c("p", "d", "q")] > 2) ||
      any(orders[c("bp", "bd", "bq")] > 1))
    stop("orders must name p, d, q in 0..2 and bp, bd, bq in 0..1.", call. = FALSE)
  .native_tramo_differencing(as.numeric(series), domain$frequency, as.integer(orders[keys]),
    native_regression_flag(mean, "mean"), native_regression_flag(log, "log"),
    native_regression_flag(seasonal, "seasonal"), native_regression_flag(initial_round, "initial_round"))
}
