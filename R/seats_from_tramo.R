#' Decompose a monthly series using a fitted TRAMO model
#'
#' Use supplied TRAMO coefficients and likelihood statistics to estimate
#' SEATS components. Remove fitted outlier effects for the decomposition
#' and restore them to the appropriate components on the original scale.
#'
#' @param series Complete finite univariate monthly \code{ts} covering the
#'   estimation span, in original units.
#' @param ssq Positive fitted residual sum of squares.
#' @param effective,estimated Effective observation count and number of
#'   estimated parameters. The parameter count must be less than the
#'   effective count, which must not exceed the series length.
#' @param d,bd Fitted regular and seasonal differencing orders.
#' @param phi,theta,bphi,btheta Fitted regular AR, regular MA, seasonal AR,
#'   and seasonal MA coefficients. See \code{\link{seats_extract_wk}} for
#'   the sign convention.
#' @param log,mean_correction Whether the fitted model uses logs and whether
#'   its differenced model includes a mean, respectively.
#' @param outliers Optional data frame in fitted model order with
#'   \code{type} (\code{"AO"} for additive outliers, \code{"LS"} for
#'   level shifts, or \code{"TC"} for temporary changes), one-based
#'   \code{position}, and fitted \code{coefficient}. Temporary changes use
#'   decay 0.7, which can also be supplied in \code{tc_rate}.
#' @param forecast,backcast Numbers of months to extend after and before
#'   the data. Negative counts denote years. Zero requests no extension.
#' @details Supply a model estimated on the full input span. The supported
#' deterministic terms are fitted outliers and mean correction. Calendar
#' effects, user regressors, and length-of-period adjustments are outside
#' the supplied-model specification.
#' @return A list with \code{status} and the dated \code{initial},
#'   \code{final}, and \code{merged} component tables. \code{initial} and
#'   \code{final} include conditional standard deviations. \code{merged}
#'   contains values, forecasts, and backcasts after restoring the fitted
#'   outlier effects.
#' @seealso \code{\link{seats_extract_wk}}, \code{\link{native_tramoseats_fit}}
seats_from_tramo <- function(series, ssq, effective, estimated, d = 0L, bd = 0L,
    phi = numeric(), theta = numeric(), bphi = numeric(), btheta = numeric(),
    log = FALSE, mean_correction = TRUE, outliers = NULL, forecast = -2L, backcast = 0L) {
  spec <- as.list(environment())[setdiff(names(formals(seats_from_tramo)), c("series", "outliers"))]
  domain <- seats_model_domain(series, spec)
  spec$effective <- native_horizon(effective, "effective")
  spec$estimated <- native_horizon(estimated, "estimated")
  if (!is.numeric(ssq) || is.complex(ssq) || !is.null(dim(ssq)) ||
      length(ssq) != 1L || !is.finite(ssq) || ssq <= 0 ||
      spec$effective > length(series) || spec$effective <= spec$estimated ||
      !is.null(dim(effective)) || !is.null(dim(estimated)))
    stop("The fitted likelihood requires positive ssq and 0 <= estimated < effective <= length(series).",
         call. = FALSE)
  if (is.null(outliers)) outliers <- data.frame(type = character(), position = integer(), coefficient = numeric())
  if (!is.data.frame(outliers) || !setequal(setdiff(names(outliers), "tc_rate"), c("type", "position", "coefficient")) ||
      anyDuplicated(names(outliers)) > 0L || !is.character(outliers$type) || !is.null(dim(outliers$type)) ||
      any(!outliers$type %in% c("AO", "LS", "TC")) ||
      !is.numeric(outliers$position) || is.complex(outliers$position) || !is.null(dim(outliers$position)) ||
      any(!is.finite(outliers$position)) || any(outliers$position != floor(outliers$position)) ||
      any(outliers$position < 1L | outliers$position > length(series)) ||
      !is.numeric(outliers$coefficient) || is.complex(outliers$coefficient) ||
      !is.null(dim(outliers$coefficient)) || any(!is.finite(outliers$coefficient)))
    stop("outliers must have AO, LS or TC type, one-based observed position and finite coefficient.", call. = FALSE)
  if ("tc_rate" %in% names(outliers)) {
    rates <- outliers$tc_rate[outliers$type == "TC"]
    if (!is.numeric(outliers$tc_rate) || is.complex(outliers$tc_rate) ||
        !is.null(dim(outliers$tc_rate)) || anyNA(rates) || any(rates != 0.7))
      stop("Monthly TC outliers require the stock tc_rate of 0.7.", call. = FALSE)
  }
  outliers$position <- as.integer(outliers$position)
  result <- .native_seats_from_tramo(spec, domain, as.numeric(series), outliers)
  seats_date_components(result, domain$frequency, c("initial", "final", "merged"))
}
