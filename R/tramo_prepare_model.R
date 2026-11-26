#' Prepare an initial TRAMO model specification
#'
#' Initialize the model orders, transformation and regression settings for
#' a standard TRAMO specification, and fill gaps before model fitting.
#'
#' @param series Univariate monthly numeric \code{ts}. At least one
#'   observation must be finite. Missing observations are filled by
#'   interpolation.
#' @param template One of \code{"TR0"}, \code{"TR1"}, \code{"TR2"},
#'   \code{"TR3"}, \code{"TR4"}, \code{"TR5"}, or \code{"TRfull"}, matched
#'   without regard to case.
#' @return A list with the template, initial SARIMA \code{orders},
#'   frequency, mean and transformation flags, whether automatic model
#'   selection is requested, and regression variable names. \code{missing}
#'   holds one-based indices into the original series and \code{interpolated}
#'   is the completed series on its original calendar.
#' @seealso \code{\link{native_tramo_preprocess}}
#' @examples
#' x <- AirPassengers
#' x[20] <- NA_real_
#' model <- tramo_prepare_model(x, "TR1")
#' model$missing
#' model$orders
tramo_prepare_model <- function(series, template) {
  domain <- native_regular_domain(series)
  if (abs(domain$start[1L]) > .Machine$integer.max || length(series) > .Machine$integer.max)
    stop("The TRAMO domain exceeds the native index range.", call. = FALSE)
  if (!is.character(template) || !is.null(dim(template)) || length(template) != 1L || is.na(template))
    stop("template must be a single TRAMO template name.", call. = FALSE)
  result <- .native_tramo_prepare_model(as.numeric(series), template, domain$start[1L],
                                        domain$start[2L], domain$frequency)
  attributes(result$interpolated) <- list(tsp = stats::tsp(series), class = "ts")
  result
}
