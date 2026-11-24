native_x13_estimated_fit <- function(series, template = "RSA0") {
  templates <- c(RSA0 = "RSA0", RSA1 = "RSA1", RSA2C = "RSA2c")
  if (!is.character(template) || length(template) != 1L || !is.null(dim(template)) ||
      is.na(template) || !toupper(template) %in% names(templates))
    stop("The supported stock estimated X-13 templates are RSA0, RSA1 and RSA2c.", call. = FALSE)
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)))
    stop("Stock estimated X-13 fits require complete finite monthly observations.", call. = FALSE)
  if (domain$start[1L] < 1L || domain$start[1L] > 9999L)
    stop("Stock estimated X-13 fits require a start year between 1 and 9999.", call. = FALSE)
  fit <- .native_x13_estimated_fit(as.numeric(series), domain$start[1L], domain$start[2L], templates[[toupper(template)]])
  fit$sigma_provenance <- "stock_template"
  fit <- date_estimated_x13_result(fit, series, domain$start, domain$frequency)
  if (length(fit$residuals$tsres)) {
    first <- domain$start_ordinal + fit$model$d + domain$frequency * fit$model$bd
    fit$residuals$tsres <- stats::ts(fit$residuals$tsres,
      start = c(first %/% domain$frequency, first %% domain$frequency + 1), frequency = domain$frequency)
  }
  structure(fit, class = c("native_x13_estimated_fit", "list"))
}
