native_tramo_preprocess <- function(series, template = "TR1") {
  if (!is.character(template) || length(template) != 1L || !is.null(dim(template)) ||
      is.na(template) || !toupper(template) %in% c("TR1", "RSA1"))
    stop("The native preprocessing template is TR1 (RSA1).", call. = FALSE)
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)) || length(series) < 24L)
    stop("TR1 requires at least 24 complete finite monthly observations.", call. = FALSE)
  if (abs(domain$start[1L]) > .Machine$integer.max)
    stop("The TR1 start year exceeds the supported native range.", call. = FALSE)
  result <- .native_tramo_tr1(as.numeric(series), domain$start[1L], domain$start[2L])
  result$template <- "TR1"
  orders <- result$model$orders[c("p", "bp", "q", "bq")]
  blocks <- rep(c("phi", "bphi", "theta", "btheta"), orders)
  for (name in c("phi", "bphi", "theta", "btheta"))
    result$model[[name]] <- result$model$parameters[blocks == name]
  first <- domain$start_ordinal + length(series) - result$likelihood$effective
  result$residuals <- list(tsres = stats::ts(result$residuals, frequency = 12,
    start = c(first %/% 12, first %% 12 + 1)))
  class(result) <- c("native_tramo_preprocess", "list")
  result
}
