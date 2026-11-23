native_tramo_preprocess <- function(series, template = "TR1") {
  templates <- c(TR1 = "TR1", RSA1 = "TR1", TR3 = "TR3", RSA3 = "TR3")
  if (!is.character(template) || length(template) != 1L || !is.null(dim(template)) ||
      is.na(template) || !toupper(template) %in% names(templates))
    stop("The native preprocessing templates are TR1 (RSA1) and TR3 (RSA3).", call. = FALSE)
  template <- unname(templates[[toupper(template)]])
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)) || length(series) < 24L)
    stop(template, " requires at least 24 complete finite monthly observations.", call. = FALSE)
  if (abs(domain$start[1L]) > .Machine$integer.max)
    stop("The ", template, " start year exceeds the supported native range.", call. = FALSE)
  fit <- if (template == "TR3") .native_tramo_tr3 else .native_tramo_tr1
  result <- fit(as.numeric(series), domain$start[1L], domain$start[2L])
  result$template <- template
  result <- format_tramo_preprocessing(result, series, domain)
  class(result) <- c("native_tramo_preprocess", "list")
  result
}

format_tramo_preprocessing <- function(result, series, domain) {
  orders <- result$model$orders[c("p", "bp", "q", "bq")]
  blocks <- rep(c("phi", "bphi", "theta", "btheta"), orders)
  for (name in c("phi", "bphi", "theta", "btheta"))
    result$model[[name]] <- result$model$parameters[blocks == name]
  first <- domain$start_ordinal + length(series) - result$likelihood$effective
  result$residuals$tsres <- stats::ts(result$residuals$tsres, frequency = 12,
    start = c(first %/% 12, first %% 12 + 1))
  result
}
