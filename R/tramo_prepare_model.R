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
