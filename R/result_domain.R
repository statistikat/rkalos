# R-only result formatting, migrated from Persephone 8450d6e R/result_schema.R.
native_result_domain <- function(series) {
  if (is.null(series)) return(NULL)
  if (!stats::is.ts(series)) stop("A normalized series must have a time domain.", call. = FALSE)
  list(start = stats::start(series), end = stats::end(series),
       frequency = stats::frequency(series), length = length(series),
       missing = is.na(series))
}

