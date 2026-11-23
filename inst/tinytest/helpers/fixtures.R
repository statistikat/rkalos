# Load frozen source bytes unchanged, adapting only the historical namespace
# lookup in one helper's lexical environment. No original fixture is rewritten.
source_fixture <- function(file, local = FALSE, ...) {
  destination <- if (isTRUE(local)) parent.frame() else if (is.environment(local)) local else .GlobalEnv
  result <- base::source(file, local = destination, ...)
  if (identical(basename(file), "arma-fit-cases.R")) {
    helper <- get("arma_fit_alternate", envir = destination, inherits = FALSE)
    environment(helper) <- list2env(list(getFromNamespace = function(x, ns) {
      utils::getFromNamespace(x, if (identical(ns, "persephone3")) "rkalos" else ns)
    }), parent = environment(helper))
    assign("arma_fit_alternate", helper, envir = destination)
  }
  invisible(result)
}

# The frozen decimal strings encode binary64 values. On platforms where long
# double == double, R's decimal conversion can round these one ULP away.
fixture_doubles <- function(values) {
  result <- rep(NA_real_, length(values))
  result[values %in% c("Inf", "+Inf")] <- Inf
  result[values %in% "-Inf"] <- -Inf
  finite <- !is.na(values) & !values %in% c("NA", "NaN", "", "Inf", "+Inf", "-Inf")
  if (any(finite)) {
    json <- paste0("[", paste(values[finite], collapse=","), "]")
    result[finite] <- jsonlite::fromJSON(json)
  }
  result
}

read_numeric_fixture <- function(path, columns = "value") {
  data <- utils::read.csv(path, colClasses=stats::setNames(rep("character", length(columns)), columns))
  for (column in columns) data[[column]] <- fixture_doubles(data[[column]])
  data
}
