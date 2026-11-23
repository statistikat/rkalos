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

# Recreate these two historical generators with their capture environment's
# trigonometry and binary64 literals. The staged kernel checks still read the
# independent JVM CSVs, and the generated inputs must match them bit for bit.
source_x11_capture_generator <- function(file, destination = parent.frame()) {
  env <- new.env(parent = destination)
  trig <- getFromNamespace(".capture_trig", "rkalos")
  env$sin <- function(x) trig(x, TRUE)
  env$cos <- function(x) trig(x, FALSE)
  code <- readLines(file)
  for (literal in c("1e300", "1e298")) {
    hex <- sprintf("%a", fixture_doubles(literal))
    code <- gsub(literal, hex, code, fixed = TRUE)
  }
  eval(parse(text = code), envir = env)
  names <- ls(env, pattern = "^x11_(sigma|replacement)_(input|cases)$")
  for (name in names) assign(name, env[[name]], envir = destination)
  invisible(NULL)
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

fixture_hex_doubles <- function(tokens) {
  b <- vapply(seq(1L, 15L, 2L), function(i) strtoi(substr(tokens, i, i + 1L), 16L),
              integer(length(tokens)))
  stream <- if (length(tokens) == 1L) b else t(b)
  readBin(as.raw(stream), double(), n = length(tokens), endian = "big")
}
