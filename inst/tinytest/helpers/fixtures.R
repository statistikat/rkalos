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
