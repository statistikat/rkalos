#!/usr/bin/env Rscript
# Reproduction driver for tools/svd-tolerance-evidence.md. Build the probe first:
#   PKG_LIBS="$(R CMD config LAPACK_LIBS) $(R CMD config BLAS_LIBS) $(R CMD config FLIBS)" \
#     R CMD SHLIB -o svd-tolerance-probe.so tools/svd-tolerance-probe.cpp
# then run with rkalos installed:
#   Rscript --vanilla tools/svd-tolerance-probe.R [path/to/svd-tolerance-probe.so]
# It changes nothing; it only measures.
so <- commandArgs(trailingOnly = TRUE)[[1]]
if (is.na(so) || !nzchar(so)) so <- "svd-tolerance-probe.so"
if (!file.exists(so)) stop("probe shared object not found: ", so,
                           " (build it with R CMD SHLIB first)", call. = FALSE)
dyn.load(normalizePath(so))
dgesvd <- function(m) as.numeric(.Call("probe_dgesvd", as.matrix(m), PACKAGE = sub("[.].*$", "", basename(so))))
dgesdd <- function(m) as.numeric(.Call("probe_dgesdd", as.matrix(m), PACKAGE = sub("[.].*$", "", basename(so))))
stopifnot(requireNamespace("rkalos", quietly = TRUE))

parameter_map <- getFromNamespace("native_arma_parameter_map", "rkalos")
identify <- getFromNamespace("native_arma_identification", "rkalos")
ulp <- function(x) ifelse(x == 0, .Machine$double.xmin, abs(x) * .Machine$double.eps)
norms <- function(m) sqrt(colSums(m^2))
jac_fixture <- function(a, difference) rbind(c(0, 0), c(1, 1), c(-a + difference, -a))

cat("provider:", format(rkalos::core_version()$LAPACK$library), "LAPACK",
    format(rkalos::core_version()$LAPACK$version), "| BLAS:", format(rkalos::core_version()$BLAS), "\n")

cat("\n== failing case: ar = -0.6, ma = 0.6 + 1e-6 ==\n")
for (period in c(4L, 12L)) {
  p <- parameter_map(ar = NA_real_, ma = NA_real_, initial = list(ar = -0.6, ma = 0.6 + 1e-6))
  r <- identify(p, period)
  J <- unname(r$numerator_jacobian)
  expected <- jac_fixture(-0.6, 1e-6)
  A <- sweep(J, 2, r$column_norms, "/")          # native matrix, native norms
  B <- sweep(expected, 2, norms(expected), "/")  # fixture matrix, R norms
  target <- svd(B, nu = 0, nv = 0)$d
  cat(sprintf("period %d: kappa2 %.4e sigma_max %.17g sigma_min %.17g\n",
              period, max(target) / min(target), max(target), min(target)))
  cat(sprintf("  native norms bitwise equal to sqrt(colSums(expected^2)): %s\n",
              identical(as.numeric(r$column_norms), as.numeric(norms(expected)))))
  cat(sprintf("  entrywise |native - fixture| in ULP: %s\n",
              paste(sprintf("%.1f", abs(as.vector(J) - as.vector(expected)) / ulp(as.vector(expected))), collapse = " ")))
  cat(sprintf("  native == DGESVD(native matrix): %s | R svd == DGESVD(fixture matrix): %s | DGESVD == DGESDD: %s\n",
              identical(dgesvd(A), as.numeric(r$singular_values)),
              identical(target, dgesvd(B)), identical(dgesvd(B), dgesdd(B))))
  cat(sprintf("  |A - B|max %.3e | ||A - B||_2 (Weyl bound on |d sigma|) %.3e | observed |d sigma_min| %.3e -> within bound: %s\n",
              max(abs(as.vector(A) - as.vector(B))), max(svd(A - B, nu = 0, nv = 0)$d),
              abs(min(dgesvd(A)) - min(dgesvd(B))),
              abs(min(dgesvd(A)) - min(dgesvd(B))) <= max(svd(A - B, nu = 0, nv = 0)$d)))
  bumps <- vapply(which(as.vector(B) != 0), function(i) {
    Bp <- B; Bp[i] <- Bp[i] + ulp(Bp[i]); min(dgesvd(Bp))
  }, numeric(1))
  cat(sprintf("  sigma_min shift from +1 ULP on a single entry: %s\n",
              paste(sprintf("%.2e", abs(bumps - min(dgesvd(B)))), collapse = " ")))
  cat(sprintf("  all.equal at 1e-14: %s | sigma-scaled bound 64*eps*sigma_max = %.3e vs gap %.3e\n",
              paste(all.equal(unname(r$singular_values), target, tolerance = 1e-14), collapse = " "),
              64 * .Machine$double.eps * max(target), max(abs(as.vector(r$singular_values) - target))))
}

cat("\n== survey of every first-loop case at tolerance 1e-14 ==\n")
for (period in c(4L, 12L)) for (a in c(-.6, 0, .4)) for (difference in c(0, 1e-6, .2)) {
  p <- parameter_map(ar = NA_real_, ma = NA_real_, initial = list(ar = a, ma = -a + difference))
  r <- identify(p, period)
  B <- sweep(jac_fixture(a, difference), 2, norms(jac_fixture(a, difference)), "/")
  target <- svd(B, nu = 0, nv = 0)$d
  gap <- max(abs(as.vector(r$singular_values) - target))
  cat(sprintf("period=%2d a=%5s diff=%-6s kappa2=%10.3e abs gap=%9.2e 1e-14 ok=%-5s 64*eps*sigma_max bound ok=%s\n",
              period, format(a), format(difference),
              if (min(target) == 0) Inf else max(target) / min(target), gap,
              isTRUE(all.equal(unname(r$singular_values), target, tolerance = 1e-14)),
              gap < 64 * .Machine$double.eps * max(target)))
}

cat("\n== structured seasonal block (well conditioned, passes) ==\n")
multiply <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) out[i + seq_along(b) - 1L] <- out[i + seq_along(b) - 1L] + a[i] * b
  out
}
source(system.file("tinytest/helpers/fixtures.R", package = "rkalos"), local = TRUE)
source_fixture(system.file("tinytest/fixtures/native/arma-fit-cases.R", package = "rkalos"), local = TRUE)
for (period in c(4L, 12L)) {
  p <- parameter_map(ar = c(NA_real_, -.2), ma = NA_real_, seasonal_ma = NA_real_,
    initial = list(ar = c(.4, -.2), ma = .25, seasonal_ma = -.3))
  r <- identify(p, period)
  base <- arma_fit_polynomials(p$initial, period)
  expected <- matrix(0, length(base$ar) + length(base$ma) - 1L, 3L)
  for (j in seq_along(c("ar", "ma", "seasonal_ma"))) {
    key <- c("ar", "ma", "seasonal_ma")[j]
    perturbed <- p$initial; perturbed[[key]][1] <- perturbed[[key]][1] + .125
    next_poly <- arma_fit_polynomials(perturbed, period)
    expected[, j] <- (multiply(base$ar, next_poly$ma - base$ma) -
      multiply(base$ma, next_poly$ar - base$ar)) / .125
  }
  B <- sweep(expected, 2, norms(expected), "/")
  target <- svd(B, nu = 0, nv = 0)$d
  cat(sprintf("period=%2d dimensions %2dx%d kappa2=%.3e abs gap=%.2e 1e-14 ok=%s\n", period, nrow(B), ncol(B),
              max(target) / min(target), max(abs(as.vector(r$singular_values) - target)),
              isTRUE(all.equal(unname(r$singular_values), target, tolerance = 1e-14))))
}
