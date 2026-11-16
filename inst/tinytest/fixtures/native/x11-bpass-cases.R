# Case inventory and deterministic inputs for the fixed additive B-pass
# fixture (item 22b.4 S1). No RNG is used anywhere; every input is a
# closed-form function of the case identity. The fixture is only consulted
# by validation/reference/capture-x11-bpass.R.
x11_bpass_cases <- function() {
  cases <- list()
  add <- function(category, frequency, n, start_period, initial, final,
                  trend = 13L, lsig = 1.5, usig = 2.5) {
    id <- paste(category, frequency, n, start_period, trend,
      format(lsig, trim = TRUE, scientific = FALSE),
      format(usig, trim = TRUE, scientific = FALSE), sep = '-')
    cases[[id]] <<- list(id = id, category = category, frequency = as.integer(frequency),
      n = as.integer(n), start_period = as.integer(start_period),
      initial = as.integer(initial), final = as.integer(final), trend = as.integer(trend),
      lsig = lsig, usig = usig)
  }
  for (f in c(12L, 4L)) {
    # Baseline: smooth series with one extreme spike and one masked row;
    # uniform fixed seasonal filters on both stages.
    add('basic', f, 9L * f, 1L, rep(5L, f), rep(9L, f))
    # Rotating first-period coverage against the same shape.
    add('rotation', f, 8L * f, if (f == 12L) 7L else 3L, rep(5L, f), rep(9L, f))
    # Stable phases on both stages: the stable phases enter the
    # normalization restoration set and the processor keeps them identity.
    add('stable', f, 9L * f, 2L, { v <- rep(5L, f); v[c(1L, 3L)] <- 0L; v },
        { v <- rep(9L, f); v[seq_len(min(f, 4L))] <- 0L; v })
    # Per-phase mixed requests on both stages including 3x15 (and, for
    # monthly, 3x1) so both normalizer span rules are exercised.
    if (f == 12L) {
      add('mixed', f, 11L * f, 4L,
          c(1L, 3L, 5L, 9L, 15L, 5L, 3L, 1L, 9L, 15L, 5L, 0L),
          c(15L, 9L, 5L, 3L, 1L, 9L, 15L, 5L, 3L, 1L, 0L, 9L))
      # Wider sigma band and a deeper spike: all three corrector rounds flag.
      add('spike2', f, 10L * f, 1L, rep(3L, f), rep(5L, f), lsig = 1, usig = 3)
    } else {
      add('mixed', f, 12L * f, 2L, c(3L, 15L, 5L, 0L), c(15L, 3L, 1L, 5L))
      add('spike2', f, 11L * f, 1L, rep(3L, f), rep(5L, f), lsig = 1, usig = 3)
    }
    # Long Henderson: wider asymmetric band and deeper endpoint musgrave.
    add('h23', f, 10L * f, 3L, rep(5L, f), rep(9L, f), trend = 23L)
    # Zero-dispersion revival: constant series keeps every final weight at
    # one, so all three correction rounds cancel and B17/B20 stay default.
    add('revive', f, 8L * f, 1L, rep(5L, f), rep(5L, f))
    # Short center span: normalization short_span rules plus the tight
    # seasonal-processor gates.
    add('short3', f, if (f == 12L) 4L * f + 2L else 6L * f, 1L,
        rep(5L, f), rep(9L, f))
    # Masked internal rows feeding the seasonal processors and detection.
    add('masked', f, 9L * f + 1L, if (f == 12L) 5L else 2L, rep(5L, f), rep(5L, f), lsig = 2, usig = 3)
  }
  cases
}
# Category shapes (all deterministic):
# basic/rotation/stable/mixed/spike2/h23 - smooth baseline with one huge row
#   near the start, one deep row near the end and (for spike2) an extra
#   mid-series spike; revive is exactly constant; short3 is the smooth shape;
# masked is smooth with a masked interior same-phase triplet.
x11_bpass_input <- function(case) {
  n <- case$n; f <- case$frequency
  t <- seq_len(n)
  values <- if (case$category == 'revive') {
    rep(37.5, n)
  } else {
    smooth <- 200 + t * .8 + 12 * sin(t * 2 * pi / f) +
      5 * cos(t * sqrt(2)) + ((t * 29L) %% 7L) / 2.5
    if (case$category == 'spike2') {
      smooth[3L] <- smooth[3L] + 4000; smooth[n - 2L] <- smooth[n - 2L] - 3500
      smooth[(n %/% 2L) + 1L] <- smooth[(n %/% 2L) + 1L] + 1800; smooth
    } else {
      smooth[4L] <- smooth[4L] + 3000; smooth[n - 1L] <- smooth[n - 1L] - 2600
      smooth
    }
  }
  missing <- rep(FALSE, n)
  if (case$category == 'masked') {
    rows <- c(f + 2L, f + 2L + f, f + 2L + 2L * f)
    missing[rows[rows <= n]] <- TRUE
  }
  list(values = as.numeric(values), missing = missing)
}
