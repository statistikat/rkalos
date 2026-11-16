# Case inventory and deterministic inputs for the fixed multiplicative B/C/D
# pass fixtures (item 22b.5 M3a). No RNG is used anywhere; every input is a
# closed-form strictly positive function of the case identity, mirroring the
# additive x11-bpass-cases.R inventory (same categories, spans, start
# periods, seasonal requests, Henderson lengths and sigma bounds) so the two
# fixtures differ only in the decomposition mode. The shapes are
# multiplicative by construction: a level with exponential drift times a
# seasonal factor, with multiplicative spikes. The deep/spike2 categories
# place extreme positive outliers inside the endpoint regions so the ratio
# tables swing by orders of magnitude across the asymmetric Henderson
# endpoints — the regime where the pinned makePositivity repair of B7/C7/D7/
# D12 guards the trend tables. The fixture is only consulted by
# validation/reference/capture-x11-bpass.R --mode=Multiplicative (and the C/
# D/selection captures through the same case names).
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
    # Baseline: smooth multiplicative shape with one huge and one deep row;
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
      # Wider sigma band and deeper outliers: all three corrector rounds flag.
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
    # Near-zero deep rows in both endpoint regions: the ratio tables swing
    # by six orders of magnitude across the asymmetric endpoints, which is
    # where the pinned trend tables can undershoot zero and the
    # makePositivity repair guards B7/C7/D7/D12.
    add('deep', f, 9L * f, if (f == 12L) 3L else 2L, rep(5L, f), rep(9L, f), lsig = 1, usig = 2)
    # Isolated huge endpoint spikes: an x60 row inside the asymmetric
    # endpoint band places an extreme value on the negative taps of the
    # mmsre/Musgrave endpoint weights, which drives the un-repaired trend
    # table negative — the regime that fires the pinned makePositivity
    # repair of B7/C7/D7/D12.
    add('spikeend', f, 9L * f, 1L, rep(5L, f), rep(9L, f), lsig = 1, usig = 3)
    # The same endpoint spike with an unusably wide sigma band (z-scores of
    # the ratio tables never reach nine): nothing is corrected, so the raw
    # spike reaches C1/D1 and the C7/D7/D12 endpoint repairs as well.
    add('wideband', f, 9L * f, 1L, rep(5L, f), rep(9L, f), lsig = 9, usig = 11)
  }
  cases
}
# Category shapes (all deterministic, all strictly positive):
# basic/rotation/stable/mixed/h23/short3 - smooth multiplicative baseline
#   with one huge (x3) row near the start and one deep (x1/3) row near the
#   end; spike2 multiplies by 5 and 3 and by 1/5 instead; deep multiplies
#   rows 4 and n-2 by 1e-6; revive is exactly constant; masked is smooth
#   with a masked interior same-phase triplet.
x11_bpass_input <- function(case) {
  n <- case$n; f <- case$frequency
  t <- seq_len(n)
  if (case$category == 'revive') {
    values <- rep(37.5, n)
    missing <- rep(FALSE, n)
    return(list(values = values, missing = missing))
  }
  level <- 200 * exp(.004 * t) *
    (1 + .06 * sin(t * 2 * pi / f) + .02 * cos(t * sqrt(2)) + ((t * 29L) %% 7L) / 40)
  values <- level
  if (case$category == 'basic' || case$category == 'rotation' ||
      case$category == 'stable' || case$category == 'mixed' ||
      case$category == 'h23' || case$category == 'short3' ||
      case$category == 'masked') {
    values[4L] <- values[4L] * 3
    values[n - 1L] <- values[n - 1L] / 3
  } else if (case$category == 'spike2') {
    values[3L] <- values[3L] * 5
    values[n - 2L] <- values[n - 2L] / 5
    values[(n %/% 2L) + 1L] <- values[(n %/% 2L) + 1L] * 3
  } else if (case$category == 'deep') {
    values[4L] <- values[4L] * 1e-6
    values[n - 2L] <- values[n - 2L] * 1e-6
  } else if (case$category == 'spikeend' || case$category == 'wideband') {
    values[n - 1L] <- values[n - 1L] * 60
  }
  missing <- rep(FALSE, n)
  if (case$category == 'masked') {
    rows <- c(f + 2L, f + 2L + f, f + 2L + 2L * f) # same phase, three cycles
    missing[rows[rows <= n]] <- TRUE
  }
  stopifnot(all(values[!missing] > 0))
  list(values = values, missing = missing)
}
