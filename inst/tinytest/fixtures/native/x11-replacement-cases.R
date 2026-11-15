# Case inventory and deterministic inputs for the same-phase replacement,
# endpoint and correction-factor fixture (item 22b.3b). No RNG is used
# anywhere; every input is a closed-form function of the case identity. The
# sigma weights feeding the corrector are produced by the pinned analyse pass;
# the fixture is only consulted by validation/reference/capture-x11-replacement.R.
x11_replacement_cases <- function() {
  cases <- list()
  add <- function(category,frequency,n,start_period,backcasts=0L,forecasts=0L,
                  exclude=FALSE,lsig=1.5,usig=2.5,multiplicative=FALSE) {
    id <- paste(category,frequency,n,start_period,backcasts,forecasts,
      if (exclude) 'ex' else 'in', if (multiplicative) 'mul' else 'add',
      format(lsig, trim = TRUE, scientific = FALSE),
      format(usig, trim = TRUE, scientific = FALSE), sep = '-')
    cases[[id]] <<- list(id=id,category=category,frequency=as.integer(frequency),
      n=as.integer(n),start_period=as.integer(start_period),backcasts=as.integer(backcasts),
      forecasts=as.integer(forecasts),exclude=exclude,lsig=lsig,usig=usig,
      multiplicative=multiplicative)
  }
  for (f in c(12L,4L)) {
    # Flagged rows at the series start/end: the neighbor walk must extend the
    # opposite side because fewer than two same-phase rows exist on one side.
    add('start',f,7L*f+1L,1L)
    add('end',f,7L*f+2L,2L)
    # Flagged same-phase rows in between: fully weighted neighbors are further
    # away and the walk must skip every non-unit weight.
    add('skip',f,9L*f,1L)
    # One calendar phase entirely flagged: no fully weighted same-phase
    # neighbor exists anywhere, so every flagged row takes the period-average
    # fallback, which equals the original value for this shape.
    add('phaseflag',f,8L*f,1L)
    # Nearly everything in the ramp band: only a handful of unit weights.
    add('allramp',f,6L*f,1L)
    # Fully masked same-phase rows act as unit-weight neighbors and poison
    # the blended replacement with NaN, so applyCorrections keeps the original.
    add('maskedphase',f,7L*f,1L)
    # Flagged rows inside the projections: captured with and without the
    # pinned projection masking of the replacement pass.
    add('projection',f,7L*f,1L,backcasts=2L,forecasts=f+1L,exclude=TRUE)
    add('projection',f,7L*f,2L,backcasts=f,forecasts=f,exclude=FALSE)
    # Legacy revival leaves every final weight at one: no replacements at all,
    # corrections cancel and the factors stay at the mode default.
    add('revive',f,7L*f,1L)
    # Multiplicative mode exercises the division-form correction factors.
    if (f == 12L) {
      add('mul',f,8L*f,1L,multiplicative=TRUE)
      add('mulramp',f,7L*f,3L,multiplicative=TRUE)
    }
    # Wider sigma band: more ramp rows, deeper neighbor walks.
    add('wide',f,8L*f+2L,3L,lsig=0.6,usig=3.0)
  }
  cases
}
# Values are deterministic closed forms. Categories control the shape:
# start/end     - smooth with two huge rows at the first/last positions;
# skip          - smooth with huge rows at rows 3, 3+f, 3+2f and a medium
#                 spike at 3+3f whose same-phase neighborhood is damaged;
# phaseflag     - zero everywhere except one calendar phase set to 40;
# allramp       - a quarter of the rows (whole calendar phases) at 40, rest zero;
# maskedphase   - smooth with a huge row at f+2 and masked rows at f+2+f(+2f);
# projection    - smooth with huge rows inside the backcast and forecast zones;
# revive        - exact zeros with isolated spikes (legacy revival no-op);
# mul/mulramp   - multiplicative around one: spike variant and a flagged-phase ramp set;
# wide          - smooth with a graduated set of large deviations.
x11_replacement_input <- function(case) {
  n <- case$n; f <- case$frequency
  t <- seq_len(n)
  smooth <- 3 * sin(t * .7) + 2 * cos(t * sqrt(2)) + ((t * 37L) %% 11L) / 3.0 - 1.7
  values <- switch(case$category,
    start = { v <- smooth; v[1L] <- 8000; v[2L] <- -6000; v },
    end   = { v <- smooth; v[n] <- 8000; v[n - 2L] <- -6000; v },
    skip  = { v <- smooth; v[3L] <- 9000; v[3L + f] <- -9000; v[3L + 2L * f] <- 9000
              v[3L + 3L * f] <- 90; v },
    phaseflag = { v <- numeric(n); v[((0L * f + 1L):(n)) %% f == 0L] <- 40; v },
    allramp   = { v <- numeric(n); v[seq_len(n) %% f %in% 1L:(f/4L)] <- 40; v },
    maskedphase = { v <- smooth; v[f + 2L] <- 9e5; v },
    projection = { v <- smooth; v[1L] <- 5e3; v[n - 1L] <- -5e3; v },
    revive  = { v <- numeric(n); v[f + 3L] <- 9; v[3L * f + 7L] <- -4; v },
    mul     = { v <- 1 + .35 * sin(t * .9) + .2 * cos(t * sqrt(3)); v[5L] <- 2.6
                v[5L + 3L * f] <- .3; v },
    mulramp = { v <- rep(1, n); v[seq_len(n) %% f %in% 1L:(f/4L)] <- 1.5; v },
    wide    = { v <- smooth; v[f + 1L] <- 60; v[2L * f + 2L] <- -45; v[3L * f + 3L] <- 30
                v[4L * f + 4L] <- -20; v[5L * f + 5L] <- 12; v },
    stop('unknown replacement case category'))
  missing <- rep(FALSE, n)
  if (case$category == 'maskedphase') missing[c(f + 2L + f, f + 2L + 2L * f)] <- TRUE
  list(values = as.numeric(values), missing = missing)
}
