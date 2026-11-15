# Case inventory and deterministic inputs for the moving five-year dispersion
# and two-pass sigma weights fixture (item 22b.3a). No RNG is used anywhere;
# every input is a closed-form function of the case identity. The pinned
# reference is only consulted by validation/reference/capture-x11-sigma.R.
x11_sigma_cases <- function() {
  cases <- list()
  add <- function(category,frequency,n,start_period,backcasts=0L,forecasts=0L,exclude=FALSE,
                  lsig=1.5,usig=2.5,multiplicative=FALSE) {
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
    # The pinned reference crashes (ArrayIndexOutOfBounds) when the first
    # calendar-year boundary exceeds the series length (n < frequency); the
    # fixture matrix therefore starts at one full period, matching the
    # qualified native input contract.
    for (n in unique(c(f, f + 1L, 2L*f, 4L*f, 5L*f - 1L, 5L*f, 5L*f + 1L,
                       6L*f, 7L*f, 8L*f + 2L, 10L*f + 5L, 19L*f)))
      for (sp in unique(c(1L, 2L, f)))
        add('smooth',f,n,sp)
    # Window geometry around partial first/last calendar years.
    for (sp in c(1L,3L,f)) {
      add('leading',f,6L*f + 3L,sp)
      add('trailing',f,7L*f + 2L,sp)
      add('both',f,9L*f + 4L,sp,backcasts=2L,forecasts=f + 1L,exclude=TRUE)
      add('protect',f,9L*f + 4L,sp,backcasts=2L,forecasts=f + 1L,exclude=FALSE)
    }
    # Projection horizons with exclude on/off and one full-span rejection case.
    add('horizon',f,7L*f,1L,backcasts=f,forecasts=f,exclude=TRUE)
    add('horizon',f,7L*f,2L,backcasts=f,forecasts=f,exclude=TRUE)
    # Sigma-limit variants allowed by upper > lower > 0.5.
    add('limits',f,7L*f + 5L,2L,lsig=0.6,usig=3.0)
    add('limits',f,7L*f + 5L,2L,lsig=2.0,usig=2.25)
    # Zero/tiny dispersion: the legacy 1e-5 rule and its unit dependence.
    add('tiny',f,8L*f,1L)
    add('zeroish',f,8L*f,3L)
    # Huge finite inputs: the pinned IEEE accumulation overflows to an infinite
    # dispersion and detection then keeps every weight at one.
    add('huge',f,7L*f+1L,2L)
    # Fully masked series: every window has no usable row, so no dispersion
    # exists and the oracle returns all-one weights after a single pass.
    add('maskall',f,6L*f,3L)
    # Spikes, ramps and the lone spike among zeros.
    add('spike',f,10L*f,1L)
    add('spike',f,10L*f,2L)
    add('revive',f,7L*f,1L)
    add('revive',f,7L*f,4L)
    add('mask',f,8L*f + 1L,2L)
    if (f == 12L) {
      add('mul',f,9L*f,1L,multiplicative=TRUE)
      add('mul',f,9L*f,7L,multiplicative=TRUE)
      add('mulzero',f,6L*f + 3L,1L,multiplicative=TRUE)
    }
  }
  cases
}
# Values are deterministic closed forms. Categories control the shape:
# smooth  - bounded pseudo-AR irregular with a moderate deterministic part;
# leading/trailing/both/protect - same smooth shape (window geometry focus);
# horizon - smooth plus one extreme late row inside the protected projections;
# limits  - smooth with one row per sigma-limit regime;
# tiny    - |values| < 1e-6 scale, triggering the legacy < 1e-5 dispersion rule;
# zeroish - exactly zero everywhere except one spike among eleven zeros;
# huge    - finite values near 1e300 so the pinned accumulation overflows;
# maskall - every row masked: no dispersion exists anywhere;
# spike   - smooth plus one huge row;
# revive  - smooth with one isolated spike among exact zeros in one year;
# mask    - smooth with explicit missing rows;
# mul     - multiplicative regime around 1 with deviations;
# mulzero - multiplicative ones except isolated deviations.
x11_sigma_input <- function(case) {
  n <- case$n; f <- case$frequency
  t <- seq_len(n)
  smooth <- 3 * sin(t * .7) + 2 * cos(t * sqrt(2)) + ((t * 37L) %% 11L) / 3.0 - 1.7
  values <- switch(case$category,
    smooth = smooth, leading = smooth, trailing = smooth, both = smooth,
    protect = smooth,
    horizon = { v <- smooth; v[n - 2L] <- 6 * 10^3; v[n - f] <- 6 * 10^3; v },
    limits  = { v <- smooth; v[f + 1L] <- 50; v[2L * f + 2L] <- 1e-9; v[3L * f + 3L] <- 25; v },
    tiny    = 1e-7 * smooth,
    zeroish = { v <- numeric(n); v[3L * f + 5L] <- 5e-3; v },
    huge    = 1e300 + 1e298 * smooth,
    maskall = smooth,
    spike   = { v <- smooth; v[4L * f + 3L] <- 1e6; v[7L * f + 9L] <- -8; v },
    revive  = { v <- numeric(n); v[f + 3L] <- 9; v[3L * f + 7L] <- -4; v },
    mask    = smooth,
    mul     = 1 + .35 * sin(t * .9) + .2 * cos(t * sqrt(3)),
    mulzero = { v <- rep(1, n); v[2L * f + 1L] <- 1.5; v[4L * f + 4L] <- .9999; v },
    stop('unknown sigma case category'))
  missing <- rep(FALSE, n)
  if (case$category == 'mask') missing[c(1L, 4L, f + 5L, 2L * f, n - 1L)] <- TRUE
  if (case$category == 'maskall') missing[] <- TRUE
  list(values = as.numeric(values), missing = missing)
}
