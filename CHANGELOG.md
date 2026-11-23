# Changelog

## Development version

Version 0.1.0.9017 expands monthly TR3 and RSA3 fits to the seasonal
controller restart and a fourth estimation pass after the mean is removed.

Version 0.1.0.9016 adds dated original-series forecast standard errors
from the fitted RSA1 preprocessing model on the level and log scales.

Version 0.1.0.9015 adds monthly TR3 automatic identification and the
RSA3 decomposition for fits accepted by the first reference comparison
and the subsequent controllers without an automatic restart.

Version 0.1.0.9014 adds the public TR1 residual standard error,
Ljung--Box independence test and Doornik--Hansen normality test to
preprocessing and complete RSA1 seasonal adjustments.

Version 0.1.0.9013 restores stock monthly transitory changes in raw RSA1
and caller-supplied TRAMO fits on both level and log scales.

Version 0.1.0.9012 restores additive AO and LS effects to fitted and raw
monthly SEATS decompositions. The raw RSA1 fit now returns the public
products for the two level-scale outlier cases.

Version 0.1.0.9011 extends the raw monthly TRAMO–SEATS product fit to
RSA1, including its selected transformation, detected AO/LS effects
and dated final components and forecasts.

Version 0.1.0.9010 adds TR1 preprocessing from the original monthly
series, including automatic transformation and outliers, model and
likelihood estimates, and dated residuals.

Version 0.1.0.9009 adds the public RSA0 Doornik–Hansen residual
normality statistic and its chi-squared upper-tail probability.

Version 0.1.0.9008 supports two monthly legacy SEATS retries, including
the fitted mixed AO/LS decomposition at default and custom horizons.

Version 0.1.0.9007 adds the public RSA0 residual Ljung–Box statistic,
upper-tail probability, degrees of freedom and lag count.

Version 0.1.0.9006 exposes the dated RSA0 preprocessing residual series
and its scalar standard error.

Version 0.1.0.9005 adds dated preprocessing standard errors for RSA0
original-series forecasts, with a separate projection status.

Version 0.1.0.9004 adds native monthly RSA0 estimation and decomposition
from the observed series, including fitted likelihood, uncertainty and
dated product forecasts and backcasts.

Added deterministic-effects restoration and merged product components for
supplied full-span fitted monthly TRAMO model coefficients.

Added the pinned TRAMO first-round differencing decision for monthly and
quarterly series with standard module thresholds.

Version 0.1.0.9003 adds standard TRAMO template model preparation and
fixed-model SEATS component extraction with dated initial and final series
and their standard deviations.

Version 0.1.0.9002 exposes dated automatic X-13 forecast and backcast
tables, including raw series, adjusted series, trend and seasonal factors.

Added R interfaces for the TRAMO critical value, monthly seasonality detector,
and non-robust log/level decision against the pinned Kalos 23-B contract.

## 0.1.0 — initial repository split

Preserved the delivered Kalos numerical scope through 22b.2.
See IMPLEMENTATION_STATUS.md for repository-specific work and qualification.
