# Changelog

## Development version

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
