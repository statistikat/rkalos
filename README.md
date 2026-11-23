# rkalos

A lean R/Rcpp interface to [Kalos](https://github.com/statistikat/kalos), with
committed core sources and an LP64 OpenBLAS CBLAS/LAPACKE provider. No Java is required.

```r
rkalos::core_version()
y <- ts(100 + sin(seq_len(120)/3), start=c(2000,1), frequency=12)
prepared <- rkalos::prepare_inputs(y)
objective <- rkalos::arma_objective_setup(y, ar=NA_real_, initial=list(ar=0.2))
fit <- rkalos::arma_optimize(objective)
rkalos::fixed_likelihood(y, ar=0.3)$likelihood
```

`capabilities()` describes delivered numerical operations. Use the exported `native_x13_fit()` and `native_x13_fit_automatic()` for the restricted and automatic monthly X-13 routes described there. Use explicit status
fields when deciding whether inference or prediction results are available.
The accurate covariance policy is the default; legacy comparison modes remain
Java-free. Existing `persephone_native_*` result class names are retained for
serialization compatibility. Persephone consumes this package as a required
runtime dependency; its application objects stay in Persephone.

For TRAMO preprocessing, use `tramo_critical_value()` for the outlier threshold,
`tramo_seasonality()` for the monthly detector, and `tramo_transform_decision()`
for the non-robust log/level decision. For a complete monthly series,
`native_tramoseats_fit()` estimates the RSA0 airline model and returns the
seasonally adjusted series, trend, seasonal and irregular components.
The `final` table includes component standard deviations and forecast
uncertainty. The `merged` table contains the final product forecasts.
For the preprocessing original-series forecast error (`y_ef`), use
`series_forecast_stdev` and its `series_forecast_stdev_status` field.
The fitted residuals and their standard error are available as
`fit$residuals$tsres` and `fit$residuals$ser`. The 24-lag Ljung–Box
test of residual independence is available as `fit$residuals$lb`. The
Doornik–Hansen test of residual normality is available as
`fit$residuals$doornikhansen`. Each test reports its statistic,
upper-tail probability and degrees of freedom.

```r
fit <- rkalos::native_tramoseats_fit(AirPassengers, template = "RSA0")
fit$finalSA
fit$merged$values[[which(fit$merged$component == "Trend" &
                         fit$merged$information == "Forecast")]]
```

`capabilities()$tramo_seats_fit$templates` lists the qualified raw-series
templates. Use `seats_from_tramo()` to decompose supported fitted TRAMO
coefficients supplied by a caller. The monthly legacy SEATS retries cover
airline seasonal MA replacement and the mean-free (2,1,0)(1,1,0) seasonal
AR-to-MA refit, including the merged product series for fitted AO/LS effects.

Use `native_tramoseats_fit(series, "RSA1")` for monthly TR1 estimation
and SEATS extraction from the original observations. It restores detected
additive outliers, level shifts and transitory changes (monthly decay 0.7) in the final components. Its per-template capability
descriptor lists available residual diagnostics and forecast errors.

`native_tramo_preprocess(series, "TR1")` fits the stock monthly TR1
preprocessing model. Its result contains the selected transformation,
fitted airline parameters, dated residuals and estimated outlier effects.
The residual standard error and the Ljung–Box and Doornik–Hansen statistics
and probabilities are available on both monthly templates.

Use `native_tramo_preprocess(series, "TR3")` or
`native_tramoseats_fit(series, "RSA3")` for monthly automatic SARIMA
order selection. Fitted models accepted at the first reference comparison
and subsequent controller sweep expose the orders, coefficients, dated
residuals and, for RSA3, SEATS components. A further identification
round, fast mean-test retry or seasonality restart raises an explicit error.

Install OpenBLAS with CBLAS and LAPACKE (qualification version 0.3.34). Set `KALOS_OPENBLAS_ROOT` to its prefix before building and installing the source archive. For custom or static providers, set `KALOS_OPENBLAS_CPPFLAGS` and `KALOS_OPENBLAS_LINK` to the full compiler and linker flags including dependencies.
C++17 and Rcpp are required. Configure probes optional OpenMP with R's actual
compiler; failed probes select serial execution. `KALOS_NUM_THREADS` takes
precedence over the legacy `PERSEPHONE_NUM_THREADS` alias. The automatic X-13 fit temporarily sets the process-wide OpenBLAS thread count to one and restores its previous value on return.

## Core synchronization

```sh
python3 tools/sync-core.py --core ../kalos --kind r --revision <commit>
python3 tools/sync-core.py --verify
```

Do not edit `src/kalos` or generated fixture snapshots. Change Kalos, commit it,
then refresh and review `core.lock.json`. Verification needs no Git/network.
Source installation never fetches core code. Tests live in `inst/tinytest`.
Numerical contracts and canonical fixtures belong to Kalos; this repository owns
R conversion, result formatting and R test qualification.

See [implementation status](IMPLEMENTATION_STATUS.md).
