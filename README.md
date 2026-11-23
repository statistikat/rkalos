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
for the non-robust log/level decision. The full TRAMO/SEATS fit is pending.

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
