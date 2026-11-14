# rkalos

A lean R/Rcpp interface to [Kalos](https://github.com/statistikat/kalos), with
committed core sources and R's configured BLAS/LAPACK. No Java is required.

```r
rkalos::core_version()
y <- ts(100 + sin(seq_len(120)/3), start=c(2000,1), frequency=12)
prepared <- rkalos::prepare_inputs(y)
objective <- rkalos::arma_objective_setup(y, ar=NA_real_, initial=list(ar=0.2))
fit <- rkalos::arma_optimize(objective)
rkalos::fixed_likelihood(y, ar=0.3)$likelihood
```

`capabilities()` describes delivered numerical operations. Full native seasonal
adjustment and automatic model selection remain unavailable. Use explicit status
fields when deciding whether inference or prediction results are available.
The accurate covariance policy is the default; legacy comparison modes remain
Java-free. Existing `persephone_native_*` result class names are retained for
serialization compatibility. Persephone consumes this package as a required
runtime dependency; its application objects stay in Persephone.

Build locally with `R CMD build .` and `R CMD INSTALL rkalos_0.1.0.tar.gz`.
C++17 and Rcpp are required. Configure probes optional OpenMP with R's actual
compiler; failed probes select serial execution. `KALOS_NUM_THREADS` takes
precedence over the legacy `PERSEPHONE_NUM_THREADS` alias. The configured BLAS
provider's thread settings are left untouched.

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
