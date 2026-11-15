# rkalos implementation status

Updated: **2026-11-14**. Initial independent R package implemented and locally
qualified. Persephone requires rkalos; this package has no Persephone or Java
dependency. Numerical roadmap and canonical reference evidence belong to
[kalos](../kalos/IMPLEMENTATION_STATUS.md).

## Completed split steps

1. Preserved 102 relevant historical commits and original-to-filtered maps.
2. Added the R-declared nine-routine BLAS/LAPACK provider outside the core snapshot.
3. Adapted every delivered native R/Rcpp entry point through 22b.2 to Kalos.
   Shared window/alignment/regression/Jacobian preparation now runs in C++.
   R retains ts conversion, naming, serialization and result display metadata.
4. Exposed lean public engine interfaces, capabilities and embedded version/provider
   metadata. Required runtime dependencies are only Rcpp and stats, plus base R.
5. Pinned **84 core/provider-notice files** and **112 fixture files** to Kalos
   `4f0625cf44d9cba3c512fdd3d97e1ec4b10132da`. Snapshots are generated and verified
   offline; source installation never fetches GitHub. Source archives include the
   verification script and lock under `inst/core.lock.json`.
6. Migrated 36 native test files and added a boundary/metadata/thread-control test.
   Frozen fixture bytes are untouched; a test-only lexical adapter handles one
   historical Persephone namespace lookup. Small legacy domain/forecast result
   formatters retain serialized compatibility without importing Persephone.

## Qualification

- Source archive `R CMD check`: **266,549 passing assertions / 37 files**,
  zero errors and zero notes. One existing R-header/Clang warning about
  `-Wfixed-enum-extension` remains; it is not a package-code diagnostic.
- All numerical tests and the boundary assertion run without loading Java/rjd3.
- Focused preparation equivalence: **4,281 assertions**.
- Archive audit: every generated core and fixture byte matches its source lock.
  Corruption checks reject changed, missing and unexpected snapshot files.
- R 4.5.2, macOS ARM64, Homebrew Clang 23.1.5, R-configured BLAS/LAPACK,
  two actual OpenMP workers. The same core passes all 37 standalone serial,
  OpenMP and sanitizer programs using OpenBLAS.
- The final source archive reports the configured BLAS library and LAPACK
  version. The serial archive passed **266,545 assertions**, and Persephone
  passed 35 backend assertions with Java-related packages absent.

## Linux (x86_64) baseline qualification — 2026-11-15

First qualification run on a new Linux machine. Nothing here restates macOS
results; only what was executed on this host is recorded.

Environment:

- Ubuntu 24.04.4 LTS, kernel 6.8.0-137-generic, x86_64; g++ 13.3.0
  (Ubuntu 13.3.0-6ubuntu2~24.04.1); R 4.5.3 (x86_64-pc-linux-gnu).
- Installed only into a scratch library (`/tmp/opencode/rlib-b`); no system
  R library was written. Note: this machine's site-wide `Rprofile.site`
  hardcodes `.libPaths()` and ignores `R_LIBS_USER`, so test sessions ran with
  `Rscript --vanilla` and an explicit `R_LIBS` listing the scratch library
  first. Machine configuration, not package code; nothing in the repository
  was changed for it.
- Numerical libraries actually used: R's `-llapack -lblas` resolve through
  Ubuntu alternatives to OpenBLAS 0.3.26 (pthread build); LAPACK is the
  netlib-derived LAPACK 3.12.0 inside `libopenblasp-r0.3.26.so`, reached via
  the `F77_CALL` provider in `src/lapack_r.cpp` (no LAPACKE). Confirmed by
  `core_version()` and `ldd` of the installed shared object.
- `configure`'s C++17 OpenMP compile/link/runtime probe passed (OpenMP 201511
  via libgomp); parallel execution is enabled in this build.

Results:

- `tools/sync-core.py --verify`: passed offline; 84 core files and 112 fixture
  files match commit `4f0625cf44d9cba3c512fdd3d97e1ec4b10132da`.
- `R CMD build` and `R CMD INSTALL` of the archive into the scratch library
  succeeded with zero compiler warnings and zero notes.
- Test suite (`tinytest::run_test_file` over all 37 files, no Java loaded):
  240,563 assertions ran; 239,905 passed; 658 failed; 14 of 37 files aborted
  partway. Every X-11 file, plus inputs, likelihood, linear algebra,
  parallel, prediction-moments, primitives and the package-boundary file ran
  fully green (235,011 assertions); the boundary test confirms rJava/rjd3
  namespaces are not loaded. The 658 failures sit in arma-residuals (288),
  fitted-models (354), identification (2) and optimizer (14). Aborted:
  arma-objective, diffuse, fixed-likelihood, forecast-arma, forecasts,
  information, interpolation, interpolation-arma, joint-information,
  legacy-residuals, regression, regression-information, tramo-information,
  tramo-legacy.
- Not run on this machine: full `R CMD check` of the archive, sanitizer and
  serial-rebuild lanes, the focused preparation-equivalence set, and the
  macOS/Windows lanes.

Blocker (pinned core, not fixed here):

- `src/kalos/src/finite_state_factor.cpp` `checked()` throws
  `underflow_error("finite state factor arithmetic underflows the finite
  numeric range")` when an `long double` intermediate is nonzero but its
  cast to `double` is zero. On x86_64 `long double` is 80-bit extended
  precision; on the macOS ARM64 baseline `long double` is identical to
  `double`, so the same products flush to zero before the guard and the
  evaluation returns `available`. A gdb backtrace for the failing case
  shows `FiniteStateFactor::observation` (finite_state_factor.cpp:135)
  accumulating an exact sum of 6.91e-334 during the diffuse filter of the
  period-12 seasonal AR likelihood (`sarima_innovations` →
  `regression_likelihood_impl`), which casts to 0 and throws.
- Effect on Linux: every likelihood/profile evaluation of a model containing
  seasonal AR/MA factors (periods 4 and 12, airline included) returns
  `numerical_failure`, producing `gradient_unavailable`/`unresolved` fits,
  NULL sensitivities and the aborted/failed test files above. Forcing
  `OMP_NUM_THREADS=1` does not change any of it, so this is a deterministic
  extended-precision semantic difference, not thread nondeterminism and not
  the LAPACK provider. Fixing it requires a Kalos core decision (guard
  semantics or accumulation representation); `src/kalos` is generated and was
  not touched here.
- Secondary divergence: two `test-native-identification.R` assertions fail
  with mean relative difference 1.29e-10 between native `dgesvd` singular
  values of structured seasonal numerator Jacobians and the R `svd()` target
  (tolerance 1e-14). OpenBLAS rounding at borderline conditioning; the
  affected statuses and ranks are correct.
- README example on this machine: `prepare_inputs`, `arma_objective_setup`
  and `fixed_likelihood` behave as documented (`fixed_likelihood(ar=0.3)`
  returns `available`, log-likelihood −680.6269, Q = 4939.992 on the sample
  series), but `arma_optimize` on that level-plus-sinusoid series ends at
  ar = 0.9999972 with status `gradient_unavailable` at the unit-root
  boundary. Not yet cross-checked against the macOS baseline; tracked as an
  open item.

## Compatibility

Accurate final-fit covariance and accurate residual/interpolation reporting stay
at their existing defaults; explicit legacy modes remain Java-free. JDemetra MSR
iteration behavior stays the default, with the published cap also available.
Historical `persephone_native_*`/forecast result classes are retained. This
experimental package does not implement full seasonal-adjustment fitting.

## Next work

1. Track integration and platform checks after the initial GitHub publication.
2. Add R conversions/tests only as Kalos delivers each next numerical slice.
3. Keep pinned source/fixture refreshes reviewable, regenerate Rcpp registrations
   when signatures change, and run archive checks and serial/no-Java lanes.
4. Expand platform qualification through the configured Linux/macOS/Windows CI.
   Remote workflow configuration is not a claim that those runs have passed.

The unresolved sigma-weight default for numerical item 22b.3a belongs to Kalos;
it is unrelated to the completed required-dependency decision.
