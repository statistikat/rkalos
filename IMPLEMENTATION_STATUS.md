# rkalos implementation status

Updated: **2026-11-15**. Initial independent R package implemented and locally
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

## Linux (x86_64) M1 groundwork — 2026-11-15

Same host and toolchain as the baseline above (Ubuntu 24.04.4, g++ 13.3.0,
R 4.5.3, OpenBLAS 0.3.26 providing LAPACK 3.12.0). The pinned core was **not**
rebuilt or edited here; `src/kalos` is byte-identical to commit
`4f0625cf44d9cba3c512fdd3d97e1ec4b10132da` and the upstream underflow-guard fix
is still pending integration. Only what was executed on this host is recorded.

### Archive `R CMD check`

- `R CMD build .` then `R CMD INSTALL` into a scratch library: clean, **zero
  compiler warnings**. The macOS `-Wfixed-enum-extension` R-header warning does not
  occur with g++ 13.3, so no warning of any kind is left on this platform.
- `R CMD check --no-tests` (default options, `_R_CHECK_FORCE_SUGGESTS_=false`, PDF
  manual built): **Status OK** — zero errors, zero warnings, zero notes; the only
  non-OK line is the informational installed-size report (58.3 Mb, of which 27.5 Mb
  are the shipped test fixtures). Documentation, NAMESPACE, Rd usage, examples,
  portability, `$(BLAS_LIBS)`/`PKG_*FLAGS` and compiled-code checks all pass. **No
  repository-level check issue had to be fixed.**
- Full `R CMD check` including tests: **1 error, 0 warnings, 0 notes**. That single
  error is the known pinned-core blocker: `tests/tinytest.R` stops on the first test
  file with `Error: finite state factor arithmetic underflows the finite numeric
  range`, so the check does not reach a per-file report. Known, not new.
- No check issue introduced by this milestone: the two archived runs differ only by
  the test phase.

### Java independence

`R CMD check --no-tests` was repeated with `JAVA_HOME` unset in a library set where
`rJava`, `rjd3*` and `persephone3` are **not installed** (only `Rcpp`, `tinytest`
and base R): **Status OK**, exit 0. In that same session the boundary, ten X-11 and
five preparation/native test files run green with `package:rJava` absent from the
search path. This machine's default `java` is OpenJDK 21.0.11 (17 and 8 also
installed); neither the check nor the tests consult any JVM.

### Serial (no-OpenMP) lane

- `configure` accepts no arguments and reads no dedicated switch: its single lever
  is the C++17 OpenMP compile/link/runtime probe, which by design honours user
  Makevars overrides. The serial lane was therefore selected the supported way — a
  `R_MAKEVARS_USER` file blanking `SHLIB_OPENMP_CXXFLAGS` — so the probe fails and
  `configure` prints "OpenMP unavailable for this C++17 toolchain; using serial
  execution". A `--disable-openmp`/`KALOS_OPENMP` style flag does not exist; adding
  one is an open usability item.
- Verification that the lane really is serial: `src/Makevars` gets empty OpenMP
  flags, no compile line carries `-fopenmp`, the installed shared object has no
  `libgomp` dependency and no `GOMP_*`/`omp_*` dynamic symbols (the OpenMP build has
  all three), and `core_version()$parallel` reports `compiled = FALSE`,
  `openmp_version = 0`. Installed separately in `rlib-b-serial`.
- Whole suite in the serial lane: **240,565 assertions, 658 failures, 14 of 37 files
  aborting**; the OpenMP lane re-measured here gives **240,563 / 658 / 14**, matching
  the baseline above. Failure and abort sets are identical, and so are the abort
  reasons — the blocker is deterministic extended-precision semantics, not threading.
  The only per-file difference is `test-native-parallel.R` (37 assertions serial,
  35 with OpenMP, from its `compiled = FALSE` branch), zero failures in both.
- X-11 lane plus preparation files (12 X-11 files, plus inputs, primitives,
  likelihood, linear algebra and prediction moments): **234,968 assertions per lane**
  (234,320 X-11 + 648 preparation), identical per-file counts, zero failures.
- Equality stronger than counts: each of those 17 files was re-run under a scratch
  harness that hashes the serialized bytes of every value an assertion compares. All
  17 file digests match between the two lanes across 408,959 recorded values, so
  serial and OpenMP return byte-identical results rather than merely equally-rounded
  ones. (The harness itself mis-scores 4 argument-validation assertions in
  `test-native-prediction-moments.R`, identically in both lanes; the official
  tinytest run reports no failure there.)

### SVD fixture tolerance

`tools/svd-tolerance-evidence.md` (with `tools/svd-tolerance-probe.cpp` and
`tools/svd-tolerance-probe.R`) analyses the two remaining
`test-native-identification.R` failures. Summary: they are line 18
(`expect_equal(result$singular_values, values, tolerance = 1e-14)`) in the
`ar = -0.6, difference = 1e-6` case at period 4 and 12 — proven by relaxing that one
line, which takes the file from 2 failures to 0 — and they concern the small
3x2 near-redundant Jacobians (`kappa_2 = 2.72e6`), not the structured seasonal ones,
whose `kappa_2` is 2.98 and which pass at 1e-14. R's `svd()`, a direct `DGESVD` and a
direct `DGESDD` agree **bitwise** on identical inputs here, and the native values are
bitwise `DGESVD` of the native-normalized Jacobian, so the LAPACK target is excluded.
The gap is 6.73e-17 absolute — 0.21 ULP of `sigma_max`, inside the Weyl bound from the
single 0.8-ULP input difference between native and fixture matrices — amplified by the
conditioning into the reported 1.29e-10 mean *relative* difference, whose theoretical
floor is `kappa * eps = 6.0e-10`. **Recommendation: change the fixture tolerance**
(compare singular values against `64 * .Machine$double.eps * max(values)` in the norm,
or derive the tolerance from `kappa`), **keep `DGESVD` as the target**. The test file
is rkalos-owned (`core.lock.json` generates only `src/kalos/**` and
`inst/tinytest/fixtures/**`), so this needs no Kalos round-trip and touches no frozen
fixture byte. The test was **not** changed in this milestone.

### `arma_optimize` unit-root boundary (README example)

Root-caused on Linux; it is **not** the pinned-core long-double issue. The README fit
ends at coordinate 6.7469984232914353 (`ar = 0.99999723157599829`, value
1.1100490271451684) with status `gradient_unavailable`, 257 evaluations, and a
rejection histogram of `parameter_boundary = 97` — no `numerical_failure` anywhere, no
underflow or overflow guard message in a 121-point scan of the profile, and identical
bytes with 1, 2 and 8 threads and in the serial build. Mechanism:
`arma_optimizer.cpp` searches in `w = tanh(reflection coordinate)`, so the optimizer's
fixed step `difference_step = 1e-4` is applied in a saturating variable where
`dc/dw = 1/(1-w^2) = 1.8e5` at that point; every proposal with `w + h >= 1` is rejected
as `parameter_boundary`, which removes the central stencil, and the surviving
one-sided stencils span `c = 4.6 … 6.4` — far too wide for the declared accuracy
budget. All 11 refinement attempts were reproduced exactly: `error` 1.9e1 … 1.0e4
against budgets of 1.4e0 … 1.4e2, all rejected. The objective itself is smooth and
available there (its own second difference falls off as `h^2` down to 1e-12, and a
stencil taken in the reflection coordinate gives `error = 7.7e-8` against a budget of
`3.9e-5`, i.e. accepted), and the profile minimum is interior: a Newton step from the
stopping point gives `c = 6.7655` (`ar = 0.99999733`), so the optimizer arrived at the
optimum but could not certify its gradient. This is boundary/coordinate handling in the
optimizer wrapper (Kalos side), double-only arithmetic, and expected to reproduce on
macOS; that cross-check still has to be run there.

### Not run on this machine

Sanitizer lanes, the focused preparation-equivalence set as its own target, Windows,
macOS, and any run against the fixed core (still upstream).

## Core integration checkpoint — 2026-11-15

The pending upstream fixes have since landed and the pin advanced to
`cd004b08bcb6cdf1e5b22b2cf34af721ccbaf06c` (via `tools/sync-core.py`, `--verify`
clean; 88 core + 132 fixture files). Two core commits moved in: `f96603f` (the
ARIMA(0,1,1) underflow guard, which removes all fourteen suite aborts and the
658 stale-slope failures) and `cd004b0` (profiled variances published as binary64
quotients of the named binary64 sum). The second was found through this package:
on `f96603f` the full suite ran with a single remaining failure class — 15 of 64
pinned forecasts disagreed in the `mean` attribute at `tolerance = 0` (1 ULP on
the profiled variance) — a scratch-width redivision in the forecast kernels that
the core's standalone suite did not cover. Fixed in Kalos (39/39 per lane, both
lanes), re-synced, re-qualified here.

- SVD comparison: the recommendation above is applied (commit `0696fb1`);
  `test-native-identification.R` is 196/196.
- Whole suite, OpenMP lane: **266,549 assertions, 0 failures, 0 skipped, 0 errors**
  across all 37 files — the first all-green Linux state.
- Whole suite, serial lane (`compiled = FALSE` recompiled): **266,545 / 0**. The
  four-assertion difference from the OpenMP lane is `test-native-tramo-information.R`'s
  six OpenMP-only thread-reproducibility checks versus `test-native-parallel.R`'s
  two serial-specific checks; both lanes report zero failures, errors and skips.
- Qualification ran from an `R CMD build` tarball because `sync-core.py --verify`
  rejects in-tree object files; the serial lane again used the `R_MAKEVARS_USER`
  probe-override and `LD_LIBRARY_PATH` was unset at test time so the pinned
  0.3.34 provider (not the system BLAS) served every LAPACK call.
- pykalos advanced to the same pin; its wheel rebuilt warning-free and smoke-verified
  (`core_version()` reports `cd004b08…`).

## Item 22b.3a/3b at the R boundary — 2026-11-16

- New native entry points `.native_x11_sigma_weights` (input, frequency, start
  period, explicit `legacy`/`accurate` policy, projection counts and exclusion,
  sigma limits, mode) and `.native_x11_replacements` (input, frequency, final
  weights, projection handling, mode), with internal R wrappers
  `native_x11_sigma_weights()`/`native_x11_replacements()` for ts-domain
  handling, 1-based window/exclusion rows and NA-padded neighbor positions.
  The policy argument has **no default**: which sigma rule becomes the product
  default remains an open Kalos maintainer decision, so nothing is promoted.
  Neither entry point is exported yet; the objects were added to the
  `Makevars(.in)/(.win.in)` object lists (the two new core translation units
  were synced but not yet compiled).
- `test-native-x11-sigma.R`: all 125 pinned JDemetra cases reproduce bit-exactly
  through R under the legacy policy (final and first-pass weights, both
  dispersion-window sets, pass counts, exclusion sets, per-row states), with
  recomputation determinism, state/weight consistency and the accurate-policy
  invariants (finite, in-range, masked rows weight one, no legacy cutoff state).
- `test-native-x11-replacement.R`: all 22 pinned cases reproduce bit-exactly
  including the 22b.3a→3b chain (the computed legacy weights equal the fixture
  weights), both projection-masking captures (`mode` 0/1), neighbor provenance
  and counts, period averages, correction factors, the corrected series and the
  apply-corrections identity, fallback counting, validation errors and
  wrapper output shapes.
- Full suite after the addition: **270,653 assertions, 0 failures** across all
  39 files on the OpenMP lane; serial lane **270,649, 0 failures** — the same
  four-assertion conditional difference as before (six OpenMP-only thread
  reproducibility checks versus two serial-specific parallel-info checks);
  both new files contribute identical counts per lane (3,390 sigma + 714
  replacement). The serial lane was installed from a fresh `R CMD build`
  tarball with the `R_MAKEVARS_USER` probe override (0 warnings,
  `compiled = FALSE`, no `omp_*` symbols).

## Core integration checkpoint and B pass at the R boundary — 2026-11-16

- The pin advanced to `16f27859896cc2c2faae69e72e7c24e088f93153` (via
  `tools/sync-core.py`, `--verify` clean offline): one core commit,
  `16f2785` "Qualify the complete fixed additive B pass (22b.4 S1)". The lock
  delta was verified from the diff: 88→90 core files (added
  `src/x11_b_pass.cpp` and `include/kalos/x11_b_pass.hpp`), 132→138 fixture
  files (`native/x11-bpass-cases.R` plus the five `native/x11-bpass/` files);
  only four pre-existing hashes changed (commit, source manifest,
  `CMakeLists.txt`, `revision.hpp`). The snapshot matches the core's own
  source manifest exactly: every `add_library` translation unit,
  `x11_b_pass.cpp` included, is present in `src/kalos`.
- New native entry point `.native_x11_b_pass` (input, frequency, start period,
  explicit `legacy`/`accurate` policy, fixed per-phase initial and final
  seasonal requests, fixed Henderson length, sigma limits) with the internal
  wrapper `native_x11_b_pass()` for ts-domain handling: the 19 pinned member
  tables come back as `ts` on their original calendars (B2..B4g on the
  cropped span via `.native_time_domain`, as in the initial-table wrapper),
  named per-phase filters align by calendar period, and the correction
  history carries 1-based analysis rows with analysis periods plus
  row-basis/available-row provenance. The policy argument has **no default**.
  Neither entry point is exported; `rcpp_x11_b_pass`/`x11_b_pass` were added
  to the `Makevars(.in)/(.win.in)` object lists and the Rcpp registration was
  regenerated (additive-only diff).
  Boundary notes: the 20th public table (B9d) is not a struct member but the
  second round's detection series, so the wrapper exposes 19 member tables
  plus three rounds carrying detection, analysis position, policy, sigma
  weights, replacements, correction factors and replacement counts — the
  final round's sigma weights and correction factors are exactly the pinned
  B17 and B20. The core's validation vocabulary (monthly/quarterly only,
  in-range first period, odd Henderson length ≥ 5, span gates, per-period
  request support) passes through unchanged.
- `test-native-x11-bpass.R`: all 18 pinned fixed additive B-pass cases
  reproduce bit-exactly through R under the legacy policy — 20 tables per
  case on their pinned spans (NaN-aware), correction-round policies, pinned
  analysis positions, replacement counts cross-checked against the B4/B9
  correction tables, round-3 flagged-row counts against B17, the detection
  and B17/B20 correction-history linkages, the additive composition
  identities recomputed in R, determinism, input immutability, the
  accurate-policy runs, the wrapper calendar/named-filter/plain-vector paths
  and the core's validation errors (1,395 assertions).
- Full suite after the addition: **272,048 assertions, 0 failures** across all
  40 files on the OpenMP lane; serial lane **272,044, 0 failures** — the same
  four-assertion conditional lane difference as before. The serial lane was
  installed from a fresh `R CMD build` tarball with the `R_MAKEVARS_USER`
  probe override (`nm` reports no `gomp` symbols).

## Core integration checkpoint and C pass at the R boundary — 2026-11-16

- The pin advanced to `7e8e38dad6d458e5a2827c6feca8dcf3c9326e80` (via
  `tools/sync-core.py`, `--verify` clean offline): one core commit, `7e8e38d`
  "Qualify the complete fixed additive C pass (22b.4 S2)". The lock delta was
  verified from the diff: 90→93 core files (added `src/x11_c_pass.cpp`,
  `include/kalos/x11_c_pass.hpp` and the new shared internal header
  `include/kalos/detail/x11_pinned.hpp`, which the S2 commit factored out of
  the B pass and which the manifest does export into the snapshot), 138→143
  fixture files (the five `native/x11-cpass/` files); six pre-existing hashes
  changed (commit, source manifest, `CMakeLists.txt`, `revision.hpp`,
  `x11_b_pass.hpp`, `x11_b_pass.cpp` — the S2 commit reworked the B-pass
  sources onto the shared header). The snapshot matches the core's own source
  manifest exactly.
- New native entry point `.native_x11_c_pass` (input, frequency, start period,
  explicit `legacy`/`accurate` policy, fixed per-phase initial and final
  seasonal requests, fixed Henderson length, **the required B-pass companion**,
  sigma limits) with the internal wrapper `native_x11_c_pass()`. Design
  decision on the companion: the wrapper takes the caller's finished
  `native_x11_b_pass()` result object and rebuilds the pinned
  `kalos::X11BPass` from it faithfully (NA marks the pinned missing flags,
  1-based analysis rows back to the 0-based core convention) instead of
  silently recomputing the B pass. Recomputing would double the cost and
  could paper over a stale or mismatched B result; passing the caller's object
  through makes the pinned core's own coupling errors fire exactly as in the
  pinned `X11Kernel` composition (completed-status, setting-echo, round- and
  sigma-policy match, B20 coverage, B1-equals-input). The 14 pinned C tables
  come back as `ts` on their original calendars, `C2`/`C4`/`C5a` on the
  cropped span (length n−2·`c2drop`, `c2drop` = frequency/2) and every other
  table on the full span, with the final-stage C13 correction round (analysis
  position, period, C13/C17/C20 linkage) and row-basis/available-row
  provenance. The policy argument has **no default**; neither entry point is
  exported; `rcpp_x11_c_pass`/`x11_c_pass` were added to the
  `Makevars(.in)/(.win.in)` object lists and the Rcpp registration was
  regenerated (additive-only diff).
- `test-native-x11-cpass.R`: all 18 pinned fixed additive C-pass cases
  reproduce bit-exactly through R under the legacy policy — 14 tables per
  case on their pinned spans (NaN-aware), echoed settings and `c2drop`, the
  pinned flagged/corrected row counts, the round policy, position and
  C13/C17/C20 linkages, the additive composition identities recomputed exactly
  in R (C1 = Y − B20, C4 = crop(C1) − C2, C6 = C1 − C5, C9 = C1 − C7,
  C11 = Y − C10, C13 = C11 − C7), determinism, input immutability, the
  accurate-policy composition, the B-coupling rejections (mismatched round
  policy, broken status, altered B1 series, non-covering B companion, malformed
  companion object), the core's validation errors and the wrapper
  calendar/named-filter/plain-vector/missing-companion paths
  (1,124 assertions).
- Full suite after the addition: **273,172 assertions, 0 failures** across all
  41 files on the OpenMP lane; serial lane **273,168, 0 failures** from the
  same `R CMD build` tarball with the `R_MAKEVARS_USER` probe override — the
  same four-assertion conditional lane difference as before.
- Upstream feedback for the Kalos maintainers (recorded here so the milestone
  report carries it; no contract files were touched):
  - (S1) The public B view lists 20 tables but `X11BPass` holds 19 members —
    B9d lives only as the second round's `detection` series (rounds pinned in
    b4d, b9d, b13 order); the contract should state that mapping explicitly.
  - (S1) B17/B20 are bit-identical duplicates of the final round's sigma
    weights and `replacements.correction_factors` — the replacement values
    themselves exist only in the round history — so the contract should say
    plainly that B20 holds correction *factors*, not replacements.
  - (S1) Gate ordering: the quarterly "three annual radii plus one cycle"
    threshold (n ≥ 10) is shadowed by the trend gate (n ≥ trend_length) for
    every quarterly request with `trend_length > 10`, including the pinned
    fixture default of 13, so that message is only reachable with
    `trend_length ≤ 9`; the monthly threshold (n ≥ 30) surfaces naturally.
  - (S2) The cpass fixture's `inputs.csv` is bit-identical to the bpass
    fixture's (verified with `cmp` against the pinned tree) — the same 18
    series and case ids — so the cpass cases replay directly on the bpass
    case list (rkalos reuses `native/x11-bpass-cases.R`); a note in the
    fixture metadata would save integrators from building a second loader.

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
