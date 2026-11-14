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
