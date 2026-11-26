# rkalos implementation changes needed

Implementation checked 2026-11-23 on macOS ARM64. Development version `0.1.0.9001` uses the existing Kalos snapshot `0e064dd`.
Refreshed 2026-11-27 at commit `7ae3949`, version `0.1.0.9032`, core `6fee68c`; all five review items remain fixed and were re-certified through every subsequent re-pin. The historical verification section is retained as recorded.

## Status refresh (2026-11-27)

Since the review the package grew to the full product surface and carried every core
advance: the estimated X-13 (RSA0/RSA1/RSA2c) and TRAMO-SEATS (RSA0/RSA1/RSA3/RSAfull)
routes, the automatic route with live-parity TD F-test numerics, and Persephone's
native-default cutover as its backend. The maintainer's roxygen help-files pass and
configure restructure merged cleanly into the re-pin line (`186e612`, build surface
re-certified inert at `f317565`). Core pins advanced `0e064dd → … → ec17543` (M11: the
vendored SLICOT closure is built from source through the Fortran probe in `configure`,
with the loud QR-fallback clause and `core_version()$slicot_stationary_initialization`
on the probe surface) `→ 6fee68c` (M12: certified post-pivot diffuse filter). The M12
re-cert executed the maintainer ruling of 2026-11-26 with measured numbers: fit-matrix
stop gate PASS (1517 checks both lanes, coefficient max 4.87286502e-07 — still below the
phase-1 level, 0.49% of its 1e-4 assert), **zero hex fixtures moved** (re-pin set empty;
measured fitted-models drift 3.34e-14 / 5.18e-16), the two ruled pins raised
`2e-8 → 1e-7` with dated in-test justification, and `field-tolerances.csv` kept
byte-identical because its governed JVM-replay path is immune to native-core changes
(the expected theta raise measured away; recorded in the ledger). Suite state:
**310,308 serial / 310,312 OMP, 69 files, zero fails**, per-file tables byte-identical
to the certified baselines; `verify-nojava --native` transcripts byte-identical to the
m8/m10/m10-s2 baselines in both lanes.

## Remaining items (2026-11-27)

1. **Hosted platform qualification**: Windows/macOS CI with a Fortran-capable toolchain
   for the vendored SLICOT closure; a Windows build without `FC` runs the documented
   QR-everywhere fallback (loud configure warning, qualified OFF state per the upstream
   ledger) — that mode needs its first hosted witness.
2. **Publication**: CRAN/GitHub release stays maintainer-side; `R CMD check` Status OK
   with vignettes was certified at the M8 review.
3. **Maintainer repo hygiene (blocks in-repo export)**: stale untracked build residue —
   `src/kalos/**/*.o` (91 files quarantined to `tmp_opencode/m10-s2/stale-objects/` plus
   ~41 under `src/`), older `src/Makevars`, `src/rkalos.so`, `config.log`, an old
   tarball — must be purged for in-repo `sync-core.py` runs to pass the export guard.
4. No open functional items: the roxygen follow-up was discharged by the maintainer's
   own help-files pass, and the standing re-pin convention (sync-core → fresh lanes →
   suites/verifies/fit-matrix) is the routine for future core advances.

## P1 — Fix the X-11 adapter build

**Fixed.** The B, C and D pass adapters call the qualified table converter. A clean source archive installed and loaded. All **4,162** B/C/D assertions passed. The build uses the core’s unfused floating-point setting, and the R test reader preserves binary64 values from the frozen decimal fixtures. [B adapter](src/rcpp_x11_b_pass.cpp) · [Fixture reader](inst/tinytest/helpers/fixtures.R)

## P1 — Give automatic-fit selection tables the correct calendar

**Fixed.** The automatic wrapper now derives observed and extended calendars from the same start and the effective core horizons. Selection metadata, B/C/D rows, cropped spans and original-domain metadata follow those calendars. In-sample finals keep the observed start. All **524** horizon assertions passed with both Homebrew and checksum-pinned OpenBLAS. They cover a year-crossing backcast, forecast-only selection and equality between vector and time-series inputs. [Calendar decoration](R/native_x13_automatic_fit.R) · [Horizon regressions](inst/tinytest/test-native-x13-horizon.R)

## P1 — Reject malformed inputs before integer conversion and allocation

**Fixed.** The R boundary validates scalars and a univariate numeric series before coercion, checks supplied dates against a time series and uses a fixed twelve-filter decoder. The native boundary checks integer ranges, ordered one-based spans, unsigned budgets and extended lengths before calling the pinned core. Quarterly and other unsupported frequencies retain the monthly-only status. [R request](R/native_x13_automatic_fit.R) · [Native request](src/rcpp_x13_automatic.cpp) · [Regression tests](inst/tinytest/test-native-x13-automatic.R)

## P1 — Provision and check the provider before package installation

**Fixed for the current build. Platform qualification remains planned for M6.** Configure compiles, links, loads and exercises an LP64 CBLAS/LAPACKE probe with the package flags. Shared OpenBLAS on macOS installs without the stale Fortran linker path. Static archives use R’s Fortran link flags. The CI workflow provisions the checksum-pinned provider before dependency installation and checks its identity. The README documents the automatic fit’s temporary process-wide thread override. [Configure](configure) · [CI workflow](.github/workflows/check.yml) · [Installation instructions](README.md)

## P2 — Provide a supported fit interface for Persephone

**Fixed.** Both qualified fits are documented exports in version `0.1.0.9001`. Installed-package fixture tests invoke each export directly. The restricted-fit suite passed all **645** assertions. [Public interface](man/native_x13_fit.Rd) · [Restricted-fit tests](inst/tinytest/test-native-x13-fit.R) · [Automatic-fit tests](inst/tinytest/test-native-x13-horizon.R)

## Verification of the current fixes

The clean source archive loaded under R 4.5.2 with LP64 OpenBLAS 0.3.34. All **6,539 assertions** across the six affected suites passed, including **1,208/1,208** automatic-fit assertions. The offline core, fixture and evidence verification passed. Hosted CI and Windows installation remain unverified here.

## Original review checks at `b90d367`

Offline core, fixture, and evidence verification passed. The horizon probe used a scratch installation with qualified converter calls and the redundant `FLIBS` removed for the local shared provider. The full R test suite was outside this review run.
