# rkalos implementation changes needed

Implementation checked 2026-11-23 on macOS ARM64. Development version `0.1.0.9001` uses the existing Kalos snapshot `0e064dd`.

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
