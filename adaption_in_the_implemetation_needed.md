# rkalos implementation changes needed

Status checked 2026-11-23 at `7d72e6c`, with Kalos pinned to `0e064dd`. `R/` and the Rcpp adapters own R input validation, result classes, and calendar metadata. The generated core and fixtures are managed by `tools/sync-core.py`.

## P1 — Fix the X-11 adapter build

**Needs work now.** The core pin advanced, but the B-pass adapter still fails the local Clang compilation check. The same naming conflict is present in C and D.

A clean `R CMD INSTALL` of the checkout failed with R 4.5.2, Clang and Rcpp 1.1.2 on macOS ARM64. In each of the B, C and D pass adapters, the `SEXP values` parameter conflicts with `using kalos_x11_rcpp::values;`. Subsequent table-conversion calls are parsed as calls to the parameter. Qualify the converter calls or rename the local arguments, preserving the R-facing argument names and result schemas. Require a clean source archive install and run the B/C/D boundary tests against that installation. [B pass](src/rcpp_x11_b_pass.cpp#L101) · [C pass](src/rcpp_x11_c_pass.cpp#L208) · [D pass](src/rcpp_x11_d_pass.cpp#L183)

## P1 — Give automatic-fit selection tables the correct calendar

**Needs work now.** The M4 horizon route already ships. Its selection serializer and R calendar decoration still use the original period and observed span.

The `multi_b3h6` horizon fixture fitted successfully in a scratch installation. Its extended `out` and B1 both had 105 rows. `out` starts in October 1989, while B1 was labeled from year 1 and its input domain and row basis described 96 rows. The R wrapper decorates nested B/C/D tables with the observed length and start. Its plain-vector ordinal uses `frequency + start_period - 1`, omitting `first_year`. The C++ projection also supplies the original first period and zero horizons to the selection serializer. [R calendar decoration](R/native_x13_automatic_fit.R#L130) · [Selection serialization](src/rcpp_x13_automatic.cpp#L399) · [Shared pass decoration](R/native_x11_pass_decorate.R#L8) · [Horizon fixture tests](inst/tinytest/test-native-x13-horizon.R#L144)

Build the observed calendar once and derive the extended calendar from the effective core horizons. Use the extended start, first period and length for selection tables, and retain the observed domain for in-sample finals. Test dates, domain lengths, horizon metadata, and row indices for B/C/D tables with a backcast crossing a year boundary, a forecast-only run, and matching `ts` and plain-vector inputs.

## P1 — Reject malformed inputs before integer conversion and allocation

**Needs work now.** The current R wrapper and native conversions still accept malformed arguments before validating their shape and range.

By code inspection, the automatic wrapper applies `as.integer()` to calendar and horizon arguments before checking exactness and range, flattens matrices with `as.vector()`, and sizes the seasonal request vector using the resulting frequency. A fractional argument can be truncated and a multivariate series can be treated as one longer series. The native adapter also converts signed outlier budgets to `size_t` and subtracts one from span endpoints before validating their native range. [R request construction](R/native_x13_automatic_fit.R#L30) · [R series and date conversion](R/native_x13_automatic_fit.R#L103) · [Native budget and span conversion](src/rcpp_x13_automatic.cpp#L298)

Reuse the applicable `native_inputs.R` checks and enforce a single univariate numeric series with one consistent calendar. For `ts` inputs with explicit dates, require agreement or use the same explicit date for both the request and the returned series. Reject fractional, missing, and out-of-range integers, malformed spans, and negative budgets before conversion. Keep core statuses for supported argument types outside the fit envelope. Adopt the checked horizon normalization from Kalos when refreshing the snapshot. [Existing R input checks](R/native_inputs.R#L1) · [Core findings](https://github.com/statistikat/kalos/blob/24e4452/adaption_in_the_implemetation_needed.md)

## P1 — Provision and check the provider before package installation

**Needs work now; platform qualification is planned for M6.** Correct the provider probe, link flags, and installation instructions for the current package. M6 owns the broader macOS and Windows installation checks. [Platform qualification record](https://github.com/statistikat/persephone3/blob/4b51d49/IMPLEMENTATION_STATUS.md#L915)

The current build uses the OpenBLAS C API. The README still tells readers it uses R’s configured BLAS/LAPACK, and the CI workflow reaches dependency setup without explicitly provisioning the qualified C API provider. `configure` link-tests two routines, while the final Makevars also adds `FLIBS`. On this host, after fixing the C++ conflicts in the scratch copy, the final link failed on a stale `gfortran` path even though the shared OpenBLAS provider link probe passed. [README](README.md#L1) · [Configure provider probe](configure#L40) · [Final link flags](src/Makevars.in#L1) · [CI installation](.github/workflows/check.yml#L29)

Provision a known LP64 CBLAS/LAPACKE build before CI package installation on each supported platform. Use the final provider and linker flags in the configure probe, including any Fortran dependencies needed by a static build. Use a shared-provider path without redundant `FLIBS` when that provider already resolves its runtime. Verify that a clean installed package loads and reports the expected provider. Update the installation instructions and clarify the scope of the automatic fit’s provider thread override. [Thread override](src/rcpp_x13_automatic.cpp#L379)

## P2 — Provide a supported fit interface for Persephone

**Needs work now.** Publish the interface used by the two existing X-13 routes. The integration record identifies the private namespace dependency without assigning its replacement to a later milestone.

Fit capabilities identify the qualified routes, but Persephone retrieves both functions by name from the private namespace. Publish a narrow exported or explicitly versioned adapter for those two fits, with the existing result schemas and scope gates. Add an installed-package contract test for both calls, and migrate Persephone once it can pin a release with the adapter. [Fit capabilities](R/capabilities.R#L12) · [Persephone callers](https://github.com/statistikat/persephone3/blob/4b51d49/R/backend_native.R#L536)

## Original review checks at `b90d367`

Offline core, fixture, and evidence verification passed. The horizon probe used a scratch installation with qualified converter calls and the redundant `FLIBS` removed for the local shared provider. The full R test suite was outside this review run.
