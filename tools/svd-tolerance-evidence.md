# Identification-fixture SVD tolerance: evidence and recommendation

Investigation date: **2026-11-15** (milestone M1 groundwork). Host: Ubuntu 24.04.4,
g++ 13.3.0, R 4.5.3, BLAS/LAPACK `libopenblasp-r0.3.26.so` (LAPACK 3.12.0,
pthread build) reached through `R CMD config` `-lblas -llapack`. rkalos built from
`rkalos_0.1.0.tar.gz` at core commit `4f0625cf44d9cba3c512fdd3d97e1ec4b10132da`.

Status: analysis only. **No test was changed.** The recommendation is at the end.

## 1. Which assertions actually fail

`inst/tinytest/test-native-identification.R` runs 196 assertions here, 2 fail.
Both are line 18, in the first loop (simple two-column transfer Jacobians):

```r
expect_equal(result$singular_values, values, tolerance = 1e-14)
```

for `initial = list(ar = -0.6, ma = 0.6 + 1e-6)` at period 4 **and** period 12,
i.e. the `near_redundant_transfer` case. Identification by relaxation:

| variant                                   | file result |
| ----------------------------------------- | ----------- |
| as committed                              | 196 / 2 fails |
| line 18 tolerance relaxed to `1e-8`       | 196 / 0 fails |
| line 19 slack relaxed to `< 1e-4`         | 196 / 2 fails |

Reported gap: `Mean relative difference: 1.294733e-10` (matches the 1.29e-10 in the
M0 note). Line 19 (`singular_ratio` vs `min(values)/max(values)`) passes: its
absolute error is 4.7e-17 because the ratio is ~3.7e-7.

Correction to the M0 wording: the failing Jacobians are **not** the structured
seasonal ones. The structured seasonal block (8x3 at period 4, 16x3 at period 12)
has `kappa_2 = 2.98` and passes the same comparison with room to spare at 1e-14.
The failing matrices are the small 3x2 near-redundant ones.

## 2. Measured quantities (period 4 and 12 give identical numbers)

Normalized fixture matrix `B = sweep(expected, 2, sqrt(colSums(expected^2)), "/")`:

| quantity | value |
| --- | --- |
| `sigma_max` | 1.4142135623729997 |
| `sigma_min` | 5.1993122734984004e-07 |
| `kappa_2` | 2.7200e+06 |
| `kappa_2 * eps` (forward-error floor of any backward-stable SVD) | 6.04e-10 |
| `ULP(sigma_max)` | 3.140e-16 |
| observed `abs(sigma_min_native - sigma_min_target)` | 6.732e-17 = 0.21 ULP(sigma_max) |
| observed relative gap on `sigma_min` (what `all.equal` reports) | 1.2947e-10 |
| requested tolerance | 1e-14 |

Per-case survey of all 18 first-loop combinations (`64*eps*sigma_max` = 2.0e-14):

```
period  a     diff    kappa2      abs gap    1e-14 ok   bound ok
  4    -0.6    0      8.06e+15   1.85e-17    TRUE       TRUE
  4    -0.6   1e-6    2.72e+06   6.73e-17    FALSE      TRUE
  4    -0.6   0.2     1.49e+01   5.55e-17    TRUE       TRUE
  4     0      0      inf        0.00e+00    TRUE       TRUE
  4     0     1e-6    2.00e+06   0.00e+00    TRUE       TRUE
  4     0     0.2     1.01e+01   2.78e-17    TRUE       TRUE
  4     0.4    0      1.67e+16   6.66e-16    TRUE       TRUE
  4     0.4   1e-6    2.32e+06   2.22e-16    TRUE       TRUE
  4     0.4   0.2     1.09e+01   4.44e-16    TRUE       TRUE
 12    (identical to period 4 for every case)
```

`a = 0, difference = 1e-6` has the same conditioning (`kappa_2 = 2.0e+06`) and a
zero gap; `a = 0.4, difference = 1e-6` has a **larger** absolute gap (2.2e-16) and
still passes. The pass/fail line is therefore not conditioning alone but whether
the last bits of the two input matrices happen to agree: it is a knife-edge, and
any rebuild, compiler or LAPACK build can move other cases across it in either
direction.

## 3. Attribution: input rounding, not the LAPACK target

Measured with `tools/svd-tolerance-probe.cpp` (direct `F77_CALL(dgesvd)` /
`F77_CALL(dgesdd)` through the same provider rkalos uses) and
`tools/svd-tolerance-probe.R`:

1. **Provider and algorithm are excluded.** For the fixture-normalized matrix,
   `svd(B, nu = 0, nv = 0)$d`, a direct `DGESVD("N","N")` call and a direct
   `DGESDD("N")` call return *bitwise identical* singular values, and rkalos'
   native values are bitwise identical to `DGESVD` of the native-normalized native
   Jacobian. R's own LAPACK module (`R_HOME/modules/lapack.so`) resolves `dgesdd_`
   from the same `/usr/lib/x86_64-linux-gnu/liblapack.so.3` -> `libopenblasp-r0.3.26.so`
   that rkalos links, so both sides of the comparison already use one LAPACK build.
   There is no "native DGESVD vs R svd()" divergence to fix.
2. **The inputs differ by less than one ULP.** Native column norms are bitwise
   equal to `sqrt(colSums(expected^2))`; of the six matrix entries, five are
   bitwise equal and one differs by 0.8 ULP, giving `max|A - B| = 1.110e-16`.
3. **The gap is fully explained by that difference.** Weyl's bound gives
   `|d sigma_min| <= ||A - B||_2 = 1.110e-16`, and the observed gap is 6.732e-17.
   Pushing any single nonzero fixture entry by +1 ULP moves `sigma_min` by up to
   8.08e-17, i.e. the same order as the observed gap. `sigma_min ~ 5.2e-07` is only
   determined to absolute accuracy `eps * sigma_max`.
4. **Not threading, not determinism.** Repeating the identification with
   `KALOS_NUM_THREADS=1,2,8` and with the serial (no-OpenMP) build reproduces the
   same bytes.

## 4. What the macOS baseline did

On the macOS ARM64 baseline (R 4.5.2, Homebrew Clang, R-configured BLAS/LAPACK) the
same two assertions passed at 1e-14, so the tolerance held there. The mechanism
that makes the assertion satisfiable is input agreement, not a more accurate LAPACK:
on AArch64 `long double` is identical to `double`, so the C++ assembly of the
Jacobian and of its column norms performs exactly the same rounding steps as the
double-only R expression in the fixture, the two matrices are bitwise equal, and
`all.equal` then compares two identical SVD inputs. On x86_64 `long double` is
80-bit extended, so a wider intermediate (`arma_identification.cpp:63` accumulates
the squared column norms in `long double`, and the ARMA parameter/polynomial paths
do the same) leaves the native matrix one fraction-of-a-ULP away from the fixture
matrix, which the conditioning of 2.7e6 amplifies into a 1.3e-10 relative error on
`sigma_min`.

Supporting measurements from this host: the native and R column norms are bitwise
equal here, and exactly one Jacobian entry differs by 0.8 ULP, which is the
signature of a different intermediate-rounding path rather than a different LAPACK.
The macOS side itself was not re-run for this note, so that half is inference from
the recorded baseline plus the platform rounding model; the recommendation below
does not depend on it.

## 5. Recommendation: change the fixture tolerance, keep the native target

Keep `DGESVD("N","N")` in `kalos::singular_values()`. It is the correct routine, it
agrees bitwise with R's `svd()` and with `DGESDD` on identical inputs on this host,
and the residual gap (0.21 ULP of `sigma_max`) is below the agreement any
backward-stable SVD can reach for inputs that differ by a ULP. A target-side change
would have to chase a relative accuracy of 1e-14 on a singular value whose
theoretical relative-error floor is `kappa * eps = 6.0e-10`; that is unreachable and
would only be "fixable" by extended-precision or one-sided orthogonal iteration,
which is out of scope for a rank-diagnostics helper.

Change the fixture instead, in `inst/tinytest/test-native-identification.R`
(that file is rkalos-owned: `core.lock.json` `generated_files` covers only
`src/kalos/**` and `inst/tinytest/fixtures/**`, verified; no frozen fixture byte is
involved, so no Kalos round-trip is needed). Suggested replacement for line 18:

```r
  # Singular values are only determined to eps * sigma_max; sigma_min of a
  # near-redundant transfer has relative sensitivity kappa. Compare in the norm.
  expect_lt(max(abs(result$singular_values - values)),
            64 * .Machine$double.eps * max(values))
```

Properties of that form on this fixture family: bound 2.0e-14, worst observed gap
6.7e-17, so ~300x margin on Linux while still failing loudly on any real regression
of the singular values (the well-conditioned block-B cases sit at 1e-16 absolute,
the rank-deficient cases at 1e-16 with `sigma_min` two orders apart). If a relative
form is preferred instead, derive it from the conditioning of the target:

```r
  expect_equal(result$singular_values, values,
               tolerance = max(1e-13, 8 * max(values) / min(values) * .Machine$double.eps))
```

Leave alone: line 17/50 (matrix comparisons, passing), line 19 (absolute ratio
slack, passing and structurally safe), and lines 21-22 (status classification:
`near_redundant_transfer`, and `rank` 2, are correct on this host, which is the
actual contract of the file). Do not widen the two `1e-13`/`1e-14` matrix
comparisons; they are unaffected by this effect.

## 6. Reproduction

```sh
PKG_LIBS="$(R CMD config LAPACK_LIBS) $(R CMD config BLAS_LIBS) $(R CMD config FLIBS)" \
  R CMD SHLIB -o svd-tolerance-probe.so tools/svd-tolerance-probe.cpp
Rscript --vanilla tools/svd-tolerance-probe.R svd-tolerance-probe.so
```

Verified with this tree: the probe compiles warning-free with g++ 13.3 and the
driver reproduces every number above.

The two failing assertions can be confirmed without the probe:

```sh
Rscript --vanilla -e 'library(rkalos); tinytest::run_test_file(
  file.path(system.file("tinytest", package = "rkalos"),
            "test-native-identification.R"), verbose = 0)'
```
