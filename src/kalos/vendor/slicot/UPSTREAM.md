# Vendored SLICOT subset — provenance and license

This directory contains an unmodified, minimal transitive closure of SLICOT
Fortran sources, vendored to build the stationary-initialization solver seam
(SB03MD) from source inside the Kalos core. No SLICOT library is downloaded or
linked at build time; the closure's BLAS/LAPACK references resolve against the
same pinned OpenBLAS provider as the rest of the core.

## Upstream

| Item | Value |
| --- | --- |
| Project | SLICOT (`SLICOT/SLICOT-Reference`), version 5.9.1-dev (CMake project 5.9.1 on the post-5.9.1/v5.9.2-dev line, CITATION.cff v5.9.1) |
| Commit | `5eeb845d84f5a60f6069754b7819dad468951933` |
| Snapshot | cloned 2026-11-25 (M9 performance track, `tmp_opencode/m9-perf/candidates/SLICOT-Reference`); vendored here 2026-11-26; every file byte-identical to the pinned tree (`cmp` clean, per-file SHA-256 in `MANIFEST.json`) |
| License | BSD-3-Clause, "Copyright (c) 1996-2026, The SLICOT Team" — full text below; every vendored file carries `SPDX-License-Identifier: BSD-3-Clause` in its header. No NOTICE-file obligations, no third-party code of a different license, no copyleft components. |

## Vendored routines and why each is needed

The closure was computed as the exact set of objects the linker extracts when
the two `SB03MD` `JOB='B'` call paths (the seam's covariance solve with
`FACT='N'` and the certificate solve reusing the Schur with `FACT='F'`, both
`DICO='D'`, `TRANA='T'`) are linked: `-Wl,-Map` member extraction from the
reference build of the pinned tree, then verified complete by extracting the
members and confirming their only remaining undefined references are the
BLAS/LAPACK/gfortran symbols listed below (full record:
`tmp_opencode/m10-slicot/closure/`).

| Routine | Source | Role |
| --- | --- | --- |
| SB03MD | `SB03MD.f` | Discrete-time Lyapunov/Sylvester solver: solves `op(A) X op(A)' - SCALE*X = -CA*C'` (here the transposed Stein equations `P - T P T' = Q` and `S - T S T' = I`, `JOB='B'` also producing the informational SEP/FERR) |
| SB03MV | `SB03MV.f` | SEP / eigenvalue-condition computation for `JOB='B'` |
| SB03MW | `SB03MW.f` | Sylvester-equation solve kernel used by SB03MD |
| SB03MX | `SB03MX.f` | FERR (forward-error estimate) computation for `JOB='B'` |
| SB03MY | `SB03MY.f` | Schur-form ordering/exchange helper used by SB03MD |
| SB04PX | `SB04PX.f` | Triangular Sylvester solver (`DTRSYL`-equivalent with SLICOT's scaling logic) used by SB03MW |
| MB01RD | `MB01RD.f` | Matrix operation `C = op(A)*op(B) ± op(C)` used to assemble the equation right-hand sides |
| SELECT | `select.f` | Eigenvalue selection function passed to `DGEES` by SB03MD |

The SLICOT-internal routine `SB03MU` and all other archive members are **not**
referenced by these call paths and are deliberately not vendored. The
SLICOT-bundled `src/lapack_aux/` fallback routines are **not** vendored: every
external symbol the closure references exists in the pinned OpenBLAS 0.3.34
provider (verified against both `libopenblas.so` and `libopenblas.a` with
`nm --defined-only`, F77 names), so no fallback is needed:

```
daxpy_  dcopy_  ddot_  dgees_  dgemm_  dlabad_  dlacn2_  dlacpy_  dlaln2_
dlamch_ dlanhs_ dlascl_ dlaset_ dlasy2_ dscal_  dswap_  dsymv_  dtrmm_
idamax_ lsame_  xerbla_
```

## Runtime dependency note

The vendored sources are legacy fixed-form FORTRAN 77 compiled by gfortran, so
every deliverable that links them additionally carries `libgfortran` (and, for
gfortran itself, `libquadmath` where present). `libgfortran` is GPL-3 **with
the GCC Runtime Library Exception** and is distributable in linked binaries; it
is recorded in provenance the same way the pinned OpenBLAS provider is.

## BSD-3-Clause license text (upstream `LICENSE`)

```
Copyright (c) 1996-2026, The SLICOT Team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
