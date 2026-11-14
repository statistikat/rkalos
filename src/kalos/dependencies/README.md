# Numerical provider baseline

OpenBLAS 0.3.34 is pinned in `openblas.lock.json`. The upstream tag archive was
downloaded and its SHA-256 independently checked on 2026-11-14; it matches the
installed Homebrew source recipe. Release information:
https://github.com/OpenMathLib/OpenBLAS/releases/tag/v0.3.34

Standalone/Python builds currently discover a preinstalled LP64 provider with
CBLAS and LAPACKE; R uses its own configured BLAS/LAPACK. No provider binaries
are bundled in this initial split. The source archive checksum is a reproducible
baseline, not a claim that arbitrary system installations were built identically.

Future distributable wheels must compile Fortran LAPACK, enable CBLAS/LAPACKE,
use 32-bit integers, a conservative architecture baseline and CPU dispatch,
private symbols/library names, and single-threaded OpenBLAS with concurrent-call
locking. Validate both NumPy import orders before publishing those wheels.
No `-march=native` is prescribed. Kalos never changes global BLAS thread settings.
The bundled notices reproduce the upstream provider/LAPACK/CBLAS licenses.
