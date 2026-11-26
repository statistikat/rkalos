# Diffuse steady-state filter probe corpus (roadmap item 29, probe slice)

Frozen probe inputs, exact-path goldens and JVM cross-validation captures
supporting `docs/item29-diffuse-filter-probe-2026-11-26.md`. This is a probe
corpus, not a replay golden: no production entry point consumes it yet. Its
role is (a) to freeze the regime windows (stationary / stationary-after
differencing / near-unit-root / rank-deficient) with bit-exact decision
patterns for any future steady-state implementation, and (b) to anchor the
JVM cross-validation numbers quoted in the dossier.

Provenance

- kalos build: `main` at `ec175430643caa293e5db24ad3f2ca62e3afc44d`, Release
  (`-O3 -DNDEBUG -ffp-contract=off -fopenmp`), single thread
  (`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=KALOS_NUM_THREADS=1`), linked to
  OpenBLAS 0.3.34 (built in scratch), GCC 13.3.0, Intel Xeon Gold 6342.
  `goldens/` bytes re-ran bit-identical (`KALOS-GOLDEN-RERUN-IDENTICAL`).
- JVM captures: `jvm/` are pass-1 files; every case re-ran in a second
  independent JVM process with byte-identical deterministic lines
  (`jvm-determinism.txt`; TIMING lines excluded by design). Temurin
  OpenJDK 21.0.11, `-XX:ActiveProcessorCount=1`, jars
  `jdplus-toolkit-base-core-3.9.0.jar`
  `1880e23dc1665be2822c30c4da11ea0f2482a8da5c3e9c80b372ff6b8375dd71` and
  `jdplus-toolkit-base-api-3.9.0.jar`
  `fed5fe554c1e64e371d40a1a55d449323cc2630a0d5ea33a3a34e47f7926c2c8` (the
  pinned R-library copies, SHA-matched against the probe cache).
- Host load was 2.7–4.0 (8 cores) for every timed phase, recorded in
  `tables/cost-phases.txt` before and after each phase; no phase was
  discarded.
- Synthetic series come from `drivers/gen_series.cpp` (fixed 64-bit LCG +
  Box-Muller, seeds in `build_corpus.sh`), byte-reproducible without host
  RNG version sensitivity. Real series are byte-copies of the pinned
  M9 gap-analysis probe inputs (`probe-12-airline-*`, digests below via
  `SHA256SUMS.csv`).

JVM coefficient convention (calibrated, not assumed): the JVM
`SarimaModel.Builder.phi/bphi` values equal the negated kalos `ar/sar`
lists (the JVM prints the polynomial verbatim); `theta/btheta` equal the
kalos `ma/sma` values verbatim. `jvm/jvm-stat-arma21.txt` (uncalibrated) vs
`jvm/jvm-stat-arma21n.txt` (calibrated, ssq bit-identical to the kalos
golden) demonstrates both branches.

Layout

- `cases.csv` — 15 probe cases (regime, orders, coefficients, series,
  design, masks); `build_corpus.sh` regenerates `series/` and `x/`.
- `goldens/<id>.txt` — exact-path output of `drivers/filter_probe.cpp` at
  variance 1: sums (`%.17g`), numerics counters, `FTRACE` (finite
  variance trace), `UKIND`, `FBAR`, and `TSTAR tau=…` transient lengths.
- `jvm/jvm-<case>.txt` — `drivers/Item29FilterProbe.java` captures:
  `MACHINE` (dim, steadyPos, C0, h0), `STRACE`, frozen gain rows,
  `LL`/`FULLCOMPUTER`/`GLSCOEF` and per-evaluation timing lines.
- `tables/cross-validation.txt` — kalos-vs-JVM agreement table
  (`drivers/summarize.py`). `tables/cost-phases.txt` — per-phase medians
  (eval/init/reg/nsweep/tail) with loadavg brackets.

Regime inventory (capture windows): `real-airline-*` + `leadgap-airline-x`
(flagship airline, gaps, leading-mask delay); `syn-airline12/4`,
`syn-mixed4` (stationary after differencing); `syn-stat-arma21`,
`syn-stat-sarma11` (intrinsically stationary); `nearur-airline`,
`nearni-airline`, `slow-seasonal` (near-boundary: AR 0.995 vs d=1, MA 0.98,
sma 0.999); `d2-arma12` (double differencing);
`real-airline-redundant-x`, `real-airline-zerocol-x` (rank-deficient
designs; both must return `status=rank_deficient` with the frozen NaN
pattern). Any future steady-state path must replay every decision bit and
every value in-band against these bytes.
