# k = 0 stationary-factorization routing probe corpus (roadmap item 30, probe slice)

Frozen probe measurements supporting
[`docs/item30-k0-routing-probe-2026-11-28.md`](../../../../docs/item30-k0-routing-probe-2026-11-28.md).
This is a probe corpus, not a replay golden: no production entry point
consumes it. Its role is to freeze (a) the current-state per-evaluation
ratios of the k = 0 shapes at `1e2e3c8` against the pinned item-29 JVM
reference, (b) the interleaved route A/B medians and per-step slopes for
the five measured k = 0 routes, (c) the correctness gates (filter-copy
fidelity, option-(ii) bit-identity, factorized-route deviation), and
(d) the cost attribution (instrumented sections + SIGPROF self-sampling).
Inputs are the frozen item-29 diffuse-steady-state-probe corpus
(`../diffuse-steady-state-probe/`), which already pins the only two k = 0
shapes (`syn-stat-sarma11`, `syn-stat-arma21`); no new series were minted.

Provenance

- kalos build: `main` at `1e2e3c896ade4f67124f057c0aa5a0fedd4fc734`, Release
  (`-O3 -DNDEBUG -ffp-contract=off`), single thread
  (`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=KALOS_NUM_THREADS=1`), OpenBLAS
  0.3.34 built in scratch (`LD_LIBRARY_PATH=.../openblas/lib`), GCC 13.3.0,
  Intel Xeon Gold 6342 (8 cores), host load 3.2-5.6 on 8 cores for every
  timed bracket (`brackets.txt`, lane brackets per lane); no round was
  discarded (zero reload retries recorded).
- `drivers/m15-item30-filter_copy.cpp` is a byte-copy of
  `src/diffuse_state_space.cpp` at the same commit, with the namespace
  renamed to `probe30` (`sed '18s/^namespace kalos {$/namespace probe30 {\nusing namespace kalos;/'`)
  plus the probe-only additions, each marked with an inline `item-30 probe`
  comment: a route
  parameter, the TSC section scopes, the k = 0 factor-arming routes (1, 2),
  and the option-(ii) structural-witness path (route 3). Copy fidelity is
  verified by the `lib == copy0` byte-identity gate in `tables/correctness.txt`.
- Drivers are compiled against `libkalos.a` of that build with
  `-O2 -DNDEBUG -ffp-contract=off -fopenmp -no-pie -g` (see
  `drivers/m15-item30-build_drivers.sh`; library build recipe in
  `drivers/m15-item30-build.sh`). The whole lane set regenerates with
  `drivers/m15-item30-capture.sh` after editing its `root` to the scratch
  checkout holding the build tree.

Routes measured (all via one process, per-rep route-order rotation)

- `lib` — HEAD `sarima_innovations` (the dense `sandwich`/`joseph` fallback).
- `copy0` — faithful namespaced copy of the HEAD path (fidelity control).
- `fac1` — k = 0 armed with `FiniteStateFactor` + zero-dimension certificate
  and frame (the k > 0 canonical call shape).
- `fac2` — factor at k = 0 without the history frame (filter guard relaxed
  in the copy; isolates the frame bookkeeping).
- `struct3` — option-(ii) skeleton: same covariance domain, same ascending-k
  term order, exact-zero structural terms skipped under a witness
  (`Z = e_1`, companion `T`, `h = 0`), persistent buffers, raw `values()`
  pointer views (the pointer views replaced the out-of-line bounds-checked
  `DenseMatrix::operator()` after sampling attributed 25.8% of the dense
  evaluation to that accessor; the bit-identity was re-verified after the
  replacement — see `correctness.txt`).
- `stat` — the existing `stationary_innovations` primitive (option-(iv)
  falsification; timed only).

Tables

- `current-state-ratios.txt` — kalos per-case eval medians (median of three
  cost-probe rounds) vs JVM warm eval=1..5 medians extracted from the frozen
  `jvm/` captures of the item-29 corpus (`ratio_table.py`).
- `eval-round*.txt`, `init-round*.txt` — raw cost-probe rounds (capture lane).
- `ab-*-round*.txt` — interleaved route A/B: `REP` lines per repetition,
  `MEDIAN` per route. Lanes: default (`ab-*`), `laneB-*` (quietest lane,
  load 3.2-4.4), `laneC-*` (adds the `stat` route), `laneA-structv1-*`
  (pre-pointer-replacement skeleton, kept for the accessor-attribution
  step). Per-phase load brackets: `brackets.txt`, `lane*-brackets*.txt`.
- `nsweep-*.txt` — per-length medians (n = 24..240) for the slope/intercept
  decomposition.
- `instr-*.txt` — TSC-calibrated section decomposition (`SECTIONS`,
  `INSTRWALL`); the instrumented run carries ~4 scope pairs per step, so
  sections are quoted as fractions of their own `INSTRWALL`, not as wall.
- `sample-syn-stat-sarma11.raw.txt` / `.symbols.txt` — SIGPROF self-sampling
  (1 ms ITIMER_PROF, PC from `ucontext`, `-no-pie` symbolization via
  `m15-item30-symbolize.py`) over 150 library evaluations; perf is blocked
  on this host (paranoid=4), so the in-process sampler replaces the item-29
  gdb-as-parent recipe with the same attribution role.
- `correctness.txt` — byte-identity gates (cmp) and, for the factorized
  route, decision-line comparison + max relative/absolute deviation
  (`m15-item30-dump_compare.py`; the dump is every named `DiffuseLikelihood`
  quantity at `%.17g`, statuses and counters included).
- `golden-rerun.txt` — the frozen item-29 goldens re-ran at HEAD: byte-
  identical on the two k = 0 cases (and the gapped case), D4 in-band
  movement on the k > 0 arming cases, as adopted by M12.

Reading guide for the two k = 0 corpus shapes: `syn-stat-sarma11` is the
m = 14 expanded SARMA(1,0,1)(1,0,1,12); `syn-stat-arma21` is m = 2. The
JVM reference medians for these shapes (0.4052 ms / 0.1602-0.1504 ms) come
from the item-29 corpus `jvm/jvm-stat-sarma11.txt` and
`jvm/jvm-stat-arma21n.txt` (calibrated twin), unchanged here.
