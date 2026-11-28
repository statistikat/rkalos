# Item-30 k = 0 structural byte-identity goldens (roadmap item 30, adoption slice)

Frozen byte-identity oracle supporting the adoption of option (ii) — the
structure-aware dense k = 0 steps shipped in `src/diffuse_state_space.cpp`.
Design and acceptance budget:
[`docs/item30-k0-routing-probe-2026-11-28.md`](../../../../docs/item30-k0-routing-probe-2026-11-28.md)
§7; maintainer adoption ruling of 2026-11-28. This is the acceptance object of
that budget made replayable: every named `DiffuseLikelihood` quantity, decision
and status, captured at `64c5804` (before the change) and replayed
line-for-line at byte level through the shipped library by
[`tests/test_diffuse_k0_goldens.cpp`](../../../../tests/test_diffuse_k0_goldens.cpp)
(ctest, 24 arms, 46 857 compared lines).

Provenance

- Capture build: `main` at `64c580496` (the probe-slice tip), Release
  (`-O3 -DNDEBUG -ffp-contract=off -fopenmp`), single thread
  (`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=KALOS_NUM_THREADS=1`), OpenBLAS 0.3.34
  from the pinned provider (`LD_LIBRARY_PATH=.../openblas/lib`), GCC 13.3.0,
  Intel Xeon Gold 6342. Capture-lane load brackets recorded by
  `drivers/m16-capture-goldens.sh` (the panel is deterministic; the brackets
  document the capture, as in the k0-routing-probe corpus).
- Replay build: `main` at `2306b1e` with the option-(ii) steps. All 24 arms are
  byte-identical to this panel (capture-vs-replay `cmp` clean on the full 46
  883-line file), so the change is value-neutral everywhere it runs.

Layout

- `k0-goldens-x86_64.txt` — the oracle panel. Each arm is an `ARM <id>` header,
  a `CFG` request line, then the dossier §7 dump protocol (the same lines as
  `../k0-routing-probe/drivers/m15-item30-route_probe.cpp dump`: STATUS,
  COUNTS, six SUM lines, NUMERICS, RES, UKIND, and every PRED/FVAR/INNO/STDI/
  FSTATE/INITCOV/INITDIFF/FINALCOV/FINALDIFF entry at `%.17g`).
- `drivers/m16-item30-goldens.cpp` — the capture driver (HEAD libkalos.a only;
  the single oracle for both lanes).
- `drivers/m16-build.sh`, `drivers/m16-capture-goldens.sh` — the
  configure/build and capture recipes with the load brackets.
- `SHA256SUMS.csv` — digests of the files in this directory.

Arms (24)

- Both frozen k = 0 corpus shapes — `syn-stat-sarma11` (m = 14 expanded
  SARMA(1,0,1)(1,0,1,12)) and `syn-stat-arma21` (m = 2) — the two shapes the
  probe measured and the only k = 0 shapes in the frozen corpus.
- The thirteen non-k0 controls: `real-airline-{noreg,x,gaps-x,redundant-x,
  zerocol-x}`, `syn-airline12`, `syn-airline4`, `syn-mixed4`, `nearur-airline`,
  `nearni-airline`, `d2-arma12`, `slow-seasonal`, `leadgap-airline-x` — the
  D4-armed flagship, the gapped and near-unit-root, period-4 and differenced
  small shapes, which must stay byte-unchanged (they run the factor path,
  gated off from the structural steps).
- The §7 coverage-row arms the budget assigns to the implementation slice:
  `k0-sarma11-gap` (interior gapped k = 0), `k0-arma21-mean` (mean ≠ 0, the
  legal stationary case), `k0-nearboundary` (AR 0.999) and
  `k0-nearboundary-seasonal` (seasonal AR 0.999) near-boundary stationary
  roots, and `k0-zerohorizon` (n = 0 empty series through the full assembly).
- The canonical k = 0 refusal arms `k0-refuse-period` (period ≠ 4, 12),
  `k0-refuse-dorder` (d = 3), `k0-refuse-dmean` (integrated with mean ≠ 0) and
  `k0-refuse-nanmean` (non-finite mean), recorded as the exact `THROW` message
  line the driver emits.

Inputs are the frozen item-29 `../diffuse-steady-state-probe/` corpus
(`cases.csv` + `series/`); no new series were minted.

Negative controls (measured at landing, see the M16 ledger entry): a one-ulp
edit of the `k0-nearboundary-seasonal` `SUM ssq` line makes the ctest assert
fire (`mismatches == 0`), and the HEAD library itself passes the panel —
pinning the claim in both directions.
