# calendar-detection corpus (roadmap 21-I stage 1: CEDM ATD step)

Frozen evidence supporting stage 1 of item 21-I
([`docs/item21-calendar-detection-probe-2026-11-28.md`](../../../../docs/item21-calendar-detection-probe-2026-11-28.md)):
the CEDM automatic-trading-days (ATD) firing, the LoP preadjust-vs-column
fork and the per-round re-fires, replayed bit-for-bit against the JVM. This
is a replay-golden corpus: `tests/test_x13_automatic_td.cpp` consumes the
generated pin table (via `tools/gen_x13_automatic_td_golden.py`) and the
production entry point `kalos::x13_automatic_fit` with
`automatic_trading_days=true` is qualified here. Staging follows the spec's
option (iv): the TD battery lands now; the easter battery and the refusal-
ladder extension are stages 2-3 and stay deferred (calendar/easter/user
flags still refuse verbatim at the narrowed gate).

Provenance

- kalos build: worktree HEAD `4b06d78` + the 21-I stage-1 changes
  (`src/x13_automatic.cpp`, `src/transformation_decision.cpp`,
  `include/kalos/x13_automatic.hpp`, `include/kalos/
  transformation_decision.hpp`, the regtest-caller fix in
  `src/regression_variable_tests.cpp`'s caller path), Release
  (`-O3 -DNDEBUG -ffp-contract=off`), single thread
  (`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=KALOS_NUM_THREADS=1`), OpenBLAS
  0.3.34 built in scratch, GCC 13.3.0, Intel Xeon Gold 6342 (host xlwp0139).
- JVM side: OpenJDK 21 with the diagnostic intrinsics switch
  (`-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp
  -XX:DisableIntrinsic=_dlog`, the qualified bridge domain), jars from the
  frozen jarcache; driver `drivers/run_td1_battery.sh`
  (`Td1Probe.java`/`Td1Corpora.java`, observer-style probes only — no
  decision rule is re-implemented in the probe).
- Corpora (`drivers/td1-corpora.txt` + the kAir array in
  `drivers/td1_replay.cpp`): `clean` (96 monthly, synthetic, degenerate
  variance), `tdeff` (156 monthly with a strong trading-days effect), `air`
  (the 144-point Box-Jenkins airline series). Battery templates
  `td2noe` (td2, auto-adjust), `td7noe` (td7, auto-adjust), `td2noadj`
  (td2, no auto-adjust); `td2lev` is captured as probe context but is not
  a replay target (the level-request variant had no native scope).

Contents

- `captures/` — the frozen first JVM capture of each battery case
  ({4 templates} x {clean,tdeff,air}: `.out` + `.err`) plus the three
  `rsa3-*` no-TD baseline captures from the m16-p7 runA battery (the
  zero-movement anchors).
- `javap/` — the JVM bytecode dumps the stage-1 decode cites (kernel calc
  order, refine/AMI strip, RegressionVariablesTest/DerivedRegressionTest/
  TRegressionTest, the ODM pair, ModelController, LogLevelModule/
  LoPTransformation, the CEDM builder, ConcentratedLikelihood,
  SingleOutlierDetector).
- `drivers/` — the capture battery, probes and corpora; `td1_replay.cpp`
  (the Tier-1 replay harness: chain records, ATD arm tables, scan tables);
  `td1_replay_head.cpp` (the HEAD-source baseline harness used to attribute
  the clean-corpus gap); `make_tables.py` (regenerates `tables/` from
  `captures/` + `replay/`).
- `replay/` — the native Tier-1 replay outputs for the nine replayed cases
  plus the two `rsa3` baselines, and `replay-head-rsa3-clean.out`: the
  HEAD-source binary reproducing the clean-corpus outcome, which proves the
  clean divergence pre-dates stage 1.
- `tables/decision-bits.csv` — per case: JVM vs native log/lop/orders/
  outlier bits and the ll/aicc ulp distances (classification column:
  `exact` <= 4 ulp, `ulp` <= 64, else `divergent`).
- `tables/atd-firings.csv` — every ATD firing, JVM events paired with
  native firings by ordinal (explicit `native-only`/`jvm-only` rows; the
  `best_match` column is "yes" on all paired rows).
- `tables/tier3-timing.txt` — loadavg-bracketed single-thread runtimes and
  the two-run byte-identity determinism lane.
- `SHA256SUMS.csv` — checksums of every file in this directory.

Capture and replay protocol

- Double JVM capture (`run1/`, `run2/` by the battery script): all twelve
  `.out` files are byte-identical modulo the `timing ... elapsed_ms` line;
  all `.err` files are byte-identical (verified in the landing session).
- Native replay determinism: two consecutive single-thread runs per case
  are byte-identical (`tables/tier3-timing.txt`, `run2_match`), and the
  frozen `replay/` outputs reproduce byte-identically from the landing
  build (`FROZEN_MATCH`).

Qualification result (the pin policy lives in
`tests/test_x13_automatic_td.cpp` and the generator: decision bits exact,
ll/aicc pinned as native bit patterns, only decisions with raw AICC margin
>= 0.1 pinned — all margins are in the pin table)

- ATD arm selections: exact JVM agreement on every paired firing of all
  nine cases (and the same arms/AICC to ulps at the pre-round firing).
- Preadjust-vs-column fork: exact JVM agreement on all nine (the JVM
  `lopTransform` bit vs the native `auto_calendar->preadjusted`).
- Final state exact-or-ulp (<= 10 ulp): `td2noe/tdeff`, `td2noe/air`,
  `td7noe/tdeff`, `td2noadj/air`, and both `rsa3` baselines (zero
  movement). `td2noadj/tdeff` reproduces the JVM orders but not the final
  likelihood; `td7noe/air` additionally loses the JVM's outlier removal.
- Documented residuals (pinned as native baselines, not JVM claims):
  1. `clean` x3 — a pre-existing 21-H baseline gap on the degenerate
     corpus (ssq ~7.5e-27), not TD-related: the HEAD-source binary
     (`replay/replay-head-rsa3-clean.out`, no ATD) reproduces the same
     native outcome. Out of 21-I stage-1 scope.
  2. `td7noe/air`, `td2noadj/tdeff` — refine/AMI ODM-boundary
     state-machine residuals: the JVM's round-1 refine path strips state
     differently (visible as the JVM's extra/narrower mid-round ATD
     re-fires in `tables/atd-firings.csv`); every ATD decision bit still
     matches. Carried forward to the 21-I stage-2 decode.
  3. `td7noe/air` firing count: JVM 1, native 2 (the native round-1
     re-fire survives the state the JVM's loop consumed) — recorded as a
     `native-only` row.

Regeneration

- JVM side: `drivers/run_td1_battery.sh` (writes `run1/`, `run2/` beside
  the drivers; the captured `run1/` in `captures/` is canonical).
- Tables: `python3 drivers/make_tables.py` from this directory.
- Pin header: `python3 tools/gen_x13_automatic_td_golden.py` (repo root).
- Native replays: build `drivers/td1_replay.cpp` against the library
  (OpenBLAS link line as in the x13 consumers) and run
  `./td1_replay <td2noe|td7noe|td2noadj|rsa3> <clean|tdeff|air>
  drivers/td1-corpora.txt` with the single-thread environment above.
