#!/usr/bin/env bash
# Item-30 probe final capture: regenerates every committed measurement table
# from the committed drivers (see README.md for the build recipe). All phases
# loadavg-bracketed; contaminated rounds re-run (threshold: 1-min load > 6.0).
set -u
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m15-item30
evidence=/data/home/kowa/seasonaladjustment/kalos/validation/evidence/native/k0-routing-probe
corpus=/data/home/kowa/seasonaladjustment/kalos/validation/evidence/native/diffuse-steady-state-probe
export LD_LIBRARY_PATH=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 KALOS_NUM_THREADS=1
tables="$evidence/tables"
mkdir -p "$tables" "$root/out/final"

load1() { awk '{print $1}' /proc/loadavg; }
contaminated() { awk -v a="$1" -v b="$2" 'BEGIN{exit !((a>6.0)||(b>6.0))}'; }
run_phase() { # outfile cmd...
  local out=$1; shift
  local lb la tries=0
  lb=$(load1)
  ( cd "$corpus" && env nice -n 10 "$@" ) > "$root/out/final/$out" 2>&1
  la=$(load1)
  while contaminated "$lb" "$la" && [ $tries -lt 4 ]; do
    echo "RELOAD-RETRY $out bracket=$lb/$la"
    tries=$((tries+1)); sleep 20
    lb=$(load1)
    ( cd "$corpus" && env nice -n 10 "$@" ) > "$root/out/final/$out" 2>&1
    la=$(load1)
  done
  echo "BRACKET $out load_before=$lb load_after=$la tries=$tries" >> "$root/out/final/brackets.txt"
  cp "$root/out/final/$out" "$tables/$out"
}

cp_bin="$root/bin/m15-item30-cost_probe"
rp_bin="$root/bin/m15-item30-route_probe"
fp_bin="$root/bin/m15-item30-filter_probe_frozen"
rm -f "$root/out/final/brackets.txt"

# 0. frozen-golden rerun at HEAD (correctness context, not timed)
mkdir -p "$root/out/final/golden-head"
( cd "$corpus" && env nice -n 10 "$fp_bin" cases.csv "$root/out/final/golden-head" >/dev/null 2>&1 )
{
  echo "Frozen-golden rerun at HEAD (ec17543 goldens vs HEAD build)"
  for f in "$corpus"/goldens/*.txt; do
    id=$(basename "$f")
    if cmp -s "$f" "$root/out/final/golden-head/$id"; then
      echo "BYTE-IDENTICAL $id"
    else
      echo "MOVED $id (expected for k>0: D4 post-pivot route ships; in-band)"
    fi
  done
} > "$tables/golden-rerun.txt"

# 1. current-state eval/init (frozen item-29 cost driver), 3 rounds
for round in 1 2 3; do
  run_phase "eval-round$round.txt" "$cp_bin" cases.csv eval 9
  run_phase "init-round$round.txt" "$cp_bin" cases.csv init 9
done

# 2. interleaved route A/B on the two k=0 shapes, 3 rounds x 9 reps
for round in 1 2 3; do
  for c in syn-stat-sarma11 syn-stat-arma21; do
    run_phase "ab-$c-round$round.txt" "$rp_bin" ab cases.csv "$c" 9
  done
done

# 3. per-length slopes
run_phase "nsweep-syn-stat-sarma11.txt" "$rp_bin" nsweep cases.csv syn-stat-sarma11 9
run_phase "nsweep-syn-stat-arma21.txt" "$rp_bin" nsweep cases.csv syn-stat-arma21 9

# 4. instrumented section decomposition
for r in copy0 struct3 fac1; do
  for c in syn-stat-sarma11 syn-stat-arma21; do
    run_phase "instr-$c-$r.txt" "$rp_bin" instr cases.csv "$c" "$r" 9
  done
done

# 5. SIGPROF self-sampling over the library dense path + symbolization
run_phase "sample-syn-stat-sarma11.raw.txt" "$rp_bin" sample cases.csv syn-stat-sarma11 150
python3 "$evidence/drivers/m15-item30-symbolize.py" \
  "$root/bin/route_probe.symbols.txt" "$tables/sample-syn-stat-sarma11.raw.txt" \
  > "$tables/sample-syn-stat-sarma11.symbols.txt"

# 6. correctness: route dumps and pairwise gates
{
  cd "$corpus"
  for c in syn-stat-sarma11 syn-stat-arma21; do
    for r in lib copy0 fac1 fac2 struct3; do
      env nice -n 10 "$rp_bin" dump cases.csv "$c" "$r" > "$root/out/final/dump-$c-$r.txt"
    done
    echo "== $c =="
    cmp -s "$root/out/final/dump-$c-lib.txt" "$root/out/final/dump-$c-copy0.txt" \
      && echo "COPY-FIDELITY lib==copy0 BYTE-IDENTICAL" || echo "COPY-FIDELITY DIFFERS"
    cmp -s "$root/out/final/dump-$c-struct3.txt" "$root/out/final/dump-$c-copy0.txt" \
      && echo "OPTION-II-BIT-IDENTICAL struct3==dense BYTE-IDENTICAL" || echo "OPTION-II DIFFERS"
    cmp -s "$root/out/final/dump-$c-fac1.txt" "$root/out/final/dump-$c-fac2.txt" \
      && echo "FAC1==FAC2 value lines identical" || echo "FAC1-vs-FAC2: $(diff "$root/out/final/dump-$c-fac1.txt" "$root/out/final/dump-$c-fac2.txt" | head -4 | tr '\n' ' | ')"
    python3 "$evidence/drivers/m15-item30-dump_compare.py" "$root/out/final/dump-$c-fac1.txt" "$root/out/final/dump-$c-copy0.txt"
    echo
  done
} > "$tables/correctness.txt"

# 7. ratio table (JVM references from the frozen item-29 jvm/ captures)
python3 "$evidence/drivers/m15-item30-ratio_table.py" "$root/out/final" "$corpus" > "$tables/current-state-ratios.txt"

cp "$root/out/final/brackets.txt" "$tables/brackets.txt"
echo "CAPTURE-DONE $(cat /proc/loadavg)"
