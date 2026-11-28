#!/usr/bin/env bash
# Item-30 probe timing lane: all rounds loadavg-bracketed, warm-up first,
# single-thread provider, nice 10. Rounds with a contaminated bracket
# (1-min load > 6.0 on 8 cores, front or back) are re-run.
set -u
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m15-item30
corpus=/data/home/kowa/seasonaladjustment/kalos/validation/evidence/native/diffuse-steady-state-probe
export LD_LIBRARY_PATH=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 KALOS_NUM_THREADS=1
out="$root/out/timing"
mkdir -p "$out"

load1() { awk '{print $1}' /proc/loadavg; }
contaminated() { awk -v a="$1" -v b="$2" 'BEGIN{exit !((a>6.0)||(b>6.0))}'; }

run_phase() { # name round cmd...
  local name=$1 round=$2; shift 2
  local lb la tries=0
  lb=$(load1)
  ( cd "$corpus" && env nice -n 10 "$@" ) > "$out/$name-$round.txt" 2> "$out/$name-$round.err"
  la=$(load1)
  while contaminated "$lb" "$la" && [ $tries -lt 4 ]; do
    echo "RELOAD-RETRY $name round=$round bracket=$lb/$la" >> "$out/brackets.log"
    tries=$((tries+1)); sleep 20
    lb=$(load1)
    ( cd "$corpus" && env nice -n 10 "$@" ) > "$out/$name-$round.txt" 2> "$out/$name-$round.err"
    la=$(load1)
  done
  echo "BRACKET $name round=$round load_before=$lb load_after=$la tries=$tries" >> "$out/brackets.log"
}

cp_bin="$root/bin/m15-item30-cost_probe"
rp_bin="$root/bin/m15-item30-route_probe"

# warm-ups (not timed)
( cd "$corpus" && env nice -n 10 "$cp_bin" cases.csv eval 3 >/dev/null 2>&1 )
( cd "$corpus" && env nice -n 10 "$rp_bin" ab cases.csv syn-stat-sarma11 1 >/dev/null 2>&1 )
( cd "$corpus" && env nice -n 10 "$rp_bin" ab cases.csv syn-stat-arma21 1 >/dev/null 2>&1 )

for round in 1 2 3; do
  run_phase "eval-$round" "$round" "$cp_bin" cases.csv eval 9
  run_phase "init-$round" "$round" "$cp_bin" cases.csv init 9
  run_phase "ab-sarma11-$round" "$round" "$rp_bin" ab cases.csv syn-stat-sarma11 9
  run_phase "ab-arma21-$round" "$round" "$rp_bin" ab cases.csv syn-stat-arma21 9
done
run_phase "nsweep-sarma11" 1 "$rp_bin" nsweep cases.csv syn-stat-sarma11 5
run_phase "nsweep-arma21" 1 "$rp_bin" nsweep cases.csv syn-stat-arma21 5
run_phase "instr-sarma11-copy0" 1 "$rp_bin" instr cases.csv syn-stat-sarma11 copy0 9
run_phase "instr-sarma11-struct3" 1 "$rp_bin" instr cases.csv syn-stat-sarma11 struct3 9
run_phase "instr-sarma11-fac1" 1 "$rp_bin" instr cases.csv syn-stat-sarma11 fac1 9
run_phase "instr-arma21-copy0" 1 "$rp_bin" instr cases.csv syn-stat-arma21 copy0 9
run_phase "instr-arma21-struct3" 1 "$rp_bin" instr cases.csv syn-stat-arma21 struct3 9
echo "TIMING-LANE-DONE $(cat /proc/loadavg)"
