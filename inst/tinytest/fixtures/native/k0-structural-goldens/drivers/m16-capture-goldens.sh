#!/usr/bin/env bash
# Item-30 adoption: capture the full byte-identity oracle panel for a given
# goldens binary into a file, bracketed by /proc/loadavg (capture lane only).
# usage: m16-capture-goldens.sh <bin-name> <out-file>
set -euo pipefail
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m16-item30adopt
corpus=/data/home/kowa/seasonaladjustment/kalos/validation/evidence/native/diffuse-steady-state-probe
export LD_LIBRARY_PATH=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib
bin="$1"; out="$2"
echo "CAPTURE-BEFORE $(cat /proc/loadavg)"
if nice -n 10 "$bin" "$corpus" > "$out" 2> "$out.stderr"; then ec=0; else ec=$?; fi
echo "CAPTURE-AFTER $(cat /proc/loadavg) exit=$ec"
[ $ec -eq 0 ] && echo "CAPTURE-OK $out $(wc -l < "$out") lines $(grep -c '^ARM' "$out") arms"
