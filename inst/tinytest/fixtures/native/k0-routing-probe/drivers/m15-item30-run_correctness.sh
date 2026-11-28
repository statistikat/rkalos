#!/usr/bin/env bash
# Item-30 probe: correctness gates. For each k=0 case, dump the full named
# likelihood through every route and diff pairwise:
#   lib vs copy0   -> filter-copy fidelity (must be byte-identical)
#   struct3 vs copy0 -> option-(ii) bit-identity claim (must be byte-identical)
#   fac1 vs fac2   -> frame bookkeeping at k=0 (values identical, observability differs)
#   fac1 vs copy0  -> factorized-route deviation (report, band vs decision audit)
set -eu
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m15-item30
corpus=/data/home/kowa/seasonaladjustment/kalos/validation/evidence/native/diffuse-steady-state-probe
export LD_LIBRARY_PATH=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib
out="$root/out/dumps"
mkdir -p "$out"
cd "$corpus"
for case_id in syn-stat-sarma11 syn-stat-arma21; do
  for route in lib copy0 fac1 fac2 struct3; do
    env nice -n 10 "$root/bin/m15-item30-route_probe" dump cases.csv "$case_id" "$route" \
      > "$out/$case_id.$route.txt" 2> "$out/$case_id.$route.stderr" || echo "RUN-FAILED $case_id $route"
  done
  echo "== $case_id: lib vs copy0 (copy fidelity) =="
  cmp "$out/$case_id.lib.txt" "$out/$case_id.copy0.txt" && echo IDENTICAL-BYTES || true
  echo "== $case_id: struct3 vs copy0 (option-(ii) bit identity) =="
  cmp "$out/$case_id.struct3.txt" "$out/$case_id.copy0.txt" && echo IDENTICAL-BYTES || true
  echo "== $case_id: fac1 vs fac2 =="
  cmp "$out/$case_id.fac1.txt" "$out/$case_id.fac2.txt" && echo IDENTICAL-BYTES || true
done
echo done
