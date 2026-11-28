#!/usr/bin/env bash
# Item-30 probe: compile the route probe (HEAD filter copy + driver) and the
# reused item-29 drivers against build-head/libkalos.a (same link line as the
# item-29/M12 drivers, plus -no-pie -g for the self-sampler symbolization).
set -euo pipefail
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m15-item30
kalos=/data/home/kowa/seasonaladjustment/kalos
blas=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib/libopenblas.so
obinc=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/include
build="$root/build-head"
mkdir -p "$root/bin"
cd "$root"
env nice -n 10 /usr/bin/c++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off -fopenmp -no-pie -g \
  -I "$kalos/include" -I "$build/generated" -I "$root" -isystem "$obinc" \
  m15-item30-route_probe.cpp m15-item30-filter_copy.cpp \
  -o "$root/bin/m15-item30-route_probe" \
  "$build/libkalos.a" "$blas" /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so \
  -lpthread > "$root/bin/route_probe.build.log" 2>&1 \
  || { echo FAILED route_probe; cat "$root/bin/route_probe.build.log"; exit 1; }
echo BUILT route_probe
env nice -n 10 /usr/bin/c++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off -fopenmp \
  -I "$kalos/include" -I "$build/generated" -isystem "$obinc" \
  m15-item30-cost_probe.cpp -o "$root/bin/m15-item30-cost_probe" \
  "$build/libkalos.a" "$blas" /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so \
  -lpthread > "$root/bin/cost_probe.build.log" 2>&1 \
  || { echo FAILED cost_probe; cat "$root/bin/cost_probe.build.log"; exit 1; }
echo BUILT cost_probe
env nice -n 10 /usr/bin/c++ -std=c++17 -O2 -DNDEBUG -ffp-contract=off -fopenmp \
  -I "$kalos/include" -I "$build/generated" -isystem "$obinc" \
  m15-item30-filter_probe_frozen.cpp -o "$root/bin/m15-item30-filter_probe_frozen" \
  "$build/libkalos.a" "$blas" /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so \
  -lpthread > "$root/bin/filter_probe_frozen.build.log" 2>&1 \
  || { echo FAILED filter_probe_frozen; cat "$root/bin/filter_probe_frozen.build.log"; exit 1; }
echo BUILT filter_probe_frozen
nm -C --defined-only "$root/bin/m15-item30-route_probe" > "$root/bin/route_probe.symbols.txt"
echo "drivers ready"
