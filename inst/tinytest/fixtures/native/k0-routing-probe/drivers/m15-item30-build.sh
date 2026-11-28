#!/usr/bin/env bash
# Item-30 probe: configure+build HEAD Release libkalos.a against the pinned
# OpenBLAS provider (same practice as item-29/M12: single-threaded lane).
set -euo pipefail
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m15-item30
src=/data/home/kowa/seasonaladjustment/kalos
export LD_LIBRARY_PATH=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib
cmake -S "$src" -B "$root/build-head" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DKALOS_REVISION=1e2e3c8 \
  -DOpenBLAS_ROOT=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas \
  > "$root/build-head.cmake.log" 2>&1
env nice -n 10 cmake --build "$root/build-head" --target kalos -j6 > "$root/build-head.build.log" 2>&1
echo "BUILD-DONE head $(cd "$src" && git rev-parse HEAD)"
