#!/usr/bin/env bash
# Item-30 adoption: configure+build a Release libkalos.a tree from a given
# source state against the pinned OpenBLAS provider (m15 probe recipe).
# usage: m16-build.sh <build-dir-name> [target...]   (default target: kalos)
set -euo pipefail
root=/data/home/kowa/seasonaladjustment/tmp_opencode/m16-item30adopt
src=/data/home/kowa/seasonaladjustment/kalos
export LD_LIBRARY_PATH=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas/lib
build="$root/$1"; shift || true
targets=("$@"); [ ${#targets[@]} -gt 0 ] || targets=(kalos)
nice -n 10 cmake -S "$src" -B "$build" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DOpenBLAS_ROOT=/data/home/kowa/seasonaladjustment/tmp_opencode/openblas \
  > "$build.cmake.log" 2>&1
nice -n 10 cmake --build "$build" --target "${targets[@]}" -j6 > "$build.build.log" 2>&1
echo "BUILD-DONE $build $(cd "$src" && git rev-parse --short HEAD) targets=${targets[*]}"
