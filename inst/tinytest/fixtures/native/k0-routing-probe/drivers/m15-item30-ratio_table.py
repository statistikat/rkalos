#!/usr/bin/env python3
"""item-30 probe: current-state ratio table.
kalos medians: EVAL phase medians from the three frozen-cost rounds (item-29
convention: median of per-round medians). JVM reference: median of warm
eval=1..5 total-ns in the frozen jvm/ captures (item-29 convention)."""
import glob, os, re, statistics, sys

out = sys.argv[1]  # timing out dir
corpus = sys.argv[2]

kalos = {}
per_case = {}
for path in sorted(glob.glob(os.path.join(out, "eval-*.txt"))):
    for line in open(path):
        if line.startswith("EVAL "):
            parts = line.split()
            per_case.setdefault(parts[1], []).append(float(parts[3].split("=")[1]))
for case, values in per_case.items():
    kalos[case] = statistics.median(values)

jvm = {}
for path in sorted(glob.glob(os.path.join(corpus, "jvm", "jvm-*.txt"))):
    base = os.path.basename(path)[len("jvm-"):-len(".txt")]
    if base.endswith("n"):  # jvm-stat-arma21n is the calibrated twin
        base = base.rstrip("n") if (path[:-4] + "n.txt").endswith("arma21n.txt") else base
    totals = []
    for line in open(path):
        m = re.search(r"eval=(\d+) .*total-ns=(\d+)", line)
        if m and 1 <= int(m.group(1)) <= 5:
            totals.append(int(m.group(2)) / 1e6)
    if totals:
        jvm[base] = statistics.median(totals)

alias = {"real-airline-noreg": "airline-noreg", "real-airline-x": "airline-x",
         "real-airline-gaps-x": "airline-gaps-x", "syn-airline12": "airline-noreg",
         "syn-stat-arma21": "stat-arma21", "syn-stat-sarma11": "stat-sarma11",
         "syn-airline4": "airline4", "syn-mixed4": "mixed4", "d2-arma12": "d2",
         "nearur-airline": "nearur", "nearni-airline": "nearni",
         "slow-seasonal": "slow", "leadgap-airline-x": "airline-gaps-x",
         "real-airline-redundant-x": None, "real-airline-zerocol-x": None}

print(f"{'case':28s} {'kalos_head_ms':>13s} {'jvm_ms':>8s} {'ratio':>7s}")
for case in sorted(kalos):
    j = alias.get(case, case)
    line = f"{case:28s} {kalos[case]:13.3f}"
    if j and j in jvm:
        print(f"{line} {jvm[j]:8.3f} {kalos[case] / jvm[j]:7.1f}")
    else:
        print(f"{line} {'-':>8s} {'-':>7s}")
print("\njvm warm eval=1..5 medians (ms):")
for case, v in sorted(jvm.items()):
    print(f"  {case:20s} {v:.4f}")
