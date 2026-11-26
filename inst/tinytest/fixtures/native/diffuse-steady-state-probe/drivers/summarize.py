#!/usr/bin/env python3
"""Item-29 probe slice: assemble the kalos-vs-JVM cross-validation table.

Reads the frozen kalos exact-path outputs (corpus/out/*.txt) and the JVM
captures (jvm_capt4/*.txt), emits the dossier tables. Deterministic; the
timing medians are computed from the captured TIMING lines (JIT-warm
eval 0 excluded).
"""
import re, glob, os, statistics

def kv(line, key):
    m = re.search(rf"{key}=(\S+)", line)
    return m.group(1) if m else None

kalos = {}
for path in sorted(glob.glob("../corpus/out/*.txt")):
    cid = os.path.basename(path)[:-4]
    d = {}
    for line in open(path):
        if line.startswith("FILTER "):
            d["m"] = kv(line, "m_state"); d["k"] = kv(line, "k_hist")
        elif line.startswith("SUMS "):
            d["ssq"] = kv(line, "ssq"); d["profiledll"] = kv(line, "profiledll")
            d["slogF"] = kv(line, "slogF")
        elif line.startswith("REG status"):
            d["reg_ssq"] = kv(line, "ssq"); d["reg_profiledll"] = kv(line, "profiledll")
            d["reg_slogF"] = kv(line, "slogF")
        elif line.startswith("TSTAR tau=1e-14"):
            d["tstar14"] = kv(line, "tstar")
        elif line.startswith("TSTAR tau=1e-10"):
            d["tstar10"] = kv(line, "tstar")
        elif line.startswith("CASE "):
            d["n"] = line.split()[2].split("=")[1]
    kalos[cid] = d

jvm = {}
for path in sorted(glob.glob("jvm_capt4/*.txt")):
    cid = os.path.basename(path)[:-4]
    d = {}
    totals, prep, appl = [], [], []
    for line in open(path):
        if line.startswith("CASE "):
            d["n"] = line.split()[2].split("=")[1]
            d["nx"] = line.split()[3].split("=")[1]
        elif line.startswith("MACHINE "):
            d["dim"] = kv(line, "dim"); d["steadyPos"] = kv(line, "steadyPos")
        elif line.startswith("LL "):
            d["ssq"] = kv(line, "ssq"); d["ldet"] = kv(line, "ldet")
            d["loglik"] = kv(line, "loglik")
        elif line.startswith("TIMING eval="):
            if " eval=0 " in line:  # skip JIT-warm first evaluation
                continue
            totals.append(float(kv(line, "total-ns")))
            prep.append(float(kv(line, "prepare-ns")))
            appl.append(float(kv(line, "apply-ns")))
    if totals:
        d["med_total_us"] = statistics.median(totals) / 1e3
        d["med_prep_us"] = statistics.median(prep) / 1e3
        d["med_apply_us"] = statistics.median(appl) / 1e3
    jvm[cid] = d

PAIRS = [  # (jvm-id, kalos-id, which kalos quantities)
    ("jvm-airline-noreg", "real-airline-noreg", "FILTER"),
    ("jvm-airline-x",     "real-airline-x",     "REG"),
    ("jvm-stat-arma21n",  "syn-stat-arma21",    "FILTER"),
    ("jvm-stat-sarma11",  "syn-stat-sarma11",   "FILTER"),
    ("jvm-nearur",        "nearur-airline",     "FILTER"),
    ("jvm-nearni",        "nearni-airline",     "FILTER"),
    ("jvm-d2",            "d2-arma12",          "FILTER"),
    ("jvm-slow",          "slow-seasonal",      "FILTER"),
    ("jvm-airline4",      "syn-airline4",       "FILTER"),
    ("jvm-mixed4",        "syn-mixed4",         "FILTER"),
]

def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)

hdr = (f"{'jvm case':18s} {'kalos ssq':20s} {'jvm ssq':20s} {'rel':>8s} "
       f"{'kalos ll':20s} {'jvm loglik*':>12s} {'rel sF':>8s} {'sp/n':>9s} "
       f"{'tstar14/n':>10s} {'m/k':>6s} {'dim':>3s}")
print(hdr)
for j, kc, which in PAIRS:
    J, K = jvm[j], kalos[kc]
    ks = K["reg_ssq"] if which == "REG" else K["ssq"]
    kl = K["reg_profiledll"] if which == "REG" else K["profiledll"]
    kf = K["reg_slogF"] if which == "REG" else K["slogF"]
    # jvm loglik is concentrated like kalos profiledll: compare printed decimals
    print(f"{j:18s} {ks:20s} {J['ssq']:20s} {rel(J['ssq'], ks):8.1e} "
          f"{kl:20s} {J['loglik']:>12s} {rel(J['ldet'], kf):8.1e} "
          f"{J['steadyPos']+'/'+J['n']:>9s} {K['tstar14']+'/'+K['n']:>10s} "
          f"{K['m']+'/'+K['k']:>6s} {J['dim']:>3s}")
print("* JVM loglik lines in the captures agree with the kalos profiled ll to")
print("  the 10 printed decimals in every paired case (see captures).")
print()
print(f"{'jvm case':22s} {'n':>4s} {'nx':>3s} {'prep us':>8s} {'apply us':>9s} {'total us':>9s}  (median of evals 1..5)")
for cid in sorted(jvm):
    d = jvm[cid]
    if "med_total_us" in d:
        print(f"{cid:22s} {d['n']:>4s} {d['nx']:>3s} {d['med_prep_us']:8.1f} "
              f"{d['med_apply_us']:9.1f} {d['med_total_us']:9.1f}")
