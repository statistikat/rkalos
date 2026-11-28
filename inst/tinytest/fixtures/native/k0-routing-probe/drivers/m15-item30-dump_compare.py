#!/usr/bin/env python3
"""item-30 probe: pairwise dump comparison.
Decision lines (STATUS/COUNTS/NUMERICS/RES/UKIND) must match bit-for-bit;
numeric lines report max relative deviation between named quantities."""
import sys, math

DECISION = ("STATUS", "COUNTS", "NUMERICS", "RES", "UKIND")
VALUE_PREFIXES = ("SUM", "PRED", "FVAR", "INNO", "STDI", "FSTATE",
                  "INITCOV", "INITDIFF", "FINALCOV", "FINALDIFF")

def load(path):
    decisions, values = [], []
    for line in open(path):
        line = line.rstrip("\n")
        fields = line.split(" ")
        if line.startswith(DECISION):
            decisions.append(line)
        elif fields[0] in VALUE_PREFIXES and (
                len(fields) == 2 or (fields[0] == "SUM" and len(fields) == 3)):
            values.append((fields[0] + (" " + fields[1] if fields[0] == "SUM" and len(fields) == 3 else ""),
                           float(fields[-1])))
        else:
            decisions.append(line)
    return decisions, values

def rel(a, b):
    if a == b:
        return 0.0
    if a == 0.0 or b == 0.0 or math.isnan(a) or math.isnan(b):
        return float("inf") if (math.isnan(a) != math.isnan(b) or a != b) else 0.0
    return abs(a - b) / max(abs(a), abs(b))

def main(pa, pb):
    da, va = load(pa)
    db, vb = load(pb)
    decision_match = da == db
    assert len(va) == len(vb), "value count mismatch"
    worst, worst_tag, worst_pair = 0.0, None, None
    worst_abs, worst_abs_tag, worst_abs_pair = 0.0, None, None
    nonzero = 0
    for (ta, a), (tb, b) in zip(va, vb):
        assert ta == tb
        d = abs(a - b)
        if d > worst_abs:
            worst_abs, worst_abs_tag, worst_abs_pair = d, ta, (a, b)
        if max(abs(a), abs(b)) < 1e-8:
            continue  # relative metric on rounding-floor entries is meaningless
        r = rel(a, b)
        if r > 0:
            nonzero += 1
        if r > worst:
            worst, worst_tag, worst_pair = r, ta, (a, b)
    print(f"A={pa}\nB={pb}")
    print(f"decision-lines-identical={decision_match} value-lines={len(va)} differing(|ref|>=1e-8)={nonzero}")
    print(f"max-rel-dev(|ref|>=1e-8)={worst:.3e} at {worst_tag} pair={worst_pair}")
    print(f"max-abs-dev(all entries)={worst_abs:.3e} at {worst_abs_tag} pair={worst_abs_pair}")
    if not decision_match:
        for x, y in zip(da, db):
            if x != y:
                print("DECISION-DIFF:\n  A:", x[:200], "\n  B:", y[:200])

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
