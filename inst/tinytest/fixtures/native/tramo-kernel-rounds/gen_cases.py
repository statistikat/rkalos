#!/usr/bin/env python3
"""Rebuild cases.txt from the pinned sources (kept next to the corpus for
regeneration): the sc_*/mr_* value lines are copied verbatim from
../tramo-outlier-rounds/cases.txt (same series, KER-form spec headers) and the
airline series from R datasets::AirPassengers (airline.txt, %.17g)."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "tramo-outlier-rounds", "cases.txt")
AIR = os.path.join(HERE, "airline.txt")
OUT = os.path.join(HERE, "cases.txt")

# ---- series from the 23-F.2 probe ----
vals = {}
cur = None
with open(SRC) as f:
    for ln in f:
        t = ln.strip()
        if t.startswith("SCAN "):
            cur = t.split()[1]
        elif t.startswith("values ") and cur:
            vals.setdefault(cur, []).extend(t.split()[1:])

# ---- airline (R datasets::AirPassengers, 144 obs) ----
with open(AIR) as f:
    air = ["%.17g" % float(x) for x in f.read().split()]
assert len(air) == 144, len(air)

CASES = [
    ("kr_noout", "K va=calc", 12, 1990, 1, "sc_noout"),
    ("kr_ao", "K va=calc", 12, 1990, 1, "sc_ao"),
    ("kr_ls", "K va=calc", 12, 1990, 1, "sc_ls"),
    ("kr_multi", "K va=calc", 12, 1990, 1, "sc_multi"),
    ("kr_scso", "K+noso va=calc", 12, 1990, 1, "sc_so"),
    ("kr_s30", "K+noso va=calc", 12, 1990, 1, "sc_s30"),
    ("kr_s30ml", "K+noso+ml va=calc", 12, 1990, 1, "sc_s30ml"),
    ("kr_spanout", "K va=calc", 12, 1990, 1, "sc_spanout"),
    ("kr_spannone", "K va=calc", 12, 1990, 1, "sc_spannone"),
    ("kr_spanfull", "K va=calc", 12, 1990, 1, "sc_spanfull"),
    ("kr_aoonly", "K+nols+notc+noso va=calc", 12, 1990, 1, "sc_aoonly"),
    ("kr_tc05", "K va=calc tcrate=0.5", 12, 1990, 1, "sc_tc05"),
    ("kr_tc09", "K va=calc tcrate=0.9", 12, 1990, 1, "sc_tc09"),
    ("kr_gs2", "K va=calc", 12, 1990, 1, "sc_gs2"),
    ("kr_shortfail", "K va=calc", 12, 1990, 1, "sc_shortfail"),
    ("kr_missing", "K va=calc", 12, 1990, 1, "sc_missing"),
    ("kr_cv4", "K va=4.0", 12, 1990, 1, "sc_cv4"),
    ("kr_maxout2", "K va=calc", 12, 1990, 1, "sc_maxout2"),
    ("kr_ml", "K+ml va=calc", 12, 1990, 1, "sc_ml"),
    ("kr_air", "K", 12, 1980, 1, None),
    ("kr_accept_air", "K+accept", 12, 1980, 1, None),
    ("kr_mean_air", "K+mean", 12, 1980, 1, None),
    ("kr_oc_air", "K+oc", 12, 1980, 1, None),
    ("kr_fct_air", "K fct=2.0", 12, 1980, 1, None),
    ("kr_tr1_air", "TR1", 12, 1980, 1, None),
    ("kr_tr1_scso", "TR1+noso", 12, 1990, 1, "sc_so"),
    ("kr_fs_none", "K+fs", 12, 1980, 1, None),
    ("kr_fs_log", "K+fs fct=1.0", 12, 1980, 1, None),
    ("kr_td_auto", "K+td", 12, 1980, 1, None),
    ("kr_q_air", "K", 4, 1980, 1, None),
    ("kr_mr_stack", "K va=calc", 12, 1990, 1, "mr_stack"),
    ("kr_mr_sostack", "K va=calc", 12, 1990, 1, "mr_sostack"),
    ("kr_mr_outcap", "K va=calc", 12, 1990, 1, "mr_outcap"),
    ("kr_mr_roundcapx", "K va=calc", 12, 1990, 1, "mr_roundcapx"),
    ("kr_mr_lsstart", "K+noso va=calc", 12, 1990, 1, "mr_lsstart"),
]

with open(OUT, "w") as f:
    for cid, spec, freq, y0, m0, src in CASES:
        v = air if src is None else vals[src]
        f.write("KER %s %s %d %d %d\n" % (cid, spec, freq, y0, m0))
        f.write("values " + " ".join(v) + "\n\n")
print("wrote", OUT, len(CASES), "cases")
