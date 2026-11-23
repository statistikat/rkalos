#!/usr/bin/env python3
"""Rebuild cases.txt for the 23-F.4 tramo_fit probe (kept next to the corpus
for regeneration). The sc_*/mr_* value lines are copied verbatim from
../tramo-outlier-rounds/cases.txt (same series, KER-form spec headers), the
airline series comes from R datasets::AirPassengers (airline.txt, %.17g),
and the two boundary series are deterministic closed forms:

  nearur: 120 obs, y_t = 0.995*y_{t-1} + 10*u_t + 6*sin(2*pi*(t mod 12)/12),
          u_t = 2*LCG_t - 1 with the MMIX 64-bit LCG seeded at
          0x00000246FCF39C1D (x = (6364136223846793005*x +
          1442695040888963407) mod 2**64, u = (x >> 11)/2**53) - a
          near-unit-root AR coefficient with a strong monthly cycle.
  const:  96 observations at exactly 100.0 (the pinned engine aborts on it
          at every TR template; the row is the capture-level failure
          witness and the golden generator skips it by case id).
"""
import math
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

# ---- near-unit-root boundary series (deterministic MMIX LCG, docstring) ----
x = 0x00000246FCF39C1D
y = 0.0
nearur = []
for t in range(120):
    x = (6364136223846793005 * x + 1442695040888963407) % (1 << 64)
    u = (x >> 11) / float(1 << 53) * 2.0 - 1.0
    y = 0.995 * y + 10.0 * u + 6.0 * math.sin(2.0 * math.pi * (t % 12) / 12.0)
    nearur.append("%.17g" % y)
const = ["%.17g" % 100.0 for _ in range(96)]

CASES = [
    # families and transform branches (airline)
    ("fit_tr0_air", "TR0", 12, 1980, 1, None),
    ("fit_tr1_air", "TR1", 12, 1980, 1, None),
    ("fit_tr2_air", "TR2", 12, 1980, 1, None),            # calendar REFUSED
    ("fit_tr3_air", "K", 12, 1980, 1, None),
    ("fit_accept_air", "K+accept", 12, 1980, 1, None),
    ("fit_mean_air", "K+mean", 12, 1980, 1, None),
    ("fit_fct_air", "K fct=2.0", 12, 1980, 1, None),
    ("fit_fsnone", "TR0-o+fs", 12, 1980, 1, None),        # fully specified
    ("fit_tr1_ao", "TR1", 12, 1990, 1, "sc_ao"),
    # planted-outlier families (the 23-F.2 series)
    ("fit_noout", "K va=calc", 12, 1990, 1, "sc_noout"),
    ("fit_ao", "K va=calc", 12, 1990, 1, "sc_ao"),
    ("fit_ls", "K va=calc", 12, 1990, 1, "sc_ls"),
    ("fit_multi", "K va=calc", 12, 1990, 1, "sc_multi"),
    ("fit_so", "K+noso va=calc", 12, 1990, 1, "sc_so"),
    ("fit_spanout", "K va=calc", 12, 1990, 1, "sc_spanout"),
    ("fit_gs2", "K va=calc", 12, 1990, 1, "sc_gs2"),
    ("fit_s30ml", "K+noso+ml va=calc", 12, 1990, 1, "sc_s30ml"),
    ("fit_aoonly", "K+nols+notc+noso va=calc", 12, 1990, 1, "sc_aoonly"),
    ("fit_cv4", "K va=4.0", 12, 1990, 1, "sc_cv4"),
    ("fit_shortfail", "K va=calc", 12, 1990, 1, "sc_shortfail"),
    ("fit_mrstack", "K va=calc", 12, 1990, 1, "mr_stack"),
    ("fit_mrls", "K+noso va=calc", 12, 1990, 1, "mr_lsstart"),
    ("fit_missing", "K va=calc", 12, 1990, 1, "sc_missing"),  # NaN REFUSED
    # regression/mean boundaries
    ("fit_nomean", "TR1-o+nomean", 12, 1990, 1, "sc_noout"),
    ("fit_pure", "TR3-o+nomean", 12, 1980, 1, None),
    # series-shape boundaries
    ("fit_nearur", "K va=calc", 12, 1990, 1, "nearur"),
    ("fit_const", "K", 12, 1990, 1, "const"),             # crash witness
    # refusal witnesses (the golden expects refused:, no stream)
    ("fit_q_air", "K", 4, 1980, 1, None),                 # D11
    ("fit_trfull_air", "TRfull", 12, 1980, 1, None),      # calendar + FTEST
    ("fit_td_air", "K+td", 12, 1980, 1, None),            # automatic TD
    ("fit_aic_air", "K amth=AIC", 12, 1980, 1, None),     # I-10
    ("fit_bic_air", "K amth=BIC", 12, 1980, 1, None),     # I-10
    ("fit_wald_air", "K amth=WALD", 12, 1980, 1, None),   # I-10
    ("fit_oc_air", "K+oc", 12, 1980, 1, None),            # +oc
    ("fit_easter_air", "K+easter", 12, 1980, 1, None),    # easter test
]

EXTRA = {"nearur": nearur, "const": const}

with open(OUT, "w") as f:
    for cid, spec, freq, y0, m0, src in CASES:
        v = air if src is None else (EXTRA[src] if src in EXTRA else vals[src])
        f.write("KER %s %s %d %d %d\n" % (cid, spec, freq, y0, m0))
        f.write("values " + " ".join(v) + "\n\n")
print("wrote", OUT, len(CASES), "cases")
