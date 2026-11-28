#!/usr/bin/env python3
"""m17-td1 stage-1: generate the frozen evidence tables from the captures
(run1/ JVM outputs) and native replays (replay/ outputs). Ordinal pairing is
explicit: JVM firings are tagged by event number, native firings by chain
stage; a count mismatch shows as a '-' row, never as a silent shift."""
import csv
import re
import struct
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
ROOT = _here.parent if _here.name == "drivers" else _here
CAP = ROOT / "captures"
REP = ROOT / "replay"
TAB = ROOT / "tables"
TAB.mkdir(exist_ok=True)

TPLS = ["td2noe", "td7noe", "td2noadj"]
KINDS = ["clean", "tdeff", "air"]


def ulp_diff(a: float, b: float) -> int:
    """Distance in ulps between two doubles (IEEE total-order trick)."""
    def order(x: float) -> int:
        i = struct.unpack("<q", struct.pack("<d", x))[0]
        return -(i & 0x7FFFFFFFFFFFFFFF) - 1 if i < 0 else i

    return abs(order(a) - order(b))


def jvm_facts(path: Path) -> dict:
    txt = path.read_text()
    out = {}

    def one(pat):
        m = re.search(pat, txt)
        return m.group(1) if m else ""

    out["log"] = one(r"fact \S+ log (\S+)")
    out["lop"] = one(r"lopTransform (\S+)")
    out["orders"] = one(r"fact \S+ orders (\d \d \d \d \d \d)")
    out["variables"] = one(r"fact \S+ variables (\d+)")
    out["ll"] = one(r"fact \S+ ll (\S+) adjll")
    out["aicc"] = one(r"adjll \S+ aicc (\S+) ssq")
    out["ssq"] = one(r"aicc \S+ ssq (\S+) npar")
    out["npar"] = one(r"ssq \S+ npar (\d+)")
    return out


def native_facts(path: Path) -> dict:
    out = {}
    txt = path.read_text()
    m = re.search(r"^fact log=(\d) mean=(\d) orders (\d \d \d \d \d \d)$",
                  txt, re.M)
    if m:
        out["log"] = "true" if m.group(1) == "1" else "false"
        out["mean"] = "true" if m.group(2) == "1" else "false"
        out["orders"] = m.group(3)
    m = re.search(r"^fact ll=(\S+) aicc=(\S+) eff=(\d+)", txt, re.M)
    if m:
        out["ll"], out["aicc"], out["eff"] = m.groups()
    m = re.search(r"^landed .*?(\btd\b)?\s*$", txt, re.M)
    out["lop"] = "td-column" if (m and m.group(1)) else "none"
    m = re.search(r"^chain final .*?n_outliers=(\d+)", txt, re.M)
    if m:
        out["variables"] = m.group(1)
    return out


def jvm_atd(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if "ATD names" not in line:
            continue
        evt = re.match(r"evt (\d+)", line).group(1)
        best = re.search(r"best=(\d+|-+)", line).group(1)
        aiccs = re.findall(r"aicc=([0-9.eE+-]+)", line)
        ssqs = re.findall(r"ssq=([0-9.eE+-]+)", line)
        names = re.search(r"names=\[([^\]]*)\]", line).group(1)
        rows.append((evt, best, names, ssqs, aiccs))
    return rows


def native_atd(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith("evt "):
            continue
        stage = line.split()[1]
        best = re.search(r"best=(\d+|-+)", line).group(1)
        aiccs = re.findall(r"aicc=([0-9.eE+-]+)", line)
        ssqs = re.findall(r"ssq=([0-9.eE+-]+)", line)
        names = re.search(r"names=\[([^\]]*)\]", line).group(1)
        rows.append((stage, best, names, ssqs, aiccs))
    return rows


# ---- decision bits + final facts ----
with open(TAB / "decision-bits.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["template", "corpus", "jvm_log", "jvm_lop", "jvm_orders",
                "jvm_variables", "native_log", "native_lop", "native_orders",
                "native_variables", "jvm_ll", "native_ll", "ll_ulp",
                "jvm_aicc", "native_aicc", "aicc_ulp", "final_match"])
    for tpl in TPLS:
        for kind in KINDS:
            j = jvm_facts(CAP / f"{tpl}-{kind}.out")
            n = native_facts(REP / f"replay-{tpl}-{kind}-native.out")
            ll_ulp = aicc_ulp = ""
            match = "state-divergent"
            try:
                ll_ulp = ulp_diff(float(j["ll"]), float(n["ll"]))
                aicc_ulp = ulp_diff(float(j["aicc"]), float(n["aicc"]))
                match = ("exact" if ll_ulp <= 4 else
                         "ulp" if ll_ulp <= 64 else "divergent")
            except (KeyError, ValueError):
                pass
            w.writerow([tpl, kind, j.get("log"), j.get("lop"),
                        j.get("orders"), j.get("variables"),
                        n.get("log"), n.get("lop"), n.get("orders"),
                        n.get("variables"), j.get("ll"), n.get("ll"),
                        ll_ulp, j.get("aicc"), n.get("aicc"), aicc_ulp,
                        match])

# ---- ATD firings, explicit ordinal alignment ----
with open(TAB / "atd-firings.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["template", "corpus", "ordinal", "jvm_evt", "jvm_best",
                "jvm_names", "jvm_aicc_no_td", "jvm_aicc_td", "native_stage",
                "native_best", "native_names", "native_aicc_no_td",
                "native_aicc_td", "best_match", "pair_status"])
    for tpl in TPLS:
        for kind in KINDS:
            jf = jvm_atd(CAP / f"{tpl}-{kind}.out")
            nf = native_atd(REP / f"replay-{tpl}-{kind}-native.out")
            for i in range(max(len(jf), len(nf))):
                j = jf[i] if i < len(jf) else None
                n = nf[i] if i < len(nf) else None
                best_match = ""
                if j and n:
                    jb, nb = j[1], n[1]
                    jb = jb.replace("−", "-")
                    nb = nb.replace("−", "-")
                    best_match = "yes" if jb == nb else "NO"
                pair = ("paired" if j and n else
                        "jvm-only" if j else "native-only")
                w.writerow([
                    tpl, kind, i + 1,
                    j[0] if j else "", j[1] if j else "",
                    j[2] if j else "",
                    j[4][0] if j and len(j[4]) > 0 else "",
                    j[4][1] if j and len(j[4]) > 1 else "",
                    n[0] if n else "", n[1] if n else "",
                    n[2] if n else "",
                    n[4][0] if n and len(n[4]) > 0 else "",
                    n[4][1] if n and len(n[4]) > 1 else "",
                    best_match, pair])

print("wrote", TAB / "decision-bits.csv", TAB / "atd-firings.csv")
