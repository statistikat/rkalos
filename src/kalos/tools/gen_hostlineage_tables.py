#!/usr/bin/env python3
"""Generate src/hostlineage_tables.cpp from pinned glibc-2.39 table sources.

MECHANICAL TRANSCRIPTION ONLY: this script performs the following
transformations and nothing else (provenance: glibc tag glibc-2.39,
sysdeps/ieee754/dbl-64/, official mirror https://github.com/bminor/glibc):

  * verify every input file against the pinned SHA256 in
    tools/hostlineage-glibc-sha256.json;
  * strip comments and glibc-internal #include plumbing;
  * evaluate the constant #if guards of the table-driven kernels (the
    x86-64 glibc-2.39 table macros EXP/LOG/POW_TABLE_BITS are fixed and
    the script REFUSES unexpected guards or entry counts);
  * keep ONLY the LITTLE_ENDI branch of the byte-swapped IBM APML table
    blocks (this transcription is validated bit-exactly against the
    little-endian host libm; everything else is dropped);
  * fold the compile-time `* N` / `/ N` scaling factors to literal hex
    floats (exact: N is a power of two);
  * convert C designated initializers to positional initializers in
    struct-declaration order (C++17 has no designated initializers);
  * rename __sincostab/__exp_data/__log_data/__pow_log_data (namespace
    kalos::detail::hostlineage_src).

Every numeric literal is copied verbatim (or exact-folded); count
assertions guard every table. tools/verify_hostlineage.cpp re-checks
element-by-element against values parsed independently from the original
glibc files AND bit-equality of every kernel against the host libm.
"""
import hashlib
import json
import re
import struct
import sys
from pathlib import Path

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "third_party/glibc-2.39"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[1] / "src/hostlineage_tables.cpp"

PINNED_SHA256 = json.loads(Path(__file__).with_name("hostlineage-glibc-sha256.json").read_text())


def load(name: str) -> str:
    text = (SRC / name).read_text()
    digest = hashlib.sha256(text.encode()).hexdigest()
    if PINNED_SHA256[name] != digest:
        raise ValueError(f"SHA256 mismatch for {name}: {digest}")
    return text


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return "\n".join(re.sub(r"//.*", "", line) for line in text.splitlines())


def little_endian_branch(text: str, name: str) -> str:
    """Return the body between '#ifdef LITTLE_ENDI' and its closing
    '#endif' (inside the '#ifdef BIG_ENDI ... #else ...' ladder)."""
    lines = text.splitlines()
    out, depth, state = [], 0, "seek"  # state: seek -> in
    started_at = None
    for i, line in enumerate(lines):
        s = line.strip()
        if state == "seek" and s == "#ifdef LITTLE_ENDI":
            state, started_at, depth = "in", i, 0
            continue
        if state == "in":
            if s.startswith("#if"):
                depth += 1
            elif s.startswith("#endif"):
                if depth == 0:
                    return "\n".join(out)
                depth -= 1
                continue
            if s == "#endif" or s.startswith("#define LITTLE") :
                pass
            out.append(line)
    raise ValueError(f"{name}: no LITTLE_ENDI branch found")


GUARD_OK = {"N == 128", "LOG_POLY1_ORDER == 12", "N == 128 && EXP_POLY_ORDER == 5",
            "N == 128 && EXP2_POLY_ORDER == 5", "N == 128 && LOG_POLY_ORDER == 6",
            "N == 128 && POW_LOG_POLY_ORDER == 8", "1"}


def drop_preproc(text: str) -> str:
    kept = []
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("#"):
            kept.append(line)
            continue
        m = re.match(r"#\s*if\b(.*)", s)
        if m and m.group(1).strip() not in GUARD_OK:
            raise ValueError(f"unexpected #if condition: {s!r}")
        # dropped: includes, defines, known-true guards, #endif/#else
    return "\n".join(kept)


def fold_scaled(text: str) -> str:
    def ev(m):
        lit, op = m.group(1), m.group(2)
        v = float.fromhex(lit)
        return float.hex(v * 128.0 if op == "*" else v / 128.0)
    return re.sub(r"(0x[0-9a-fA-F.]+p[-+0-9]+)\s*([*/])\s*N\b", ev, text)


def wrap_int_words(text: str) -> str:
    """Wrap integer hex literals >= 2**31 in IW() (C++ forbids the silent
    unsigned->int narrowing glibc's C initializers rely on)."""
    def ev(m):
        return f"IW({m.group(0)})" if int(m.group(1), 16) >= 0x80000000 else m.group(0)
    return re.sub(r"0[xX]([0-9A-Fa-f]{1,8})\b", ev, text)


def find_block(text: str, start: int) -> str:
    brace = text.index("{", start)
    depth = 0
    for pos in range(brace, len(text)):
        depth += text[pos] == "{"
        depth -= text[pos] == "}"
        if depth == 0:
            return text[brace:pos + 1]
    raise ValueError("unbalanced braces")


def split_top(block: str) -> list[str]:
    inner = block[1:-1]
    parts, depth, cur = [], 0, ""
    for ch in inner:
        if ch in "{(":
            depth += 1
        elif ch in "})":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur.strip())
            cur = ""
        elif ch not in "\r\n" or depth > 0 or cur.strip():
            cur += ch
    if cur.strip():
        parts.append(cur.strip())
    return parts


def designated_to_positional(text: str, fields: list[tuple[str, int]], label: str) -> str:
    block = find_block(text, text.index("="))
    entries = {}
    for part in split_top(block):
        m = re.match(r"\.\s*([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.*)$", part, flags=re.S)
        if not m:
            raise ValueError(f"{label}: not a designated initializer: {part[:60]!r}")
        entries[m.group(1)] = re.sub(r"\s+", " ", m.group(2).strip())
    if set(entries) != {name for name, _ in fields}:
        raise ValueError(f"{label}: fields differ: {sorted(entries)} vs {[f[0] for f in fields]}")
    out = []
    for name, count in fields:
        value = entries[name]
        if count:
            if not value.startswith("{"):
                raise ValueError(f"{label}: {name} expected brace block")
            n = len(split_top(find_block(value, 0)))
            if n != count:
                raise ValueError(f"{label}: {name} has {n} entries, expected {count}")
        out.append(f"  /* .{name} */ {value},")
    return "\n".join(out)


def main() -> int:
    # ---- sincostab.c (sin/cos table) -----------------------------------------
    sc = drop_preproc(little_endian_branch(strip_comments(load("sincostab.c")), "sincostab.c"))
    m = re.search(r"union\s*\{[^{}]*\}\s*__sincostab\s*=", sc)
    sincos_items = split_top(find_block(re.sub(r"\.\s*i\s*=", "", find_block(sc, m.end() - 1))[1:], 0))
    if len(sincos_items) != 880:
        raise ValueError(f"sincostab: {len(sincos_items)} != 880")

    # ---- asincos.tbl (asin/acos table) ----------------------------------------
    ac = drop_preproc(little_endian_branch(strip_comments(load("asincos.tbl")), "asincos.tbl"))
    m = re.search(r"union\s*\{[^{}]*\}\s*asncs\s*=", ac)
    asncs_items = split_top(find_block(re.sub(r"\.\s*i\s*=", "", find_block(ac, m.end() - 1))[1:], 0))
    if len(asncs_items) != 5136:
        raise ValueError(f"asncs: {len(asncs_items)} != 5136")

    # ---- root.tbl / powtwo.tbl -------------------------------------------------
    root = strip_comments(load("root.tbl"))
    inroot = find_block(root, root.index("inroot[128]"))
    if len(split_top(inroot)) != 128:
        raise ValueError("inroot count")
    pt = strip_comments(load("powtwo.tbl"))
    powtwo = find_block(pt, pt.index("powtwo[]"))
    if len(split_top(powtwo)) != 28:
        raise ValueError("powtwo count")

    # ---- uatan.tbl (atan2 table cij[241][7]) -----------------------------------
    ua = drop_preproc(little_endian_branch(strip_comments(load("uatan.tbl")), "uatan.tbl"))
    cij = find_block(ua, ua.index("cij[241][7]"))
    rows = split_top(cij)
    if len(rows) != 241 or any(len(split_top(find_block(r, 0))) != 7 for r in rows):
        raise ValueError("cij shape")

    # ---- table-driven kernel data ----------------------------------------------
    def expand_pow_tab(text: str) -> str:
        # e_pow_log_data.c emits its table via A(invc,logc,logctail) -> {invc,0,logc,logctail}
        out, n = [], 0
        for line in text.splitlines():
            m = re.match(r"A\((.*)\)\s*$", line.strip())
            if m:
                a, b, c = [x.strip() for x in m.group(1).split(",")]
                out.append(f"{{{a}, 0.0, {b}, {c}}},")
                n += 1
            else:
                out.append(line)
        if n != 128:
            raise ValueError(f"pow_log_data: expanded {n} A() rows, expected 128")
        return "\n".join(out)

    exp_body = designated_to_positional(fold_scaled(drop_preproc(strip_comments(load("e_exp_data.c")))), [
        ("invln2N", 0), ("shift", 0), ("negln2hiN", 0), ("negln2loN", 0),
        ("poly", 4), ("exp2_shift", 0), ("exp2_poly", 5),
        ("invlog10_2N", 0), ("neglog10_2hiN", 0), ("neglog10_2loN", 0),
        ("exp10_poly", 5), ("tab", 256)], "exp_data")
    log_body = designated_to_positional(fold_scaled(drop_preproc(strip_comments(load("e_log_data.c")))), [
        ("ln2hi", 0), ("ln2lo", 0), ("poly", 5), ("poly1", 11),
        ("tab", 128), ("tab2", 128)], "log_data")
    pow_body = designated_to_positional(fold_scaled(drop_preproc(strip_comments(expand_pow_tab(load("e_pow_log_data.c"))))), [
        ("ln2hi", 0), ("ln2lo", 0), ("poly", 7), ("tab", 128)], "pow_log_data")
    for label, body in (("exp", exp_body), ("log", log_body), ("pow", pow_body)):
        if re.search(r"\bN\b", body):
            raise ValueError(f"{label}: unfolded N token remains")

    header = """// GENERATED by tools/gen_hostlineage_tables.py — DO NOT EDIT.
// Vendored glibc-2.39 (tag glibc-2.39, sysdeps/ieee754/dbl-64) TABLE DATA for
// the host-lineage kernels in src/hostlineage.cpp. Literal-for-literal
// transcription; the generator docstring lists every transformation applied
// and tools/hostlineage-glibc-sha256.json pins the input digests.
// Licensing: these table sources are under the GNU LGPL-2.1-or-later; see
// docs/contracts/x13-automatic-pipeline-contract.md (CI-FDLIBM section,
// hostlineage provenance note) for the full note.

#include "hostlineage_internal.hpp"

namespace kalos::detail::hostlineage_src {

"""
    body = f"""// ---- sincostab.c: __sincostab -> sincostab (880 int4, LITTLE_ENDI) ----------
const sincostab_t sincostab = {{
    {", ".join(wrap_int_words(s) for s in sincos_items)}}};

// ---- asincos.tbl: asncs (5136 int4, LITTLE_ENDI) ------------------------------
const asncs_t asncs = {{
    {", ".join(wrap_int_words(s) for s in asncs_items)}}};

// ---- root.tbl / powtwo.tbl -----------------------------------------------------
const double inroot[128] = {inroot};
const double powtwo[28] = {powtwo};

// ---- uatan.tbl: cij[241][7] (LITTLE_ENDI) ---------------------------------------
const mynumber cij[241][7] = {wrap_int_words(cij)};

// ---- e_exp_data.c: __exp_data -> exp_data ----------------------------------------
const struct exp_data exp_data = {{
{exp_body}
}};

// ---- e_log_data.c: __log_data -> log_data -----------------------------------------
const struct log_data log_data = {{
{log_body}
}};

// ---- e_pow_log_data.c: __pow_log_data -> pow_log_data ------------------------------
const struct pow_log_data pow_log_data = {{
{pow_body}
}};

}} // namespace kalos::detail::hostlineage_src
"""
    OUT.write_text(header + body)
    print(f"wrote {OUT} ({len((header + body).splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
