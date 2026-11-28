#!/usr/bin/env python3
"""item-30 probe: symbolize SIGPROF PC histogram against the -no-pie binary."""
import bisect, sys, collections

symbols = []
for line in open(sys.argv[1]):  # nm -C --defined-only output
    parts = line.split(" ", 2)
    if len(parts) != 3:
        continue
    addr, kind, name = parts
    if kind not in "TtWwVvIiRr":
        continue
    try:
        symbols.append((int(addr, 16), name.strip()))
    except ValueError:
        continue
symbols.sort()
addrs = [a for a, _ in symbols]

hist = collections.Counter()
total = 0
for line in open(sys.argv[2]):
    fields = line.split()
    if len(fields) != 2:
        continue
    try:
        pc = int(fields[0], 16)
    except ValueError:
        continue
    count = int(fields[1])
    total += count
    i = bisect.bisect_right(addrs, pc) - 1
    hist[symbols[i][1] if i >= 0 else "?"] += count

print(f"total samples {total}")
for name, count in hist.most_common(25):
    print(f"{100.0 * count / total:6.2f}%  {count:6d}  {name[:110]}")
