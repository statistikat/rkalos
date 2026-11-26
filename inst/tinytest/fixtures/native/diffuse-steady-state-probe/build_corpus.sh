#!/bin/bash
# Item-29 probe slice — frozen corpora builder.
#
# Deterministic: synthetic series come from gen_series (fixed 64-bit LCG +
# Box-Muller; byte-reproducible on any IEEE-754 binary64 platform). Real
# series are byte-copies of the pinned gapanalysis probe inputs (SHA256 in
# SHA256SUMS.txt; jar/commit provenance in the dossier). Re-running this
# script from the kalos probe scratch tree must reproduce every listed file
# byte-for-byte.
set -eu
cd "$(dirname "$0")"
GS=../cpp/gen_series

# --- real series: byte-copies of the M9 gap-analysis probe inputs ----------
GP=/data/home/kowa/seasonaladjustment/tmp_opencode/m9-perf/gapanalysis
cp "$GP/probe-12-airline-complete-none-yadj.csv" series/real-airline-y.csv
cp "$GP/probe-12-airline-gaps-none-yadj.csv"    series/real-airline-gaps.csv
cp "$GP/probe-12-airline-complete-none-x.csv"   x/real-x.csv
cp "$GP/probe-12-airline-gaps-none-x.csv"       x/real-gaps-x.csv

# leading-gap diffuse-delay variant: first five observations masked
awk -F, 'NR==1{print;next} NR<=6{print "NA,1";next}{print}' \
    series/real-airline-y.csv > series/real-airline-leadgap.csv

# redundant-factor design (col2 = 2*col1 exactly) and zero-column design
awk -F, 'NR==1{print "x1,x2";next}{printf "%s,%.17g\n", $1, 2*$1}' \
    x/real-x.csv > x/redundant2.csv
awk -F, 'NR==1{print "x1,x2";next}{printf "%s,%.17g\n", $1, 0}' \
    x/real-x.csv > x/zerocol2.csv

# --- synthetic series -------------------------------------------------------
$GS series/syn-airline12.csv   12 1 1 ";-0.45;;-0.55"        156 20261101
$GS series/syn-stat-arma21.csv 12 0 0 "0.6,-0.25;-0.4;;"     240 20261102
$GS series/syn-stat-sarma.csv  12 0 0 "0.5;-0.3;0.25;-0.2"   240 20261103
$GS series/syn-airline4.csv    4  1 1 ";-0.35;;-0.6"         96  20261104
$GS series/syn-mixed4.csv      4  1 1 ";;0.3;-0.5"           96  20261105
$GS series/syn-d2.csv          12 2 0 ";-0.5;;"              144 20261106
$GS series/syn-slowseasonal.csv 12 1 1 ";-0.2;;0.999"        156 20261107
$GS series/syn-nearur.csv      12 1 1 "0.995;-0.4;;-0.55"    144 20261108
$GS series/syn-nearni.csv      12 1 1 ";0.98;;0.95"          144 20261109
