#!/bin/bash
# Item-29 JVM probe runner: pinned model/series, deterministic captures x2.
# Usage: run_jvm.sh OUTDIR
set -eu
OUT=$1
mkdir -p "$OUT"
J=/usr/lib/jvm/java-21-openjdk-amd64/bin/java
CP=.:/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3jars/java/jdplus-toolkit-base-core-3.9.0.jar:/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3jars/java/jdplus-toolkit-base-api-3.9.0.jar
S=../corpus/series
X=../corpus/x
run() { # name y x orders coefs evals
  nice -n 10 "$J" -cp "$CP" -XX:ActiveProcessorCount=1 Item29FilterProbe \
      "$1" "$2" "$3" "$4" "$5" "$6" > "$OUT/$1.txt" 2>&1 || true
}
for pass in 1 2; do
  if [ "$pass" = 2 ]; then OUT="$OUT.b"; mkdir -p "$OUT"; fi
  echo "load_before_pass$pass=$(cat /proc/loadavg)"
  run jvm-airline-x       "$S/real-airline-y.csv" "$X/real-x.csv"   0,1,1,0,1,1,12 ";-0.45;;-0.55"            6
  run jvm-airline-noreg   "$S/real-airline-y.csv" NULL              0,1,1,0,1,1,12 ";-0.45;;-0.55"            6
  run jvm-stat-arma21     "$S/syn-stat-arma21.csv" NULL             2,0,1,0,0,0,12 "0.6,-0.25;-0.4;;"         6
  run jvm-stat-arma21n    "$S/syn-stat-arma21.csv" NULL             2,0,1,0,0,0,12 "-0.6,0.25;-0.4;;"        6
  run jvm-stat-sarma11    "$S/syn-stat-sarma.csv" NULL              1,0,1,1,0,1,12 "-0.5;-0.3;-0.25;-0.2"     6
  run jvm-nearur          "$S/real-airline-y.csv" NULL              1,1,1,0,1,1,12 "-0.995;-0.4;;-0.55"        6
  run jvm-nearni          "$S/real-airline-y.csv" NULL              0,1,1,0,1,1,12 ";0.98;;0.95"              6
  run jvm-d2              "$S/syn-d2.csv" NULL                      0,2,1,0,0,0,12 ";-0.5;;"                  6
  run jvm-slow            "$S/syn-slowseasonal.csv" NULL            0,1,1,0,1,1,12 ";-0.2;;0.999"             6
  run jvm-airline4        "$S/syn-airline4.csv" NULL                0,1,1,0,1,1,4  ";-0.35;;-0.6"             6
  run jvm-mixed4          "$S/syn-mixed4.csv" NULL                  0,1,0,1,1,1,4  ";;-0.3;-0.5"              6
  run jvm-airline-gaps-x  "$S/real-airline-gaps.csv" "$X/real-gaps-x.csv" 0,1,1,0,1,1,12 ";-0.45;;-0.55"      6
done
echo "load_after=$(cat /proc/loadavg)"
