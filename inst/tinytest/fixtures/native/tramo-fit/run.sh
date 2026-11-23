#!/bin/sh
# 23-F F.4 probe runner (Scratch23f, TramoKernel round driver + fit capture
# surface): one JVM per
# case (per-file JVM capture hygiene), G2 intrinsic-disable recipe.
# Usage: ./run.sh <CaseFile> <outdir-tag> [witness] [mode]
set -e
cd "$(dirname "$0")"
JARS=/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3jars/java
TK=/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3toolkit/java
TS=/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3tramoseats/java
CP="$JARS/jdplus-sa-base-api-3.9.0.jar:$JARS/jdplus-sa-base-core-3.9.0.jar:$JARS/jdplus-toolkit-base-core-3.9.0.jar:$JARS/jdplus-toolkit-base-api-3.9.0.jar:$JARS/protobuf-java-4.36.1.jar:$TK/jdplus-toolkit-base-protobuf-3.9.0.jar:$TK/jdplus-sa-base-r-3.9.0.jar:$TK/jdplus-toolkit-base-r-3.9.0.jar:$TK/jdplus-sa-base-protobuf-3.9.0.jar:$TS/jdplus-tramoseats-base-protobuf-3.9.0.jar:$TS/jdplus-tramoseats-base-api-3.9.0.jar:$TS/jdplus-tramoseats-base-r-3.9.0.jar:$TS/jdplus-tramoseats-base-core-3.9.0.jar"
[ -d classes ] || mkdir -p classes
/usr/lib/jvm/java-21-openjdk-amd64/bin/javac --release 21 -cp "$CP" -d classes Scratch23f.java
cf=$1
tag=$2
witness=$3
mode=$4
mkdir -p "out/$tag"
splitdir="out/$tag/cases"
rm -rf "$splitdir"
mkdir -p "$splitdir"
awk -v d="$splitdir" 'BEGIN{n=0; f=sprintf("%s/%03d.txt", d, n)} /^\S/ && $1=="KER" {n++; f=sprintf("%s/%03d.txt", d, n)} {print > f}' "$cf"
: > "out/$tag/$tag.out"
: > "out/$tag/$tag.err"
XLOG=""
[ "$witness" = witness ] && XLOG=on
: > "out/$tag/$tag.classload.log"
for c in "$splitdir"/*.txt; do
  if [ "$XLOG" = on ]; then
    # The class+load log MUST go to its own file: -Xlog on stdout interleaves
    # with large System.out writes and splices data lines.
    env _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog" \
      /usr/lib/jvm/java-21-openjdk-amd64/bin/java -Xss16m \
      "-Xlog:class+load=info:file=$c.load.log" -cp "classes:$CP" \
      jdplus.tramoseats.base.core.tramo.Scratch23f "$c" ${mode:+"$mode"} >> "out/$tag/$tag.out" 2>> "out/$tag/$tag.err"
    sed 's/^\[[^]]*\]\[[^]]*\]//' "$c.load.log" >> "out/$tag/$tag.classload.log"
  else
    env _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog" \
      /usr/lib/jvm/java-21-openjdk-amd64/bin/java -Xss16m -cp "classes:$CP" \
      jdplus.tramoseats.base.core.tramo.Scratch23f "$c" ${mode:+"$mode"} >> "out/$tag/$tag.out" 2>> "out/$tag/$tag.err"
  fi
done
cp "out/$tag/$tag.out" "out/$tag/$tag.clean.out"
echo "wrote out/$tag/$tag.clean.out ($(wc -l < "out/$tag/$tag.clean.out") lines)"
