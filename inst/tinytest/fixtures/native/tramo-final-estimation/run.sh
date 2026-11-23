#!/bin/sh
# 23-F F.1 probe runner (Scratch23f.java, split package tramo) against the
# pinned rjd3 3.9.0 jars; G2 intrinsic-disable recipe (21-D.2 precedent).
# Usage: ./run.sh <CaseFile> <outdir-tag> [extra args]
set -e
cd "$(dirname "$0")"
JARS=/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3jars/java
TK=/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3toolkit/java
TS=/data/home/kowa/R/x86_64-pc-linux-gnu-library/4.5/rjd3tramoseats/java
CP="$JARS/jdplus-sa-base-api-3.9.0.jar:$JARS/jdplus-sa-base-core-3.9.0.jar:$JARS/jdplus-toolkit-base-core-3.9.0.jar:$JARS/jdplus-toolkit-base-api-3.9.0.jar:$JARS/protobuf-java-4.36.1.jar:$TK/jdplus-toolkit-base-protobuf-3.9.0.jar:$TK/jdplus-sa-base-r-3.9.0.jar:$TK/jdplus-toolkit-base-r-3.9.0.jar:$TK/jdplus-sa-base-protobuf-3.9.0.jar:$TS/jdplus-tramoseats-base-protobuf-3.9.0.jar:$TS/jdplus-tramoseats-base-api-3.9.0.jar:$TS/jdplus-tramoseats-base-r-3.9.0.jar:$TS/jdplus-tramoseats-base-core-3.9.0.jar"
[ -d classes ] || mkdir -p classes
/usr/lib/jvm/java-21-openjdk-amd64/bin/javac --release 21 -cp "$CP" -d classes Scratch23f.java
tag=$2
mkdir -p "out/$tag"
shift 2
env _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog" \
  /usr/lib/jvm/java-21-openjdk-amd64/bin/java -Xss16m -cp "classes:$CP" \
  jdplus.tramoseats.base.core.tramo.Scratch23f "$@" > "out/$tag/$tag.out" 2> "out/$tag/$tag.err"
echo "wrote out/$tag/$tag.out ($(wc -l < "out/$tag/$tag.out") lines)"
