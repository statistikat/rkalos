#!/usr/bin/env bash
# m17-td1 stage-1 TD-only battery: {td2noe,td7noe,td2lev,td2noadj} x
# {clean,tdeff,air}, bridge domain (the committed capture domain), double
# capture (run1/run2) for the byte-identity protocol modulo the timing line.
set -u
cd /data/home/kowa/seasonaladjustment/tmp_opencode/m17-td1
CP='classes:/data/home/kowa/seasonaladjustment/tmp_opencode/jarcache/*'
JAVA=/usr/lib/jvm/java-21-openjdk-amd64/bin/java
mkdir -p run1 run2
for tpl in td2noe td7noe td2lev td2noadj; do
  for kind in clean tdeff air; do
    env _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog" nice -n 10 "$JAVA" -cp "$CP" Td1Probe "$tpl" "$kind" > "run1/$tpl-$kind.out" 2> "run1/$tpl-$kind.err"
    env _JAVA_OPTIONS="-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog" nice -n 10 "$JAVA" -cp "$CP" Td1Probe "$tpl" "$kind" > "run2/$tpl-$kind.out" 2> "run2/$tpl-$kind.err"
  done
done
echo TD1_BATTERY_DONE
