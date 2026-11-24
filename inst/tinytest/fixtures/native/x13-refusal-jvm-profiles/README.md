# X-13 refusal-demand JVM profiles

The nine monthly demand records cover the three series and the RSA0, RSA1
and RSA2c templates used by Persephone's M4 refusal test. The original
records live under `inst/tinytest/fixtures/native/refusals-m4gap` at
Persephone revision `b625852e5ea1ff8cfa6f64b9fe6a7c1fdb940c09`.
The three CSVs contain the start year/month and the exact round-trip inputs.
Their hashes, the full jar classpath hashes and the original revision are
recorded in [reference.json](reference.json).

[demand.tsv](demand.tsv) records the projected theta, seasonal theta and
regression coefficients, with the SHA-256 of each **entire** 20-line demand
record. All nine x86-64 default records reproduce the historical bytes. Eight
portable records have differing coefficients. All other fields and
the ninth complete record agree. The native M4 refusal boundary is unchanged.

## Reproduction

Use Temurin 21.0.11 and the pinned jars shipped by rjd3jars 0.0.6,
rjd3x13 3.9.0 and rjd3toolkit 3.9.1. Make these R packages and
RProtoBuf 0.4.28 available to Rscript, including a host-compatible rJava.
Run from the Kalos repository root with the paths to those installations.

```sh
python3 validation/reference/capture-x13-refusal-profiles.py \
  --jdk "$REFERENCE_JDK" \
  --jar-dir "$R_LIBRARY/rjd3jars/java" \
  --jar-dir "$R_LIBRARY/rjd3x13/java" \
  --jar-dir "$R_LIBRARY/rjd3toolkit/java"
```

The probe executes the actual X13 kernel, serializes each result and reads
the forecast length from the kernel. R projects the protobuf to the original
text format. The runner checks the full digest and the displayed coefficient
fields. On ARM64 the default JVM reproduces the portable records. On x86-64
the default JVM reproduces the historical records. Pass `--portable-x86`
to independently reproduce the portable records on x86-64 by disabling
the exp, log, pow, cos, sin and tan JVM intrinsics.

Capture and replay succeeded on macOS ARM64 with both native ARM64 and
translated x86-64 Temurin builds. This panel describes reference decisions
used to test refusal. It does not qualify execution of the templates by Kalos.
