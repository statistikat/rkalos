# SEATS WK factorizer tiers

The refused-profile JVM oracle for the Wiener-Kolmogorov factorizer chain. The
input is [seas.csv](../x13-refusal-jvm-profiles/seas.csv) (digest pinned in
[reference.json](reference.json)): the monthly profile whose SEATS WK
extraction the native mirror refused because the total error variance filter
of the Burman variance estimates resolved on no factorizer tier. Its
one-sided weights are recorded as `tevar_filter` and appear in the call ledger
below as calls 3 and 4.

## Captures

All four captures replay the same pinned jdplus 3.9.0 stack and print every
value as a raw binary64 bit pattern. Each was produced twice, in two
independent JVM invocations, byte-identical, before freezing. They use the
pinned x86-64 numerical profile and the same JVM math flags as the raw
TRAMO-SEATS preset stages under `../tramo-seats-bridge/`.

| File | Probe | Content |
| --- | --- | --- |
| [wk-x86_64.txt](wk-x86_64.txt) | `SeatsWkFactorizerProbe` | RSA3 and RSAfull on seas.csv: fit facts and residuals, the ucarima model and its compact form with the per-component symmetric MA and all four tier results, then every final-component and merged series |
| [factorize-calls-x86_64.txt](factorize-calls-x86_64.txt) | `SeatsWkFactorizeCallsProbe` | every `SymmetricFilter.factorize()` of the RSA3 run (twelve calls) with caller frames, the standard-chain result and one row per installed tier |
| [robust-stages-x86_64.txt](robust-stages-x86_64.txt) | `SeatsWkRobustStagesProbe` | stage walk of the robust tier on the refused filter: two-sided weights, w0 normalization, `UnitRootsSolver(0)`, the square-root and least-squares route, the symmetric Muller-Newton solver, the tier result |
| [muller-adversarial-x86_64.txt](muller-adversarial-x86_64.txt) | `SeatsWkMullerAdversarialProbe` | the two adversarial filters of the M6 factorizer fuzz corpus (indices 144 and 304) through all five tiers, with the corpus input row |

The ledger fixes the JVM's tier contract: the refused filter (calls 3 and 4)
resolves on `robust` alone with scaling `3f01dca02266dac2`; every other filter
of the run resolves at `ev2`. The wk capture freezes the extraction the native
mirror must now reproduce: 29 component rows per preset (58 in total) where it
previously raised `ArimaException` with zero finals.

The adversarial capture freezes the two fuzz corpus filters that found the
`+-0.0` sensitivity of the robust tier's `SymmetricMullerNewtonSolver`: its
first `initialize()` division computes `cart(-2.0, -0.0)` and `api.math.Complex`
answers with the `NEG_TWO` singleton `(-2.0, +0.0)`; the normalized sign of the
zero imaginary part decides the parabola-root square-root branch two steps
later and both filters resolve on `robust ok`. A port that builds its complex
values without the `cart` normalization loses that sign and throws
`ArrayIndexOutOfBoundsException` from the solver iteration guard instead.

## Reproduction

Use Temurin 21.0.11 and the pinned jars shipped by rjd3jars 0.0.6,
rjd3x13 3.9.0 and rjd3tramoseats 3.9.1. Run from the Kalos repository root:

```sh
python3 validation/reference/capture-seats-wk-factorizer-tiers.py \
  --jdk "$REFERENCE_JDK" \
  --jar-dir "$R_LIBRARY/rjd3jars/java" \
  --jar-dir "$R_LIBRARY/rjd3toolkit/java" \
  --jar-dir "$R_LIBRARY/rjd3x13/java" \
  --jar-dir "$R_LIBRARY/rjd3tramoseats/java"
```

The script re-checks the input and jar digests, compiles the probes with
`javac --release 21`, runs each probe twice and requires byte-identical
outputs and the frozen capture digests.
