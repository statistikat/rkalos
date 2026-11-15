list(schema_version = 1L, captured_utc = "2026-11-15 17:31:30 UTC",
    source_sha256 = c(`capture-x11-sigma.R` = "f5e25963f28cb65e5a4da80b26472f4b37e207ea6423311dff61fdd24272e8a0",
    `x11-sigma-cases.R` = "73c7877256fcb750d6a4eaad5709970cfe5cf8c936481376c5a3ac9c3ec48ed4"
    ), java = list(version = "21.0.11", vendor = "Ubuntu", runtime_version = "21.0.11+10-1-24.04.2-Ubuntu",
        vm_name = "OpenJDK 64-Bit Server VM", architecture = "amd64",
        oracle = "rjd3 3.9.0/3.9.1 pinned jars via rjd3jars classpath; DefaultExtremeValuesCorrector (the pinned default rules, equivalent to X11SigmaPolicy::legacy) driven through its public analyse/getObservationWeights API, with its calcStdev/extremeValuesDetection pass primitives and private removeExtremes exclusion helper re-run for per-pass provenance"),
    R = "R version 4.5.3 (2026-03-11)", cases = 125L, public_scope = "Pinned corrector dispersion windows, two-pass weights and exclusions only; captured artifacts reflect the pinned default rules (native X11SigmaPolicy::legacy); the native accurate policy is a deliberate alternative qualified by contract analysis, not by this oracle; no replacement values, correction factors, endpoint policy or decomposition claim")
