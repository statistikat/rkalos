# rkalos implementation status

Updated: 2026-11-14. Repository extraction in progress.

## Completed

- Preserved relevant Git history from Persephone commit `8450d6ec8c205f52715fa36d8c6c6d0826e78c99`.
- Recorded original-to-filtered commit mappings and extraction path rules.
- Adopted the `kalos` C++ namespace and include paths; numerical algorithms are unchanged.
- Preserved all original reference fixtures byte for byte.

## Next steps

1. Separate R and OpenBLAS providers and qualify the standalone core.
2. Move shared preparation into C++ and connect the lean R package.
3. Pin generated core snapshots, add minimal Python bindings, validate artifacts.
4. Connect Persephone through a **required** `rkalos` dependency and publish the repositories.

The inherited numerical implementation extends through fixed additive initial X-11 tables
(22b.2). Full native seasonal-adjustment fitting remains unavailable. The historical
baseline and numerical contracts belong to Kalos; R boundary tests belong to rkalos.
