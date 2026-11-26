# Secondary baseline: native-AFH P4 natural-xenon halo fold

This archive imports the completed 2026-11-22 P4 native-AFH result from the source repository without rerunning the calculation.

## Source identity

- Source repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Result/archive commit: `55e06b9403f7ac3e2913d1dd5d192167bf80617d`
- P4 run source/workspace commit: `fada259460ac42460baf933021dc16dc7a9cfb9e`
- Original result directory: `research/three_body_full_current_nucleon_match/results/2026-11-22/p4_native_afh_regeneration/`
- Native AFH dataset: `nuclear/afh/data/native_v6_2026-11-22/` (already imported separately)
- AFH package: DMFormFactor v6.0; SHA-256 `dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372`. The package itself is not redistributed.

## Contents

- Five-recoil natural-xenon result and full 81-point recoil spectrum for two halo models.
- Per-isotope contributions and natural-xenon recoil-spectrum CSVs.
- Matched provisional/native AFH comparison, P4 audit, and quadrature convergence records.
- Original runtime provenance and original SHA-256 manifest.
- Original plots (PNG/PDF).
- Seven AFH default-density OBDME extraction tables and extraction metadata, plus natural-xenon isotope weights.

The OBDME tables here are the AFH-B&D/default-density inputs, not SN100PN data.

## Numerical status

The source archive reports:
- 16,807 native AFH response rows passed pointwise comparison.
- Five-recoil and 81-point P4 audits passed.
- Speed quadrature 24-to-48 relative difference: (8.9978914\times10^{-6}), against (10^{-5}).
- Three-body phase-space quadrature 8x6x4-to-10x8x6 relative difference: (2.5525175\times10^{-9}), against (10^{-8}).
- Matched provisional/native rate differences: at most (1.53\times10^{-10}) on the continuous grid.

## Physical interpretation

This is a **secondary leading-static one-body vector-current candidate**, not the paper Fig. 3 event-spectrum reproduction and not a complete bound-nucleus VV current prediction. The source records retain `full_current_complete=false` and `P3_tensor_validation_complete=false`. The P3-C bound-current amplitude gate remains separate.

The native AFH package and run were executed locally on Windows 11 with Wolfram 14.2 and Python 3.14.7; the exact runtime, CPU/GPU, packages, source hashes, input hashes and output hashes are preserved in `provenance.json`.

## Integrity

- `SHA256SUMS.txt` is the original source-result checksum list; it has not been regenerated or edited.
- `IMPORT_MANIFEST.json` records source commits, original Git blob IDs, and imported target blob IDs.
- All imported result payloads were transferred as Git blobs with matching source/target Git blob SHA-1 values.
