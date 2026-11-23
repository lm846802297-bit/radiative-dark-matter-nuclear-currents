# External baselines and source pins

This repository is independent of the original model repository. The following immutable source revisions define the initial reuse boundary.

## Original VV/P3/P4 baseline

- Repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Branch at the time of pinning: research/p4-provisional-halo-fold-20261119
- Exact commit: dd9a7387814cd157d9cde739fa80181e4cd7d428
- Role: source anchor for the original radiative VV benchmark, validated AFH v6.0 inputs, and the P3/P4 implementation.
- Recorded scope at that commit: native seven-isotope/eight-response AFH production and P4 halo-fold outputs are validated; the full chiral-EFT nuclear current remains incomplete.
- Import status: PR #2 imported 15 pinned VV/P3 source files, the frozen Fig. 3 benchmark definition, and curated Native AFH response outputs. The upstream AFH package itself is not redistributed. No Fig. 3 event spectrum or new cross section has been produced in this repository; see docs/VV_BASELINE_MIGRATION.md.

Do not replace this pin with a moving branch name. Any later import must record the exact file path, source commit, Git blob ID, SHA-256, license/redistribution status, and independent reproduction result.

## Migrated general-current theory note

- Source repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Source branch: research/general-radiative-nuclear-currents-20261123
- Source commit: de39bf15e08870af82e5217ac30092e0d8e89878
- Source path: research/radiative_dark_matter_nuclear_current_theory/GENERAL_2TO3_CURRENT_FRAMEWORK.md
- Source Git blob SHA-1: e86874ab7d9bba447909b86364f915980e23004b
- Source file SHA-256: ba5bd93752244f8627504a9af868d813fbc864cd4b06372f2f5cd80d9746a04b

The source note is retained as a historical record in the original repository. This repository contains a corrected migration at docs/GENERAL_2TO3_CURRENT_FRAMEWORK.md; see docs/MIGRATION_NOTES.md and results/provenance/INITIALIZATION_MANIFEST.json.

## Reuse policy

The original repository remains the authoritative home for the first paper's code and historical calculations. This repository may import selected, pinned components after checking their physical scope, source identity, input hashes, dependencies, and redistribution rights. The new project must not silently inherit the original repository's assumptions or mark the original leading one-body impulse result as a complete chiral-EFT current calculation.

## Universal vector-current literature audit

- Repository-local theory note: `docs/VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md`.
- Provenance record: `results/provenance/VV_VECTOR_CURRENT_ADDENDUM_2026-11-23.json`.
- The note uses the original LZ P4 baseline pin above and the P3-C bound-current audit at source commit `d02e5e009f32a1db1b2520a789148bbd0160ccfd`.
- This particular audit is a literature/matching import only. The separate curated VV/P3 baseline import below contains selected LZ source code and raw Native AFH response outputs; no xenon OBDME wave functions or new numerical rates were imported.
- The narrow q_A/q_N interface PR (#111) was merged into the source project's P4 research branch at commit `f63c58249cdf40f4551f0b0b475332e1314e2af7`; PR #112 then corrected its regression assertion only, at P4 HEAD `55e06b9403f7ac3e2913d1dd5d192167bf80617d`. Neither change has been merged to the source project's `main`; the follow-up changes no physics kernel or numerical input.

## Curated VV/P3 baseline import (2026-11-23)

- Core source repository: \`lm846802297-bit/LZ-3-body-for-hdm\`.
- Exact source commit: \`dd9a7387814cd157d9cde739fa80181e4cd7d428\`.
- Paper Fig. 3 benchmark-definition commit:
  \`5ea4d360fe03ce58079220bb0b7624b3faf03fb8\`.
- Imported Native AFH dataset: 7 xenon isotopes, 16,807 eight-response rows,
  9,604 rank-resolved spatial multipole rows; Native imaginary residual zero.
- AFH DMFormFactor v6.0 SHA-256:
  \`dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372\`.
- Curated import record:
  \`results/provenance/VV_BASELINE_IMPORT_2026-11-23.json\`.

The upstream AFH package source is not redistributed. The raw Native response
tables retain AFH isoscalar normalization; apply the equal-\(p=n\) factor
four exactly once downstream. The P3-C current-completeness gate remains
open. The source repository's exploratory eigenmode analysis and secondary
P4 halo spectrum are not promoted as the paper's main result.
