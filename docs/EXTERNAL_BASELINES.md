# External baselines and source pins

This repository is independent of the original model repository. The
following immutable source revisions define the initial reuse boundary.

## Original VV/P3/P4 baseline

- Repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Branch at the time of pinning: \`research/p4-provisional-halo-fold-20261119\`
- Exact source commit: \`dd9a7387814cd157d9cde739fa80181e4cd7d428\`
- Role: source anchor for the original radiative VV benchmark, the
  three-body spin-resolved kernel, and the P3/P4 implementation.
- Curated imported files: see
  \`benchmarks/vv_recovery/source_snapshot/\` and
  \`nuclear/afh/data/native_v6_2026-11-22/\`.
- Native AFH evidence: 16,807 response rows and 9,604 rank-resolved spatial
  amplitude rows passed the archived Native comparison.
- AFH package SHA-256:
  \`dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372\`.
- Import record:
  \`results/provenance/VV_BASELINE_IMPORT_2026-11-23.json\`.

The upstream \`dmformfactor-V6.m\` source package is not copied. The native
Wolfram exporter scripts and generated Native data are included; the
external package must be installed separately and match the pinned hash.

## Paper Fig. 3 benchmark definition

- Source repository: same as above.
- Source commit: \`5ea4d360fe03ce58079220bb0b7624b3faf03fb8\`.
- Source path:
  \`research/three_body_full_current_nucleon_match/benchmarks/paper_fig3_vector_natural_xenon.json\`.
- Source Git blob SHA-1: \`af3731952da389872de21358050750235f9a688d\`.
- Imported SHA-256:
  \`940fad45c3d408a94ed755adc29e11e966d899317faa49956527e221e2035308\`.

## Migrated general-current theory note

- Source repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Source branch: \`research/general-radiative-nuclear-currents-20261123\`
- Source commit: \`de39bf15e08870af82e5217ac30092e0d8e89878\`
- Source path: \`research/radiative_dark_matter_nuclear_current_theory/GENERAL_2TO3_CURRENT_FRAMEWORK.md\`
- Source Git blob SHA-1: \`e86874ab7d9bba447909b86364f915980e23004b\`
- Source file SHA-256: \`ba5bd93752244f8627504a9af868d813fbc864cd4b06372f2f5cd80d9746a04b\`

The source note is retained as a historical record in the original repository.
This repository contains a corrected migration at
\`docs/GENERAL_2TO3_CURRENT_FRAMEWORK.md\`; see
\`docs/MIGRATION_NOTES.md\` and
\`results/provenance/INITIALIZATION_MANIFEST.json\`.

## Universal vector-current literature audit

- Repository-local theory note: \`docs/VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md\`.
- Provenance record: \`results/provenance/VV_VECTOR_CURRENT_ADDENDUM_2026-11-23.json\`.
- The note uses the original LZ P4 baseline pin above and the P3-C
  bound-current audit at source commit
  \`d02e5e009f32a1db1b2520a789148bbd0160ccfd\`.
- This is a literature/matching import only. The numerical Native AFH data
  are imported separately in this PR.

## Reuse policy

The original repository remains the historical calculation archive. This
repository imports selected, pinned components after checking physical
scope, source identity, input hashes, dependencies, and redistribution
rights. It does not silently promote the original leading one-body impulse
result to a complete chiral-EFT current calculation.
