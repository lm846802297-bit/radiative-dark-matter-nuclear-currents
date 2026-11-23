# General Radiative Dark-Matter Nuclear-Current Theory

A standalone research project for

\[
\chi_1+A\rightarrow\chi_1+A+s,
\]

where a neutral scalar \(s\) is radiated from the dark-matter line and the
target nucleus remains in its elastic ground state.

The project develops reusable three-body response kernels and matches them
to vector, axial, and scalar nuclear currents, including one- and two-body
nuclear operators.

## Research scope

- Exact three-body kinematics and coherent radiation from both dark-line orderings.
- Explicit channels \(VV,VA,AV,AA,SS\), distinguishing dark-matter and nucleon currents.
- Relativistic matching and nonrelativistic nuclear-current reduction.
- AFH multipoles and chiral-EFT one-/two-body currents.
- Reproducible xenon isotope responses, recoil spectra, halo folding, and detector predictions.

## Current status

| Work package | Status | Evidence boundary |
|---|---|---|
| Literature audit and theory architecture | Initial framework established | See docs/LITERATURE_REVIEW.md and docs/GENERAL_2TO3_CURRENT_FRAMEWORK.md |
| VV recovery benchmark | Source-pinned kernel and Native AFH data imported; Fig. 3 reproduction pending | See benchmarks/vv_recovery/ and nuclear/afh/data/native_v6_2026-11-22/ |
| VA / axial nuclear current | Theory task opened | Requires explicit \(G_A,G_P\), PCAC, and xenon axial-response normalization |
| SS / scalar nuclear current | Theory task opened | Requires scalar Wilson-source convention and one-/two-body matching |
| New numerical V/A/S predictions | Not yet performed | No new cross sections or rates are claimed by this migration |

## Start here

1. Research conventions — docs/CONVENTIONS.md
2. Literature review — docs/LITERATURE_REVIEW.md
3. Universal \(u+d+s\) VV current audit — docs/VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md
4. VV source migration and scope — docs/VV_BASELINE_MIGRATION.md
5. General 2-to-3 current framework — docs/GENERAL_2TO3_CURRENT_FRAMEWORK.md
6. Research plan and gates — docs/RESEARCH_PLAN.md
7. External baseline pins — docs/EXTERNAL_BASELINES.md
8. Research and provenance protocol — AGENTS.md

## Relationship to the original project

The original model and historical calculations remain in
[LZ-3-body-for-hdm](https://github.com/lm846802297-bit/LZ-3-body-for-hdm).
This repository imports only selected, immutable source snapshots and Native
AFH data with exact commits, Git blob IDs, SHA-256 hashes, and validation
boundaries. The upstream AFH package itself is not redistributed.

The initial general-current framework note was migrated from original
repository commit de39bf15e08870af82e5217ac30092e0d8e89878. The original
VV/P3 source snapshot is pinned to dd9a7387814cd157d9cde739fa80181e4cd7d428.
See results/provenance/VV_BASELINE_IMPORT_2026-11-23.json.

## Scientific status policy

A derivation, code test, AFH database validation, and complete nuclear-current
prediction are different milestones. No result is promoted to complete
without source identity, input hashes, environment/hardware metadata, raw
outputs, convergence/error analysis, and independent physical validation.
All research changes enter main through reviewed pull requests.
