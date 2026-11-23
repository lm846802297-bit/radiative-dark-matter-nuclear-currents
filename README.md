# General Radiative Dark-Matter Nuclear-Current Theory

A standalone research project for

\[
\chi_1+A\rightarrow\chi_1+A+s,
\]

where a neutral scalar \(s\) is radiated from the dark-matter line and the target nucleus remains in its elastic ground state.

The project develops reusable three-body response kernels and matches them to vector, axial, and scalar nuclear currents, including one- and two-body nuclear operators.

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
| VV recovery benchmark | Reusable external baseline | Existing calculations remain in the original repository; no numerical data have yet been copied here |
| VA / axial nuclear current | Theory task opened | Requires explicit \(G_A,G_P\), PCAC, and xenon axial-response normalization |
| SS / scalar nuclear current | Theory task opened | Requires scalar Wilson-source convention and one-/two-body matching |
| New numerical V/A/S predictions | Not yet performed | No new cross sections or rates are claimed by this bootstrap |

## Start here

1. Research conventions — docs/CONVENTIONS.md
2. Literature review — docs/LITERATURE_REVIEW.md
3. General 2-to-3 current framework — docs/GENERAL_2TO3_CURRENT_FRAMEWORK.md
4. Research plan and gates — docs/RESEARCH_PLAN.md
5. External baseline pins — docs/EXTERNAL_BASELINES.md
6. Research and provenance protocol — AGENTS.md

## Relationship to the original project

The original model and validated VV calculations remain in [LZ-3-body-for-hdm](https://github.com/lm846802297-bit/LZ-3-body-for-hdm). This repository is independent. It will reuse only explicitly pinned source revisions and data snapshots with recorded SHA-256 hashes; it will not track a moving branch of the original repository.

The initial framework note was migrated from original repository commit de39bf15e08870af82e5217ac30092e0d8e89878. The migration record and file hashes are in results/provenance/.

## Scientific status policy

A derivation, code test, AFH database validation, and complete nuclear-current prediction are different milestones. No result is promoted to complete without source identity, input hashes, environment metadata, native outputs, convergence/error analysis, and independent physical validation. All research changes enter main through reviewed pull requests.
