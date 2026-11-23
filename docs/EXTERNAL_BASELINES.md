# External baselines and source pins

This repository is independent of the original model repository. The following immutable source revisions define the initial reuse boundary.

## Original VV/P3/P4 baseline

- Repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Branch at the time of pinning: research/p4-provisional-halo-fold-20261119
- Exact commit: dd9a7387814cd157d9cde739fa80181e4cd7d428
- Role: source anchor for the original radiative VV benchmark, validated AFH v6.0 inputs, and the P3/P4 implementation.
- Recorded scope at that commit: native seven-isotope/eight-response AFH production and P4 halo-fold outputs are validated; the full chiral-EFT nuclear current remains incomplete.
- Import status: no raw AFH tables, native Wolfram outputs, source code, or numerical spectra have been copied into this repository yet.

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
