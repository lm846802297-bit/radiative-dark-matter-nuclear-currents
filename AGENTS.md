# Research and repository operating rules

This repository adopts the shared research standards pinned in [docs/standards-lock.yaml](docs/standards-lock.yaml). Start substantive work from the [central standards index](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/CENTRAL_STANDARDS_INDEX.md). Shared governance, compute, task handoff, quality gates, archival, RCRPS, document review, and Git integration are governed by the pinned central documents; this file retains the nuclear-current-specific physics requirements.

Before substantive research, derivation, numerical, code, or validation work, follow the [central Research Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and use the [central Research Task Brief template](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/templates/research-task-brief.md). Trivial editorial/maintenance changes may omit a brief when the PR explains why it is not applicable.

## Scientific scope

The fixed process is dark-matter-line scalar radiation, \(\chi_1+A\to\chi_1+A+s\), with the target in its elastic ground state. Use two-sided channel labels: VV, VA, AV, AA, SS. The first letter labels the dark-matter transition bilinear; the second labels the nucleon-side current. Do not use an unqualified “A channel”.

## Shared provenance, compute, and archival requirements

Follow the pinned central [Research Compute Policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/RESEARCH_COMPUTE_POLICY.md), [Calculation-result Archival Policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/CALCULATION_RESULT_ARCHIVAL_POLICY.md), and [RCRPS v1.0.1](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/STRUCTURED_RESULT_DATA_STANDARD_v1.0.1.md).

The local [nuclear-current RCRPS extension](docs/NUCLEAR_CURRENT_RCRPS_EXTENSION.md) defines the project-specific `extensions.nuclear_currents` fields. The existing [repository provenance validator](tools/validate_provenance.py) remains responsible for checking the top-level `results/provenance/SHA256SUMS.txt` inventory; it is separate from RCRPS record validation.

The existing Native AFH manifest is a project-specific historical record. Preserve it and its exact source/data fingerprints; do not retrofit or overwrite it as an RCRPS record.

## Physical validation gates

As applicable, test direct covariant amplitudes against tensor contractions; both radiation diagrams and their coherent interference; Lorentz-index and metric conventions; exact phase-space boundaries; vector Ward/continuity identities only when justified by the selected Hamiltonian/current; axial PCAC and induced-pseudoscalar relations; scalar-source normalization and trace anomaly; AFH multipole reconstruction, signed interference, positivity and Cauchy-Schwarz bounds; analytic limits; independent quadrature and precision refinements.

A passing unit test proves only the tested software property. It does not prove a physical model complete. Apply the central [Research Quality Gate](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/RESEARCH_QUALITY_GATE.md) to classify the task-level evidence.

## Claims and status language

Follow the central [Research Governance Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_GOVERNANCE_STANDARD.md). Use “derived”, “implemented”, “tested”, “numerically converged”, and “independently validated” only for the exact scope supported by archived evidence. Never promote a leading one-body result to a full chiral-EFT current result. Never label a response-level contraction as an amplitude-level nuclear validation.
