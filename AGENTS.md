# Research and repository operating rules

These rules apply to every calculation and manuscript claim in this repository.

## Scientific scope

The fixed process is dark-matter-line scalar radiation, \(\chi_1+A\to\chi_1+A+s\), with the target in its elastic ground state. Use two-sided channel labels: VV, VA, AV, AA, SS. The first letter labels the dark-matter transition bilinear; the second labels the nucleon-side current. Do not use an unqualified “A channel”.

## Required provenance for numerical work

Every production calculation must archive:

1. Exact source Git commit SHA and clean/dirty working-tree state.
2. SHA-256 for every input file, external data table, parameter file, and generated source snapshot.
3. Actual software environment: OS, architecture, Python/Mathematica/compiler versions, package versions, BLAS/LAPACK implementation, and relevant environment variables.
4. Actual hardware: CPU model/core count, RAM, GPU model/driver when used, and accelerator settings.
5. Full computational parameters: units, conventions, precision, quadrature rules/nodes, tolerances, interpolation method, cutoffs, and width prescription.
6. Native, unrounded calculation outputs. Plots and tables are derived artifacts, not substitutes for raw results.
7. Numerical convergence and error analysis, separating phase-space, velocity, recoil-grid, interpolation, precision, chiral-truncation, and nuclear many-body uncertainties.
8. Independent physical validation, preferably using a distinct implementation or representation.
9. A run manifest, checksum list, concise result report, and reproducible command line.

A missing item is marked pending; it must not be inferred or fabricated.

## Physical validation gates

As applicable, test direct covariant amplitudes against tensor contractions; both radiation diagrams and their coherent interference; Lorentz-index and metric conventions; exact phase-space boundaries; vector Ward/continuity identities only when justified by the selected Hamiltonian/current; axial PCAC and induced-pseudoscalar relations; scalar-source normalization and trace anomaly; AFH multipole reconstruction, signed interference, positivity and Cauchy-Schwarz bounds; analytic limits; independent quadrature and precision refinements.

A passing unit test proves only the tested software property. It does not prove a physical model complete.

## Local-first computation

Run scientific production calculations in the active local research environment first. GitHub Actions is an independent reproduction/regression check, not the default production computer. Record actual environment and hardware.

## Immutable outputs

Keep accepted raw outputs immutable. Store new analyses in dated result directories. Record parent source hashes and derived-file hashes. Never overwrite an accepted baseline to make a new model appear consistent. Keep provisional, validated, and production results distinct.

## Git workflow

Work on a descriptive research branch. Do not commit research changes directly to main. Open a pull request with scope, exact base/head commits, validation evidence, and open gates. Review the diff and run relevant local tests before merging. Do not bundle unrelated projects.

## Claims and status language

Use “derived”, “implemented”, “tested”, “numerically converged”, and “independently validated” only for the exact scope supported by archived evidence. Never promote a leading one-body result to a full chiral-EFT current result. Never label a response-level contraction as an amplitude-level nuclear validation.
