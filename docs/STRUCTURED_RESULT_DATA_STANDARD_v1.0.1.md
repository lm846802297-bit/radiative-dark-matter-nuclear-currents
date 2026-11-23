# RCRPS v1.0.1 — nuclear-current project extension

**Normative common standard:** [Central RCRPS v1.0.1](https://github.com/lm846802297-bit/high-energy-theory-research/blob/db503f0eb579bdbcac702e44237a6610c34690e7/docs/STRUCTURED_RESULT_DATA_STANDARD_v1.0.1.md)  
**Pinned central revision:** `db503f0eb579bdbcac702e44237a6610c34690e7`  
**Scope:** New substantial numerical runs in this repository. Historical migration and AFH records remain unchanged.

This file defines only the project-specific extension. It does not redefine common RCRPS fields or replace the central schema.

## 1. Required record bundle

For each new substantial numerical run, use the central record roles:

- `run_manifest.json`: what actually ran, exact source and input identities, configuration, environment, command, outputs, checksums, and archive status;
- `result_summary.json`: benchmark, scope, numeric observables with separate units, scientific status, and interpretation;
- `validation_report.json`: the analytic, numerical, convergence, independent, physical, and integrity checks actually performed;
- native unrounded CSV/JSON/HDF5/ROOT or other directly usable outputs;
- `SHA256SUMS.txt` for the archived native outputs where appropriate;
- a concise `README.md` with reproduction instructions and status.

A record must describe the execution that actually occurred. Do not create a successful run manifest before execution, fill unknown environment fields from assumptions, or mark unrun checks as passed.

## 2. Project extension namespace

Store project-specific metadata under:

```json
{
  "extensions": {
    "radiative_dark_matter_nuclear_currents": {
      "process": "chi1 + A -> chi1 + A + s",
      "channel": "VV | VA | AV | AA | SS",
      "dark_sector": {},
      "nucleon_current": {},
      "nuclear_current": {},
      "kinematics": {},
      "nuclear_structure": {},
      "numerics": {}
    }
  }
}
```

This is a namespace example, not a complete valid run manifest. Populate only fields relevant to the actual run.

### 2.1 Physics identity

Record, as applicable:

- process and two-sided channel label (VV, VA, AV, AA, SS); the first letter identifies the dark-matter transition bilinear and the second the nucleon-side current;
- particle nature, masses, couplings, mediator and scalar parameters, propagator/width prescription, and coupling normalization;
- exact amplitude/source implementation and whether both radiation orderings and their interference are included;
- current type and matching level: relativistic, nonrelativistic, one-body, two-body, or combined;
- vector/axial/scalar form-factor inputs, source references, parameterization, momentum argument, and normalization;
- (q_A), (q_N), (t_A), (t_N), (E_R), frame and sign conventions as separate quantities where applicable.

### 2.2 Nuclear structure identity

Record, as applicable:

- isotope and isotope-specific nuclear mass;
- nuclear Hamiltonian/wave-function source and version;
- AFH/DMFormFactor version, OBDME source and hashes, oscillator-length prescription;
- nuclear-current operator, one-/two-body truncation, chiral order and regulator;
- response names, multipole/rank range, isospin indices, signed interference convention, and normalization;
- external nuclear input provenance and redistribution/access status.

### 2.3 Numerical and environment identity

Record the actual:

- execution target, OS, architecture, CPU, core count, RAM, GPU/driver if used;
- Python, Mathematica/Wolfram, compiler, BLAS/LAPACK and dependency versions;
- precision, quadrature variables/orders, grids, interpolation, tolerances, cutoffs, random seeds, and exact command;
- input file paths/immutable URIs, byte counts, versions, and SHA-256;
- output paths, formats, byte counts, SHA-256, and CSV row/column/unit metadata.

## 3. Scientific validation

Select checks from the actual calculation. Record methods, observed residuals, thresholds, units, and evidence paths separately. Depending on the claim, relevant checks include:

- direct covariant amplitude versus tensor contraction;
- coherent sum of both radiation diagrams and interference;
- exact three-body kinematics, on-shell conditions, four-momentum conservation, and phase-space Jacobian;
- vector Ward/continuity identities only for a conserved current and consistent Hamiltonian;
- axial PCAC/induced-pseudoscalar relations;
- scalar-source normalization and trace-anomaly matching;
- AFH multipole reconstruction, signed interference, positivity and Cauchy-Schwarz bounds;
- precision, quadrature, phase-space, recoil-grid, interpolation, and nuclear many-body convergence;
- independent implementation or representation.

Keep the statuses distinct: computed, archived, numerically validated, and physically validated. A software test or a matching checksum does not by itself establish a physical result.

## 4. Historical records

The existing `nuclear/afh/data/native_v6_2026-11-22/run_manifest.json`, initialization/migration manifests under `results/provenance/`, and their checksum inventories are historical records with their own established formats. Preserve them byte-for-byte. Do not retroactively convert them or claim they satisfy RCRPS v1.0.1.

New substantial numerical runs use RCRPS v1.0.1. If an old record is later supplemented, add a separately identified addendum rather than silently rewriting its original bytes.
