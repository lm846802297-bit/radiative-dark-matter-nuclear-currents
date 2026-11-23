# Nuclear-Current Domain Extension to RCRPS v1.0.1

**Status:** Project-specific extension to the centrally maintained RCRPS v1.0.1 standard.  
**Central source commit:** 142267bc803c90aa681110293b66b73947bbcaa1

The common JSON record roles, source-control provenance, input/output integrity, CSV contract, status semantics, and validation requirements are defined by the [central RCRPS v1.0.1 standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/STRUCTURED_RESULT_DATA_STANDARD_v1.0.1.md). This file defines only the reusable dark-matter nuclear-current metadata required by this project.

The local v1.0.0/v1.0.1 schemas, validator, dependency pin, and regression tests are synchronized byte-for-byte with the central files at the pinned source commit. Update them only through a reviewed downstream PR that also updates `docs/standards-lock.yaml`.

## 1. Local implementation

- Historical v1.0.0 schema: `docs/schemas/research-result-record-v1.0.0.schema.json`
- Adopted v1.0.1 schema: `docs/schemas/research-result-record-v1.0.1.schema.json`
- Validator: `docs/schemas/validate_result_records.py`
- Regression tests: `docs/schemas/tests/`
- Pinned dependency: `docs/schemas/requirements-result-validation.txt`

Install the pinned dependency and run the local regression suite:

```bash
python -m pip install -r docs/schemas/requirements-result-validation.txt
python -m unittest discover -s docs/schemas/tests -v
python -m compileall -q docs/schemas
```

Validate a v1.0.1 run manifest with the explicit schema:

```bash
python docs/schemas/validate_result_records.py --schema docs/schemas/research-result-record-v1.0.1.schema.json --repo-root . <run_manifest.json> <result_summary.json> <validation_report.json>
```

The v1.0.0 schema remains available only for historical compatibility. New substantial runs use v1.0.1 and must record the complete `execution.source_control` object required by the central schema. Do not rewrite historical custom manifests or fabricate missing provenance.

## 2. Nuclear-current domain metadata

Place project-specific fields under `extensions.nuclear_currents`. Do not add them to common top-level RCRPS fields.

### 2.1 Process and channel identity

Record:

- the process `chi_1 + A -> chi_1 + A + s`, including the target state and whether the nuclear transition is elastic;
- the two-sided channel label `VV`, `VA`, `AV`, `AA`, or `SS`;
- the first channel letter as the dark-matter transition bilinear and the second as the nucleon-side current;
- particle nature, scalar/mediator content, coupling conventions, widths, and benchmark parameters.

Do not use an unqualified `A channel`.

### 2.2 Kinematics and conventions

Define each quantity independently, including its frame, direction, external legs, and units. As applicable, record:

- `q_A`: spatial momentum transferred to the nucleus;
- `q_N`: momentum transferred to the nucleon;
- `q_chi`: dark-matter-line momentum transfer;
- `k_s`: emitted-scalar four-momentum;
- `t_A` and `t_N`: separately defined Lorentz invariants;
- `E_R`: nuclear recoil energy;
- phase-space boundaries, angular variables, reference frames, and integration variables.

Never substitute a spatial momentum magnitude for a Lorentz invariant without an explicit derivation.

### 2.3 Nucleon and nuclear structure

Record the isotope, nuclear mass convention, abundance, nuclear state, source/version of the nuclear matrix elements, density-matrix provenance, response channels, and normalization.

For AFH inputs, record:

- AFH/DMFormFactor version and source SHA-256;
- isotope and density-matrix source/version/hash;
- the spatial argument (for example, (y_A=(b_A q_A/2)^2));
- response names and multipole ranks;
- raw package normalization and any downstream conversion, including the equal-proton/equal-neutron factor-four conversion where applicable;
- signed interference terms and whether the output is a raw response, reconstructed amplitude, or contracted observable.

Keep raw AFH package responses distinct from physically rescaled nuclear responses and from point-nucleus approximations.

### 2.4 Chiral-EFT and one-/two-body currents

For chiral-EFT calculations, record:

- chiral order and truncation;
- regulator, cutoff, and low-energy constants;
- one-body and two-body operator content;
- matching and renormalization conventions;
- nucleon form factors and their momentum arguments;
- current conservation/PCAC assumptions and the Hamiltonian/current pair to which they apply;
- omitted operators and the exact claim scope.

### 2.5 Observable and validation level

Distinguish the calculation level explicitly: amplitude, spin-summed squared amplitude, nuclear-response contraction, (d\sigma_A/dE_R), fixed-speed rate, halo-folded rate, natural-isotope rate, detector-response prediction, or likelihood-level observable.

Keep phase-space, velocity, recoil-grid, interpolation, precision, chiral-truncation, and nuclear many-body uncertainties separate. The status of a nuclear-response input does not automatically promote the complete particle–nuclear current or event prediction to a validated result.

## 3. Historical AFH manifest

The existing `native-afh-eight-response-run-v1` manifest under `nuclear/afh/data/native_v6_2026-11-22/` is a project-specific historical record, not an RCRPS v1.0.0/v1.0.1 record. Its schema name, source provenance, native Wolfram metadata, response comparison, and validation boundary remain unchanged.

Do not rewrite that manifest to make it appear to conform to RCRPS. A future RCRPS sidecar may reference the historical manifest and native outputs only when the exact source, inputs, outputs, and checksums are verifiable.
