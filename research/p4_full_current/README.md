# P4-FC production archive

This directory contains the curated, immutable P4-FC Stage A/B/C1/C2 production package imported from `lm846802297-bit/LZ-3-body-for-hdm` at source commit `14ad781ee15db5c043ba29aefa5b23e48a301ce2`.

## Included material

- **Stage A:** native seven-isotope, 81-recoil-point response, multipole, and OBDME-term exports, with run metadata, provenance, and validation.
- **Stage B:** finite-`t_N` direct Dirac–Pauli current outputs, convergence records, run manifest, and validation.
- **Stage C1:** the accepted full Lorentz-tensor grid, tensor index, convergence report, run manifest, and PASS report.
- **Stage C2:** the production one-body current contraction, isotope/halo contributions, natural-xenon rates, validation samples, run manifest, and 12/12 unit-test record.
- **Independent checks:** AFH A13 coefficient validation and continuity-projection stress/replay records.
- **Implementation:** source scripts, audit scripts, task brief, and regression tests.

## Scientific status boundary

The production record states:

- `stage_c2_complete=true`
- `one_body_impulse_current_complete=true`
- `P3_tensor_validation_complete=true`
- `independent_nuclear_amplitude_validation_complete=false`
- `full_current_complete=false`

The result is the completed **one-body impulse-current contraction**. It is not a completed exchange-current or nuclear-Hamiltonian calculation. No detector efficiency or exposure is applied, and no LZ expected-count prediction is claimed.

The C2 production archive covers 7 isotopes × 81 recoil points × 2 halo models × 24 speed nodes (27,216 contractions), 1,134 isotope/recoil/halo points, and 162 natural-xenon rate rows. The archived C2 unit-test suite passed 12/12.

For final C2 production status, use `output/stage-c2-full-current-contraction/C2_PASS_REPORT.md` and the final 2026-11-26 disposition at the end of `C2_IMPLEMENTATION_STATUS.md`. Earlier dated status paragraphs in that implementation-history file are historical snapshots, not the final production disposition.

## Integrity and provenance

- The 78 imported source artifacts have exact source/target Git blob SHA matches and exact byte-size matches.
- [`IMPORT_MANIFEST.json`](IMPORT_MANIFEST.json) records source and target paths, source and target Git blob IDs, byte sizes, SHA-256 values, stages, dispositions, and redistribution status.
- [`SHA256SUMS.txt`](SHA256SUMS.txt) contains exactly 78 SHA-256 entries, one for each imported source artifact.
- Original stage-level checksum files and run manifests are preserved unchanged.
- The four previously transferred large files retain their locally verified SHA-256 values, matching the source-stage checksum records.

## Exclusions and transfer boundary

The superseded C1 `nw8-fail` trial is retained as provenance only; its failure report, checksum record, log, convergence report, and run manifest are kept, while its failed data payloads are excluded. The separate `pilot_248` archive is not part of this production import. See [`FAILURE_TRIAL_EXCLUSIONS.md`](FAILURE_TRIAL_EXCLUSIONS.md) for exact excluded paths and hashes.

No nuclear calculation was run or rerun during migration. Native outputs, source files, run metadata, and validation records were transferred as existing Git objects.

## Redistribution status

No explicit LICENSE, COPYING, or NOTICE file was found at the pinned source repository root. This archive transfer is between repositories controlled by the same owner and does not grant an open-source reuse license.
