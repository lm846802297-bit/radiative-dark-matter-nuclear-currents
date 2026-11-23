# Provenance archive

Every accepted numerical run gets a dated directory with a run manifest, SHA-256 list, exact source commit, clean/dirty state, input hashes, environment/hardware, physics and numerical parameters, native unrounded outputs, convergence/error tables, independent physical checks, and a concise report.

The initialization manifest records the source commit and file hash of the imported theory note. It is a repository-migration record, not a numerical physics result. Never infer SHA-256 from a Git SHA-1 blob identifier.

## Theory and literature audits

A source-matching or literature audit is not a numerical production run. Record the exact papers/versions, equations or tables used, source-project commit pins, imported document SHA-256 values, and explicit claims that were not established. Do not invent hardware/runtime fields for a documentation-only review. Numerical production still requires the full run manifest and native outputs described above.

The universal-vector VV audit is recorded in `VV_VECTOR_CURRENT_ADDENDUM_2026-11-23.json`; its readable theory note is `../../docs/VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md`.

