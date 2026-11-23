# Task-level standards provenance

For each substantive research task, store a machine-readable record at:

`docs/standards-provenance/tasks/<task-id>.json`

Use the central standards-adoption schema copied to `schemas/standards-adoption-v1.0.0.schema.json`. The record must bind the exact `docs/standards-lock.yaml` and finalized task brief bytes by SHA-256, identify the exact pinned central commit, list only the profiles actually applied, reference project-local extensions, and state the predeclared acceptance criteria.

For each numerical run, repeat the task-level provenance under `run_manifest.extensions.standards_provenance`. The standards-adoption workflow validates task records and run manifests when present.

The lock file records project-wide adoption; it does not prove that a particular task applied those standards. Do not create a task/run record before the corresponding task or execution exists, and do not use placeholder hashes.
