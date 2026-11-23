# Provenance archive

Every accepted numerical run gets a dated directory with a run manifest, SHA-256 list, exact source commit, clean/dirty state, input hashes, environment/hardware, physics and numerical parameters, native unrounded outputs, convergence/error tables, independent physical checks, and a concise report.

The initialization manifest records the source commit and file hash of the imported theory note. It is a repository-migration record, not a numerical physics result. Never infer SHA-256 from a Git SHA-1 blob identifier.
