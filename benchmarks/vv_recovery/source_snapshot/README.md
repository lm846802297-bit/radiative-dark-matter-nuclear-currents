# Immutable VV/P3 source snapshot

This directory preserves a curated, byte-identical source snapshot from
\`lm846802297-bit/LZ-3-body-for-hdm\` at commit
\`dd9a7387814cd157d9cde739fa80181e4cd7d428\`.

The original relative directory layout is retained so the historical Python
modules keep their local import and path assumptions. The snapshot contains
the Eq. (6) kinematic kernel, spin-resolved two-diagram dark-matter current,
four-vector response audit, phase-space-integrated Gram tensor, convergence
audit, one-body Dirac--Pauli form factors, AFH particle-side coefficients,
and fixed-speed response contraction.

No source file in this snapshot was modified during migration. Git blob IDs
and SHA-256 hashes are recorded in \`SHA256SUMS.txt\` and the migration
provenance record.

The AFH \`DMFormFactor v6.0\` package itself is not copied here. Its exact
SHA-256 is pinned in the Native AFH dataset manifest. Reproduction requires
the upstream package to be installed separately and supplied through the
documented \`DMFormFactorPath\` / environment interface.

This is an immutable historical reference, not yet the new project's
production API. New reusable modules should be implemented under \`src/\`
with tests and explicit conventions, while this snapshot remains unchanged.

The imported one-body Dirac--Pauli current and AFH Eq. (107) coefficients
are a leading one-body impulse baseline. They do not constitute a complete
chiral-EFT nuclear current or a validated full-current event prediction.
