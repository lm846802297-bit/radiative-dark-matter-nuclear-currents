# VV baseline migration and PR decision

## Decision

The original-project PR #67 is **not merged wholesale**. At the initial migration
audit it remained open with 200 changed files and targeted
`research/p3b-p3c-p4-gates-20261119`, not `main`. Its scope combines the
P3/P4 integration branch, native AFH work, recoil scans, interpolation
audits, and exploratory Gram eigenmode analysis. At the latest migration review, PR #67 had advanced to head `4e0d6628d6446ee3caf91c779ef59e80145c595c` with 224 changed files; the P3-C native cross-archive and Gram-eigenmode checks were failing while other P4 jobs remained pending or in progress. Those changes should not be treated as one indivisible contribution to the new general-current project.

Instead, this PR imports a curated, source-pinned subset. The original
repository and PR remain available as the historical calculation archive.
The new repository owns the paper-facing VV recovery benchmark and the
future reusable V/A/S current framework.

## Exact source pins

### Core VV/P3 source snapshot

- Repository: `lm846802297-bit/LZ-3-body-for-hdm`
- Commit: `dd9a7387814cd157d9cde739fa80181e4cd7d428`
- Branch at source time: `research/p4-provisional-halo-fold-20261119`
- Imported code: 15 files under
  `benchmarks/vv_recovery/source_snapshot/research/`.
- The source files are byte-identical to their source Git blobs. Both Git
  blob IDs and SHA-256 hashes are retained.

### Paper Fig. 3 benchmark

- Source commit: `5ea4d360fe03ce58079220bb0b7624b3faf03fb8`
- Source path:
  `research/three_body_full_current_nucleon_match/benchmarks/paper_fig3_vector_natural_xenon.json`
- Source Git blob SHA-1: `af3731952da389872de21358050750235f9a688d`
- Imported SHA-256:
  `940fad45c3d408a94ed755adc29e11e966d899317faa49956527e221e2035308`.

### Native AFH data

- AFH package: DMFormFactor v6.0.
- Package SHA-256:
  `dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372`.
- Native Wolfram version: 14.2.0 for Windows.
- Native response coverage: 7 isotopes x 2,401 momentum points = 16,807 rows.
- Responses: all eight AFH v6.0 channels.
- Native spatial multipoles: 9,604 rank-resolved Xe-129/Xe-131 rows.
- Maximum native imaginary residual: zero.
- Pointwise comparison: \(r_{\rm tol}=5\times10^{-7}\),
  \(a_{\rm tol}=2\times10^{-12}\); maximum combined-tolerance utilization
  0.05879.

The AFH package source itself is not copied. The original Native output,
metadata, comparison, and WLS exporter sources are imported. The package
must be supplied separately with the pinned hash.

## What is imported

1. The original Eq. (6) three-body kinematic kernel and spin-resolved
   two-diagram dark-matter current.
2. The four-vector current audit and phase-space-integrated Gram tensor.
3. Component-wise quadrature convergence code and tests.
4. Universal \(u+d+s\) vector-nucleon form-factor and AFH Eq. (107)
   particle-side matching code.
5. The fixed-speed AFH contraction source.
6. Native AFH eight-response tables for all seven xenon isotopes and
   Native rank-resolved spatial multipoles for Xe-129/Xe-131.
7. P3-C bound-current and spatial-helicity matching notes as reference
   documents, with their original incomplete-status flags preserved.

## What is deliberately not promoted

- The imported AFH response validation is a nuclear-input validation,
  not a complete particle--nuclear current validation.
- The P3-C report retains
  `physically_validated_full_current=false` and
  `full_current_complete=false`.
- The source repository's P4 81-point halo spectrum is a different
  secondary benchmark, a leading-static one-body candidate, and does not
  include the paper's LZ efficiency or 2.84 tonne-year exposure. It is not
  labelled as a Fig. 3 reproduction.
- The 121-point dense-q spatial interpolation/native audit and the
  Gram-eigenmode study are not part of this primary migration. They remain
  available in the source repository for later targeted use.

## Immediate research objective

The next scientific task is not another eigenmode or interpolation study.
It is an independent reproduction of the paper's Fig. 3 VV event spectrum:

1. Reproduce the exact two-diagram \(2\to3\) amplitude and phase-space
   integral at the frozen paper benchmark.
2. Rebuild the \(2\to2\) comparison with its stated integrated-rate
   normalization.
3. Rebuild the loop-induced elastic curve using the paper's \(F_3\).
4. Import the paper's LZ nuclear-recoil efficiency and apply the
   2.84 tonne-year exposure.
5. Compare the three event spectra against the published figure and
   archive residuals, figures, and provenance.

The imported code/data are a source baseline, not a claim that these
five reproduction steps have already passed.


## Provenance redaction

The imported Native AFH metadata and Wolfram log preserve the original source-file hashes in `SOURCE_SHA256SUMS.txt`. Absolute local Windows workspace/output paths have been redacted from the public copies; the AFH package SHA-256, Wolfram version, hardware, isotope coverage, grid, normalization, and numerical validation data are retained. `SHA256SUMS.txt` verifies the public, sanitized files.
