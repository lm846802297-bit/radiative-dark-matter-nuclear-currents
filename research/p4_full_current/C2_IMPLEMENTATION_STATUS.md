# C2 implementation status

## Conclusion

The C2 contraction runner has now been executed on the full pinned C1 grid
and passed its production acceptance. The implementation preserves the frozen
C1 and Stage B inputs and does not perform any phase-space or Stage B
recomputation.

## Completed in this preparation pass

- Generalized multipole matrix construction to use the actual target spin `J_A`, including the `J_A=0`, `K=0`, `L=1` longitudinal current.
- Implemented the full four-component Lorentz contraction with charge/spatial and longitudinal/transverse interference bookkeeping.
- Preserved the signed before/after/diagram-interference decomposition.
- Encoded the physical normalization `4*C0(A)*K_raw`; the raw AFH matching check remains `C0(A)*K_raw`.
- Encoded the frozen B benchmark and the natural-xenon speed-node fold.
- Added exact upstream SHA-256 checks and a clean-Git-tree/external-output requirement.
- Added independent explicit-loop versus vectorized contraction tests and target-spin tests.

## Validation of the original C2 implementation

- Python syntax compilation: PASS.
- Original synthetic unit tests: 6/6 PASS.
- These results predate the continuity-projection integration below; they are not a syntax/test pass for the updated C2 runner. The original preparation-only JSON record is preserved as [C2_PREPARATION_UNIT_TEST_REPORT.json](C2_PREPARATION_UNIT_TEST_REPORT.json); the production test result is recorded separately in [C2_UNIT_TEST_REPORT.json](C2_UNIT_TEST_REPORT.json) and the full output log at `output/stage-c2-full-current-contraction/unit_test_output.txt`.

## Production input location and blocker (2026-11-25)

The frozen C1 production archive has now been located. It is committed on the repository branch:

- Branch: `research/p4-fc-stage-c1-full-lorentz-grid-20261124`
- Source commit: `6ff1f357539adde80392f6aa587d695949050edb`
- C1 status: `PASS`, `stage_c1_complete=true`
- Production coverage: 7 isotopes x 81 recoil points x 2 halo models x 24 speed nodes = 27,216 open tensor rows.
- C1 production quadrature: (14,6,4); comparison quadrature: (16,14,12).
- Native run: Windows local research environment, Python 3.13.18, NumPy 2.3.5, SciPy 1.17.0, 8 workers, wall time 825.30 s. No GitHub Actions runner was used.

The exact archive directory is:

`research/three_body_full_current_nucleon_match/p4_full_current/output/stage-c1-full-lorentz-grid/`

Pinned file identities:

| File | Size (bytes) | SHA-256 |
|---|---:|---|
| `full_lorentz_tensor_grid.npz` | 16,067,358 | `98898fbdac7a08a81453972d013ad027e46679a02e67d41a3d59a529f1671c76` |
| `tensor_index.csv` | 7,417,583 | `3647a0b4261be85b146ffe040f4e0241fd1508d85cc0d5f1f9b502dc1b33c06c` |
| `run_manifest.json` | 17,097 | `18d36a643205ad0f8cde581cb04f1ff6dbd3f2aa6e86e1877837b6b8678fc8e8` |
| `convergence_report.json` | 10,344 | `1cb9df212451374652e0713a477f86e8954e3f9d8aee9aa71123e2910b3a0878` |

The C1 archive is not present on `main`, and the inspected C1 workflow attempts produced no downloadable workflow artifact. The GitHub connector can inspect the private branch, file sizes, manifests, and checksums, but the current ChatGPT Linux runtime cannot materialize the raw binary NPZ/index from the connector. Consequently, the production C1 tensors are present in the repository but are not yet available as local files for C2.

Stage A native responses and Stage B production CSV/validation report are already archived on `main` and remain reusable read-only. Stage B production inputs are not to be recomputed.

## Continuity-projection update (2026-11-25)

The C2 runner now compares the AFH six-response basis against the continuity-projected leading-NR current, while retaining the actual longitudinal current and its analytic nonconserved residual in the full contraction. Three focused regression tests have been added for conserved-current reduction, the nonconserved residual formula, and full-kernel reconstruction.

The standalone algebraic stress test is `validation/continuity_projection/test_c2_continuity_projection_stress.py` (SHA-256 `17a793d74df9562c62d6da5ffcd95fbc95bfa846da177015078d8937b9bbf67f`). The original archived run passed 4/4 tests, including 20,000 random Hermitian-PSD tensor/current samples.

An independent local replay has also passed 4/4 in the current ChatGPT Linux runtime; its exact environment, output, and checksums are submitted in [PR #179](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/pull/179). That replay is an algebraic software test, not a production C2 result. The three new tests in the C2 runner test module are added but have not yet been run against the full production checkout.

## State boundary

`stage_c2_complete`, `P3_tensor_validation_complete`, and `one_body_impulse_current_complete` remain **NOT CLAIMED** until the runner is executed against the pinned production C1 inputs and all acceptance checks pass. The generic `full_current_complete` flag remains false pending a consistent nuclear Hamiltonian and exchange-current closure. No C2 production output or production manifest is claimed.


## Full Proca numerator implementation update (2026-11-25)

The earlier C2 implementation-status record correctly identified that the
runner did not yet apply the Proca numerator. A C2 v2 source update now
addresses that concrete issue.

For \(q_A^\mu=(E_R,0,0,+Q_A)\), metric \((+---)\), the implementation forms

\[
M^\mu{}_\nu=\delta^\mu{}_\nu-\frac{q_A^\mu q_{A,\nu}}{m_V^2},
\qquad
L_{\rm Proca}=M^T L M.
\]

The Proca-dressed tensor is now used consistently for:
- exact and leading-NR four-current contractions;
- continuity-projected AFH comparisons and nonconserved residuals;
- before/after scalar-emission diagrams and their signed interference;
- AFH reference coefficients and deterministic explicit-loop comparisons.

The scalar denominator \((m_V^2-t_A)^{-2}\) remains a separate rate factor
and is not duplicated in the numerator. The C2 source records Proca-dressed
Hermiticity, PSD, and diagram-decomposition diagnostics.

Three new synthetic tests check the direct covariant-amplitude identity,
zero-transfer limit, and preservation of the before/after/interference tensor
sum. A local algebraic smoke check on 1,000 random complex current pairs gave
a maximum relative residual \(3.37\times10^{-15}\). This is not a full
repository test-suite result.

**Current disposition (2026-11-26):** Proca numerator implementation is
present on the C2 production branch. The complete C2 unit-test suite and the
full production grid have now been executed on a clean checkout; the run
passed and is archived under
`output/stage-c2-full-current-contraction/`. The frozen C1 archive and its
SHA-256 pins remain unchanged; no C1, Stage A, or Stage B calculation was
repeated. `stage_c2_complete`, `one_body_impulse_current_complete`, and
`P3_tensor_validation_complete` are true; `full_current_complete` remains
false pending a consistent nuclear Hamiltonian and exchange-current closure.
