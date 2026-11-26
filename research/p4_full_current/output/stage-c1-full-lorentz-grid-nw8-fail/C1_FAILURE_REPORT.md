# Stage C1 production run: FAIL_CONVERGENCE

Status: **FAIL_CONVERGENCE**, not C1 PASS.

> **Archive note**: this is the superseded \((8,6,4)\) trial. The C1 production
> specification was subsequently revised to \((14,6,4)\) with a
> \((16,14,12)\) comparison; the canonical PASS archive lives in
> `../stage-c1-full-lorentz-grid/`. This directory is retained as the
> failed-preflight provenance record and is not used for C2.

This is the actual local production run of the frozen C1 grid on commit
`e3cb7b81339b16336dc203182192623123a81ab0` (branch
`research/p4-fc-stage-c1-full-lorentz-grid-20261124`). It was executed on
Windows 11 with Python 3.13.18, NumPy 2.3.5, SciPy 1.17.0, 8 workers, and
Wolfram was not involved. The run used no GitHub Actions runner.

## What passed

- Coverage: 7 isotopes, 81 recoil points, 2 halo domains, 24 speed nodes:
  1,134 isotope/recoil/halo points and 27,216 tensor rows; all speed intervals
  open.
- Every stored `L` is finite and `4x4`; Hermiticity residual `0.0`; negative
  PSD eigenvalue relative measure `0.0`; Gram partial trace `0.0`; before +
  after + signed interference decomposition `3.778e-16`; q-frame rotation
  residual `1.018e-15`.
- Speed convergence `24 -> 48`: `9.1975603181786e-06`, below `1e-5`.
- Independent Dirac-trace residual: `6.875444658049901e-11`, below `1e-8`.
- Runtime: `470.11 s`; `run_manifest.json` records the clean source tree and
  all input/output SHA-256 values.

## What failed

Phase-space tensor convergence at the frozen production order `(8,6,4)`:

- `(8,6,4) -> (10,8,6)`: maximum relative Frobenius difference
  `3.6785399008114782e-06`, above the required `1e-6`.
- The four failures are the `mb_untruncated` anchors:
  Xe-131/268 keV `3.6785e-6`, Xe-129/268 keV `3.6100e-6`,
  Xe-131/201.5 keV `3.1646e-6`, Xe-129/201.5 keV `3.1281e-6`.
- The four `shm_cutoff` anchors pass at roughly `6.6e-7`.
- This is a real under-resolution of the `u` integration at high speed, not a
  NaN, non-finite tensor, Hermiticity, PSD, Gram, or decomposition failure.

At the worst anchor (Xe-131, 268 keV, `mb_untruncated`, speed 1126.029 km/s),
the production order is not close to the high-order reference:

| comparison | relative Frobenius |
|---|---:|
| `(8,6,4)` vs `(10,8,6)` | `3.6785399008114782e-06` |
| `(8,6,4)` vs `(12,10,8)` | `5.256875076892477e-06` |
| `(8,6,4)` vs `(14,12,10)` | `6.063310952475214e-06` |
| `(8,6,4)` vs `(16,14,12)` | `6.518610762050962e-06` |
| `(10,8,6)` vs `(16,14,12)` | `2.8591612828473853e-06` |
| `(12,10,8)` vs `(16,14,12)` | `1.2679670997438545e-06` |
| `(14,12,10)` vs `(16,14,12)` | `4.572537915221440e-07` |

The limiting quadrature is the `nw` (u) order; changing `nz` or `nphi` alone
does not remove the discrepancy.

## Superseding revision

The adopted revision raises the production order to `(14,6,4)` and compares
against `(16,14,12)`. At all eight anchors, the maximum `(14,6,4)` vs
`(16,14,12)` relative difference is `4.572538091907602e-07`, below `1e-6`.
The original `(12,6,4)` vs `(14,6,4)` diagnostic gave
`8.107133197075289e-07`, but the higher `(14,6,4)` production order was chosen
for margin. The original `(8,6,4)` run remains archived here as failed
provenance and is not used for C2.

## Boundary

`stage_c1_complete=false` in `run_manifest.json`. No Stage C2 contraction, no
nuclear-current contraction, no AFH factor-four or C0 application, and no
full-current or halo-fold claim are made here. Stage A and Stage B archives
were not regenerated.
