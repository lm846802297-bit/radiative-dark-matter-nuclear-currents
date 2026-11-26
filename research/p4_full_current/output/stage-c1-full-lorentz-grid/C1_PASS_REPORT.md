# Stage C1 production PASS archive

Status: **PASS**, `stage_c1_complete=true`.

This is the canonical C1 full Lorentz-tensor archive. It was produced on the
Windows local research environment from the clean source tree at commit
`6ff1f357539adde80392f6aa587d695949050edb` (branch
`research/p4-fc-stage-c1-full-lorentz-grid-20261124`). The run used 8 workers,
Python 3.13.18, NumPy 2.3.5, SciPy 1.17.0; no GitHub Actions runner was used.

## Production grid

- 7 AFH isotopes, 81 recoil points, 2 halo domains, 24 speed nodes;
- 1,134 isotope/recoil/halo points and 27,216 tensor rows, all open;
- phase-space production order `(14,6,4)`, comparison order `(16,14,12)`,
  mapping `threshold_u2`;
- B benchmark unchanged: `m_chi1=1000 GeV`, `m_s=200 keV`,
  `delta_m=1000 keV`, `g2=1`, asymmetric Pole width.

## Numerical acceptance

- Hermiticity `0.0`; negative PSD eigenvalue relative measure `0.0`;
  Gram partial-trace `0.0`.
- Before + after + interference decomposition residual
  `3.682924815528239e-16`.
- q-frame rotation residual `1.0182568436090762e-15`.
- Speed convergence `24 -> 48`: `9.197564387522144e-06` (limit `1e-5`).
- Phase-space convergence `(14,6,4) -> (16,14,12)`:
  `4.572538091907602e-07` (limit `1e-6`).
- Independent Dirac-trace residual: `6.63189503313788e-11` (limit `1e-8`).
- Wall time `825.30 s`.

The superseded `(8,6,4)` production trial is archived separately in
`../stage-c1-full-lorentz-grid-nw8-fail/`. It remains `FAIL_CONVERGENCE` and
is not used for C2.

## Boundary

The archive contains uncontracted q-aligned `L^{mu nu}(v,E_R)` tensors only.
It does not contain the Stage B nuclear-current contraction, `C_0(A)`, the AFH
factor four, mediator/coupling prefactors, or an event rate. Those belong to
C2. `full_current_complete=false`; Stage A and Stage B archives were reused
read-only and not regenerated.
