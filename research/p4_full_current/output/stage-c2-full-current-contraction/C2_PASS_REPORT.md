# Stage C2 production PASS archive

Status: **PASS**.

This archive is the first complete one-body impulse-current C2 production
contraction. It reuses the frozen C1 tensor archive and the audited Stage A/B
inputs read-only; no C1 phase-space tensor, Stage A multipole, or Stage B
nuclear-current calculation was regenerated.

## Provenance

- C2 branch: `research/p4-fc-stage-c2-production-20261126`
- Clean C2 source commit during the run:
  `90d36811107df51120042db752c5bb3a98d322ce`
- Parent C1 archive: commit `6ff1f357539adde80392f6aa587d695949050edb`
- C2 branch is a merge of the C1 archive branch and the current `main` C2
  implementation; the C1 and C2 inputs coexist in one clean checkout.
- Run environment: Windows local research host, Python 3.13.18, NumPy 2.3.5,
  SciPy 1.17.0, SymPy 1.14.0, no GitHub Actions runner.
- Wall time: `21.04 s`.

## Coverage and physics

- 7 AFH xenon isotopes x 81 recoil points x 2 halo models x 24 speed nodes.
- 1,134 isotope/recoil/halo points and 27,216 speed-node contractions.
- 162 natural-xenon rate rows (81 recoil points x 2 halos).
- C1 production quadrature `(14,6,4)`; convergence `(16,14,12)` inherited.
- Full Proca numerator applied once:
  `L_Proca = M^T L M`,
  `M^mu{}_nu = delta^mu_nu - q_A^mu q_A,nu / m_V^2`,
  `q_A^mu=(E_R,0,0,+Q_A)`, metric `(+---)`.
- The scalar denominator `(m_V^2 - t_A)^2` remains a separate rate factor.
- AFH physical normalization `4*C0(A)*K_raw` applied once; the raw
  leading-NR/AFH comparison uses `C0(A)*K_AFH6` without the factor four.
- The actual Stage B `J^3` is retained; the AFH six-response comparison uses
  the continuity-projected `J^3_proj=(E_R/Q_A)J^0`, with the nonconserved
  longitudinal residual recorded separately.

## Numerical acceptance

- Proca-dressed tensor Hermiticity: `2.78477603892925e-18`.
- Proca-dressed negative PSD eigenvalue relative measure: `0.0`.
- Proca-dressed diagram decomposition: `7.2378758358977e-16`.
- Contracted component decomposition: `1.70162980834016e-15`.
- Contracted diagram decomposition: `2.05546450884744e-15`.
- Longitudinal/transverse decomposition: `5.7667218704072e-16`.
- AFH projection decomposition: `5.7667218704072e-16`.
- Nonconserved longitudinal formula, conditioned against the large `K_0z`
  scale: `3.11253155197955e-15`; the unconditioned relative-to-residual value
  is `7.20255224620872e-10` and is dominated by cancellation in the direct
  difference of two `~10^7` kernels.
- Spinful leading-NR/AFH six-response difference: `1.55702919021284e-12`.
- Spin-zero leading-NR charge/AFH LC difference: `8.62584020097197e-12`.
- Explicit-loop vs vectorized contraction: `1.85804981412264e-16`.
- Minimum exact raw kernel: `0.00317871639930551` (positive).
- C2 unit tests: `12/12` PASS in the production run.

## Boundary

`stage_c2_complete=true`, `one_body_impulse_current_complete=true`,
`P3_tensor_validation_complete=true`.
`independent_nuclear_amplitude_validation_complete=false` and the generic
`full_current_complete=false` remain unchanged: the one-body impulse current
is completely contracted, but a consistent nuclear Hamiltonian and
exchange-current closure are not claimed. No detector efficiency or exposure
is applied; no LZ expected-count claim is made.
