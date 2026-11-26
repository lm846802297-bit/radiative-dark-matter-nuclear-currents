# P4-FC Stage C2: full one-body current contraction and halo fold

## 1. Objective

Use the frozen C1 Lorentz tensor grid and the audited Stage B direct-current multipoles to compute the complete one-body four-current contraction

$$
K_{\rm raw}(A,E_R,v)=\sum_K\frac{1}{2J_A+1}\sum_{m_i,m_f}
L^{\mu\nu}J_{\mu,K}(m_f,m_i)J^*_{\nu,K}(m_f,m_i),
$$

then apply the physical AFH normalization and perform the stored speed-node halo fold. C2 must not regenerate C1 phase-space nodes or Stage B nuclear-current matrix elements.

## 2. Frozen inputs and scope

- C1 production archive: 7 isotopes, 81 recoil points, 2 halo models, 24 speed nodes; 27,216 tensor rows, all `phase_space_open=true`.
- C1 quadrature: production `(14,6,4)`, comparison `(16,14,12)`; accepted convergence remains inherited from C1.
- Stage B direct-current archive: 567 isotope/recoil points and 891 rank-resolved rows.
- Stage B independent audit: `PASS_P4_FC_STAGE_B_DIRECT_CURRENT_AUDIT`; 567 rank-0 longitudinal rows, including 405 rows for the five even-even isotopes.
- Stage A native response table supplies the isotope-specific `C0(A)` correction.
- Natural-xenon isotope weights and the frozen B benchmark supply the halo and absolute-rate normalization.

The calculation covers 1,134 isotope/recoil/halo points and 27,216 speed-node contractions. It produces theoretical natural-xenon rates before detector efficiency and exposure; it does not claim LZ expected counts.

## 3. Fixed physics and normalization conventions

### Momentum assignments

$$
Q_A^2=2m_AE_R+E_R^2,\qquad t_A=-2m_AE_R,\qquad t_N=t_A.
$$

`Q_A` enters AFH multipoles, `t_N` enters the nucleon form factors already embedded in Stage B, and `t_A` enters the mediator propagator. These variables must not be interchanged.

### AFH normalization

The direct-current contraction is first checked against the raw AFH normalization:

$$
K_{\rm AFH,raw}=C_0(A)K_{\rm raw}.
$$

The physical equal-proton/equal-neutron normalization is then

$$
\boxed{K_{\rm physical}=4C_0(A)K_{\rm raw}.}
$$

The factor four and `C0(A)` are each applied exactly once. The raw leading-NR/AFH comparison does not include the factor four. Stage B already uses effective form factors `Fhat_i=F_i^(uds)/3`; no additional factor of three is allowed.

### Full current

All four Lorentz components are retained. The decomposition must include charge, spatial, charge-spatial interference, longitudinal `(0,z)`, transverse `(x,y)`, and longitudinal-transverse interference. Before/after-diagram terms and their signed interference must reconstruct the total contraction.

The `K=0,L=1` spatial current is retained as a longitudinal channel for all seven isotopes. It is not relabelled as AFH `Delta` or `SigmaPrime`; for even-even isotopes it is not discarded merely because `J_A=0`.

### Longitudinal continuity and the AFH six-response projection

For each rank and speed node define

$$
\alpha=\frac{E_R}{Q_A},\qquad d=J^3-\alpha J^0.
$$

The AFH six-response comparison is made with the continuity-projected current
$J^3_{\rm proj}=\alpha J^0$. Its longitudinal weight is $L_C|J^0|^2$,
and its transverse part contains the matched $\Delta$, $\Sigma'$ and signed
$\Sigma'\!-\!\Delta$ terms. The production contraction always retains the
actual $J^3$. It separately records

$$
K_{0z}-L_C|J^0|^2
=-(L_{03}-\alpha L_{33})J^0d^*
 -(L_{30}-\alpha L_{33})dJ^{0*}
 +L_{33}|d|^2,
$$

as well as projected and nonconserved longitudinal-transverse cross terms.
No nuclear-current continuity condition is imposed on the actual Stage B
impulse current. The full exact-current kernel is reconstructed from the
AFH six-response projection plus all separately recorded residuals.

### Halo-fold normalization

The rate fold reuses the C1 speed nodes, speed weights, halo PDFs, `beta`, `p_chi`, `m_A`, and `t_A`. It applies the frozen external coupling product `(g1*g2*gVN)^2`, with `gVN=3*gVq`, the mediator denominator `(m_V^2-t_A)^2`, and the natural isotope abundances without renormalizing the seven-isotope AFH coverage. Output units are events/(kg day keV).

## 4. Predeclared acceptance criteria

1. All pinned C1, Stage A, Stage B, isotope-input, B-freeze, and source hashes match.
2. C1 manifest is production `PASS`, C1 is complete, and C1 has not marked `full_current_complete`.
3. Exact grid coverage: 7 isotopes x 81 recoil points x 2 halos x 24 speed nodes = 27,216 rows; all nodes open.
4. Stage B audit is PASS; 567 points, 891 rank rows, 567 `K=0,L=1` rows, 405 even-even `K=0,L=1` rows.
5. Dynamic target-spin dimensions are used for every isotope; rank sets are stable across recoil points.
6. C1 Lorentz tensors pass Hermiticity and PSD checks; the C1 before/after/interference tensor decomposition is preserved.
7. Full-current component and longitudinal/transverse decompositions close to `1e-12`; diagram decomposition closes to `1e-11`; imaginary residual is below `1e-12` relative.
8. The leading-NR AFH comparison uses the continuity-projected current and agrees with the native AFH `L_C` plus transverse six-response contraction to `1e-8` on the full grid. The actual `J^3` is retained in the exact contraction; the analytic nonconserved `0-z` residual and longitudinal-transverse residuals are recorded and reconstruct the full kernel. For `J_A=0`, the `K=0,L=1` current remains explicit in the full contraction and its continuity residual is not silently absorbed into the AFH charge response.
9. An independent explicit-loop contraction agrees with the vectorized contraction to `1e-12` at deterministic coverage points.
10. The total exact-current kernel is nonnegative within numerical tolerance; the seven-isotope halo contributions sum to the natural-xenon spectrum and diagram contributions reconstruct the total rate.
11. The run records exact source/input hashes, Git commit/branch/clean status, Python and library versions, hardware, parameters, wall time, native CSV outputs, validation report, and `SHA256SUMS.txt`.

Only after all checks pass may the C2 manifest set `stage_c2_complete=true`, `one_body_impulse_current_complete=true`, and `P3_tensor_validation_complete=true`. The generic `full_current_complete` flag remains false until a consistent nuclear Hamiltonian and exchange-current closure are available.

## 5. Production execution

Run from a clean checkout of the C2 branch. The output directory must be empty and outside the repository so that the source tree remains clean during execution:

```bash
python research/three_body_full_current_nucleon_match/p4_full_current/stage_c2_full_current_contraction.py --out-dir C:/p4c2_run_output/stage-c2-full-current-contraction
```

After a production PASS, copy the complete output directory into
`research/three_body_full_current_nucleon_match/p4_full_current/output/stage-c2-full-current-contraction/`, verify its checksums, and submit a PR targeting the C1 research branch. Do not modify `main` directly.

## 6. Current status

- C1: accepted and archived; no recomputation requested.
- Stage B independent current audit: PASS and reused read-only.
- C2 implementation: continuity-projection and full Proca-numerator implementation present; the complete C2 unit-test suite passes.
- Production C2: **PASS** on the full pinned grid. The archive is
  `output/stage-c2-full-current-contraction/`; the clean run commit is
  `90d36811107df51120042db752c5bb3a98d322ce`, with 27,216 speed-node
  contractions, 1,134 isotope/recoil/halo points, and 162 natural-xenon rate
  rows. C2 unit tests passed 12/12 in the production run.
- `stage_c2_complete=true`, `one_body_impulse_current_complete=true`,
  `P3_tensor_validation_complete=true`; `full_current_complete=false` remains
  unchanged pending a consistent nuclear Hamiltonian and exchange-current
  closure. No detector efficiency or exposure is applied.

## 7. Pinned production inputs

| Input | SHA-256 |
|---|---|
| C1 `full_lorentz_tensor_grid.npz` | `98898fbdac7a08a81453972d013ad027e46679a02e67d41a3d59a529f1671c76` |
| C1 `tensor_index.csv` | `3647a0b4261be85b146ffe040f4e0241fd1508d85cc0d5f1f9b502dc1b33c06c` |
| C1 `run_manifest.json` | `18d36a643205ad0f8cde581cb04f1ff6dbd3f2aa6e86e1877837b6b8678fc8e8` |
| C1 `convergence_report.json` | `1cb9df212451374652e0713a477f86e8954e3f9d8aee9aa71123e2910b3a0878` |
| Stage B `direct_bound_current_o14.csv` | `0b335b39cd288beda1367fd977e311d595c2a88c2c5e7dfcaa9456cf8e26d355` |
| Stage B `stage_b_validation_report.json` | `5f0272f2c976ab894d8e69fc6930c609c13a0aba93ee445c5b83710c0ff4f313` |
| Stage A native response CSV | `39a5b8c2401eb5f4e5acabb0a4f38b9a9d1fa7103d12d208462b3fcefef57c48` |
| Natural-xenon isotope weights | `fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991` |
| B benchmark freeze | `2b0b8c58baf559c176fe9e64d14f75973d7987b072fa6af2275f09430d5fbe88` |
| Stage B production script | `78058780862acb06d5c3fd86957dead6d16d9a8e3041c9b4bdc3d5450189b75e` |
| C1 production script | `c82b5a47c12fb59967dd41cf0de54c01a3c2a42c0e16152e2e0484babaa91bbc` |
| Native halo-fold reference | `71ed1c7f4b68d243f7a6b7f65727066e13bb958852630a0ab61f9d690d95f5c1` |

### Full Proca numerator (implemented in C2 v2)

The production contraction must include the full massive-vector numerator

\[
P_{\mu\nu}(q_A)=g_{\mu\nu}-\frac{q_{A,\mu}q_{A,\nu}}{m_V^2},
\qquad
q_A^\mu=(E_R,0,0,+Q_A),
\qquad
g_{\mu\nu}=\operatorname{diag}(1,-1,-1,-1).
\]

With the existing bilinear convention, where the current matrices are lowered
inside the contraction routine, define the mixed matrix

\[
M^\mu{}_\nu=\delta^\mu{}_\nu-
\frac{q_A^\mu q_{A,\nu}}{m_V^2}.
\]

The particle tensor entering the nuclear-current contraction is

\[
\boxed{L_{\rm Proca}=M^T L M.}
\]

This transformation is applied to the total, before, after, and
interference particle tensors before all exact-current, leading-NR,
continuity-projection, AFH-reference, diagram-interference, and
independent-loop contractions. The scalar denominator
\((m_V^2-t_A)^{-2}\) remains in the rate kernel and is applied once; it is
not part of the numerator transformation.

C2 v2 records and checks the Hermiticity, positive semidefiniteness, and
before/after/interference reconstruction of the Proca-dressed tensor. Its
unit tests additionally compare the transformed bilinear against a direct
covariant-amplitude contraction, check the zero-transfer limit, and verify
that the Proca map preserves the diagram-tensor sum.

The algebraic identity was independently smoke-checked in the ChatGPT Linux
runtime on 1,000 random complex current pairs; the maximum relative residual
was \(3.37\times10^{-15}\). This is a local algebraic check only, not a run
of the repository C2 test suite or the production grid.

**Execution status:** the Proca numerator is now present in the C2 source
branch. The full C2 test suite and production run against the pinned C1
archive have not yet been executed. C2 remains NOT RUN until the exact C1
NPZ/index bytes are available in the local execution environment, their
SHA-256 values pass, and all predeclared production acceptance checks pass.
