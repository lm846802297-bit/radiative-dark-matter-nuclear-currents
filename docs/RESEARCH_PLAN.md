# Research plan and promotion gates

## Stage 0 — repository and theory baseline

- [x] Separate public repository created.
- [x] General-current framework note migrated with source provenance.
- [x] Kinematic, Lorentz, AFH, and channel conventions established.
- [x] Prior-work boundary and reusable references recorded.
- [x] Universal \(u+d+s\) VV source, one-nucleon matching, and channel-specific two-body-current literature audited.
- [ ] Import source-pinned VV baseline data and an independent reproduction report.
- [ ] Freeze machine-readable benchmark configuration and environment lock.

## Stage 1 — dark-matter-side three-body kernels

1. VV recovery: reproduce the inherited two-diagram vector kernel from an exact source commit; compare direct covariant amplitude and integrated spin-resolved tensor at fixed \(E_R,v\).
2. VA: keep dark transition vector and replace only the nucleon current by axial; implement \(G_A,G_P\), PCAC, axial charge, and spatial multipoles.
3. AV and AA: derive a new axial dark transition kernel; do not reuse the VV tensor by changing a response label.
4. SS: define the scalar mediator/source and Wilson coefficients at quark/gluon level; derive the scalar dark-line kernel and match scalar one-/two-body nuclear source.

## Stage 2 — nuclear-current matching

- Vector: reproduce the inherited universal \(u+d+s\) VV/AFH benchmark with explicit \(g_{VN}=3g_{Vq}\) normalization and finite-\(t_N\) form factors; compare the historical, 2020, and 2026 strange-form-factor input sets. Add a two-body vector operator only after source-specific matching shows it contributes to VV at the selected order. The published AV two-body vector current is not a VV correction.
- Axial: match \(G_A,G_P\), axial charge and spatial multipoles; compare xenon axial structure factors and uncertainty bands.
- Scalar: match quark/gluon Wilson coefficients, nucleon scalar form factors, pion-exchange two-body current, and trace-anomaly terms.
- Preserve amplitudes and signed interference until after nuclear contraction and squaring.

## Stage 3 — xenon nuclear responses

- Vector coherent channel: AFH seven-isotope set.
- Axial elastic spin channel: \(^{129}\mathrm{Xe}\) and \(^{131}\mathrm{Xe}\) first.
- Scalar coherent channel: published xenon scalar structure factors as cross-checks; map normalization before use.
- Keep AFH, shell-model, and chiral-EFT uncertainties separate.

## Stage 4 — physical recoil predictions

1. Fixed-speed \(d\sigma/dE_R\).
2. Full three-body phase-space integration.
3. Isotope-by-isotope recoil spectra.
4. Natural-xenon abundance weighting.
5. Halo folding under declared halo models.
6. Detector efficiency, resolution, exposure, and event likelihood.

## Promotion gates

| Gate | Required evidence |
|---|---|
| Theory derivation | Explicit source/Lagrangian, amplitude, index conventions, derivation |
| Matching | Relativistic-to-NR coefficients and independent algebraic checks |
| Nuclear response | Source-pinned OBDME/structure factors, normalization and interference audit |
| Numerical production | Exact source/input hashes, environment/hardware, raw outputs, convergence |
| Physical validation | Independent representation, Ward/PCAC or symmetry checks, analytic limits |
| Published result | PR-reviewed code/data and archived manuscript figures/tables |

A completed VV recovery test does not imply VA, AA, or SS is complete.
