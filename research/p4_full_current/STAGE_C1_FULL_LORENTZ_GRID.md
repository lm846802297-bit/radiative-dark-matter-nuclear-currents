# P4-FC Stage C1 — full Lorentz tensor grid

## Purpose and boundary

C1 generates the missing uncontracted, q-aligned three-body Lorentz tensors
\[
L_A^{\mu\nu}(v,E_R),\qquad \mu,\nu=0,1,2,3,
\]
for the seven Native-AFH xenon isotopes, the 81-point recoil grid, and both
halo speed domains. It is a tensor-data stage only. It does not contract the
Stage B nuclear current, apply \(C_0(A)\), apply the AFH factor four, multiply
the mediator/coupling prefactors, or produce an event rate. Those operations
belong to C2.

## Immutable inputs reused

- Formal B freeze:
  \`research/nucleon_form_factor_B_baseline/B_DATA_VERSION_FREEZE_2026-11-24.json\`.
- Stage A Native AFH multipoles and responses:
  \`p4_full_current/output/native-seven-isotope-81-point/\`.
- Stage B finite-\(t_N\) direct-current production:
  \`p4_full_current/output/stage-b-direct-current/direct_bound_current_o14.csv\`.
- AME2020 isotope masses:
  \`research/afh_default_density_natural_xenon_halo/input/natural_xenon_isotope_weights.csv\`.
- Three-body integrator:
  \`three_body_phase_space_tensor.integrate_spin_resolved_tensor\`.

The runner verifies the frozen input SHA-256 values and the B-freeze version
before starting. Stage A and Stage B are read-only references; neither is
recomputed.

## Frozen kinematics and quadrature

The B benchmark is read from the frozen manifest and checked against:

- \(m_{\chi_1}=1000\) GeV;
- \(m_s=200\) keV;
- \(\delta_m=1000\) keV;
- \(g_2=1\);
- asymmetric Pole width (the tensor integrator's \`pole\` option).

Production grid:

- isotopes: Xe-128, Xe-129, Xe-130, Xe-131, Xe-132, Xe-134, Xe-136;
- recoil: 40–800 keV, 81 points, 9.5-keV spacing;
- halo domains: SHM cutoff \(v_{\max}=776\) km/s and untruncated MB \(v_{\max}=1600\) km/s;
- speed Gauss-Legendre order: \(n_v=24\);
- phase-space orders: \((n_w,n_z,n_\phi)=(14,6,4)\);
- phase-space convergence comparison order: \((16,14,12)\);
- phase-space mapping: \`threshold_u2\`.

The initial C1 trial used \((8,6,4)\) with a \((10,8,6)\) comparison. That
trial was executed on the complete grid and is archived separately as
\`output/stage-c1-full-lorentz-grid-nw8-fail/\`; its manifest is
\`FAIL_CONVERGENCE\` because the untruncated-MB anchors reached
\(3.6785\times10^{-6}\) in relative Frobenius norm, above \(10^{-6}\). The
failure is an under-resolution of the \(u\) direction, not a matrix or
non-finite-tensor failure. The revised production order above changes only
the quadrature; the B benchmark, AFH inputs, isotope grid, and physics
boundary are unchanged.

At each isotope/recoil/halo point, the speed interval is
\([v_{\min}^{3\to2},v_{\max}]\). The speed nodes and weights are retained in
the index so C2 can apply the halo PDF, \(\beta\), and \(p_\chi^{-2}\) without
reconstructing the quadrature.

## Stored products

The runner writes:

- \`full_lorentz_tensor_grid.npz\`: complex \(4\times4\) \(L\), diagram-before,
  diagram-after, and signed diagram-interference tensors at every speed node;
- \`tensor_index.csv\`: row-to-array index, isotope, recoil, halo, speed,
  quadrature weight, PDF, \(\beta\), \(p_\chi\), \(Q_A\), \(t_A=t_N\), and
  per-node diagnostics;
- \`point_coverage.csv\`: all isotope/recoil/halo points, including closed
  speed intervals;
- \`convergence_report.json\`;
- \`run_manifest.json\` and \`SHA256SUMS.txt\`.

The Lorentz tensors are stored in the right-handed q-aligned frame. The DM
spin trace is already taken in \(L^{\mu\nu}\); the four-component Lorentz
indices and all complex off-diagonal entries are retained.

## Acceptance criteria

All conditions are required for a production C1 PASS:

1. Complete seven-isotope × 81-recoil × two-halo coverage; every open speed
   interval has 24 indexed speed nodes.
2. Every stored \(L\) and diagram component is finite and has shape \(4\times4\).
3. Per-node Hermiticity residual \(\le 10^{-10}\).
4. Per-node negative PSD eigenvalue, normalized by \(\|L\|\), is no larger
   than \(10^{-10}\).
5. Independent partial-trace reconstruction from the returned spin Gram
   tensor agrees with stored \(L\) to \(10^{-12}\).
6. Before + after + signed interference reconstructs total \(L\) to \(10^{-12}\).
7. q-frame rotation orthogonality/determinant residual \(\le 10^{-12}\).
8. On Xe-129/Xe-131 at 201.5/268.0 keV and both halo domains, the
   speed-integrated full tensor agrees between \(n_v=24\) and \(48\) to
   \(10^{-5}\) in relative Frobenius norm.
9. At the median speed of those anchors, the full tensor agrees between
   \((14,6,4)\) and \((16,14,12)\) phase-space orders to \(10^{-6}\) in relative
   Frobenius norm.
10. The selected anchor tensors independently reproduce the direct
    Dirac-trace contractions for the charge, spatial, interference, and full
    tensor sectors to (10^{-8}).
11. Source/input hashes, actual local runtime/hardware, raw NPZ arrays,
    point index, coverage, convergence report, and SHA-256 inventory are
    archived together.

The \(10^{-5}\) speed criterion is consistent with the previously accepted
rate-level \(24\to48\) check; the \(10^{-6}\) phase-space criterion is a
separate full-tensor check and is not inferred from the old scalar-rate
convergence result.

## Execution

From the repository root, on the local research compute environment:

\`\`\`powershell
python research/three_body_full_current_nucleon_match/p4_full_current/stage_c1_full_lorentz_tensor_grid.py --workers 8
python -m unittest discover -s research/three_body_full_current_nucleon_match/p4_full_current -p 'test_stage_c1_*.py' -v
\`\`\`

A reduced smoke run must use \`--pilot\`; it is explicitly labelled
\`PILOT_NOT_PRODUCTION\` and cannot be used to close C1. The full production runner does not dispatch GitHub Actions. Run the unit tests
locally before executing the production grid.

## Stop boundary

C1 PASS means the full per-speed Lorentz tensor archive exists and passes its
own numerical checks. It does not mean the Stage B current has been contracted
with the tensors. Keep \`P3_tensor_validation_complete=false\` and
\`full_current_complete=false\` until C2 and the independent nuclear-amplitude
validation are complete.
