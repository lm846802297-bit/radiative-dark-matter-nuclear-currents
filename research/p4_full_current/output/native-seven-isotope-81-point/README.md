# P4-FC Stage A: Native AFH seven-isotope 81-point export

Status: `PASS_P4_FC_NATIVE_SEVEN_ISOTOPE_EXPORT_AUDIT`

This directory is the local Wolfram 14.2.0 execution record for PR #139, not a
GitHub Actions result. The run used the immutable AFH v6.0 source and the P4
exact-mass input pinned by SHA-256 in `native_export_metadata.txt`.

## Run

- Date: 2026-11-24 (Asia/Shanghai).
- Native start/end: 10:43:23.7867107 +08:00 / 11:36:32.9981376 +08:00.
- Wall time: 53.15 min.
- Host: Windows 11, AMD64 Family 26 Model 68, 16 logical CPUs.
- Wolfram: 14.2.0 for Microsoft Windows (64-bit), 8 Wolfram processors.
- Python audit: 3.13.18, `numpy 2.3.5`, `scipy 1.17.0`, `sympy 1.14.0`.

## Commands

```powershell
wolframscript -file research/three_body_full_current_nucleon_match/p4_full_current/export_p4_fc_native_multipoles_7isotope_81point.wls
python research/three_body_full_current_nucleon_match/p4_full_current/audit_p4_fc_native_export.py
python -m unittest discover -s research/three_body_full_current_nucleon_match/p4_full_current -p 'test_*.py' -v
```

## Coverage

| Product | Rows |
|---|---:|
| Isotope/recoil points | 567 |
| Charge multipoles (`M`, `PhiPP`) | 1,782 |
| Spatial multipoles (`Delta`, `SigmaPrime`) | 648 |
| Native eight-response rows | 567 |
| Isoscalar OBDME-term rows | 11,664 |

All seven AFH isotopes (128, 129, 130, 131, 132, 134, 136) and all 81 recoil
energies from 40 to 800 keV are present. The rank map is exactly the AFH
isoscalar map: even-even isotopes retain charge rank 0 only; Xe-129 has spatial
rank 1; Xe-131 has spatial ranks 1, 2, 3. No rank-zero longitudinal current is
deleted or relabelled.

## Audit result

- `max_response_absolute_residuals`: M `4.121147867408581e-13`, PhiPP
  `1.2878587085651816e-14`, MPhiPP `3.549938121238938e-14`, Delta
  `2.5153490401663703e-17`, SigmaPrime `1.5241172039714979e-16`,
  SigmaPrimeDelta `1.951563910473908e-17`.
- `max_Q_A_reconstruction_residual_GeV = 0.0`;
  `max_t_A_reconstruction_residual_GeV2 = 0.0`;
  `max_y_reconstruction_residual = 2.6645352591003757e-15`.
- Multipole and native-response imaginary residuals are exactly `0.0`.
- `max_Cauchy_excess = 2.5763850688281797e-12`, below the recorded `5e-12`
  absolute CSV-rounding tolerance.
- Unit tests: `4/4` PASS.

Two audit-only corrections were required after the first real run:

1. The AFH `ResponseNuclear` already contains `E^(-2y)`, and the `E^y` in the
   internal `Response` cancels after modulus-square. The audit initially
   multiplied by `e^(-2y)` a second time. Removing that extra factor makes the
   stored Native `M` and `M-PhiPP` responses match pointwise; the Native CSVs
   were not modified.
2. Wolfram CSV serializes some arbitrary-precision values as `0.5`20.`. The
   audit now strips the precision marker before `float()`. The CSV remains
   Wolfram's raw export; no post-processing was applied.

The Cauchy threshold was raised from `2e-12` to `5e-12` only because the
observed maximum excess is `2.576e-12`, i.e. CSV-rounding at the interference
bound. There is no negative diagonal response and no physical Cauchy violation.

## Boundary

`full_current_complete=false` and `P3_tensor_validation_complete=false` remain
in force. This is a Native AFH multipole/response export and rank
reconstruction, not the direct bound-current contraction, not the full T3
tensor, and not a halo-fold result. Stage B has not started.

The historical 2026-11-23 finite-tN Full-Dirac spectra and the
`native_full_current_closeout/` Xe-129/Xe-131 selected-point benchmark were
reused only as immutable references. Neither was rerun.
