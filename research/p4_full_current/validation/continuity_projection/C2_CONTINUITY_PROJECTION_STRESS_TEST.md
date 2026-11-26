# C2 longitudinal-continuity projection algebra stress test

Date: 2026-11-25  
Scope: algebraic/software validation of the Stage C2 continuity-projection helper only. This is not a production C2 run and contains no physical recoil-rate result.

## Reproduction

From the repository root, run:

~~~bash
python research/three_body_full_current_nucleon_match/p4_full_current/validation/continuity_projection/test_c2_continuity_projection_stress.py
~~~

The test uses deterministic random seeds. It checks:
1. the conserved-current limit reduces to the LC contraction;
2. the nonconserved 0-z residual agrees with an independent direct 0,z contraction and the analytic formula;
3. the AFH six-response projection plus all retained cross/residual terms reconstructs the full kernel;
4. a 20,000-sample stress test over random Hermitian positive-semidefinite 4x4 tensors and complex currents.

## Results

- Tests: 4/4 PASS.
- Random stress samples: 20,000.
- Maximum analytic longitudinal-residual relative difference: 5.87315574842313e-14.
- Maximum full-kernel decomposition relative residual: 1.9559898855147555e-15.
- Maximum conserved-current residual: exactly zero.

## Provenance

- Test source SHA-256: 17a793d74df9562c62d6da5ffcd95fbc95bfa846da177015078d8937b9bbf67f.
- Python: 3.13.5.
- NumPy: 2.3.5.
- OS/kernel: Linux x86_64, kernel 6.18.57.
- CPU: Intel Xeon Platinum 8370C @ 2.80 GHz, 5 logical CPUs.
- Memory visible to the runtime: 5.8 GiB.
- GPU: none available (nvidia-smi unavailable).

The production C1 tensor archive and Stage B production inputs are not consumed by this test. Stage A/Stage B production artifacts and P0-B1 phase-space tensors are not recomputed.
