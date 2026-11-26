# Independent local replay: C2 continuity-projection stress test

## Result

Status: PASS_C2_CONTINUITY_PROJECTION_STRESS_REPLAY.

The exact standalone test source from the merged C2 continuity-projection change was executed locally on 2026-11-25. All four tests passed, including 20,000 deterministic random Hermitian-positive-semidefinite 4x4 tensor/current samples.

## Numerical results

- Maximum analytic longitudinal-residual relative difference: 5.87315574842313e-14.
- Maximum full-kernel decomposition relative residual: 1.9559898855147555e-15.
- Maximum conserved-current residual: exactly zero.

## Provenance

- Test source SHA-256: 17a793d74df9562c62d6da5ffcd95fbc95bfa846da177015078d8937b9bbf67f.
- Python 3.13.5, NumPy 2.3.5.
- Linux 6.18.57 x86_64; Intel Xeon Platinum 8573C; 5 logical CPUs; 6,236,844,032 bytes RAM; no GPU.
- Wall time: 4.55 s.

## Scope boundary

This is a synthetic algebraic/software replay. It does not consume the C1 production tensor grid or Stage B production currents and is not a production C2 result or detector-rate calculation.
