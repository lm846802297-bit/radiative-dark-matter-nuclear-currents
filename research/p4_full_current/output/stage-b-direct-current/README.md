# P4-FC Stage B: direct finite-transfer bound current

This main-branch archive imports the exact 19 Stage B paths from [PR #141](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/pull/141) at `280c46b443a0695e231d03d0d3f8b03e303ff962` plus its exact Borah/Sufian adapter and seven pinned AFH density CSVs. The original run is not rerun or relabelled. [Original manifest](stage_b_run_manifest.json), [native SHA-256 list](SHA256SUMS_stage_b.txt), [independent validation report](stage_b_validation_report.json), and [order-convergence report](convergence_report.json) all retain original bytes.

The archived AME2020 P4 calculation is seven isotopes × 81 points, with 891 rank-resolved rows. It explicitly retains K=0,L=1 for the five spin-zero isotopes. The external energy assignment is q_N^0=E_R; a nonzero divergence in this bound-current diagnostic is not treated as an on-shell Ward-identity failure. There is no DM tensor contraction, halo fold, detector prediction, or exchange-current calculation in this stage. The separate optional old `au` diagnostic still belongs to the original research stack and was not promoted.

## Main compatibility and historical byte replay

Main's shared `rate_kinematics.py` is retained without modification. Its SHA-256 is `ced8ad5261b6a3f4847f4d068a2dbab65ad8073623dc8b6d407fa21db38f3fe2` at main baseline `66c2ccca5baafe78e93b4cf99d41936caa6b8538`; the original run used `daf409e8246c04c290511dfb78f60992131057f75044cacaeebb159d3f9239a2`. A text comparison found only the module docstring and whitespace changed. The entire remaining Python AST, including functions and constants, is identical, and both variants passed all 7 focused tests. [The exact historical transfer-helper bytes](rate_kinematics_run_snapshot.py) are provided for the strict source-hash test, not for normal main imports.

Use a disposable worktree for strict historical replay, to avoid replacing the tracked shared helper or overwriting immutable reports:

```bash
audit_dir=$(mktemp -d)
git worktree add --detach "$audit_dir/replay" HEAD
(
  cd "$audit_dir/replay"
  p=research/three_body_full_current_nucleon_match/p4_full_current
  o=$p/output/stage-b-direct-current
  (cd "$o" && sha256sum -c SHA256SUMS_stage_b.txt)
  cp "$o/rate_kinematics_run_snapshot.py" research/three_body_calculation_and_validation/rate_kinematics.py
  python -B "$p/audit_stage_b_direct_current.py" --json "$o/direct_bound_current_o14.json" --out "$audit_dir/stage-b-replay.json"
  python -B "$p/audit_stage_b_convergence.py" --o10 "$o/convergence_o10.json" --o14 "$o/direct_bound_current_o14.json" --o18 "$o/convergence_o18.json" --o22 "$o/convergence_o22.json" --out "$audit_dir/stage-b-convergence-replay.json"
  python -B -m unittest discover -s "$p" -p 'test_stage_b_direct_current.py' -v
)
```

Integration verification on exact archived inputs: 13/13 output checksums, direct-current archive audit PASS (order-14 maximum fit residual `8.05e-15`, leading-NR-to-AFH maximum residual `1.39e-13`), convergence replay PASS (18-to-22 maximum absolute difference `6.04e-14`), and focused tests 7/7 with each transfer-helper version. No production-grid or new Wolfram run was required, no CI was run, and no manuscript or release flags were changed.
