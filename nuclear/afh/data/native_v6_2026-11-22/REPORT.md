# Native AFH Mathematica calculation — 2026-11-22

## Outcome

The native Wolfram calculation completed for seven xenon isotopes and all eight AFH v6.0 nuclear responses. The pointwise comparison passed for all 16,807 isotope/q records. **This closes the native AFH eight-response response-database gate for the pinned source, isotope set, q grid, and tolerances.**

It does not close the separate direct covariant bound-nucleus P3-C amplitude gate. The full-current status remains `full_current_complete=false`, and the P4 full-current halo result remains provisional.

## Kinematic and normalization conventions

- Invariant transfer: \(t_A=-2m_AE_R\).
- Spatial recoil momentum: \(Q_A^2=2m_AE_R+E_R^2\), with \(q_A=\sqrt{Q_A^2}\).
- AFH response argument: \(y_A=(b_Aq_A/2)^2\); the response table is evaluated as `ResponseNuclear[y,index,0,0]`.
- The propagator factor \((m_V^2-t_A)^2\) belongs to the particle/scattering assembly and is not part of this AFH nuclear-response run.
- Compared values are raw isoscalar AFH responses. The equal-proton/equal-neutron factor-four conversion is not applied to these raw columns. Signed \(M\Phi''\) and \(\Sigma'\Delta\) interference values are retained.

## Direct comparison

| Response | Maximum absolute difference | Maximum relative difference |
|---|---:|---:|
| M | 2.04636e-12 | 6.27818e-06 |
| SigmaDoublePrime | 1.40513e-16 | 1.5629e-13 |
| SigmaPrime | 2.94903e-16 | 2.19798e-07 |
| PhiPP | 8.52651e-14 | 6.00127e-07 |
| PhiTildePrime | 1.59595e-16 | 9.14563e-10 |
| Delta | 2.08167e-16 | 4.20566e-09 |
| MPhiPP | 3.12639e-13 | 1.29148e-09 |
| SigmaPrimeDelta | 1.63064e-16 | 1.88568e-10 |


Acceptance: `rtol=5e-07`, `atol=2e-12`, and native imaginary residual `<= 2e-12`. All 16,807 rows passed and the failure list is empty. Across all eight responses the largest use of the combined pointwise tolerance \( |\Delta W|/(2\times10^{-12}+5\times10^{-7}\max|W|) \) is 0.05879 (about 5.88% of the limit). A maximum relative error above `rtol` can occur at near-zero values; all records still pass the combined absolute-plus-relative criterion.

## Run environment and source

- Local Windows 10 Pro (build 26300), AMD Ryzen 7 9800X3D, 31.1 GiB RAM.
- Wolfram 14.2.0 for Microsoft Windows (64-bit) (December 26, 2024); Python 3.14.7, NumPy 2.5.4, SciPy 1.18.1, SymPy 1.14.0.
- Workspace branch/commit: `fix/pr67-halo-audit-finite` / `aac9b33d66f6a7fe1ca5d5d8cc648f5afe52a81b`.
- AFH package SHA-256: `dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372`.
- The native output was cleanly rerun; all seven per-isotope native CSV hashes match the first native export byte-for-byte.
- The restored historical runner used the AFH package's native built-in default density matrices. The source-translated Python reference used the seven archived density-matrix inputs included here.

## P3-B spatial multipole amplitude comparison

| Amplitude channel | Rows | Maximum absolute difference | Failed points |
|---|---:|---:|---:|
| `DeltaJ` | 9,604 combined rows | 8.88178e-16 | 0 |
| `SigmaPJ` / \(\Sigma'_J\) | 9,604 combined rows | 2.05391e-15 | 0 |

The 9,604 rows cover Xe-129 \(J=1\) and Xe-131 \(J=1,2,3\), each on the same 2,401-point grid. Both channels passed `rtol=5e-07`, `atol=2e-12` with zero imaginary residual. The largest relative error occurs for a channel that is zero in the Python result; its absolute difference is below \(10^{-15}\) and passes the combined tolerance.

At strictly `q=0`, the direct AFH multipole formula evaluates to `Indeterminate` because of a removable zero-over-zero. The final native run evaluates that one row at `q_eval=1e-9 GeV` but keeps the archived coordinates `q=0, y=0`. Convergence against the P3-B Python zero-momentum amplitudes was checked at both `1e-8` and `1e-9 GeV`: maximum absolute errors were 4.55e-15 and 3.89e-16, respectively; the fixed acceptance absolute tolerance is `2e-12`. The final corrected Wolfram run exited 0 with a finite CSV and zero imaginary residual.

## Comparison with the P4 workflow status

| Layer | P4 Actions evidence | Local native computation | Status after this run |
|---|---|---|---|
| Halo-fold numerical workflow | The two linked PR #67 runs passed the seven-isotope, velocity, three-body, provenance, and upload jobs | Not rerun here | P4 numerical workflow remains PASS |
| AFH nuclear responses | Previously provisional Python reference; native gate was open | Mathematica `ResponseNuclear` 1–8 compared against 16,807 Python records | **Native AFH gate PASS** |
| P3-B spatial amplitudes | Source-translated Python `DeltaJ` / `SigmaPJ` amplitudes | Native Wolfram comparison over 9,604 rank/q rows | **PASS for these multipole amplitudes** |
| Full particle–nuclear current match | Separate from the P4 integration jobs | Not computed in this task | **Still open**; `full_current_complete=false` |

The P4 job links supplied with the task: [seven-isotope Halo fold](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/actions/runs/44962224752) and [continuous natural-xenon spectrum](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/actions/runs/44962224780).

## Files and reproduction

- `native/`: seven native Wolfram CSV exports, native metadata, logs, and pointwise comparison JSON.
- `python_reference/`: local source-translated Python AFH-8 data, metadata, and internal validation record.
- `native/spatial_multipoles/`: native Xe-129/131 rank-resolved Delta and Sigma-prime amplitudes plus their P3-B comparison.
- `p3_python_reference/`: source-translated P3-B rank amplitudes used for the independent amplitude comparison.
- `diagnostics/`: q=0 limit convergence comparison and the initial exact-zero indeterminate diagnostic log.
- `inputs/density_matrices/`: seven input matrix CSVs and extraction metadata.
- `source/`: pinned AFH package plus native/Python runners and comparators.
- `run_manifest.json`: source, environment, grid, conventions, tolerances, and statuses.
- `SHA256SUMS.txt`: hashes of all archived files except the checksum list itself.

Reproduce from the repository root:

```powershell
$env:DMFormFactorPath = (Resolve-Path .\v6\dmformfactor-V6.m).Path
$env:AFH_FULL_RESPONSE_OUTDIR = 'path-to-output\native'
wolframscript -file .\research\afh_full_response_database\run_afh_full_response_database.wls
python .\research\afh_full_response_database\run_afh_full_response_database.py --outdir path-to-output\python_reference
python .\research\afh_full_response_database\validate_afh_full_response_database.py --csv path-to-output\python_reference\afh_v6_full_eight_response_7isotope.csv --out-json path-to-output\python_reference\validation.json
python .\research\afh_full_response_database\compare_native_afh_full_response.py --python-csv path-to-output\python_reference\afh_v6_full_eight_response_7isotope.csv --native-dir path-to-output\native --out-json path-to-output\native\native_python_comparison.json
```

The native response calculation itself depends on the package's `SetIsotope` default tables. The Python reference generator additionally reads the archived density-matrix CSV files.
