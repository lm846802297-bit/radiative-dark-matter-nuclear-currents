# P4 native AFH regeneration

## Scope and provenance

These are local Windows computations against the hash-verified native AFH v6.0
archive at `results/2026-11-22/native_afh_eight_response/`. The P4 source was
checked out at `fada259460ac42460baf933021dc16dc7a9cfb9e` when the calculations
were run; `provenance.json` records source-file hashes and the runtime. The
source AFH package SHA-256 is
`dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372`.

The frozen inputs are (m_\chi=1000\) GeV, (m_s=200\) keV,
\(\delta m=1000\) keV, \(m_V=10\) GeV, \(g_2=1\), squared coupling
product \(10^{-8}\), local density \(0.3\) GeV/cm\(^3\),
`bishara_linear` nucleon form factors, Pole-width prescription, and the seven
AFH-supported xenon isotopes. Both halo models are `shm_cutoff` and
`mb_untruncated`.

The kinematic implementation is fixed to
\(t_A=-2m_AE_R\),
\(Q_A^2=2m_AE_R+E_R^2\), and mediator denominator
\((m_V^2-t_A)^2\). AFH uses \(Q_A^2\); the propagator uses \(t_A\).
The current assignment is the declared one-body impulse convention
\(q_N^\mu=q_A^\mu\), not a bound-nucleon off-shell reconstruction.

## Reproduction commands

Run from the repository root in PowerShell. Set the directory variables first:

```powershell
$b = 'research/three_body_full_current_nucleon_match'
$native = "$b/results/2026-11-22/native_afh_eight_response"
$out = "$b/results/2026-11-22/p4_native_afh_regeneration"
```

Five-point baseline, both halos, \(n_v=8\):

```powershell
python "$b/native_charge_sector_halo_fold.py" --native-dir $native --out "$out/baseline-five-recoil.json" --mchi-GeV 1000 --ms-keV 200 --delta-m-keV 1000 --mV-GeV 10 --g2 1 --coupling-product-squared 1e-8 --rho-GeV-cm3 0.3 --ER-keV 150 200 267 350 400 --nv 8 --nw 8 --nz 6 --nphi 4 --form-factor bishara_linear --width-scheme pole --include-leading-static-vector-candidate
```

Speed quadrature at 200 and 267 keV, holding the three-body quadrature at
\(8\times6\times4\):

```powershell
python "$b/native_charge_sector_halo_fold.py" --native-dir $native --out "$out/speed-refined.json" --mchi-GeV 1000 --ms-keV 200 --delta-m-keV 1000 --mV-GeV 10 --g2 1 --coupling-product-squared 1e-8 --rho-GeV-cm3 0.3 --ER-keV 200 267 --nv 12 --nw 8 --nz 6 --nphi 4 --form-factor bishara_linear --width-scheme pole --include-leading-static-vector-candidate
python "$b/native_charge_sector_halo_fold.py" --native-dir $native --out "$out/speed-ultra-refined.json" --mchi-GeV 1000 --ms-keV 200 --delta-m-keV 1000 --mV-GeV 10 --g2 1 --coupling-product-squared 1e-8 --rho-GeV-cm3 0.3 --ER-keV 200 267 --nv 24 --nw 8 --nz 6 --nphi 4 --form-factor bishara_linear --width-scheme pole --include-leading-static-vector-candidate
python "$b/native_charge_sector_halo_fold.py" --native-dir $native --out "$out/speed-48.json" --mchi-GeV 1000 --ms-keV 200 --delta-m-keV 1000 --mV-GeV 10 --g2 1 --coupling-product-squared 1e-8 --rho-GeV-cm3 0.3 --ER-keV 200 267 --nv 48 --nw 8 --nz 6 --nphi 4 --form-factor bishara_linear --width-scheme pole --include-leading-static-vector-candidate
```

Three-body phase-space refinement at 200 and 267 keV, holding \(n_v=8\):

```powershell
python "$b/native_charge_sector_halo_fold.py" --native-dir $native --out "$out/kernel-refined.json" --mchi-GeV 1000 --ms-keV 200 --delta-m-keV 1000 --mV-GeV 10 --g2 1 --coupling-product-squared 1e-8 --rho-GeV-cm3 0.3 --ER-keV 200 267 --nv 8 --nw 10 --nz 8 --nphi 6 --form-factor bishara_linear --width-scheme pole --include-leading-static-vector-candidate
```

Continuous spectrum: 81 recoil energies from 40 to 800 keV in 9.5 keV steps,
both halos, \(n_v=24\), three-body quadrature \(8\times6\times4\), split
into six independent recoil shards with no reduction in per-point quadrature:

```powershell
python "$b/run_p4_native_continuous_parallel.py" --native-dir $native --out "$out/continuous-81-point-spectrum.json" --workers 6
```

The old provisional outputs used for the matched comparison are the verified
Actions artifacts from run 44962224752 / artifact 14730173313 and run
44962224780 / artifact 14731493024. Both old and new folds use common linear
response interpolation. The comparison therefore measures the change from
the Python AFH input to the native AFH input at otherwise matched P4 settings;
it does not test a complete exchange-current model. The comparison script pins
the downloaded JSON payload SHA-256 values (`76b09aee…b41bb18` for the five
recoil points and `fdacf233…1e626c` for the continuous spectrum) and records
those hashes alongside the Actions artifact identifiers and artifact digests.
Substituted JSON files are rejected before comparison.

The generated audit, comparison, plots, and checksums are summarized in
`p4_native_regeneration_audit.json`, `matched-provisional-vs-native.json`,
`quadrature-comparison.json`, `plot_manifest.json`, `provenance.json`, and
`SHA256SUMS.txt`.

## Observed results

The native AFH gate passed 16,807 response rows and verified 50 archive files.
The five-point fold contains 10 spectrum rows / 70 isotope rows. The
continuous fold contains 162 spectrum rows / 1,134 isotope rows. Both report
`native_afh_validated=true`, `full_current_complete=false`.

The maximum relative rate differences against the old provisional Python-AFH
outputs are `2.6458134e-11` on the five-point grid and `1.5274858e-10` on the
continuous grid. Both use common linear interpolation. Speed 24→48-node
convergence is `8.9978914e-6` against `1e-5`; phase-space 8×6×4→10×8×6
convergence is `2.5525175e-9` against `1e-8`.

All sampled P4 response Cauchy excesses are zero. On the full raw native
table, 95 `M–Phi''` points exceed the strict diagnostic threshold
`max(1e-12, 2e-10*scale)`; maximum absolute and relative excesses are
`3.8635540e-10` and `3.1387149e-6`. The raw values are retained without
clipping. The continuous CSV audit passes with isotope-fraction closure
`2.22e-16`, D versus C plus spatial closure `3.79e-16`, and maximum total
spatial/Coulomb ratio `2.3901e-9`.
