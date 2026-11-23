# AFH response inputs

The imported Native AFH v6.0 dataset is stored at
\`data/native_v6_2026-11-22/\`. It contains all eight response functions for
seven xenon isotopes on the 2,401-point \(q_A\in[0,0.6]\) GeV grid, plus
Native rank-resolved \(\Delta_J,\Sigma'_J\) amplitudes for Xe-129 and Xe-131.

The Native comparison passed 16,807 isotope/q rows at
\(r_{\rm tol}=5\times10^{-7}\), \(a_{\rm tol}=2\times10^{-12}\). The spatial
amplitude comparison passed 9,604 rank/q rows. See the dataset's REPORT.md,
run_manifest.json, and SHA256SUMS.txt.

The imported CSVs retain raw AFH isoscalar \(\tau_1=\tau_2=0\) normalization.
No equal-\(p=n\) factor-four conversion is applied in the raw columns.
The raw zero-momentum coherent response is \(A^2/4\); the physical equal
proton/neutron coupling convention gives \(A^2\). Apply the factor four
exactly once in the downstream conversion. Signed interference is retained.
\(W_M\) is the finite-size nuclear response, not a point-nucleus response.

The AFH DMFormFactor v6.0 package itself is an upstream source asset and is
not redistributed here. Its exact SHA-256 is pinned in the imported
run_manifest.json. Reproduction requires a separately obtained package
whose bytes match that hash.

Native AFH validation is a nuclear-input result only. It does not imply
that the complete P3-C bound-nucleus current or a full chiral-EFT current
has passed.
