# AFH response inputs

AFH v6.0 supplies elastic nuclear multipoles including \(M,\Sigma'',\Sigma',\Phi'',\widetilde\Phi',\Delta\) and associated interference responses.

For each import, record AFH source version/hash, isotope, spin, oscillator length, OBDME source, q-grid, and Fourier convention. Keep raw AFH package normalization separate from physical proton/neutron coupling normalization. In the inherited equal-\(p=n\) convention, the raw isoscalar zero-momentum response is \(A^2/4\), while the physical coupling convention gives \(A^2\); apply the factor of four exactly once. Apply any \(C_0(A)\) correction exactly once.

Validate native output against an independent implementation and retain signed interference. \(W_M\) is finite-size nuclear response, not a point-nucleus response.

No AFH response CSVs are copied during bootstrap. Import only with source license/provenance and a reproducible comparison.

## Imported Native AFH v6.0 xenon dataset

The curated dataset is stored in
\`data/native_v6_2026-11-22/\`: seven isotopes, 2,401 \(q_A\) points per
isotope, all eight AFH responses, and 9,604 rank-resolved Xe-129/Xe-131
spatial multipole rows. The Native comparison passed with
\(r_{\rm tol}=5\times10^{-7}\), \(a_{\rm tol}=2\times10^{-12}\), and zero
maximum imaginary residual. Raw AFH isoscalar normalization is retained;
the factor four is applied only once downstream. See the dataset README
and SHA-256 manifest. This validates nuclear response inputs, not the
complete P3-C current.
