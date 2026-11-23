# AFH response inputs

AFH v6.0 supplies elastic nuclear multipoles including \(M,\Sigma'',\Sigma',\Phi'',\widetilde\Phi',\Delta\) and associated interference responses.

For each import, record AFH source version/hash, isotope, spin, oscillator length, OBDME source, q-grid, and Fourier convention. Keep raw AFH package normalization separate from physical proton/neutron coupling normalization. In the inherited equal-\(p=n\) convention, the raw isoscalar zero-momentum response is \(A^2/4\), while the physical coupling convention gives \(A^2\); apply the factor of four exactly once. Apply any \(C_0(A)\) correction exactly once.

Validate native output against an independent implementation and retain signed interference. \(W_M\) is finite-size nuclear response, not a point-nucleus response.

No AFH response CSVs are copied during bootstrap. Import only with source license/provenance and a reproducible comparison.
