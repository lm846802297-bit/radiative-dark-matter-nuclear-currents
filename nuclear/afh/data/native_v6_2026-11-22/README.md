# Native AFH v6.0 xenon response dataset

This is a curated import of the Native Wolfram AFH v6.0 outputs generated
on 2026-11-22 from the pinned original-project source commit
dd9a7387814cd157d9cde739fa80181e4cd7d428.

## Coverage

- Isotopes: Xe-128, Xe-129, Xe-130, Xe-131, Xe-132, Xe-134, Xe-136.
- Momentum grid: 2,401 points per isotope, \(q_A\in[0,0.6]\) GeV,
  \(\Delta q_A=0.00025\) GeV.
- Total response rows: 16,807.
- Responses: \(M,\Sigma'',\Sigma',\Phi'',\widetilde\Phi',\Delta,
  M\Phi'',\Sigma'\Delta\).
- Rank-resolved Native \(\Delta_J,\Sigma'_J\): 9,604 Xe-129/Xe-131 rows.
- Native imaginary residual: zero.
- Wolfram: 14.2.0 for Windows.

## Normalization

The CSV columns use the raw AFH isoscalar package convention
\(\tau_1=\tau_2=0\). No equal-\(p=n\) factor-four conversion is applied
to these raw columns. The raw zero-momentum coherent response is \(A^2/4\);
the physical equal-proton/equal-neutron convention gives \(A^2\). Apply
the factor four exactly once downstream. Signed \(M\)-\(\Phi''\) and
\(\Sigma'\)-\(\Delta\) interference is retained.

The AFH argument is \(y_A=(b_Aq_A/2)^2\), with \(q_A\) the spatial nuclear
recoil momentum. The invariant transfer \(t_A=-2m_AE_R\) is separate and
is not used as the AFH spatial argument.

## Validation boundary

The Native response database passed all 16,807 isotope/q rows against the
source-translated Python reference at \(r_{\rm tol}=5\times10^{-7}\) and
\(a_{\rm tol}=2\times10^{-12}\). Maximum combined-tolerance utilization
was 0.05879. The 9,604 rank-resolved spatial amplitudes also passed their
Native comparison.

These checks validate the AFH nuclear-response inputs and the listed
spatial multipole amplitudes. They do not validate the complete P3-C
particle--nuclear current contraction. The original run manifest retains
\`full_current_complete=false\`.

The upstream AFH package source is not redistributed. Its exact version
and SHA-256 are recorded in run_manifest.json. The two Native Wolfram
exporter scripts are retained under source/.

Run \`sha256sum --check SHA256SUMS.txt\` from this directory to verify
this curated dataset import.
