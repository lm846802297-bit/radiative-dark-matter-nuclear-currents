# Pinned author SN100PN particle-count audit

This audit diagnoses the **input export**, not the WIMP response rank-2
reconstruction, a new shell-model diagonalization, or cross-section
data. Author machine-readable density files are precise-decimal GNU
GPL inputs and are linked, not vendored. They are not the original
NuShellX proton–neutron density files.

## Meaning and calculation

[Abdel Khaleq et al., arXiv:2311.15764v1, App. B](https://arxiv.org/html/2311.15764v1#A2)
gives the transformation from NuShellX proton–neutron reduced densities
to the isospin-reduced density. At J=0, write G=2j+1, D=(2Ji+1)(2T+1);
the diagonal combinations sum = n_p+n_n and diff = n_n−n_p are
obtained by multiplying Ψ^0 by sqrt(2G/D), Ψ^1 by
sqrt(6GT/[D(T+1)]), respectively, at T=|M_T|. The isovector
sign corresponds to the neutron excess and the published table.
Closed-core orbitals have a common proton and neutron occupation G.
For any normalized many-body eigenvector of the published fixed-Z/N
valence space, the number operator is exactly constant; restricting
configurations within that space does not change the trace.

Decimals are parsed exactly and weights are evaluated at 55 significant
digits. Per-field uncertainty is bounded twice: half of its last printed
place assuming rounding, or a full last place as a conservative bound
for truncation. These bounds are independent of radial integration.
Core and valence sums, both isospin-channel traces, printed versus
precise occupations and the downstream number-squared F_M(q=0)
cross-check are archived. No values have been corrected or clipped.

## Key findings

| A | Author-file proton trace | Author-file neutron trace | Maximum nearest-decimal error per species (whole file) | Maximum nearest-decimal error (valence only) |
|---:|---:|---:|---:|---:|
| 129 | 53.99704734745 | 74.97265716228 | 1.2035e-5 | 1.7664e-8 |
| 131 | 53.99665534160 | 76.96873764684 | 8.1432e-6 | 7.4015e-9 |

The expected pairs are (54,75) and (54,77). Restoring the closed core
analytically would remove only 1.16e−5 and 8.70e−6 per species,
respectively, and cannot explain the deficit. The valence-only
(S=n+p,D=n−p) deviations for Xe129/Xe131 are (-.03027226,-.02439019)
and (-.03458961,-.02791769). Separate factors needed to force the two
traces to exact values are different within each isotope by much
more than last-digit uncertainty, and differ between isotopes too.
None of these proxy factors is used. This rules out a common scalar
renormalization, an isovector sign/CG change alone, and decimal endpoint
rounding as explanations of both discrepancies. Other code or export
normalization issues upstream remain possible.

The paper's §IV.2.4 states that its Xe131 calculation was unrestricted
in the modeled shell, whereas Xe129 used the original iodine restriction
from §III.2, including a filled neutron 2d5/2 orbit. The latter trace
from the machine-readable file gives 5.993318324, not 6; the half-last-place
bound for that orbit is only 1.33e−8. Therefore neither restriction
nor any physical “missing neutron count” in the described **fixed-number**
model is identified by these numbers. A density-projection bookkeeping
difference cannot be narrowed to its cause without the author's native
pn NuShellX densities and reproducible settings. Isospin derived from
the published mapping is not a claim that the original wavefunction
has no isospin mixing. An invertible pn conversion by itself cannot
remove particle number when applied and inverted consistently.

Cross-checks:
- The [published supplemental tables](https://journals.aps.org/prd/supplemental/10.1103/PhysRevD.109.075036/Supplementary.pdf)
  cannot always be rounded versions of the pinned higher-precision
  files: Xe131 neutron 1g7/2 is printed as 7.71 but the reconstructed
  value is 7.6999612674 (rounds to 7.70). This does not identify which
  stage of data generation created the difference.
- For F_M at q=0 higher multipoles vanish. Xe131 n,n is printed 5920
  (nominal three significant figures). The precise file predicts
  5924.18657495, inside the plausible nearest rounding interval 5915–5925;
  exact 77²=5929 is outside. The raw neutron–neutron effect here is
  −8.1184e−4 relative. Xe129 checks 5620.89932197 versus printed 5620,
  with exact 75²=5625 at the rounding boundary. This is a consistency
  cross-check of the displayed outputs under the stated print-precision
  interpretation, not an authentication of native NuShellX.
- Seven isotope totals are recorded as controls, but only Xe129/131
  are focus of this qualification. A small excess for Xe136 neutron
  trace is not evidence for a global correction.

## Inputs and replay

- [Pinned author tree](https://github.com/Raghda-AK/DM-Form-Factor/tree/1395f0c161905341fc89a2978148b1ce7c9ee1dd).
  Provide its seven Xe shell-model '.dat' files in the input directory.
  Exact Git blobs and SHA-256 for each are checked by the program.
- [arXiv TeX source](https://arxiv.org/src/2311.15764v1), extract only
  member 'tek/bodysupplementary.tex'. Its SHA-256 is checked.
- Python standard library only; no network is used during replay.

```sh
python -m unittest -v test_audit_count
python audit_count.py /path/to/author_data /path/to/bodysupplementary.tex number_trace.json
sha256sum -c SHA256SUMS.txt
```
