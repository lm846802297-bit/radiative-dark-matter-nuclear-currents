# Xe-129 / Xe-131 SN100PN signed monopole qualification

**Provisional one-body nuclear responses, not a cross section and not a baseline replacement.** We reconstruct the rank-zero \(M,\Phi'',M\Phi''\) channels from published SN100PN one-body density matrix elements (OBDMEs) on the unchanged AFH v6.0/Baldridge–Dalton (B&D) native grid. This is the complete even-multipole response in these channels for Xe-129 (\(J_i=1/2\)). Xe-131 (\(J_i=3/2\)) also allows \(J=2\), so the two Xe-131 model columns are **not a like-for-like full comparison**.

## Method and normalization

- Using [Abdel Khaleq et al., supplement §I, Eq. (1)](https://journals.aps.org/prd/supplemental/10.1103/PhysRevD.109.075036/Supplementary.pdf), set \(d=(2J_i+1)(2T+1)\), \(T=(N-Z)/2=|M_T|\). For diagonal orbital \(\alpha\), \(n_p+n_n=\Psi_\alpha^0\sqrt{2(2j+1)/d}\) and \(n_n-n_p=\Psi_\alpha^1\sqrt{6(2j+1)T/[d(T+1)]}\). The isovector sign is checked against the published species occupations. The full machine-readable files include isoscalar closed-core orbitals not in the supplement's valence tables. The reconstructed valence occupations are within 0.011 nucleons of §III.G. Raw reconstructed counts are (53.99705, 74.97266) and (53.99666, 76.96874), rather than exactly (54, 75/77). This is not radial integration error. Check the density provenance or native convention before promotion; here no rescaling or clipping is applied.
- Define \(s_\alpha=j(j+1)-\ell(\ell+1)-3/4\), and \(n_r=(N_{\rm HO}-\ell)/2\). With \(\mathbf L=-i\mathbf r\times\nabla\), the gradient expression in [AFH, Eq. (39)](https://arxiv.org/abs/1308.6288) gives the single-orbital factors \(\sqrt{4\pi}M_{00}=j_0(qr)\) and \(\sqrt{4\pi}\Phi''_{00}=-s_\alpha j_1(qr)/(qr)\), since \(\hat{\mathbf r}\cdot(\boldsymbol\sigma\times\nabla)=-i\,\boldsymbol\sigma\cdot\mathbf L/r\). Fix the mixed-response sign against the negative low-q cross term in SN100PN supplement §II.G; do not take the printed static sign in AFH Eq. (42) alone as a finite-q convention check.
- Average these factors with normalized HO radial density \(x^{\ell+1/2}e^{-x}[L_{n_r}^{\ell+1/2}(x)]^2\), \(x=(r/b)^2\), \(y=(q_A b/2)^2\). Use \(b^2=41.467/(45A^{-1/3}-25A^{-2/3})\ \mathrm{fm}^2\) and \(1\,\mathrm{fm}=5.0677\,\mathrm{GeV}^{-1}\) as in the baseline. Sum the factors against the species occupations; the signed amplitudes' outer product gives \(F^{NN'}\) in supplement §II, Eq. (2). This rank-zero Gram matrix is positive semidefinite and its pairwise Cauchy–Schwarz bounds saturate, up to roundoff. Keep both orders \(F_{M_p\Phi_n}\) and \(F_{M_n\Phi_p}\). The CSV gives one mixed product per pair, not twice that product.
- The frozen native AFH \(W^{00}\) retains \((2J_i+1)/(4\pi)\). To compare physical equal p/n coupling, multiply by \(16\pi/(2J_i+1)\) once, not just four. This gives AFH \(F_M(0)=16639.1643,17159.5838\), respectively, close to \(A^2=16641,17161\). Do not apply an additional normalization correction to the immutable benchmark.

## Results

The calculation covers 0–0.6 GeV at 0.25 MeV intervals, 2,401 nodes per isotope. The 32- and 80-point generalized Gauss–Laguerre and independent analytic Laguerre-polynomial/\({}_1F_1\) representations agree within \(4.3\times10^{-14}\) absolute in signed amplitudes. This checks integration, not nuclear many-body accuracy. At \(q_A=0.10\) GeV, the full Xe-129 comparisons in physical equal-p/n normalization are:

| Response | SN100PN | Native AFH/B&D |
| --- | ---: | ---: |
| \(F_M\) | 1734.84 | 1764.04 |
| \(F_{\Phi''}\) | 120.15 | 86.97 |
| \(F_{M\Phi''}\) | −456.55 | −391.69 |

This does not quantify a change in \(d\sigma/dE_R\); that requires the matched three-body current. Ratios are unstable near diffraction zeros. The reconstructed equal-p/n monopole \(M_{00}\) amplitude changes sign in 0.15800–0.15825 GeV for Xe-129 and 0.15650–0.15675 GeV for Xe-131.

The rounded published coefficients are a literature cross-check, **not a positivity-preserving response backend**. The source expressions agree with monopole species entries at \(q=0\) within 0.5% numerically. For Xe-131 that is only an endpoint compatibility check, not a separation of its \(J=2\) spin-orbit contribution, which can survive at zero momentum. At 0.15 GeV the printed Xe-129 equal-p/n coherent response is **−4.308**, while the signed-amplitude reconstruction gives **+8.940**. This is a cancellation in rounded coefficients, not a negative physical prediction. Do not subtract printed polynomials to infer the Xe-131 \(J=2\) remainder.

The next gate is a separately checked Xe-131 \(J=2\) reconstruction or unrounded native reference and resolution of the small particle-number residual. Then perform a fixed-kinematics three-body one-body comparison before extending to more isotopes or a consistent chiral-EFT calculation. The original-paper Fig. 3 baseline gate is not addressed here.

## Reproduce

Run locally from this directory (Python 3.12; network only for an empty input cache):

~~~bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python audit.py --inputs /tmp/xe-sn100pn-inputs --output /tmp/xe-sn100pn-replay
OPENBLAS_NUM_THREADS=1 python -m unittest -v
~~~

The upstream inputs are the [authors' full density files](https://github.com/Raghda-AK/DM-Form-Factor/tree/1395f0c161905341fc89a2978148b1ce7c9ee1dd) (GPL-3.0, not vendored), the [arXiv supplement source](https://arxiv.org/src/2311.15764v1), and the [native AFH tables pinned to PR #2 at 5204f1a](https://github.com/lm846802297-bit/radiative-dark-matter-nuclear-currents/tree/5204f1a78b32b6ff5b423cc5d278d781e88add5e/nuclear/afh/data/native_v6_2026-11-22/native). The script checks input SHA-256. Each CSV gives the unrounded binary64 values. The summary landmark vectors are in order \(M,\Phi'',M\Phi''\); the printed and native AFH columns include all allowed multipoles.

The two Git-tracked `.csv.xz` files are lossless archives of the native CSV outputs (one per isotope). Verify the tracked artifacts and decode the tables in a disposable directory:

~~~bash
sha256sum -c SHA256SUMS.txt
staging=$(mktemp -d)
mkdir "$staging/native"
for f in native/*.csv.xz; do xz -dc "$f" > "$staging/${f%.xz}"; done
cp RAW_CSV_SHA256SUMS.txt "$staging/"
(cd "$staging" && sha256sum -c RAW_CSV_SHA256SUMS.txt)
~~~

[Separate audit of the frozen AFH native charge normalization for seven isotopes](native_afh_audit/NATIVE_AFH_NORMALIZATION.md).
