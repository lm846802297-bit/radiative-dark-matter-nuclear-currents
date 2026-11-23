# Xe-131 SN100PN signed quadrupole reconstruction

**Scope:** signed one-body elastic \(M,\Phi'',M\Phi''\) responses at \(J=2\), together with their \(J=0+2\) totals. Those are all the allowed multipoles for these channels in the \(J_i=3/2\) ground state. This is not a three-body cross section, a refit of SN100PN, or a consistent chiral-EFT current calculation. The AFH v6.0/Baldridge–Dalton benchmark remains unchanged.

## Inputs and calculation

Use all 14 ordered \(J=2\) entries in each isospin channel of the [authors' pinned Xe-131 OBDME file](https://github.com/Raghda-AK/DM-Form-Factor/blob/1395f0c161905341fc89a2978148b1ce7c9ee1dd/131Xe_OBDME_SN100PN.dat), with its full published machine-readable precision. The reduction and normalization follow [Abdel Khaleq et al., Phys. Rev. D 109, 075036 (2024)](https://doi.org/10.1103/PhysRevD.109.075036) and [its supplement](https://journals.aps.org/prd/supplemental/10.1103/PhysRevD.109.075036/Supplementary.pdf). No wavefunction was generated and no response was constructed by subtraction from the rounded supplement. The two nearly antisymmetric pairs of published input entries give a maximum Hermiticity residual of \(3.46\times10^{-5}\) in the density units. Inputs are retained without symmetrization or rescaling. The small particle-number residual of the monopole remains as recorded in the [pinned rank-zero audit](https://github.com/lm846802297-bit/radiative-dark-matter-nuclear-currents/tree/51da16ceb636413c0a0c5e0b98f64224fd4a216d/nuclear/sn100pn/xe129-xe131-rank0).

The operators are \(M_{J0}=j_J(qr)Y_{J0}\) and \(\Phi''_{J0}=i(\nabla M_{J0})\cdot(\boldsymbol\sigma\times\nabla)/q^2\), from [AFH, Eq. (39)](https://arxiv.org/abs/1308.6288). Use Condon–Shortley spinor spherical harmonics and positive-at-origin HO radial functions. For each ordered pair of orbitals, project at every allowed magnetic sublevel and recover the reduced element by the Wigner–Eckart theorem. With \(z=qr\), the three spin-orbit radial kernels multiply the angular components as follows:

| Reduced angular component | Radial kernel (inside the radial integral) |
| --- | --- |
| \(iY_{J0}[\sigma_\theta(\partial_\phi/\sin\theta)-\sigma_\phi\partial_\theta]\) | \(R_aR_b j'_J(z)/z\) |
| \(i(\partial_\theta Y_{J0})\sigma_\phi\) | \(R_a(r\partial_r R_b)j_J(z)/z^2\) |
| \(-i(\partial_\theta Y_{J0})\sigma_r(\partial_\phi/\sin\theta)\) | \(R_aR_b j_J(z)/z^2\) |

At \(q=0\), \(M_{20}=0\), but \(\Phi''_{20}\) can be finite because \(j_2(z)/z^2\to1/15\) and \(j'_2(z)/z\to2/15\). The \(p,n\) amplitudes include \(\sqrt{4\pi/(2J_i+1)}\). Their outer products are added **over** \(J\), without interference **between** \(J=0\) and \(J=2\). Each order of the \(p,n\) mixed response is retained. Raw native AFH \(W^{00}\) is multiplied by \(16\pi/(2J_i+1)=4\pi\) for physical equal p/n coupling; no arbitrary renormalization is applied.

## Numerical checks and results

The native grid is 0–0.60 GeV, spacing 0.25 MeV, 2,401 nodes. Generalized Gauss–Laguerre and an independent analytic Laguerre-polynomial/\({}_1F_1\) radial representation differ by at most \(3.75\times10^{-16}\) in the signed \(J=2\) amplitudes. A simultaneous angular/radial node refinement gives \(9.30\times10^{-16}\); the generalized algorithm reproduces the pinned \(J=0\) signed amplitudes to \(8.53\times10^{-14}\). Hermitian reversal holds for the newly constructed \(M\) and \(\Phi''\) single-particle operators within \(1.6\times10^{-14}\); the total Gram matrix is positive semidefinite and rank at most two up to roundoff. These checks bound numerical integration and convention errors, not shell-model uncertainties.

At \(q_A=0\), \((M_{2,p},M_{2,n},\Phi''_{2,p},\Phi''_{2,n})=(0,0,0.01490036,0.10464768)\), giving an equal-p/n quadrupole spin-orbit response of 0.01429174. At \(q_A=0.10\) GeV:

| Equal p/n response | SN100PN \(J=2\) | SN100PN full | Native AFH/B&D full |
| --- | ---: | ---: | ---: |
| \(F_M\) | \(4.6269\times10^{-5}\) | 1730.7866 | 1760.1854 |
| \(F_{\Phi''}\) | \(8.1214\times10^{-4}\) | 129.2692 | 103.2012 |
| \(F_{M\Phi''}\) | \(-1.9385\times10^{-4}\) | −473.0076 | −425.7807 |

At the first equal-p/n \(J=0\) charge zero, \(q_A\simeq0.15656997\) GeV, the \(J=2\) coherent charge leaves \(F_M=2.7235\times10^{-6}\) instead of a full zero. This is an exact-model calculation at the isolated monopole zero, not the location of an experimental minimum. Rounded published response polynomials remain useful as a coarse literature check when discrepancies are normalized to their zero-momentum scale, but cannot resolve the quadrupole or the sign of an interference close to a cancellation. No cross-section conclusion follows from these numbers alone.

## Reproduce

Keep this directory adjacent to the pinned **xe129-xe131-rank0** directory (the replay checks the SHA-256 of its helper). Use a new output directory rather than overwriting the published archive. First run requires network access for pinned inputs.

~~~bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python rank2.py --inputs /tmp/xe131-inputs --output /tmp/xe131-replay
OPENBLAS_NUM_THREADS=1 python -m unittest -v
~~~

The native data are binary64 printed at round-trip precision. The CSV contains signed amplitudes \(J=0,2\), all ten upper-triangle species/operator products for both \(J=2\) and their total, and the equal-p/n responses for reconstructed SN100PN, native AFH/B&D, and the rounded published supplement. The JSON summarizes checkpoints, input pins and validation. Vectors in its landmarks are \(M_p,M_n,\Phi''_p,\Phi''_n\) for amplitudes and \(M,\Phi'',M\Phi''\) for equal-p/n responses. Upstream GPL-3.0 densities are accessed by hash and not vendored.

## Independent checks and tracked output

- [AFH-source angular recoupling and native B&D closure](independent_recoupling/README.md): all four SN100PN signed rank-2 amplitudes agree over the entire grid to \(1.20\times10^{-15}\) absolute. The old (J=0+2) responses replay against the frozen AFH Wolfram outputs; the JSON also supplies like-for-like B&D rank-2 landmarks.
- [Published OBDME number trace](number_trace/README.md): the Xe-131 neutron trace is 76.96873765 under the published transformation. This departure from 77 cannot be accounted for by the last printed digit or the closed core. Its upstream cause is not identified; no density was rescaled.

Verify the tracked bytes and reconstruct the lossless native table in a temporary directory:

~~~sh
sha256sum -c SHA256SUMS.txt
replay=$(mktemp -d)
mkdir -p "$replay/native"
cat native/xe131_rank2.csv.xz.part-* | xz -dc > "$replay/native/xe131_rank2.csv"
cp RAW_CSV_SHA256SUMS.txt "$replay/"
(cd "$replay" && sha256sum -c RAW_CSV_SHA256SUMS.txt)
~~~
