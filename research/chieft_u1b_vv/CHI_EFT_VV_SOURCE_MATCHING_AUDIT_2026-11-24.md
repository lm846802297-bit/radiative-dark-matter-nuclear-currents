# χEFT source matching for the universal VV interaction

**Date:** 2026-11-24  
**Research branch:** vector–vector (VV), universal light-flavor source  
**Status:** source-level matching audit complete; numerical χEFT current not yet constructed  
**Promotion:** \`full_current_complete=false\`

## 1. Exact quark-level source

The P3-C/P4 production convention uses the universal three-flavor vector source
\[
J^\mu_{uds}=\bar u\gamma^\mu u+\bar d\gamma^\mu d+\bar s\gamma^\mu s.
\]
In a three-flavor low-energy theory this is three times the baryon-number current,
\[
J^\mu_{uds}=3J^\mu_B,\qquad
J^\mu_B=\frac13\sum_{q=u,d,s}\bar q\gamma^\mu q.
\]
Thus the source to be matched is a flavor-singlet \(U(1)_B\) vector source, not the electromagnetic current. Pions carry zero baryon number; nucleons carry baryon number one.

At the nucleon level the production normalization is already fixed by the external \(g_{VN}=3g_{Vq}\) convention. The Borah2020+Sufian2017 form factors remain the current one-body phenomenological baseline. They must not be counted again as an independent factor of three when constructing the χEFT source.

## 2. Channel-specific power counting

Hoferichter, Klos and Schwenk derive the one- and two-body currents for dark-matter interactions through \(\mathcal O(p^3)\). Their explicit vector-current result states that the only two-body vector current through this order occurs in the axial-DM/vector-nucleon (AV) channel. It is not a two-body correction to the present vector-DM/vector-nucleon (VV) channel.

Table 1 of that paper makes the next step explicit: the VV two-body time component first appears at \(\mathcal O(p^4)\). The spatial vector component is \(\mathcal O(p^2)\) in the nucleon current but carries the additional two-order suppression of the WIMP spatial bilinear, so it also enters the VV amplitude at effective \(\mathcal O(p^4)\). The \(\mathcal O(p^2)\) two-body vector current and \(\mathcal O(p^3)\) NLO correction belong to AV, not VV.

Consequently, the AV pion-exchange current must not be imported into this project’s VV amplitude. The existing exact Dirac–Pauli P3-C result remains a one-body impulse calculation; it is not missing the published AV operator.

## 3. Higher-order charge operators and the matching caveat

The electromagnetic two-nucleon charge operator has been derived beyond the leading one-body term. Pastore et al. derive the two-nucleon electromagnetic charge operator through \(eQ\) (N4LO in their counting), with off-shell/unitary-prescription dependence in individual N3LO/N4LO pieces; physical predictions are invariant when the potential and charge operator are transformed consistently.

Krebs, Epelbaum and Meißner derive electromagnetic one-, two- and three-nucleon four-current contributions through fourth order and verify the corresponding modified continuity equation. Their analysis makes the accompanying unitary transformation and current regularization part of the operator definition.

These results are relevant construction references, not ready-made operators for \(J^\mu_{uds}\). Electromagnetic charge is a different flavor source. The required calculation is to introduce the \(U(1)_B\) external source into the same χEFT Lagrangian/regularization used for the NN+3N interaction, derive its charge and current operators, and only then determine which electromagnetic isoscalar structures map onto the baryon source.

The continuity relation must be checked in the selected convention,
\[
\mathbf q\cdot\mathbf J(\mathbf q)=[H,\rho(\mathbf q)],
\]
with the appropriate energy-transfer terms when the chosen operator convention requires them. A continuity-reduced longitudinal current may be used only after this relation is verified for the same Hamiltonian and regulator.

## 4. Nuclear-structure route already demonstrated for xenon

Hu et al. demonstrate a VS-IMSRG route for \(^{129,131}\mathrm{Xe}\) using χEFT NN+3N interactions, natural-orbital bases and consistently transformed response operators. Their benchmark settings include the \(^{100}\mathrm{Sn}\)-core valence space
\[
\{0g_{7/2},1d_{5/2},1d_{3/2},2s_{1/2},0h_{11/2}\},
\]
with \(e^{\mathrm{NAT}}_{\max}=14\), \(e_{\max}=8\), \(\hbar\omega=16\,\mathrm{MeV}\), and \(E_{3\max}=22\) for the heavy nuclei. They compare 1.8/2.0(EM), N4LO+3Nlnl and \(\Delta\)NNLOGO interactions.

Their public [AbInitio4DarkMatter repository](https://github.com/BaishanHu/AbInitio4DarkMatter) contains a spin-dependent notebook and fitted \(\Sigma'\)/\(\Sigma''\) response data for xenon. It does not provide the vector-charge/current operators or the required VV response. It is reusable as a method and regression reference, not as a numerical substitute for the current project.

The paper also reports that heavy-nucleus spectra are not uniformly reproduced, likely reflecting limitations in deformation. This makes a direct comparison with shell-model benchmarks important when selecting the pilot Hamiltonian.

## 5. Explicit SU(2) external-source realization

The SU(2) heavy-baryon chiral Lagrangian already has the required singlet source. Kölling et al. define
\[
D_\mu N_v=\left(\partial_\mu+\Gamma_\mu-i v_\mu^{(s)}\right)N_v,
\qquad
v_{\mu\nu}^{(s)}=\partial_\mu v_\nu^{(s)}-\partial_\nu v_\mu^{(s)}.
\]
For the project’s baryon-number-normalized source, choose
\[
r_\mu=\ell_\mu=0,\qquad v_\mu^{(s)}=B_\mu,
\]
where \(B_\mu\) is the external source coupled to \(J_B^\mu\). Since \(v_\mu^{(s)}\) is proportional to the identity in isospin space, the pion covariant derivative has no \(B_\mu\) coupling: pions are baryon-neutral. The source enters the nucleon covariant derivatives and the isoscalar field-strength structures \(f_+^{\mu\nu}\), with \(f_-^{\mu\nu}=0\) for this pure singlet source.

For comparison, Kölling et al. use the electromagnetic source
\[
r_\mu=\ell_\mu=\frac{e\tau_3}{2}A_\mu,\qquad
v_\mu^{(s)}=\frac{e}{2}A_\mu.
\]
Therefore the baryon source is obtained by retaining the \(v_\mu^{(s)}\) part alone and normalizing its nucleon charge to one. In this source convention it is twice the isolated isoscalar electromagnetic source; it is **not** twice the full electromagnetic current. The universal \(u+d+s\) source in the project is \(3J_B^\mu\), while the production factor \(g_{VN}=3g_{Vq}\) and \(F_i^{uds}/3\) leave the nuclear operator normalized to \(J_B^\mu\).

This gives a concrete extraction rule for existing χEFT Lagrangians: set the isovector sources \(r_\mu,\ell_\mu\) to zero, keep \(v_\mu^{(s)}\), and functionally differentiate the same gauged NN+3N Lagrangian with respect to \(B_\mu\). In the electromagnetic contact-current expressions, this corresponds to projecting the charge matrix onto its isoscalar identity part before applying the factor of two. Isovector terms proportional to \(\tau_3\) or \([\boldsymbol\tau_1\times\boldsymbol\tau_2]_3\) are not part of the baryon source. The isoscalar short-range terms generated by the covariant derivatives must be retained; their LECs are those of the same NN interaction, not new electromagnetic fit parameters.

Pastore et al. explicitly note that the earlier two-nucleon charge operator was restricted to its isoscalar component, while the two-pion-exchange loop charge corrections are isovector. This is the relevant source decomposition to investigate for \(U(1)_B\). Kölling et al. provide the covariant derivative and contact-current realization; Krebs et al. provide the fourth-order continuity/unitary-transformation framework. The actual operator extraction and regulator-consistent reduction remain to be carried out.


## 6. Concrete operator projection to carry into the derivation

Kölling et al.'s contact Lagrangian (their Eq. 3.11) writes the electromagnetic charge matrix as
[
hat e=rac12+rac{	au_3}{2}.
]
For the baryon source, the source matrix is the identity, (hat e	o mathbf 1), because every nucleon has baryon number one. The short-range spatial-current terms proportional to (C_2,C_4,C_5) therefore survive with this replacement. In the magnetic contact term, the (L_1) structure is isovector and drops out, while the (L_2) structure is isoscalar and remains. The (C_i) and (L_2) values must come from the same NN interaction/LEC fit used in the Hamiltonian.

Their one-pion-exchange charge density (Eq. 4.26) contains both an isovector piece proportional to (	au_3/2) and an isoscalar piece proportional to (oldsymbol	au_1cdotoldsymbol	au_2), together with terms proportional to ([oldsymbol	au_1	imesoldsymbol	au_2]_3). The baryon-source projection keeps the (oldsymbol	au_1cdotoldsymbol	au_2) pieces and removes the (	au_3) and cross-product pieces. The surviving terms depend on the unitary phases (areta_8,areta_9). Those phases cannot be chosen independently of the NN potential: the same unitary convention must be used in the Hamiltonian and charge operator.

This identifies a concrete starting operator set for the (U(1)_B) derivation:
- isoscalar one-pion-exchange charge terms from the (oldsymbol	au_1cdotoldsymbol	au_2) projection;
- short-range vector-current terms obtained from the (C_2,C_4,C_5) contact structures with (hat e=mathbf 1);
- the isoscalar magnetic contact term proportional to (L_2);
- any higher-order source-insertion and three-nucleon charge terms required by the chosen truncation.

The electromagnetic formulas are not copied wholesale. The source projection, normalization, unitary phases, and regulator must be applied at the Lagrangian/operator level, then the continuity relation must be checked against the chosen NN+3N Hamiltonian.


## 7. Reusable many-body software and the exact implementation boundary

The public code review identifies a practical software stack:

- **IMSRG++** is the preferred VS-IMSRG engine. Its public C++ implementation supports Hartree–Fock, single-reference and valence-space IMSRG, natural-orbital transformations, ensemble normal ordering, and optional IMSRG(3) commutators. The public `GetNormalOrderedH.py` example explicitly reads NN and 3N matrix elements and demonstrates an `EM1.8_2.0` 3N HDF5 filename convention. Repository: [ragnarstroberg/imsrg](https://github.com/ragnarstroberg/imsrg); see the [NAT implementation](https://github.com/ragnarstroberg/imsrg/blob/d8357bd08865cb7087bdb148d8e68e3be2349266/src/HFMBPT.cc), [VS-IMSRG driver](https://github.com/ragnarstroberg/imsrg/blob/d8357bd08865cb7087bdb148d8e68e3be2349266/src/imsrg++.cc), and [NN+3N input example](https://github.com/ragnarstroberg/imsrg/blob/d8357bd08865cb7087bdb148d8e68e3be2349266/work/scripts/python/GetNormalOrderedH.py).
- **NuHamil** is a complementary matrix-element and operator generator, not the VS-IMSRG solver. It generates NN+3N matrix elements and has an existing rank-1 magnetic vector two-body multipole interface, `Tmag_2B_J1_Tz0_Q100-NLO`. Despite the `Tz0` label, the source code contracts the rank-1 isospin tensor `tau1_tau2_tensor_iso(...,1)` (equivalently the isovector cross product in the proton-neutron representation). It is an isovector operator, not the isoscalar baryon-number current. Repository: [Takayuki-Miyagi/NuHamil-public](https://github.com/Takayuki-Miyagi/NuHamil-public); the isovector structure is explicit in [`M1_2b_current.F90`](https://github.com/Takayuki-Miyagi/NuHamil-public/blob/eb65165f5ecf11fde6e188a65385c98e7526ccb6/src/TwoBody/M1_2b_current.F90), and the multipole-order limits are in [`TwoBodyMultiPoleOperator.F90`](https://github.com/Takayuki-Miyagi/NuHamil-public/blob/eb65165f5ecf11fde6e188a65385c98e7526ccb6/src/TwoBody/TwoBodyMultiPoleOperator.F90).
- **Hu et al.'s public notebook/data** are a regression reference for the nuclear solver and axial responses. They contain fitted (Sigma') and (Sigma'') response functions for Xe-129/Xe-131 and the published basis/truncation choices, but not the raw NN+3N matrix elements, transformed wavefunctions, one-/two-body transition densities, or vector VV operators.

The source review establishes two limitations of the existing NuHamil vector-current implementation. First, its current `Tmag_2B` multipole is isovector (rank-one in isospin), so it does not contribute to the pure U(1)_B source and must not be inserted into the VV calculation. Second, its `TwoBodyMultiPoleOperator` supports only the implemented NLO vector magnetic/electric two-body multipoles, explicitly reports that N3LO vector currents are not implemented, and does not implement the two-body charge branch. The multipole-generation and HO transformation infrastructure may be reused, but the U(1)_B isoscalar O(p4) charge/current operators must be derived and implemented independently.

The public Hu repository does not include the exact `EM1.8_2.0` NN+3N input files. IMSRG++'s example script shows how such files are consumed, but that is not a reproducible input archive. NuHamil's bundled interaction files include unevolved N3LO/N2LO NN inputs rather than a ready-made Hu-compatible `1.8/2.0(EM)` Hamiltonian. Its README also warns that non-locally regulated 3N interactions are not fully supported. We therefore cannot claim an exact reproduction of Hu's Hamiltonian from these public repositories alone.

**Recommended implementation split:** use IMSRG++ for the NAT/VS-IMSRG evolution and operator transformation; use NuHamil as a matrix-element generator and as a reference implementation for multipole generation/HO transformation; its currently exposed NLO Tmag two-body operator is isovector and is not a U(1)_B current. First obtain or reconstruct and fingerprint the exact NN+3N Hamiltonian used for the chosen benchmark. Then extend the operator generator with the source-projected U(1)_B O(p4) charge/current pieces, and validate the new operator against the continuity relation before any Xe production run.

## 8. Consequence for the current project

The P3-C direct current and the full-complex T3 contraction have already been exercised for the declared one-body impulse model. The remaining χEFT closure is not achieved by adding the published AV two-body current. It requires:

1. deriving the \(U(1)_B\) source operators for the actual VV interaction;
2. fixing the χEFT order, regulator, NN+3N interaction and LEC set;
3. generating Xe-129/Xe-131 states and consistently transformed one-/two-body operators/densities;
4. evaluating the complex current multipoles before squaring;
5. contracting those multipoles with the existing exact-mass Pole-width T3 tensor, then independently replaying the full complex contraction;
6. checking the continuity relation, regulator consistency, chiral truncation, basis convergence and many-body uncertainties.

The 2022 VS-IMSRG work provides a credible starting point for the nuclear solver and Xe valence space. It does not supply the missing VV current or its many-body matrix elements.

## 9. References

1. M. Hoferichter, P. Klos and A. Schwenk, “Chiral power counting of one- and two-body currents in direct detection of dark matter,” *Phys. Lett. B* **746** (2015) 410, [arXiv:1503.04811](https://arxiv.org/abs/1503.04811), DOI: [10.1016/j.physletb.2015.05.041](https://doi.org/10.1016/j.physletb.2015.05.041).
2. S. Pastore, L. Girlanda, R. Schiavilla and M. Viviani, “The two-nucleon electromagnetic charge operator in chiral effective field theory up to one loop,” *Phys. Rev. C* **84** (2011) 024001, [arXiv:1106.4539](https://arxiv.org/abs/1106.4539), DOI: [10.1103/PhysRevC.84.024001](https://doi.org/10.1103/PhysRevC.84.024001).
3. S. Kölling, E. Epelbaum, H. Krebs and U.-G. Meißner, “Two-nucleon electromagnetic current in chiral effective field theory: One-pion exchange and short-range contributions,” *Phys. Rev. C* **84** (2011) 054008, DOI: [10.1103/PhysRevC.84.054008](https://doi.org/10.1103/PhysRevC.84.054008).
4. H. Krebs, E. Epelbaum and U.-G. Meißner, “Nuclear electromagnetic currents to fourth order in chiral effective field theory,” *Few-Body Syst.* **60** (2019) 46, [arXiv:1902.06839](https://arxiv.org/abs/1902.06839), DOI: [10.1007/s00601-019-1500-5](https://doi.org/10.1007/s00601-019-1500-5).
5. B.-S. Hu et al., “Ab Initio Structure Factors for Spin-Dependent Dark Matter Direct Detection,” *Phys. Rev. Lett.* **128** (2022) 072502, [arXiv:2109.00193](https://arxiv.org/abs/2109.00193), DOI: [10.1103/PhysRevLett.128.072502](https://doi.org/10.1103/PhysRevLett.128.072502).

## 10. Decision and next gate

**Decision:** do not add the published AV two-body vector current to the VV calculation. Retain the exact-mass Dirac–Pauli P3-C result as the validated one-body impulse benchmark. The next physics task is source-specific \(U(1)_B\) χEFT matching and an explicit order-by-order operator derivation.

**Next acceptance gate:** produce the \(U(1)_B\) one- and two-nucleon vector charge/current operators for one declared χEFT Hamiltonian/regulator, demonstrate the matching to the universal \(u+d+s\) source, and verify the continuity relation before any xenon many-body calculation is promoted.
