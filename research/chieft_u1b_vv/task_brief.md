# Research Task Brief — χEFT matching of the universal VV source

**Task ID:** 20261124-chiral-vv-source-matching  
**Status:** completed (literature/source audit); operator derivation is the next task  
**Prepared on:** 2026-11-24  
**Prepared by:** ChatGPT  
**Assigned to:** ChatGPT  
**Research topic:** LZ three-body full-current nucleon matching / χEFT  
**Baseline repository:** lm846802297-bit/LZ-3-body-for-hdm  
**Baseline branch and commit:** research/xe131-native-spatial-preflight-20261123 @ 92fa04dc7beeffd3a53732edbf8814eca40ab7a4

## 1. Research question

Determine which χEFT one- and two-nucleon vector charge/current operators are relevant to the project’s universal source
\[
J^\mu_{uds}=\bar u\gamma^\mu u+\bar d\gamma^\mu d+\bar s\gamma^\mu s,
\]
and identify which published direct-detection or electromagnetic operators can and cannot be reused.

## 2. Primary literature

- Hoferichter, Klos & Schwenk (2015), arXiv:1503.04811, DOI 10.1016/j.physletb.2015.05.041: DM-channel-dependent one-/two-body chiral power counting through O(p3).
- Pastore et al. (2011), arXiv:1106.4539, DOI 10.1103/PhysRevC.84.024001: two-nucleon electromagnetic charge operator through one loop, including unitary/off-shell prescription dependence.
- Kölling, Epelbaum, Krebs & Meißner (2011), DOI 10.1103/PhysRevC.84.054008: isoscalar external vector source in the nucleon covariant derivative and short-range two-nucleon current construction.
- Krebs, Epelbaum & Meißner (2019), arXiv:1902.06839, DOI 10.1007/s00601-019-1500-5: electromagnetic four-current through fourth order and continuity equation.
- Hu et al. (2022), arXiv:2109.00193, DOI 10.1103/PhysRevLett.128.072502: VS-IMSRG/NAT xenon nuclear-structure method and public spin-dependent response archive.

## 3. Fixed source and conventions

- Quark source: universal three-flavor vector current J_uds^mu.
- Three-flavor low-energy matching: J_uds^mu = 3 J_B^mu, with J_B^mu = (1/3) sum_{q=u,d,s} qbar gamma^mu q.
- External interaction: dark-matter vector × quark/nucleon vector (VV).
- Effective nucleon normalization: g_VN = 3 g_Vq, with the existing production form-factor convention. After F_i^{uds}/3, the nuclear operator is normalized to J_B^mu.
- SU(2) source realization: r_mu = l_mu = 0 and v_mu^(s) = B_mu, where B_mu couples to the unit-baryon-charge current.
- Nuclear kinematics: q_N = q_A, t_N = -2 m_A E_R, Q_A^2 = 2 m_A E_R + E_R^2, exact AME2020 nuclear masses.
- Existing AFH/P3-C results remain the one-body impulse benchmark; no redefinition of those archived results.

## 4. Required analysis

1. Identify the chiral representation of J_uds^mu in the three-flavor low-energy theory and its U(1)_B interpretation.
2. Separate VV, AV, and electromagnetic source operators. Do not transfer an AV two-body current into VV.
3. Reconcile the power-counting conventions used in direct-detection χEFT and electromagnetic-current derivations.
4. Use the O(p4) VV two-body order in Hoferichter–Klos–Schwenk Table 1 as the target, and derive its time and spatial pieces for the matched U(1)_B source, including regulator and unitary-transformation conventions.
5. Project the SU(2) electromagnetic source formalism onto v_mu^(s) alone; isolate the isoscalar operator before applying the normalization relative to the electromagnetic isoscalar source.
6. Specify the continuity relation with the same NN+3N Hamiltonian and current.
7. Assess the Hu et al. VS-IMSRG xenon method as a pilot nuclear-structure route; identify what data/code are reusable and what is absent for VV.

## 5. Predeclared acceptance

- The u+d+s source is matched explicitly, including its normalization.
- Every cited two-body operator is assigned to the correct DM Lorentz channel and source.
- Chiral-order statements use the source paper’s stated counting; no order labels are mixed without conversion.
- Any electromagnetic operator proposed for reuse has an explicit flavor/isospin projection to the baryon-number source.
- The isoscalar external-source normalization is derived from the covariant derivative, not obtained by scaling the full electromagnetic current.
- A specific NN+3N interaction, regulator, LEC set, many-body solver and required density/operator inputs are identified as the next numerical prerequisites.
- No numerical χEFT completion flag is set by this literature audit.

## 6. Result

The u+d+s source is the three-flavor baryon-number vector source, J_uds^mu = 3 J_B^mu. In Hoferichter–Klos–Schwenk’s direct-detection counting, the AV two-body vector current appears at O(p2), with AV NLO corrections at O(p3). The first VV two-body time component appears at O(p4); the spatial component is O(p2) in the nucleon current plus the O(p2) WIMP-spatial suppression, also effective O(p4). The AV operator must not be imported into this VV calculation.

Kölling et al. provide the SU(2) isoscalar external-source slot v_mu^(s) in the nucleon covariant derivative. The baryon source is obtained by setting r_mu = l_mu = 0 and normalizing v_mu^(s) to unit nucleon baryon charge. Relative to the isolated electromagnetic isoscalar source v_mu^(s) = e A_mu/2, this is a factor of two; it is not a rescaling of the full electromagnetic current.

Hu et al. provide a VS-IMSRG/NAT method benchmark for Xe-129/Xe-131, but their public repository supplies spin-dependent response fits, not the required vector charge/current operators or VV matrix elements. Kölling et al.'s source formalism gives a concrete projection: keep the tau1 dot tau2 one-pion-exchange charge terms, set the contact-current charge matrix to the identity (retaining C2/C4/C5/C7 terms), drop isovector L1 and retain isoscalar L2, with beta8/beta9 fixed consistently with the NN potential.

## 7. Exact scope

This task closes the literature/source-channel and external-source audit. It does not derive the complete U(1)_B χEFT operator, generate a χEFT Hamiltonian or xenon states, or complete the full-current/T3 numerical calculation.

## 8. Next task

Derive the source-specific U(1)_B χEFT charge/current operators through the declared O(p4) VV order for one pinned NN+3N interaction and regulator, isolate the isoscalar source terms from the published electromagnetic construction, and verify the corresponding continuity relation before beginning the Xe-129/Xe-131 VS-IMSRG calculation.


## 9. Reusable implementation audit

- IMSRG++ is the selected VS-IMSRG code candidate: it contains NAT-basis transformations, ensemble normal ordering, valence-space evolution and optional IMSRG(3) machinery.
- NuHamil is the candidate NN+3N matrix-element and operator generator. Its existing `Tmag_2B_J1_Tz0_Q100-NLO` interface is a rank-one isovector magnetic vector-current implementation (the `Tz0` suffix specifies the isospin projection, not isospin rank); it is not a U(1)_B current. Its multipole/HO transformation infrastructure may still be reused.
- NuHamil's current vector implementation is isovector, explicitly does not implement N3LO vector currents, and does not implement the two-body charge operator. It therefore does not supply the required isoscalar U(1)_B O(p4) VV current.
- Hu et al.'s public notebook contains fitted axial `Sigma'`/`Sigma''` responses for Xe-129/Xe-131, not the raw NN+3N Hamiltonian, transformed wavefunctions, transition densities or VV operators.
- The exact Hu `1.8/2.0(EM)` NN+3N matrix-element inputs are not included in the inspected public repository. NuHamil's bundled inputs do not directly provide a ready-made matching Hamiltonian, and its README notes incomplete support for non-locally regulated 3N interactions.

The implementation split is IMSRG++ for NAT/VS-IMSRG evolution, NuHamil for interaction/operator matrix elements, and a source-specific extension for the missing U(1)_B O(p4) charge/current. Exact interaction files and provenance must be obtained or regenerated before claiming reproduction of the Hu benchmark.
