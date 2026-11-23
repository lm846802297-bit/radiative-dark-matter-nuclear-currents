# Literature review and reusable nuclear inputs

This is a curated starting bibliography, not an exhaustive search.

## Direct 2→3 scalar radiation

**Project benchmark manuscript:** Rundong Fang, Ji-Heng Guo, Ming Li, Jia Liu, and Xiao-Ping Wang, “High-Energy Nuclear Recoils from 2 → 3 Dark Matter Scattering.” The supplied manuscript studies \(\chi_1N\to\chi_1N+s\), with a neutral scalar radiated from the dark-matter line before or after vector exchange. It establishes the extra \(\chi_1s\) invariant-mass continuum and the matched recoil-endpoint comparison to inelastic \(2\to2\) scattering when \(\delta=m_s\). The source PDF is not copied into this repository.

David Curtin and Yuhsin Tsai, “The Double-Dark Portal,” JHEP 11 (2014) 136. [arXiv:1405.1034](https://arxiv.org/abs/1405.1034), [DOI](https://doi.org/10.1007/JHEP11(2014)136). A direct $2\to3$ scalar-emission precedent, $\chi N\to\chi N\phi$, in which the light scalar is emitted at a quark–two-scalar portal (Eq. (2.1) and Fig. 2). The topology is established; the interaction does not give the present VV two-diagram amplitude or its AFH nuclear contraction.

## Related processes with different radiation or initial states

Chris Kouvaris and Josef Pradler, “Probing Sub-GeV Dark Matter with Conventional Detectors,” Phys. Rev. Lett. 118, 031803 (2017). [arXiv:1607.01789](https://arxiv.org/abs/1607.01789), [DOI](https://doi.org/10.1103/PhysRevLett.118.031803). Uses photons emitted from the recoiling target to access sub-GeV dark matter. It is a comparison for nuclear-recoil bremsstrahlung, not for scalar radiation from the dark-matter line.

Nicole F. Bell et al., “Migdal effect and photon bremsstrahlung in effective field theories of dark matter direct detection and coherent elastic neutrino-nucleus scattering,” Phys. Rev. D 101, 015012 (2020). [arXiv:1905.00046](https://arxiv.org/abs/1905.00046), [DOI](https://doi.org/10.1103/PhysRevD.101.015012). Useful for radiation-plus-EFT methodology, but its photon is associated with the recoiling atom/nucleus rather than dark-matter-line scalar emission.

S. Åminneborg, L. Bergström, and B.-Å. Lindholm, “Radiative scattering of neutralino dark matter,” Phys. Rev. D 43, 2527 (1991). [DOI](https://doi.org/10.1103/PhysRevD.43.2527). An early example of photino–nucleus scattering accompanied by direct photon emission, distinct from the scalar process here.

Wei Chao, Mingjie Jin, and Ying-Quan Peng, “Directly detecting sub-MeV dark matter via 3-body inelastic scattering process,” Phys. Rev. D 107, 093009 (2023). [DOI](https://doi.org/10.1103/PhysRevD.107.093009). This is $3\to2$ with two incoming dark-matter particles and a target, not $2\to3$ radiation.

## Dirac–Pauli and AFH precedent: magnetic upscattering followed by decay

Joshua Eby, Patrick J. Fox, and Graham D. Kribs, “Earth-catalyzed detection of magnetic inelastic dark matter with photons in large underground detectors,” JHEP 06 (2024) 165. [arXiv:2312.08478v2, Secs. 2.1–2.2 and App. B](https://arxiv.org/html/2312.08478v2), [DOI](https://doi.org/10.1007/JHEP06(2024)165). Magnetic dipole scattering produces an excited dark-matter state in $\chi_1+A\to\chi_2+A$; its subsequent decay $\chi_2\to\chi_1+\gamma$ supplies the photon. The photon is not emitted by the recoiling nucleus, and the paper does not calculate coherent scalar emission before and after vector exchange.

Equation (8) uses the electromagnetic Dirac–Pauli current with $F_1(0)=Q_N$ and $F_2(0)=\kappa_N$. Equations (12)–(20) describe their magnetic-transition NREFT matching; Appendix B builds the AFH nuclear responses from BIGSTICK one-body density matrices. The operator and density-matrix workflow is useful as an independent check. The electromagnetic charges, form factors, and $2\to2$ dark-matter coefficients do not transfer to the project's universal $u+d+s$ vector source and radiative tensor.

The authors' [IsotopeResponses table](https://github.com/joshaeby/IsotopeResponses/blob/5bb6cc9b702e88b609913a053c411a5949c3b67a/Nuclear_Response_Data.dat) supplies nuclear-response polynomial coefficients for 20 isotopes, eight responses including interference, and four ordered isospin pairs. The response table does not include xenon, in particular $^{129}\mathrm{Xe}$ or $^{131}\mathrm{Xe}$; it is a candidate cross-check on other nuclei, not a replacement for the project's xenon inputs.

## AFH nuclear responses

A. Liam Fitzpatrick, Wick Haxton, Emanuel Katz, Nicholas Lubbers, and Yiming Xu, “The Effective Field Theory of Dark Matter Direct Detection,” JCAP 02 (2013) 004. [arXiv:1203.3542v3](https://arxiv.org/abs/1203.3542v3), [DOI](https://doi.org/10.1088/1475-7516/2013/02/004). The original Galilean-invariant nonrelativistic operators and their nuclear responses provide a basis for matching. The $2\to2$ scattering coefficients are not a radiative $2\to3$ kernel.

Nikhil Anand, A. Liam Fitzpatrick, and W. C. Haxton, “Weakly interacting massive particle-nucleus elastic scattering response,” Phys. Rev. C 89, 065501 (2014). [arXiv:1308.6288](https://arxiv.org/abs/1308.6288), [DOI](https://doi.org/10.1103/PhysRevC.89.065501).

This is the AFH nuclear-response framework used for elastic nuclear multipoles. Reuse requires the exact AFH v6.0 normalization, Fourier convention, OBDME input, and interference signs. \(W_M\) is finite-size nuclear response, not point nucleus; \(\Phi''\) is not the Lorentz-scalar nucleon source.

## Chiral EFT one- and two-body currents

Martin Hoferichter, Philipp Klos, and Achim Schwenk, “Chiral power counting of one- and two-body currents in direct detection of dark matter,” Phys. Lett. B 746, 410–416 (2015). [arXiv:1503.04811](https://arxiv.org/abs/1503.04811), [DOI](https://doi.org/10.1016/j.physletb.2015.05.041). Establishes common chiral power counting for vector, axial-vector, scalar, and pseudoscalar WIMP-nucleon interactions and derives one-/two-body currents through third order. Section 5.2 states that the only two-body vector current through \(\mathcal O(p^3)\) occurs in the AV channel (axial dark-matter bilinear, vector nucleon current), with an \(i(\boldsymbol\tau_1\times\boldsymbol\tau_2)_3\) structure. This is not a VV current; the statement is limited to the paper's order and counting. It does not exclude higher-order VV two-body terms or settle bound-state continuity for the radiative process.

Martin Hoferichter, Philipp Klos, Javier Menéndez, and Achim Schwenk, “Analysis strategies for general spin-independent WIMP-nucleus scattering,” Phys. Rev. D 94, 063505 (2016). [arXiv:1605.08043](https://arxiv.org/abs/1605.08043), [DOI](https://doi.org/10.1103/PhysRevD.94.063505). Provides coherent scalar/vector response formalism and xenon structure-factor fits, including scalar two-body contributions. These are not AFH \(W_M\) columns and require normalization, isospin, oscillator-length, and nuclear-model matching.

Martin Hoferichter, Philipp Klos, Javier Menéndez, and Achim Schwenk, “Nuclear structure factors for general spin-independent WIMP-nucleus scattering,” Phys. Rev. D 99, 055031 (2019). [arXiv:1812.05617](https://arxiv.org/abs/1812.05617), [DOI](https://doi.org/10.1103/PhysRevD.99.055031). Extends the coherent-response programme to additional nuclear targets and operators.

## Xenon axial responses

Javier Menéndez, D. Gazit, and A. Schwenk, “Spin-dependent WIMP scattering off nuclei,” Phys. Rev. D 86, 103511 (2012). [arXiv:1208.1094](https://arxiv.org/abs/1208.1094), [DOI](https://doi.org/10.1103/PhysRevD.86.103511). Derives one-body and leading long-range two-body axial currents and applies them to \(^{129,131}\mathrm{Xe}\).

P. Klos et al., “Large-scale nuclear structure calculations for spin-dependent WIMP scattering with chiral effective field theory currents,” Phys. Rev. D 88, 083516 (2013). [arXiv:1304.7684](https://arxiv.org/abs/1304.7684), [DOI](https://doi.org/10.1103/PhysRevD.88.083516). Provides xenon spin-dependent structure factors and nuclear-structure calculations with chiral-EFT current effects.

B.-S. Hu et al., “Ab Initio Structure Factors for Spin-Dependent Dark Matter Direct Detection,” Phys. Rev. Lett. 128, 072502 (2022). [DOI](https://doi.org/10.1103/PhysRevLett.128.072502). Provides ab initio spin-dependent structure factors for \(^{129,131}\mathrm{Xe}\) using chiral interactions and consistent axial-current effects.

## Vector-current source matching and two-body-current hierarchy

S. Kölling, E. Epelbaum, and D. R. Phillips, “The magnetic form factor of the deuteron in chiral effective field theory,” Phys. Rev. C 86, 047001 (2012). [arXiv:1209.0837](https://arxiv.org/abs/1209.0837), [DOI](https://doi.org/10.1103/PhysRevC.86.047001). Derives isoscalar electromagnetic two-nucleon charge/current contributions at higher chiral order. Its \(\bar d_9\) and \(L_2\) coefficients belong to the electromagnetic source and its chosen NN/regulator scheme; they are not automatically LECs for the universal dark-matter baryon-number source.

Hermann Krebs, “Nuclear currents in chiral effective field theory,” Eur. Phys. J. A 56, 234 (2020). [DOI](https://doi.org/10.1140/epja/s10050-020-00230-9). Review of chiral nuclear currents, including the relation between current operators, nuclear forces, regulator choices, and unitary-transformation/off-shell conventions.

Martin Hoferichter, Javier Menéndez, and Achim Schwenk, “Coherent elastic neutrino-nucleus scattering: EFT analysis and nuclear responses,” Phys. Rev. D 102, 074018 (2020). [arXiv:2007.08529](https://arxiv.org/abs/2007.08529), [DOI](https://doi.org/10.1103/PhysRevD.102.074018). Provides xenon nuclear responses and vector-current form-factor inputs. Its Table I lists \(\kappa_s=-0.017(4)\) and \(\langle r_{E,s}^2\rangle=-0.0048(6)\,\mathrm{fm}^2\). Useful for nucleon input and response methodology; its weak-current coefficients are not the universal dark-matter source.

Gang Li, Chuan-Qiang Song, Feng-Jie Tang, and Jiang-Hao Yu, “Comprehensive effective field theory framework for coherent elastic neutrino-nucleus scattering,” Phys. Rev. D 114, 055032 (2026). [arXiv:2601.19883v2](https://arxiv.org/abs/2601.19883v2), [DOI](https://doi.org/10.1103/klg7-36x7). A method reference for matching the low-energy effective field theory to chiral operators and carrying the power counting through to nuclear responses in $\mathrm{CE}\nu\mathrm{NS}$. It does not provide the dark-matter source, the radiative kernel, or a matched xenon current for this project.

Raza Sabbir Sufian et al. (\(\chi\)QCD Collaboration), “Strange Quark Magnetic Moment of the Nucleon at the Physical Point,” Phys. Rev. Lett. 118, 042001 (2017). [DOI](https://doi.org/10.1103/PhysRevLett.118.042001). Historical physical-point lattice input, \(\kappa_s=-0.064(17)\), retained as an alternative to the later lattice-informed values.

Constantia Alexandrou et al., “Nucleon strange electromagnetic form factors from \(N_f=2+1+1\) lattice QCD,” [arXiv:2603.26591](https://arxiv.org/abs/2603.26591) (2026 preprint). Table XIV reports \(G_M^s(0)=-0.0179(20)\) and the electric z-expansion coefficient \(c_{1,0}=0.0897(86)\) for \(Q_{\rm cut}^2=0.85\,\mathrm{GeV}^2\). Treat this as a separately labelled preprint sensitivity input, not a silent replacement of the frozen baseline.

For the special universal \(u+d+s\) vector source and its VV-specific scope, see [VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md](VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md).

## Few-body chiral-EFT method cross-check

E. Filandri and M. Viviani, “Dark matter scattering off \(^{2}\mathrm H\) and \(^{4}\mathrm{He}\) nuclei within chiral effective field theory,” Phys. Rev. C 110, 034002 (2024). [DOI](https://doi.org/10.1103/PhysRevC.110.034002). A useful independent method reference, but not a heavy-xenon response database.

## Novelty boundary

In the works reviewed here, direct scalar radiation, the elastic AFH multipoles, and Dirac–Pauli matching for magnetic upscattering are established separately. None of these references supplies this project's coherent two-diagram scalar-radiation tensor contracted with the elastic nuclear-current tensor for the same interaction, including the signed interference terms. This is a boundary of the reviewed sources, not a claim of exhaustive priority. Reuse the nuclear operator and density-matrix machinery where its conventions agree, and derive and validate the radiative coefficients for the chosen source.

Related project records are [LZ #51](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/pull/51), [LZ #90](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/pull/90), and [LZ #100](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/pull/100). These are links to project history, not independent literature evidence or a change to this repository's validation status.

Before numerical reuse, record paper version, table/equation, normalization, isospin convention, nuclear Hamiltonian/model, regulator, oscillator length, and source provenance.
