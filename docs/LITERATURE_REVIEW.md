# Literature review and reusable nuclear inputs

This is a curated starting bibliography, not an exhaustive search.

## Radiative dark-matter scattering

**Project benchmark manuscript:** Rundong Fang, Ji-Heng Guo, Ming Li, Jia Liu, and Xiao-Ping Wang, “High-Energy Nuclear Recoils from 2 → 3 Dark Matter Scattering.” The supplied manuscript studies \(\chi_1N\to\chi_1N+s\), with a neutral scalar radiated from the dark-matter line before or after vector exchange. It establishes the extra \(\chi_1s\) invariant-mass continuum and the matched recoil-endpoint comparison to inelastic \(2\to2\) scattering when \(\delta=m_s\). The source PDF is not copied into this repository.

David Curtin and Yuhsin Tsai, “The Double-Dark Portal,” JHEP 11 (2014) 136. [DOI](https://doi.org/10.1007/JHEP11(2014)136). An earlier radiative dark-matter–nucleus process with a very light emitted scalar; it establishes that the \(2\to3\) topology itself is not new.

Nicole F. Bell et al., “The Migdal Effect and Photon Bremsstrahlung in effective field theories of dark matter direct detection,” Phys. Rev. D 101, 015012 (2020). [arXiv:1905.00046](https://arxiv.org/abs/1905.00046), [DOI](https://doi.org/10.1103/PhysRevD.101.015012). Useful for radiation-plus-EFT methodology, but its photon is associated with the recoiling atom/nucleus rather than dark-matter-line scalar emission.

## AFH nuclear responses

Nikhil Anand, A. Liam Fitzpatrick, and W. C. Haxton, “The WIMP-Nucleus Elastic Scattering Response,” Phys. Rev. C 89, 065501 (2014). [arXiv:1308.6288](https://arxiv.org/abs/1308.6288), [DOI](https://doi.org/10.1103/PhysRevC.89.065501).

This is the AFH nuclear-response framework used for elastic nuclear multipoles. Reuse requires the exact AFH v6.0 normalization, Fourier convention, OBDME input, and interference signs. \(W_M\) is finite-size nuclear response, not point nucleus; \(\Phi''\) is not the Lorentz-scalar nucleon source.

## Chiral EFT one- and two-body currents

Martin Hoferichter, Philipp Klos, and Achim Schwenk, “Chiral power counting of one- and two-body currents in direct detection of dark matter,” Phys. Lett. B 746, 410–416 (2015). [arXiv:1503.04811](https://arxiv.org/abs/1503.04811), [DOI](https://doi.org/10.1016/j.physletb.2015.05.041). Establishes common chiral power counting for vector, axial-vector, scalar, and pseudoscalar WIMP-nucleon interactions and derives one-/two-body currents through third order. Section 5.2 states that the only two-body vector current through \(\mathcal O(p^3)\) occurs in the AV channel (axial dark-matter bilinear, vector nucleon current), with an \(i(\boldsymbol\tau_1\times\boldsymbol\tau_2)_3\) structure. This is not a VV current; the statement is limited to the paper's order and counting.

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

Raza Sabbir Sufian et al. (\(\chi\)QCD Collaboration), “Strange Quark Magnetic Moment of the Nucleon at the Physical Point,” Phys. Rev. Lett. 118, 042001 (2017). [DOI](https://doi.org/10.1103/PhysRevLett.118.042001). Historical physical-point lattice input, \(\kappa_s=-0.064(17)\), retained as an alternative to the later lattice-informed values.

Constantia Alexandrou et al., “Nucleon strange electromagnetic form factors from \(N_f=2+1+1\) lattice QCD,” [arXiv:2603.26591](https://arxiv.org/abs/2603.26591) (2026 preprint). Table XIV reports \(G_M^s(0)=-0.0179(20)\) and the electric z-expansion coefficient \(c_{1,0}=0.0897(86)\) for \(Q_{\rm cut}^2=0.85\,\mathrm{GeV}^2\). Treat this as a separately labelled preprint sensitivity input, not a silent replacement of the frozen baseline.

For the special universal \(u+d+s\) vector source and its VV-specific scope, see [VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md](VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md).

## Few-body chiral-EFT method cross-check

E. Filandri and M. Viviani, “Dark matter scattering off \(^{2}\mathrm H\) and \(^{4}\mathrm{He}\) nuclei within chiral effective field theory,” Phys. Rev. C 110, 034002 (2024). [DOI](https://doi.org/10.1103/PhysRevC.110.034002). A useful independent method reference, but not a heavy-xenon response database.

## Novelty boundary

Radiative \(2\to3\) scattering, AFH multipoles, chiral-EFT one-/two-body currents, and xenon scalar/axial responses are established ingredients. The target contribution is the channel-specific three-body dark-matter response kernel and its consistent contraction with the corresponding elastic nuclear-current tensor, followed by reproducible xenon predictions.

Before numerical reuse, record paper version, table/equation, normalization, isospin convention, nuclear Hamiltonian/model, regulator, oscillator length, and source provenance.
