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

Martin Hoferichter, Philipp Klos, and Achim Schwenk, “Chiral power counting of one- and two-body currents in direct detection of dark matter,” Phys. Lett. B 746, 410–416 (2015). [arXiv:1503.04811](https://arxiv.org/abs/1503.04811), [DOI](https://doi.org/10.1016/j.physletb.2015.05.041). Establishes common chiral power counting for vector, axial-vector, scalar, and pseudoscalar WIMP-nucleon interactions and derives one-/two-body currents through third order. Exchange-current content depends on both dark-side and nucleon-side channels.

Martin Hoferichter, Philipp Klos, Javier Menéndez, and Achim Schwenk, “Analysis strategies for general spin-independent WIMP-nucleus scattering,” Phys. Rev. D 94, 063505 (2016). [arXiv:1605.08043](https://arxiv.org/abs/1605.08043), [DOI](https://doi.org/10.1103/PhysRevD.94.063505). Provides coherent scalar/vector response formalism and xenon structure-factor fits, including scalar two-body contributions. These are not AFH \(W_M\) columns and require normalization, isospin, oscillator-length, and nuclear-model matching.

Martin Hoferichter, Philipp Klos, Javier Menéndez, and Achim Schwenk, “Nuclear structure factors for general spin-independent WIMP-nucleus scattering,” Phys. Rev. D 99, 055031 (2019). [arXiv:1812.05617](https://arxiv.org/abs/1812.05617), [DOI](https://doi.org/10.1103/PhysRevD.99.055031). Extends the coherent-response programme to additional nuclear targets and operators.

## Xenon axial responses

Javier Menéndez, D. Gazit, and A. Schwenk, “Spin-dependent WIMP scattering off nuclei,” Phys. Rev. D 86, 103511 (2012). [arXiv:1208.1094](https://arxiv.org/abs/1208.1094), [DOI](https://doi.org/10.1103/PhysRevD.86.103511). Derives one-body and leading long-range two-body axial currents and applies them to \(^{129,131}\mathrm{Xe}\).

P. Klos et al., “Large-scale nuclear structure calculations for spin-dependent WIMP scattering with chiral effective field theory currents,” Phys. Rev. D 88, 083516 (2013). [arXiv:1304.7684](https://arxiv.org/abs/1304.7684), [DOI](https://doi.org/10.1103/PhysRevD.88.083516). Provides xenon spin-dependent structure factors and nuclear-structure calculations with chiral-EFT current effects.

B.-S. Hu et al., “Ab Initio Structure Factors for Spin-Dependent Dark Matter Direct Detection,” Phys. Rev. Lett. 128, 072502 (2022). [DOI](https://doi.org/10.1103/PhysRevLett.128.072502). Provides ab initio spin-dependent structure factors for \(^{129,131}\mathrm{Xe}\) using chiral interactions and consistent axial-current effects.

## Few-body chiral-EFT method cross-check

E. Filandri and M. Viviani, “Dark matter scattering off \(^{2}\mathrm H\) and \(^{4}\mathrm{He}\) nuclei within chiral effective field theory,” Phys. Rev. C 110, 034002 (2024). [DOI](https://doi.org/10.1103/PhysRevC.110.034002). A useful independent method reference, but not a heavy-xenon response database.

## Novelty boundary

Radiative \(2\to3\) scattering, AFH multipoles, chiral-EFT one-/two-body currents, and xenon scalar/axial responses are established ingredients. The target contribution is the channel-specific three-body dark-matter response kernel and its consistent contraction with the corresponding elastic nuclear-current tensor, followed by reproducible xenon predictions.

Before numerical reuse, record paper version, table/equation, normalization, isospin convention, nuclear Hamiltonian/model, regulator, oscillator length, and source provenance.
