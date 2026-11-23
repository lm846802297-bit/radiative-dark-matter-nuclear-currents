# Universal u+d+s vector source and VV nuclear-current audit

**Status:** literature and matching audit; no new numerical cross section or xenon prediction is claimed.  
**Purpose:** define the reusable one-nucleon VV baseline and the limits of importing published two-body vector currents into the radiative 2-to-3 framework.

## 1. Source and channel definition

Consider the nucleon-side quark vector source

\[
J_q^\mu=\sum_{q=u,d,s} C_q^V\,\bar q\gamma^\mu q .
\]

The original radiative benchmark has equal light-flavor vector couplings,

\[
C_u^V=C_d^V=C_s^V=C_V .
\]

Its hadronic source is therefore

\[
J_{uds}^\mu
=\bar u\gamma^\mu u+\bar d\gamma^\mu d+\bar s\gamma^\mu s
=3J_B^\mu ,
\qquad
J_B^\mu\equiv\frac{1}{3}\sum_{q=u,d,s}\bar q\gamma^\mu q .
\]

Thus this special benchmark is a flavor-singlet vector source proportional to the three-flavor baryon-number current. In the two-flavor decomposition,

\[
C_V^+=\frac{C_u^V+C_d^V}{2},\qquad
C_V^-=\frac{C_u^V-C_d^V}{2},
\]

it has \(C_V^-=0\). The strange-quark matrix element remains a separate contribution to the nucleon form factors; it must not be confused with the SU(2) singlet-source notation.

The channel is **VV** in this repository's convention: the first V is the dark-matter transition bilinear and the second V is the nucleon-side vector current. An AV two-body current is not a VV current.

A useful notation warning comes from Filandri and Viviani (2024). Their external source \(v_\mu^{(s)}\) is the SU(2) isoscalar singlet source coupled to the \(u,d\) doublet in their Eq. (1); the superscript \(s\) means *singlet*, not strange quark. Their SU(2) treatment is a method reference for the singlet source, while the strange-quark form factor in the present \(u+d+s\) source still requires explicit matching.

## 2. One-nucleon Dirac-Pauli matching

Write the flavor-current matrix element as

\[
\langle N'|J_{uds}^\mu|N\rangle
=
\bar u_N(p_N')
\left[
F_1^{uds,N}(t_N)\gamma^\mu
+\frac{iF_2^{uds,N}(t_N)}{2m_N}
\sigma^{\mu\nu}q_{N,\nu}
\right]
u_N(p_N).
\]

With charge symmetry, the flavor decomposition used by Hoferichter, Klos, and Schwenk (2015), Eq. (21), gives

\[
F_i^{uds,N}(t)
=
3\left[
F_i^{\mathrm{EM},p}(t)
+F_i^{\mathrm{EM},n}(t)
+F_i^{s,N}(t)
\right],
\qquad i=1,2.
\]

If the external rate prefactor is already written using \(g_{VN}=3g_{Vq}\), define

\[
\widehat F_i(t_N)\equiv \frac{F_i^{uds,N}(t_N)}{3}
=
F_i^{\mathrm{EM},p}(t_N)
+F_i^{\mathrm{EM},n}(t_N)
+F_i^{s,N}(t_N).
\]

The factor of three is then included in the coupling convention exactly once. At zero transfer,

\[
\widehat F_1(0)=1,\qquad
\widehat F_2(0)=\kappa_p+\kappa_n+\kappa_s .
\]

The Pauli-term sign must be checked when moving formulae between references: Hoferichter et al. write the nucleon matrix element with a minus sign in front of \(i\sigma^{\mu\nu}q_\nu F_2/(2m_N)\) in their Eq. (19), while other conventions use a plus sign together with a different momentum/sign definition. Do not copy a sign without translating the convention.

### 2.1 Strange form-factor input sets

The following values are useful as distinct input choices, not as interchangeable measurements:

| Source | Strange magnetic moment | Strange electric radius / shape input | Use |
|---|---:|---:|---|
| Sufian et al. (2017) | \(\kappa_s=-0.064(17)\) | \(\langle r_{E,s}^2\rangle=-0.0043(21)\,\mathrm{fm}^2\) | Historical lattice baseline |
| Hoferichter, Menéndez, Schwenk (2020), Table I | \(\kappa_s=-0.017(4)\) | \(\langle r_{E,s}^2\rangle=-0.0048(6)\,\mathrm{fm}^2\) | Updated CEvNS input set |
| Alexandrou et al. (2026), arXiv:2603.26591, Table XIV | \(G_M^s(0)=-0.0179(20)\) | \(G_E^s\) z-expansion coefficient \(c_{1,0}=0.0897(86)\), \(Q_{\mathrm{cut}}^2=0.85\,\mathrm{GeV}^2\) | Recent preprint sensitivity input |

For the 2026 z-expansion, the paper chooses \(t_0=0\) and \(t_{\mathrm{cut}}=4m_K^2\), so its Eq. (28) gives

\[
\langle r_{E,s}^2\rangle=-\frac{3c_{1,0}}{2t_{\mathrm{cut}}}.
\]

The 2026 result is a preprint and should remain a separately labelled sensitivity input until its status and covariance are reviewed. Do not silently replace the historical benchmark. When propagating its uncertainty, retain the reported fit covariance rather than treating correlated z-expansion coefficients as independent.

## 3. Transfer variables in the radiative process

For an initially stationary nucleus,

\[
Q_A^2\equiv |\mathbf q_A|^2=2m_AE_R+E_R^2,\qquad
t_A=q_A^2=-2m_AE_R.
\]

The inherited formal one-body impulse prescription assigns the external transfer to the one-body nucleon-current vertex,

\[
q_N^\mu\equiv q_A^\mu,\qquad t_N=t_A.
\]

For a bound nucleus this is an **external vertex-transfer assignment**, not a reconstruction of each bound nucleon's microscopic off-shell initial and final four-momenta. Only for a free on-shell nucleon can \(q_N^\mu\) be identified literally with \(k_N'{}^\mu-k_N^\mu\).

Keep the arguments separate:

- nucleon form factors: \(F_i(t_N)\), with \(-t_N=2m_AE_R\);
- AFH momentum variable: \(y_A=(Q_Ab_A/2)^2\), with \(Q_A^2=2m_AE_R+E_R^2\);
- mediator propagator: the invariant \(t_A\), not \(Q_A^2\).

A static diagnostic with \(q_N^0=0\) and \(t_N=-Q_A^2\) is a different prescription and must not be mixed with the formal external-transfer production convention.

## 4. One-body AFH projection

For the inherited leading one-body vector current, the nonrelativistic charge and spatial pieces project onto the familiar AFH structures:

- charge: \(M\) and \(\Phi''\);
- spatial current: \(\Delta\) and \(\Sigma'\).

Retain the signed \(M\)-\(\Phi''\) and \(\Delta\)-\(\Sigma'\) interference terms. The three-body radiation changes the dark-matter-side response coefficient and phase-space weighting; it does not by itself change the elastic nuclear state or the AFH operator definitions.

The AFH responses are finite-size nuclear matrix elements. They are not nucleon form factors, and \(W_M\) is not a point-nucleus response. The exact AFH v6.0 normalization, Fourier convention, OBDME normalization, and response interference signs must be pinned before importing numerical tables.

## 5. What the two-body-current literature does and does not establish

### 5.1 Direct dark-matter chiral counting

Hoferichter, Klos, and Schwenk (2015), especially Table 1 and Sec. 5.2, derive WIMP/nucleon-channel-dependent currents through \(\mathcal O(p^3)\). They state that the only two-body vector current at that order appears in the **AV** channel, where A is the dark-matter axial bilinear and V is the nucleon vector current. Its isospin structure contains

\[
i(\boldsymbol\tau_1\times\boldsymbol\tau_2)_3
\]

and its coefficient is the AV Wilson coefficient. It is not an operator that can be inserted into the VV benchmark by changing a label.

This is an order- and convention-specific result. It does not prove that every higher-order isoscalar VV two-body current vanishes.

### 5.2 Singlet-source and few-body methodology

Filandri and Viviani (2024) formulate external SU(2) vector sources as an isovector source \(v_\mu\) and an isoscalar singlet source \(v_\mu^{(s)}\), and study WIMP scattering from deuteron and \(^{4}\mathrm{He}\) through N2LO. Their Eqs. (1), (14), and (15) provide a useful source-level organization of isoscalar and isovector vector couplings. Their few-body matrix elements and nuclear interactions are not a heavy-xenon response database; their order labels should not be identified with the \(\mathcal O(p^3)\) counting of Hoferichter et al. without an explicit translation.

### 5.3 Electromagnetic isoscalar currents are not automatically baryon currents

Kölling, Epelbaum, and Phillips (2012) derive an isoscalar electromagnetic two-body magnetic current at \(\mathcal O(eP^4)\), involving \(\bar d_9\) and the short-range LEC \(L_2\). In their calculation \(L_2\) is constrained by two-nucleon electromagnetic data, and the quoted LECs depend on the regulator and NN interaction.

This is a valuable operator-level reference, not a ready-made coefficient set for the universal \(u+d+s\) baryon-number source. A separate external-source matching is required before using these operators or LECs for the dark-matter current.

Krebs (2020) reviews the construction of chiral nuclear currents and emphasizes consistency with the nuclear interaction, regulator, and unitary-transformation/off-shell convention. Any higher-order current calculation must pin these choices together with the NN/3N Hamiltonian and fitted LECs.

## 6. Bound-state continuity is a separate validation

A free on-shell nucleon Ward identity does not establish conservation of a current after it is embedded in a bound-nucleus wave function. The current and the nuclear Hamiltonian must form a consistent pair; exchange-current terms may be required by the chosen scheme.

The original LZ project's P3-C audit of its positive-energy free-spinor one-body impulse embedding reported a maximum normalized external nuclear-current divergence of \(1.040\times10^{-3}\) after bound-state folding. This is a model-specific diagnostic, not a universal estimate for xenon or for all impulse currents. It is sufficient to prevent an unsupported claim that the truncated bound current is exactly conserved. See the [pinned P3-C audit record](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/blob/d02e5e009f32a1db1b2520a789148bbd0160ccfd/research/three_body_full_current_nucleon_match/results/2026-11-23/p3c_bound_current_audit/AUDIT_REPORT.md).

## 7. Reuse boundary for this repository

1. Recover the VV benchmark with exact source and data pins before adding other channels.
2. Use the universal-source matching and finite-\(t_N\) conventions above, retaining the external coupling normalization exactly once.
3. Treat the HKS 2015 AV two-body current as an AV result, not a VV contribution.
4. Use Filandri-Viviani as a singlet-source/few-body method reference and Krebs/Kölling as higher-order current references; derive the baryon-source matching before importing electromagnetic LECs.
5. Keep the strange form-factor input sets separate and compare them at the amplitude/kernel level.
6. Do not import the LZ-specific P3-C impulse code, xenon OBDMEs, or numerical rates as general nuclear-current results. Any later numerical import requires exact source/input hashes, a compatible nuclear model, and an independent reproduction.

This note is a theory/literature audit. It does not claim a completed higher-order VV two-body current, a conserved full bound-nucleus current, or a new xenon prediction.
