# Nuclear-current modules

## Vector

\[
J_V^\mu=\bar u_N\left[F_1^V(t_N)\gamma^\mu+
\frac{iF_2^V(t_N)}{2m_N}\sigma^{\mu\nu}q_{N,\nu}\right]u_N.
\]

Recover the inherited VV benchmark first. Keep the Dirac/Pauli sign, \(q_N\) scheme, and form-factor argument explicit. Use vector-current continuity only when the selected Hamiltonian/current pair supports it.

### Universal \(u+d+s\) VV benchmark

For equal light-flavor vector couplings, \(J_{uds}^\mu=3J_B^\mu\). Under charge symmetry,

\[
F_i^{uds,N}=3(F_i^{\mathrm{EM},p}+F_i^{\mathrm{EM},n}+F_i^{s,N}).
\]

When the external coupling is \(g_{VN}=3g_{Vq}\), use \(\widehat F_i=F_i^{uds,N}/3\), with the division by three applied once. Keep the external one-body transfer assignment \(q_N=q_A\), \(t_N=t_A\), distinct from the spatial AFH argument \(Q_A^2\).

The form-factor sign convention must be translated explicitly when using a reference: Hoferichter et al. (2015) write the Pauli term with a minus sign in their Eq. (19). Their only two-body vector current through \(\mathcal O(p^3)\) is in the AV channel, not VV. Higher-order electromagnetic isoscalar currents require a separate match to the baryon-number source. See [the full VV source audit](../../docs/VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md).

## Axial

\[
J_A^{\mu,a}=\bar u_N\left[G_A(t_N)\gamma^\mu\gamma_5+
\frac{G_P(t_N)}{2m_N}q_N^\mu\gamma_5\right]\frac{\tau^a}{2}u_N.
\]

The induced-pseudoscalar term is tied to PCAC and modifies the longitudinal spin current. Project transverse/longitudinal components onto AFH \(\Sigma'\) and \(\Sigma''\). Do not impose \(q_\mu J_A^\mu=0\). Compare Xe-129/Xe-131 against Menéndez–Gazit–Schwenk and Klos et al.

## Scalar

Define the quark/gluon source and Wilson-coefficient convention before matching, for example

\[
J_S=\sum_{q=u,d,s}C_q^S m_q\bar q q+
C_g^S\frac{\alpha_s}{\pi}G_{\mu\nu}^aG^{a\mu\nu}.
\]

The one-body density has density-type \(M\) operator structure but is not automatically the vector \(W_M\). Include pion-exchange and trace-anomaly two-body terms at the chosen order. \(J_S\) is not AFH \(\Phi''\).
