# Kinematic, Lorentz, and nuclear-response conventions

Use \(g^{\mu\nu}=\operatorname{diag}(+1,-1,-1,-1)\) unless explicitly declared otherwise.

The process is \(\chi_1(p_\chi)+A(P_A)\to\chi_1(p_\chi')+A(P_A')+s(k_s)\). Define

\[
q_A^\mu=P_A'^\mu-P_A^\mu,\quad
q_\chi^\mu=p_\chi'^\mu-p_\chi^\mu,\quad
q_N^\mu=\text{transfer through the one-body nucleon-current vertex}.
\]

Momentum conservation gives

\[
p_\chi-p_\chi'=q_A+k_s,\qquad q_\chi=-(q_A+k_s).
\]

The two internal dark-state momenta are \(r_{\rm before}=p_\chi-k_s\) and \(r_{\rm after}=p_\chi'+k_s\). The scalar momentum is not part of the AFH nuclear Fourier transfer.

For an initially stationary nucleus,

\[
Q_A^2=|\mathbf q_A|^2=2m_AE_R+E_R^2,\qquad
t_A=q_A^2=-2m_AE_R,\qquad
y_A=(Q_Ab_A/2)^2.
\]

For a free on-shell nucleon, \(q_N=k_N'-k_N\). The inherited one-body impulse production prescription assigns \(q_N=q_A\), hence \(t_N=t_A\); this does not reconstruct the microscopic off-shell four-momenta of a bound nucleon. A separate static diagnostic may use \(q_{N,\rm static}^0=0\) and \(t_{N,\rm static}=-Q_A^2\). Never mix these prescriptions.

The inherited AFH Fourier convention is \(\mathbf q_{\rm AFH}=-\mathbf q_A\). The AFH responses are finite-size nuclear matrix elements. \(W_M\) is not a point-nucleus response. Nucleon form factors \(F_i(t_N)\) and nuclear responses \(W_k(y_A)\) are different objects.

The inherited AFH v6.0 raw isoscalar zero-momentum normalization is \(A^2/4\); the physical equal-proton/equal-neutron coupling convention gives \(A^2\). Apply the factor of four exactly once. Apply any \(C_0(A)\) OBDME correction exactly once and record it.

Channel labels: VV = vector DM/vector nucleon; VA = vector DM/axial nucleon; AV = axial DM/vector nucleon; AA = axial DM/axial nucleon; SS = scalar DM/scalar nucleon. The emitted particle remains the same neutral scalar in this first stage.

The Lorentz-scalar source \(J_S\) is not AFH \(\Phi''\), which is a spin-orbit-related multipole. A leading scalar density has a density-type \(M\) operator structure, but its numerical response is not automatically vector \(W_M\).


## Universal light-flavor vector benchmark

For the special benchmark with \(C_u^V=C_d^V=C_s^V=C_V\),

\[
J_{uds}^\mu=\sum_{q=u,d,s}\bar q\gamma^\mu q=3J_B^\mu,
\qquad
J_B^\mu=\frac13\sum_{q=u,d,s}\bar q\gamma^\mu q.
\]

It is a flavor-singlet vector source. In the \(u,d\) isospin decomposition \(C_V^-=(C_u^V-C_d^V)/2\), this benchmark has \(C_V^-=0\); the strange-quark matrix element is still included explicitly in the nucleon form factors. Do not identify the SU(2) singlet-source notation \(v_\mu^{(s)}\) with a strange-quark current.

Under charge symmetry,

\[
F_i^{uds,N}=3\left(F_i^{\mathrm{EM},p}+F_i^{\mathrm{EM},n}+F_i^{s,N}\right).
\]

If the external rate convention already uses \(g_{VN}=3g_{Vq}\), the current uses \(\widehat F_i=F_i^{uds,N}/3\). Apply this factor of three exactly once. At zero transfer, \(\widehat F_1(0)=1\) and \(\widehat F_2(0)=\kappa_p+\kappa_n+\kappa_s\).

For the inherited formal one-body impulse assignment, \(q_N^\mu=q_A^\mu\) is the external four-transfer at the current vertex, not the microscopic off-shell momentum change of each bound nucleon. Therefore \(t_N=t_A=-2m_AE_R\), while \(Q_A^2=2m_AE_R+E_R^2\) remains the AFH spatial momentum squared. Use \(F_i(t_N)\), \(y_A=(Q_Ab_A/2)^2\), and the invariant \(t_A\) in their respective places.

The channel-dependent two-body-current hierarchy is recorded in [the VV universal-source audit](VV_UNIVERSAL_VECTOR_CURRENT_AUDIT.md). In particular, an AV two-body vector current is not a VV current, and an electromagnetic isoscalar LEC is not automatically a baryon-number-current LEC.
