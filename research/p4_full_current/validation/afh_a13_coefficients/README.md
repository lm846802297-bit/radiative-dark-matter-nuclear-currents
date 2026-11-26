# Independent AFH A13 coefficient verification

## Result

The exact script verify_afh_three_body_coefficients.py was executed locally in the ChatGPT Linux runtime on 2026-11-25. It evaluates the leading-NR one-body Dirac-Pauli current projections in AFH Appendix A, Eq. (A13), without importing the Stage C2 coefficient implementation.

Status: PASS_INDEPENDENT_A13_COEFFICIENT_IDENTITY.

## Conventions

- Metric: (+---).
- Nuclear transfer: q_N = +Q_A zhat.
- AFH Fourier momentum: q_AFH = -Q_A zhat.
- External coupling: g_VN = 3 g_Vq.
- Effective form factors in the particle-side coefficients: f_i = F_i^(uds)/3.

## Verified identities

$$
R_M^{(3)}=|c_M|^2L_{00},\quad
R_{\Phi''}^{(3)}=|c_\Phi|^2L_{00},\quad
R_{M\Phi''}^{(3)}=2\operatorname{Re}(c_Mc_\Phi^*)L_{00}.
$$

$$
R_\Delta^{(3)}=\frac{f_1^2Q_A^2}{2m_N^2}(L_{11}+L_{22}),
$$

$$
R_{\Sigma'}^{(3)}=\frac{(f_1+f_2)^2Q_A^2}{8m_N^2}(L_{11}+L_{22}),
$$

$$
R_{\Sigma'\Delta}^{(3)}=-\frac{f_1(f_1+f_2)Q_A^2}{2m_N^2}(L_{11}+L_{22}).
$$

At this leading NR order, the two projections are zero:

$$
R_{\Sigma''}^{(3)}=R_{\widetilde\Phi'}^{(3)}=0.
$$

## Numerical checks

- 20,000 random complex current vectors.
- 2,000 generic complex Hermitian positive-semidefinite 4x4 Gram tensors.
- Explicit charge-sector reconstruction and the factor of two in the M-Phi'' interference.
- Explicit signed Sigma'-Delta interference and its Cauchy identity.
- Explicit zero projections for Sigma'' and Phi-tilde-prime.

Maximum relative residuals are in validation_report.json. The largest is 1.1879386363489175e-14 for charge-sector reconstruction. The general complex-Gram reduction residual is 8.265594928561099e-16.

## Provenance

- Source definition: AFH v6.0, v6/dmformfactor-V6.m, SHA-256 dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372.
- The test is algebraic/numerical and uses a fixed PRNG seed; it does not read Stage A, Stage B, or C1 production data.
- It does not validate the complete bound-nucleon current contraction, the additional longitudinal multipoles, C2 production, or a detector rate.

## Reproduction

~~~bash
python verify_afh_three_body_coefficients.py > validation_report.json
sha256sum verify_afh_three_body_coefficients.py validation_report.json README.md provenance.json
~~~
