# Finite-invariant-transfer nucleon inputs for the fixed VV charge comparison

**Scope:** a reweighting of the four pinned, exact-momentum, one-body radiative VV **charge** kernels for Xe-129/Xe-131. No new radiative tensor integration, no transverse channels, and no differential cross section or event rate. The two model responses and T3 tensors from the [anchor stage, commit 1125554](https://github.com/lm846802297-bit/radiative-dark-matter-nuclear-currents/tree/112555411a63afb5e7a900f3748a03edfdbf3a10/nuclear/sn100pn/vv-charge-matching) are reused and verified by SHA-256. The old or alternative form-factor inputs do not overwrite the baseline.

## Nucleon basis and source controls

The probe couples equally to \(u,d,s\), not to the electromagnetic current. With charge symmetry \(u_p=d_n\), \(d_p=u_n\), \(s_p=s_n\), the two species have the same universal flavor sum, despite different electromagnetic form factors. In the inherited external \(g_{VN}=3g_{Vq}\) convention,
\[
 \widehat F_i^p=\widehat F_i^n=\frac{F_i^{u/p}+F_i^{d/p}+F_i^{s/p}}{3}.
\]
No independent proton/neutron parameter or isovector coupling is fitted. Equivalently, the effective charge-symmetric input is \(F_i^{p,\mathrm{em}}+F_i^{n,\mathrm{em}}+F_i^s\).

The adopted **Bishara linear** input uses [Bishara et al., Appendix A, Eqs. (A3), (A6)–(A8)](https://arxiv.org/abs/1707.06998):
\[
 \widehat F_1(t_N)=1+\frac{5.57+2.84-0.018}{3}\ {\rm GeV}^{-2}\,t_N,\qquad
 \widehat F_2(t_N)=\widehat F_2(0)=0.880-1-0.064=-0.184.
\]
This is a first-order Dirac slope and a **constant Pauli anchor**. It supplies neither a Pauli shape nor a precision estimate of missing curvature. The separately labelled **Borah isoscalar \(z\)** sensitivity uses the [Borah et al. isoscalar Sachs coefficients, Eq. (16)](https://arxiv.org/abs/2003.13640) and the same pinned strange Dirac slope and constant strange Pauli input. It has a finite electromagnetic isoscalar Pauli shape but no new strange Pauli shape. The two results are central prescriptions, **not an error band** or a covariance-propagated uncertainty. The source does not provide a validated finite-\(t\) strange Pauli shape for this replay; no newer lattice coefficients or isospin-breaking correction were assumed. The original constant anchors \((1,-0.184)\) are the comparison reference.

Use the archived external-nucleon-transfer impulse convention \(q_N=q_A\). Evaluate form factors at \(-t_N=2m_AE_R\) in functions taking a positive spacelike argument, while the HO responses and the explicit spin-orbit factor retain the exact spatial \(Q_A^2=2m_AE_R+E_R^2\). With \(D=\widehat F_1+2\widehat F_2\), the matched coefficients are \(c_M=\widehat F_1+t_ND/(8m_N^2)\), \(c_{\Phi''}=Q_A^2D/(4m_N^2)\). The signed charge pieces are \(L(c_M^2F_M,c_{\Phi''}^2F_{\Phi''},2c_Mc_{\Phi''}F_{M\Phi''})\), with the factor two counted once. The common \(L=L_{00}\) gives the temporal charge diagnostic; \(L_C\) is preserved only as a conditional nuclear-continuity comparison. No continuity is imposed on the actual bound-state current in this work.

## Results

Each pair of entries at a given recoil shares the same tensor weight, mass, speed, width, momentum transfer and nucleon input. \(K_{N,f}\) denotes the fixed charge kernel for model \(N\) and nucleon input \(f\). For each model, \(K_{N,\mathrm{Bishara}}/K_{N,\mathrm{anchor}}\) measures just the nucleon input; \(K_{\mathrm{SN},f}/K_{\mathrm{BD},f}\) measures just the nuclear-model swap at that input. The ratio of the latter contrasts is archived to show how those effects interact.

| Isotope | \(E_R\) (keV) | \(K_{\rm BD,Bishara}/K_{\rm BD,anchor}\) | \(K_{\rm SN,Bishara}/K_{\rm SN,anchor}\) | SN / B&D at anchors | SN / B&D with Bishara |
| --- | ---: | ---: | ---: | ---: | ---: |
| Xe-129 | 200 | 0.74993 | 0.74997 | 1.122690 | 1.122746 |
| Xe-129 | 267 | 0.67640 | 0.67667 | 1.180749 | 1.181229 |
| Xe-131 | 200 | 0.74640 | 0.74645 | 1.134819 | 1.134894 |
| Xe-131 | 267 | 0.67971 | 0.67948 | 1.872106 | 1.871476 |

At these four points, the adopted finite-\(t\) input lowers each model's fixed charge kernel by **25.0–32.4%**, but changes the same-input SN/B&D contrast by at most **0.041% relative**. The alternate isoscalar \(z\) input lowers the kernels by about **21–27%** and likewise leaves the same-input model contrast nearly unchanged, within **0.009% relative** of the anchor contrast. These changes are specific to the four existing kinematic points and do not establish a bound elsewhere.

For example, at Xe-131, 267 keV the Bishara \(L_{00}\)-weighted \(M,\Phi'',M\Phi''\) pieces are \((1.61135{\times}10^6,4.00836{\times}10^3,-1.52420{\times}10^5)\) in B&D and \((3.00825{\times}10^6,6.37400{\times}10^3,-2.76775{\times}10^5)\) in SN100PN. The net spin-orbit correction relative to \(M\) is **−9.21%** and **−8.99%**, respectively, not negligible at this small coherent form factor. These are charge-kernel numbers in the archived T3 normalization, not cross-section units.

## Replay and verification

Use a clean checkout of the pinned anchor-stage commit `112555411a63afb5e7a900f3748a03edfdbf3a10`:

~~~bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
sha256sum -c SHA256SUMS.txt
python -m unittest -v
python finite_t.py --repo /path/to/clean-anchor-commit --output /tmp/finite_t_replay.json
cmp native/finite_t_summary.json /tmp/finite_t_replay.json
~~~

An independently written form-factor implementation agrees with the pinned source to \(2.22\times10^{-16}\) absolute; the matched coefficients agree to \(1.74\times10^{-18}\). Explicit signed terms close to a direct positive-semidefinite \(2\times2\) Gram contraction to \(2.98\times10^{-16}\) relative. The old exact-momentum anchor kernels are reproduced bit-for-bit. The archive retains all signed terms, both weight conventions and both factorization paths at each point.
