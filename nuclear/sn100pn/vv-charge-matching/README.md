# Fixed radiative VV charge comparison: SN100PN vs AFH/B&D

This is a fixed-speed **charge-sector** calculation for \(\chi_1+A\to\chi_1+A+s\), not elastic dark-matter matching and not \(d\sigma/dE_R\). The comparison changes only the one-body \(M,\Phi'',M\Phi''\) nuclear responses. The \(J=0\) Xe-129 and \(J=0+2\) Xe-131 responses are complete for these three operators. Nuclear states, form factors or current terms outside these channels are not inferred.

## Matching and independent verification

Transfer \(q_A^\mu=(E_R,0,0,Q_A)\) is into the nucleus. For elastic nuclear recoil, \(Q_A^2=2m_AE_R+E_R^2\) goes into the HO responses and the spin-orbit coefficient. The impulse convention \(q_N=q_A\) gives \(t_N=t_A=-2m_AE_R\), the invariant for the scalar recoil correction. These are held distinct; the AFH momentum-direction convention in the historical spatial interface is \(q_{\rm AFH}=-q_A\).

For the normalized positive-energy Dirac spinors \(u/\sqrt{2E_N}\), the independently expanded temporal Dirac–Pauli current through second order in three-momenta is
\[
 J^0 =
 \left[F_1-\frac{Q_A^2}{8m_N^2}(F_1+2F_2)\right]\mathbf 1
 +\frac{i(F_1+2F_2)}{4m_N^2}
       \boldsymbol{\sigma}\cdot(\mathbf q_A\times\mathbf k)
 +O(p^4/m_N^4),\quad
 \mathbf k=(\mathbf p_{N,f}+\mathbf p_{N,i})/2 .
\]
Using the source-prescribed \(t_N\) instead of \(-Q_A^2\) in the scalar factor is the invariant impulse continuation at this order, not an exact bound-nucleon energy assignment. When the common momentum is adiabatically scaled down, direct four-component spinor evaluation verifies the quartic residual and rules out the reversed spin-orbit sign. In a plane wave \(e^{i\mathbf q_A\cdot\mathbf r}\), the gradient in [AFH Eq. (39)](https://arxiv.org/abs/1308.6288) converts the Pauli spin term into \(+Q_A^2(F_1+2F_2)\Phi''/(4m_N^2)\), fixing the relative sign without an elastic projectile coefficient.

The matched coefficients agree with the immutable [VV/P3 implementation](https://github.com/lm846802297-bit/radiative-dark-matter-nuclear-currents/tree/e045901c684419ba70af9d39aacc434dbbbfbb33/benchmarks/vv_recovery/source_snapshot) at every point:
\[
 c_M=F_1+\frac{t_N(F_1+2F_2)}{8m_N^2},\qquad
 c_{\Phi''}=\frac{Q_A^2(F_1+2F_2)}{4m_N^2},\qquad
 K_{\rm ch}/L=c_M^2F_M+c_{\Phi''}^2F_{\Phi''}+2c_Mc_{\Phi''}F_{M\Phi''}.
\]
Here \(F_1=1,\ F_2=-0.184,\ m_N=0.938918\) GeV are the frozen \(q=0\) effective anchors in the external \(g_{VN}=3g_{Vq}\) basis. They are not evaluated with a finite-\(t\) fit. The raw native AFH \(W^{00}\) is converted once using \(F_{\rm eqpn}=16\pi W^{00}/(2J_i+1)\). No separate \(C_0\) or particle-count rescaling is applied to either model.

The spin-resolved **two-diagram radiative** tensor is regenerated with the source-pinned spinor implementation, including their interference. \(L=L_{00}\) gives the temporal charge-only projection. The additional \(L_C=L_{00}-(E_R/Q_A)\operatorname{Re}(L_{03}+L_{30})+(E_R/Q_A)^2L_{33}\) column is a separately labelled diagnostic **conditional on nuclear vector-current continuity**. Actual nuclear Ward continuity of this truncated impulse current has not been assumed or tested by this calculation. Both \(L\) choices cancel from the same-point two-model ratio; they are archived for absolute charge-kernel bookkeeping.

## Results

\(m_\chi=1000\) GeV, \(m_s=300\) keV, laboratory speed \(776\) km/s, constant pole width with \(g_2=1\), and isotope-specific \(m_A=A(0.93149410242)\) GeV. The two inherited benchmarks are green (\(\delta m=450\) keV, \(E_R=200\) keV) and near-pole (\(\delta m=362\) keV, \(E_R=267\) keV). Phase-space, couplings and width convention are identical for the two nuclear models **within each** isotope and recoil. The ratios shown use exact \(Q_A\) source-level recoupling for both models, not extrapolation or rounded supplement polynomials.

| Isotope | \(E_R\) (keV) | \(Q_A\) (GeV) | SN100PN / AFH-B&D, matched charge kernel |
| --- | ---: | ---: | ---: |
| Xe-129 | 200 | 0.219238 | 1.12269 |
| Xe-129 | 267 | 0.253312 | 1.18075 |
| Xe-131 | 200 | 0.220931 | 1.13482 |
| Xe-131 | 267 | 0.255268 | 1.87211 |

The two-grid numerical replay changes \(L_{00}\) and \(L_C\) by at most \(1.11\times10^{-10}\) relative, and its independent direct Dirac-trace check agrees to \(4.14\times10^{-10}\) relative.

The last ratio is large because the AFH-B&D coherent form factor is small at that recoil. Per common radiative weight \(L\), the kernels are 0.00494675 (AFH-B&D) and 0.00926085 (SN100PN), not enhancements to an unweighted total rate. The signed spin-orbit plus interference correction is about −10% of \(M\)-only in this row, for either model. The maximum ratio difference using common PSD-preserving linear interpolation on the 0.25-MeV response grid is 0.1521% relative. Both adjacent native nodes agree with the independent source-level recoupling before interpolation is judged.

No target mediator factor, transverse or additional longitudinal multipoles, flux normalization, abundances, velocity distribution, detector response, bound-current correction or chiral two-body current is included. This response comparison is useful for deciding where a subsequent fully matched calculation needs additional nuclear structure checks; it does not alter the archived AFH baseline or the production status.

## Reproduce

Requires a clean checkout of commit `e045901c684419ba70af9d39aacc434dbbbfbb33` of this project and the Python package pins in this directory. The numerical replay runs the source archives locally; the first exact-\(Q_A\) audit additionally downloads three hash-pinned inputs. Nothing upstream is vendored.

~~~bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python -m unittest -v
python fixed_charge.py --repo /path/to/pinned-clean-repo --output /tmp/vv-charge-replay
python check_interpolation.py --repo /path/to/pinned-clean-repo --inputs /tmp/vv-charge-inputs --json /tmp/vv-charge-replay/off_grid.json --charge-summary /tmp/vv-charge-replay/summary.json
~~~

`native/summary.json` gives each matched term on the common native response grid; `native/radiative_tensors.json.xz` archives the numerical radiative integrals at both quadrature settings. `native/off_grid.json` gives exact-\(Q_A\) response/kernel checks. The preserved vacuum-nucleon matches and direct two-diagram trace check do not validate a complete off-shell many-body current.


### Verify the archived tensor bytes

The lossless XZ archive is tracked instead of the approximately 1 MiB raw JSON. From this directory, verify both representations:

~~~bash
sha256sum -c SHA256SUMS.txt
expected=$(awk '{print $1}' RAW_TENSOR_SHA256SUMS.txt)
actual=$(xz -dc native/radiative_tensors.json.xz | sha256sum | awk '{print $1}')
test "$expected" = "$actual"
~~~

The uncompressed SHA-256 also remains in the run manifest under the uncompressed JSON path.
