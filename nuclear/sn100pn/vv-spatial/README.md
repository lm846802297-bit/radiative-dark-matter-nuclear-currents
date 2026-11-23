# Signed SN100PN transverse responses in the fixed VV xenon comparison

**Scope:** four previously validated radiative fixed points, leading/static one-body transverse \(\Delta,\Sigma',\Sigma'\Delta\) only. This is not a full covariant nuclear current, proof of bound-state continuity, a spectrum, differential cross section, or a chiral-EFT prediction. The signed B&D native archive and the preceding charge-stage results are not modified.

## Calculation and normalization

For each isotope and exact inherited \(Q_A=|\mathbf q_A|\), directly recouple the [published SN100PN OBDMEs](https://github.com/Raghda-AK/DM-Form-Factor/tree/1395f0c161905341fc89a2978148b1ce7c9ee1dd) and the [AFH v6.0 source OBDMEs](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/blob/aac9b33d66f6a7fe1ca5d5d8cc648f5afe52a81b/v6/dmformfactor-V6.m) with the AFH DeltaJ and SigmaPJ operators. A hash-pinned [independent 3j/6j/9j + harmonic-oscillator recoupling helper](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/blob/ebe6e8cbd2e67d79e2a34305651d6b53ac565795/research/records/3body-hdm/xe131-dense-spatial-source-recoupling-2026-11-23/check_spatial.py) supplies the single-particle matrix elements. The LZ calculation itself is not rerun, and no LZ response or benchmark result is imported. The orbital one-body operator has abnormal parity: \(J=1\) contributes for Xe-129 and \(J=1,3\) for Xe-131; the archived \(J=2\) entries are zero by parity, not dropped by a small-value cutoff.

Use the paper's reduced-density convention with isoscalar and isovector factors \(\sqrt{2/(2T+1)}\) and \(\sqrt{6T/[(T+1)(2T+1)]}\) for \(T=(A-108)/2\). The individual SN100PN \(p,n\) amplitudes are \(\sqrt{4\pi/(2J_A+1)}(a_0\mp a_1)/2\). The physical universal total is their sum. The B&D raw isoscalar \(a_0\) is converted by \(\sqrt{4\pi/(2J_A+1)}\) once so that its squared amplitude is consistent with the \(16\pi/(2J_A+1)\) factor on the raw AFH response. No \(C_0\) or other independent rescaling is applied, as in the radiative charge archive. The mixed response is \(\sum_J A_{\Sigma',J} A_{\Delta,J}\); in species-separated supplement entries the first species labels \(\Sigma'\) and the second labels \(\Delta\).

Use the same archived \(L_T=L_{11}+L_{22}\), exact kinematics and nucleon inputs as the previous [finite-\(t_N\) charge comparison](https://github.com/lm846802297-bit/radiative-dark-matter-nuclear-currents/tree/7f82851b777e21a8618a8bd5767448dcc58f4a70/nuclear/sn100pn/vv-finite-t). The [AFH Eq. (107)](https://arxiv.org/abs/1308.6288) radiative VV one-body transverse projection closes as a sum of squares:
\[
 K_T=\frac{Q_A^2 L_T}{2m_N^2}\sum_J
 \left(\widehat F_1 A_{\Delta,J}
       -\frac{\widehat F_1+\widehat F_2}{2}A_{\Sigma',J}\right)^2.
\]
This fixes one signed mixed term; there is no extra factor of two on the nuclear response itself. The \(t_N=t_A=-2m_AE_R\) argument for the nucleon form factors is inherited and is not set to \(-Q_A^2\). Any combined charge-plus-spatial comparison below uses the separately labeled continuity-conditional \(L_C\) charge kernel; it is not a claim that this approximate nucleus obeys the exact bound-current Ward identity.

## Results

| Isotope | \(E_R\) (keV) | adopted Bishara \(K_T/K_{\rm ch}\), B&D | adopted Bishara \(K_T/K_{\rm ch}\), SN100PN |
| --- | ---: | ---: | ---: |
| Xe-129 | 200 | \(1.54\times10^{-13}\) | \(1.80\times10^{-12}\) |
| Xe-129 | 267 | \(1.10\times10^{-11}\) | \(1.36\times10^{-10}\) |
| Xe-131 | 200 | \(6.23\times10^{-12}\) | \(8.73\times10^{-12}\) |
| Xe-131 | 267 | \(6.84\times10^{-9}\) | \(5.81\times10^{-9}\) |

The signed mixed kernel term is destructive for Xe-129 and constructive for Xe-131 at these four points. Spatial responses can vary substantially between interactions: for Xe-129 at 267 keV, the SN100PN \(\Sigma'\) response is about 5.8 times B&D and the net spatial kernel is about 14.6 times B&D. These ratios concern a very small transverse contribution. Across the frozen anchors, the adopted Bishara prescription, and the alternative isoscalar \(z\) sensitivity, the addition of this transverse sector changes the SN100PN/B&D matched charge contrast by less than \(1.09\times10^{-9}\) relative. This finding is only for these four radiative points and cannot be extrapolated at a diffraction zero or to other kinematics.

### Numerical checks

- Signs of all 20 spatial species entries at \(Q_A=0\) agree with the [published rounded supplement](https://arxiv.org/abs/2311.15764); worst relative difference \(0.351\%\). The published rounded coefficients do not enter the finite-momentum calculation.
- Direct native B&D signed anchors at zero and at both bracketing native grid nodes agree within \(3.09\times10^{-16}\) absolute. The four target momenta are evaluated directly, not by interpolation.
- Independent tensor and complex-amplitude P3 projectors agree with the analytic coefficient vector within \(7.78\times10^{-15}\) relative; the three signed terms close to the rank sum of squares; each response obeys the corresponding Cauchy bound.
- Repeating the SN100PN direct recoupling at 50 and 80 decimal working digits changes signed amplitudes by at most \(4.48\times10^{-52}\) absolute. This checks arithmetic, not the physical accuracy or available precision of the OBDMEs.

## Replay

Use a clean checkout at 7f82851b777e21a8618a8bd5767448dcc58f4a70. Save the linked helper and AFH package, both published Xe OBDME files, and the [pinned arXiv source tarball](https://arxiv.org/src/2311.15764v1) locally. The replay validates all consumed inputs by SHA-256, including the individual TeX member of the source archive. Third-party OBDME files are not bundled.

~~~bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
sha256sum -c SHA256SUMS.txt
python -m unittest -v
python compute.py --repo /path/to/pinned-clean-repo --afh /path/to/dmformfactor-V6.m \
  --helper /path/to/check_spatial.py --sn-inputs /path/to/obdmes \
  --paper /path/to/2311.15764v1 --out /tmp/vv-spatial-replay.json
cmp native/summary.json /tmp/vv-spatial-replay.json
~~~

The JSON retains per-rank signed species amplitudes, each of the three transverse terms, and exact fixed benchmark comparisons for all existing nucleon prescriptions. No Woods–Saxon replacement, rate normalization, additional nuclear current or experimental projection is implied.
