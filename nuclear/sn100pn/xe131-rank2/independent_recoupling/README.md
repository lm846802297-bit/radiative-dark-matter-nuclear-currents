# Xe-131 independent recoupling and native AFH closure

This reference evaluates **signed one-nucleon reduced matrix elements** for charge (M_J) and spin-orbit (Phi''_J) by Wigner–Eckart recoupling (3j, 6j, 9j; analytic HO 1F1 radial integrals), independently of the main spherical/angular-projection integrator. It follows the [AFH v6.0 source, MJ and PhiPPJ](https://github.com/lm846802297-bit/LZ-3-body-for-hdm/blob/aac9b33d66f6a7fe1ca5d5d8cc648f5afe52a81b/v6/dmformfactor-V6.m), SHA-256 `dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372`, the exact version named in the frozen native table's manifest. AFH definitions and nuclear embedding: [Anand, Fitzpatrick, Haxton, Eq. (39)](https://arxiv.org/abs/1308.6288). SN100PN densities are the [Abdel Khaleq et al. data](https://arxiv.org/abs/2311.15764) at their pinned public release, retaining all published digits. Third-party source and raw densities are not redistributed here.

The recoupling is in the positive-at-origin HO radial convention and tests the signed static diagonal (Phi''_{00}=-(sigmacdot L)M_{00}/3), the quadrupole zero-charge limit, and the off-diagonal Hermiticity phase. For SN100PN, this uses exactly the original 14 rows in each of the two rank-2 isospin blocks, including both orbital orders; no symmetrization or rescaling. Physical p and n amplitudes are ((S-V)/2,(S+V)/2), where (S=sqrt{8pi/d}sumPsi^{J;0}R), (V=sqrt{24pi T/[d(T+1)]}sumPsi^{J;1}R), (R) is the angular-momentum reduced single-particle operator, (d=(2J_i+1)(2T+1)), and (T=11.5). Equal-p/n responses are products of *one-body nuclear amplitudes*, not two-body nuclear currents. Sum products, not amplitudes, over nuclear multipole ranks.

A second use of the same recoupling is a strict non-SN100PN convention control: reconstruct (J=0,2) directly from the original Baldridge–Dalton density entries embedded in AFH v6.0; their summed physical responses reproduce the frozen native Wolfram (M,Phi'',MPhi'') at all 2,401 nodes with absolute errors (1.64	imes10^{-11},5.69	imes10^{-13},2.28	imes10^{-12}), respectively. The native Xe-131 raw (W^{00}) receives (16pi/(2J_i+1)=4pi) once; no fit to the rounded paper responses. This certifies the numerical convention under the two published density models, not the physics accuracy of either and not any radiative differential rate.

~~~~sh
OPENBLAS_NUM_THREADS=1 python -m unittest -v
OPENBLAS_NUM_THREADS=1 python verify_reference.py \
  --inputs /tmp/xe131-recoupling-inputs --json /tmp/xe131-recoupling.json \
  --signed-csv /tmp/xe131-recoupling-signed.csv \
  --candidate-csv /tmp/xe131_rank2.csv
cmp /tmp/xe131-recoupling.json recoupling_reference.json
~~~~

The candidate must be the decompressed official output of the adjacent rank-2 audit. All four signed amplitudes agree at all 2,401 nodes with maximum absolute difference \(1.20\times10^{-15}\), including the limits at zero momentum.

Requires Python, NumPy, SciPy and SymPy; the replay checks all input SHA-256. The published JSON gives representative signed amplitudes and like-for-like rank-2 SN100PN and source-verified AFH/B&D response landmarks in order (M,Phi'',MPhi''). No independent Xe-131 many-body wavefunction calculation is claimed.
