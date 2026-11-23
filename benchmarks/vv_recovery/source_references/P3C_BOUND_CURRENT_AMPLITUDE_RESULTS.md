# P3-C Direct Covariant Bound-Current Amplitude Comparison

**Run date:** 2026-11-22
**Scope status:** COMPLETED_MODEL_SPECIFIC_DIRECT_COVARIANT_VS_AFH_COMPARISON
**Formal production flags:** native_afh_validated=false; full_current_complete=false

## Conclusion

The selected-point P3-C comparison is complete for one explicit impulse embedding of the AFH v6.0 Xe-129 and Xe-131 OBDMEs into positive-energy free Dirac spinors. The leading nonrelativistic (NR) reduction reconstructs the native Wolfram AFH M, PhiPP, Delta, and SigmaPrime amplitudes to a maximum absolute residual of \(3.575\times10^{-14}\). The direct Dirac-Pauli calculation is numerically converged through order-22 tensor-product Gauss-Hermite quadrature. It gives nonzero, model-specific covariant corrections to the NR amplitudes and to the rank-summed T3 kernels.

This closes the requested P3-C comparison for the stated model. It does not establish a unique covariant bound-nucleus current: the bound OBDMEs do not fix an off-shell energy-sharing prescription, and the calculation omits exchange currents and relativistic bound-state wavefunctions. The separate P4 production flags therefore remain false.

## Definitions and model

For isotope mass \(m_A=A(0.93149410242\ \mathrm{GeV})\),
\[
t_A=-2m_AE_R,\qquad Q_A^2=2m_AE_R+E_R^2,\qquad q_A^\mu=(E_R,0,0,Q_A).
\]
The AFH spatial form factors are evaluated at \(Q_A=\sqrt{Q_A^2}\); the three-body vector propagator continues to use \((m_V^2-t_A)^{-2}\). The impulse current is
\[
J^\mu_N(p_f,p_i)=\bar u_N(p_f)
\left[F_1\gamma^\mu+\frac{iF_2}{2m_N}\sigma^{\mu\nu}q_{A,\nu}\right]u_N(p_i),
\quad p_i=\mathbf k-\frac{\mathbf q_A}{2},\quad
p_f=\mathbf k+\frac{\mathbf q_A}{2},
\]
with the reported external-transfer prescription \(q_A^0=E_R\), \(F_1=1\), \(F_2=-0.184\), \(m_N=0.938918\ \mathrm{GeV}\), metric \((+---)\), and spinors normalized as \(u/\sqrt{2E_N}\). AFH harmonic-oscillator momentum wavefunctions use the AFH v6.0 isotope-specific oscillator length and OBDME tables. The direct calculation integrates the spin-angular HO amplitudes and the covariant current over \(\mathbf k\), then performs the same Wigner-Eckart and magnetic-substate reductions as the AFH native multipoles. All complex per-rank amplitudes are compared before squaring.

The nucleon current is an explicit positive-energy free-spinor impulse embedding, not a unique off-shell nuclear current. No exchange current is added. The alternative internal-transfer diagnostic uses \(q_N^0=E_N(p_f)-E_N(p_i)\) at each momentum node while retaining the same external AFH momentum transfer.

The native AFH export was run with Wolfram 14.2.0 for Xe-129/Xe-131 at \(E_R=200,267\ \mathrm{keV}\). It contains 40 native multipole rows and 128 OBDME-term rows. The native raw package convention is retained in the amplitude comparison. Isospin handling distinguishes the density-matrix rank \(T=0\) from the target isospin \(T_{\rm Iso}=|A/2-Z|\); the AFH target factor \(\sqrt{2/(2T_{\rm Iso}+1)}\) is applied. The inherited \(C_0\) calibration is applied only in the physical kernel comparison: Xe-129 \(C_0=1.0001103230777202\), Xe-131 \(C_0=1.0000825298069376\). The AFH package's response factor \(1/4\) cancels the physical equal-\(p=n\) factor four; raw amplitudes are checked before that bookkeeping.

## Native-input and amplitude checks

| Check | Result |
|---|---:|
| Kinematic identity, maximum absolute residual in GeV² | \(6.94\times10^{-18}\) |
| Native export: rows / OBDME terms | 40 / 128 |
| Native multipole reconstruction by summing OBDME terms | \(4.44\times10^{-16}\) maximum absolute difference |
| Native sigma-prime term reconstruction | 0 |
| Native maximum imaginary residual | 0 |
| Direct multipole-fit residual, order 22 | \(1.1873\times10^{-14}\) |
| Leading-NR amplitude vs native AFH maximum absolute residual | \(3.5749\times10^{-14}\) |
| Direct order-22 maximum imaginary amplitude | \(6.88\times10^{-16}\) |

The channel map is M (charge identity), PhiPP (spin-orbit charge), Delta (longitudinal vector identity), and SigmaPrime (transverse vector spin). For the direct current, spatial \(L=J-1,J,J+1\) components are retained in the Cartesian current and T3 contraction. The AFH Delta and SigmaPrime response convention selects its standard multipoles; it is not used to silently drop the additional direct-current \(L\ne J\) pieces.

## Convergence and exact-Dirac vs leading-NR amplitudes

The direct result has 12 isotope/recoil/rank rows: Xe-129 has charge ranks \(J=0,1\) and spatial \(J=1\); Xe-131 has charge ranks \(J=0,1,2,3\) and spatial \(J=1,2,3\), at both recoil energies.

| Gauss-Hermite order per axis | Largest absolute amplitude change from previous order | Largest relative change among amplitudes above \(10^{-8}\) |
|---|---:|---:|
| 10 → 14 | \(9.23\times10^{-12}\) | \(5.50\times10^{-9}\) |
| 14 → 18 | \(3.84\times10^{-14}\) | \(1.24\times10^{-9}\) |
| 18 → 22 | \(4.37\times10^{-14}\) | \(3.17\times10^{-9}\) |

For a transparent channel-by-channel comparison, the following table reports the largest exact-Dirac minus leading-NR amplitude difference and the largest relative difference over rows with \(|A_{\rm NR}|>10^{-10}\). The location is isotope, recoil, and rank.

| Channel | Largest absolute amplitude difference | Largest relative difference | Location of relative maximum |
|---|---:|---:|---|
| M | \(4.2765\times10^{-4}\) | 0.5667% | Xe-131, 267 keV, \(J=0\) |
| PhiPP | \(1.4344\times10^{-4}\) | 6.0170% | Xe-131, 267 keV, \(J=0\) |
| Delta | \(1.3504\times10^{-3}\) | 4.4584% | Xe-131, 267 keV, \(J=3\) |
| SigmaPrime | \(6.1444\times10^{-4}\) | 13.3547% | Xe-131, 267 keV, \(J=1\) |

The relative threshold excludes amplitudes that are numerically zero by angular-momentum or parity selection. These percentages quantify the chosen Dirac embedding relative to its leading-NR limit; they are not discrepancies in the leading-NR/native-AFH mapping.

## Pole-width T3 kernel comparison

The spin-resolved dark-matter T3 tensor archive contains eight tensors (the four selected pole-width parameter points re-evaluated for each isotope); the contraction output contains 24 isotope/recoil/rank rows. The phase-space quadrature is \(n_w=24,n_z=16,n_\phi=8\), or 3,072 points per tensor, with the threshold-\(u^2\) mapping. Below are rank-summed kernels for representative \(g_2=1\) cases. The common kernel scale includes the isotope \(C_0\) correction and target spin average, but excludes external mediator/coupling, halo-flux, detector, and cross-section normalizations. These are not rate or cross-section units.

The independent tensor acceptance checks pass: maximum Dirac-trace relative residual \(1.588\times10^{-9}\) against tolerance \(2\times10^{-9}\), maximum Gram Hermiticity residual 0, minimum Gram eigenvalue \(-3.19\times10^{-14}\) against the \(-10^{-10}\) floor, and positive minimum tensor eigenvalue \(4.446\times10^{-4}\).

| Isotope | \(E_R\) | Pole-width benchmark | AFH standard four-channel (continuity-reduced) | AFH L00 + spatial (unreduced) | Direct leading NR | Direct exact Dirac | Exact vs NR |
|---|---:|---|---:|---:|---:|---:|---:|
| Xe-129 | 200 keV | green, \(\delta m=450\) keV | \(1.33719480\times10^8\) | \(1.33720067\times10^8\) | \(1.33720067\times10^8\) | \(1.33789803\times10^8\) | +0.05215% |
| Xe-129 | 267 keV | near-pole 362 keV | \(2.88480254\times10^7\) | \(2.88481719\times10^7\) | \(2.88481719\times10^7\) | \(2.89723483\times10^7\) | +0.43045% |
| Xe-131 | 200 keV | green, \(\delta m=450\) keV | \(1.46268078\times10^8\) | \(1.46268714\times10^8\) | \(1.46268714\times10^8\) | \(1.46358502\times10^8\) | +0.06139% |
| Xe-131 | 267 keV | near-pole 362 keV | \(2.15245661\times10^6\) | \(2.15246745\times10^6\) | \(2.15246745\times10^6\) | \(2.18777001\times10^6\) | +1.64010% |

Across all rank/benchmark rows, the leading-NR direct kernel agrees with the native-AFH L00-plus-spatial contraction to a maximum relative difference of \(1.4312\times10^{-12}\). The direct kernel decomposition closes to \(8.15\times10^{-16}\); its maximum relative imaginary remainder is \(3.35\times10^{-31}\). The smallest exact-kernel-to-largest-kernel ratio is \(6.76\times10^{-15}\), so the sampled kernels are nonnegative up to the recorded numerical floor.

The model-specific exact spatial current contains small \(L\ne J\) vector multipoles; their largest effective reduced amplitudes are \(3.91\times10^{-6}\) in the identity sector and \(1.51\times10^{-6}\) in the spin sector. They remain included in the full current result.

## Signed interference

The entries below are rank-summed, raw signed T3-kernel interference terms for \(g_2=1\), before the common target \(C_0\)/factor-four presentation choice. AFH and leading-NR values preserve the signed terms to relative residuals \(6.76\times10^{-13}\) (M–PhiPP) and \(5.53\times10^{-14}\) (Delta–SigmaPrime) at worst over all rows. The exact-current correction changes their magnitudes while preserving the signs shown.

| Isotope, recoil | M–PhiPP: AFH | leading NR | exact Dirac | Delta–SigmaPrime: AFH | leading NR | exact Dirac |
|---|---:|---:|---:|---:|---:|---:|
| Xe-129, 200 keV | \(-3.38940\times10^5\) | \(-3.38904\times10^5\) | \(-3.29739\times10^5\) | \(-2.93748\times10^{-4}\) | \(-2.93716\times10^{-4}\) | \(-2.75366\times10^{-4}\) |
| Xe-129, 267 keV | \(-7.36835\times10^5\) | \(-7.36758\times10^5\) | \(-6.97391\times10^5\) | \(-4.35981\times10^{-3}\) | \(-4.35933\times10^{-3}\) | \(-4.05849\times10^{-3}\) |
| Xe-131, 200 keV | \(-4.94721\times10^5\) | \(-4.94683\times10^5\) | \(-4.77143\times10^5\) | \(+1.53905\times10^{-4}\) | \(+1.53892\times10^{-4}\) | \(+1.40734\times10^{-4}\) |
| Xe-131, 267 keV | \(-2.61766\times10^5\) | \(-2.61745\times10^5\) | \(-2.47329\times10^5\) | \(+2.93657\times10^{-3}\) | \(+2.93633\times10^{-3}\) | \(+2.62743\times10^{-3}\) |

The Delta–SigmaPrime entries are much smaller than the charge interference and total kernels; their fractional changes should therefore not be read as a comparable fractional change to the total rate. The full amplitudes and all signs are retained in the CSV/JSON rows.

## Continuity and current-component comparison

The exact impulse current does not satisfy the external nuclear continuity equation exactly for these bound HO orbitals with the free-spinor embedding. Its maximum normalized divergence is
\[
\max \frac{\|E_R J^0-Q_AJ^z\|}{(E_R+Q_A)\|J\|}=1.0423\times10^{-3}.
\]
At the nucleon-node level, the on-shell Dirac-Pauli Ward identity is satisfied when the vertex uses \(q_N^0=E_N(p_f)-E_N(p_i)\): the maximum relative \(q_\mu J^\mu\) residual over all \(10{,}648\) order-22 nodes is \(1.15\times10^{-15}\). The external \(q_N^0=E_R\) prescription deliberately uses a different energy transfer and therefore does not inherit that free-nucleon identity node by node. The q-frame tensor rotation basis has a maximum orthogonality residual of \(3.17\times10^{-16}\).
Consequently, comparing the exact \(0,z\) components to AFH's continuity-reduced longitudinal charge gives a maximum 2.0358% relative difference (Xe-131, \(267\) keV, \(J=0\)). The formal leading-NR AFH standard longitudinal-current reduction differs by \(5.1276\times10^{-6}\) at kernel level. In contrast, when the same unreduced L00 plus spatial contraction is used on both sides, the leading-NR direct and AFH kernels agree at \(1.4312\times10^{-12}\). This separates the AFH continuity convention from the covariant higher-order correction.

## Energy-prescription sensitivity

The order-18 internal-\(q_N^0=E_f-E_i\) diagnostic was compared with the external \(q_N^0=E_R\) result. M and PhiPP amplitudes are unchanged. Delta changes by at most \(6.71\times10^{-8}\) absolute; SigmaPrime changes by at most \(2.985\times10^{-5}\) absolute, or 1.575% relative in the Xe-131, 267 keV, \(J=1\) row. The rank-summed T3 kernel changes by at most \(1.76\times10^{-7}\) relative. The base-weight RMS internal energy transfer is 14.0–16.2 MeV for the tabulated cases; the largest quadrature-node values are Gaussian-tail support, not a claimed physical momentum-transfer interval.

## Reproducibility and archive

The native AFH source SHA-256 is dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372. Exact command lines, runtime versions, acceptance thresholds, input/output hashes, row counts, and environment are recorded in results/2026-11-22/p3c_bound_current_amplitude/run_manifest.json; SHA256SUMS.txt records the archived files. The native Wolfram exporter is source/export_p3c_afh_multipoles_native.wls. Python sources export direct amplitudes, generate the pole-width T3 tensors, contract both currents, and regenerate p3c_current_identity_checks.json, which records the nodewise Ward, charge-only, kinematic, q-frame, and tensor checks. The comparison data are in the same result directory as CSV and JSON files.

## Final status boundary

This completed result establishes a reproducible, native-AFH-referenced P3-C amplitude comparison for the selected Xe-129/Xe-131 points within the stated one-body impulse embedding. The direct NR limit matches native AFH multipoles and the corresponding unreduced L00-plus-spatial kernel. The exact Dirac current gives the corrections tabulated above.

native_afh_validated=false was the flag on the earlier P4 production outputs because those outputs had not yet consumed the separately validated native AFH archive. The local regeneration recorded below supersedes that P4 status only; it does not change the scope of this P3-C comparison. `full_current_complete=false` remains correct because the bound-nucleon off-shell current, exchange-current order, and nuclear states are not specified by this impulse calculation.

## 2026-11-22 native P4 regeneration and χEFT status

The P4 fold was regenerated from the validated local AFH v6.0 eight-response
archive. The five-point/two-halo output contains 10 recoil/halo rows and 70
isotope rows; the continuous spectrum contains 162 recoil/halo rows and 1,134
isotope rows. Both report `native_afh_validated=true` and
`full_current_complete=false`.

On matched grids, the maximum relative rate difference between the old
Python-AFH and new native-AFH calculations is (2.6459\times10^{-11}) on the
five-point grid and (1.5275\times10^{-10}) on the 81-point grid. Both use
the same common linear PSD-preserving response interpolation. The speed
24→48-node maximum difference is (8.9979\times10^{-6}) (gate
(10^{-5})); the phase-space refinement maximum is (2.5525\times10^{-9})
(gate (10^{-8})). The full record and checksums are under
`results/2026-11-22/p4_native_afh_regeneration/`.

The user selected a consistent χEFT force/current framework for the broader
current. The model-input audit found no matching NN+3N Hamiltonian or
consistent seven-isotope nuclear states/operators in this repository. In
addition, the published vector two-body term through \(\mathcal O(p^3)\) is
for an axial-DM/vector-nucleon channel, whereas this branch is vector–vector;
that operator cannot be imported as the missing `VV` exchange current.
`FULL_CURRENT_MODEL_CLOSURE_AUDIT.md` records the current model and required
inputs. The direct P3-C conclusion remains limited to the selected one-body
impulse embedding.
