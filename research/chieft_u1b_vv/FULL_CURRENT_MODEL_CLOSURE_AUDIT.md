# Full nuclear-current closure audit

**Date:** 2026-11-22
**User-selected framework:** mutually consistent chiral EFT (χEFT) force and
current
**Formal state:** `full_current_complete=false`

## Decision

The selected framework is a direction, not yet a numerically specified
Hamiltonian/current pair. The current repository cannot produce a complete
χEFT P3-C/P4 prediction from its AFH response archive alone. The native AFH
P4 fold can be regenerated and reported independently; that does not promote
the nuclear-current model to complete.

The existing direct covariant P3-C calculation is complete only for the
declared positive-energy free-spinor Dirac–Pauli one-body impulse embedding
and its selected Xe-129/Xe-131 recoil points. Its leading-NR/native-AFH
amplitude residual is (3.575\times10^{-14}). It does not determine a unique
off-shell bound-nucleon energy assignment, exchange current, or χEFT nuclear
state, so it is not relabeled as a χEFT calculation.

## Channel-specific χEFT finding

This repository's interaction branch is vector–vector (`VV`). In the
published dark-matter chiral power counting through \(\mathcal O(p^3)\), the
only two-body **vector** current appears in the axial-DM/vector-nucleon (`AV`)
channel. It cannot be copied into the `VV` calculation as its exchange
current. A `VV` two-body contribution requires a declared higher-order
operator set and matching convention. This channel distinction is established
in [Hoferichter, Klos & Schwenk (2015)](https://arxiv.org/abs/1503.04811).

Even if the model were deliberately truncated at an order where the `VV`
two-body operator is absent, that would not make the current calculation
complete here: the P4/P3-C matrix elements still use AFH default shell-model
inputs rather than nuclear states generated from the selected χEFT force.
Using the published `AV` exchange current, or electromagnetic exchange
currents, would change the interaction channel/source and would not satisfy
the user's selected consistent `VV` model.

## Inputs present and missing

| Requirement | Available evidence | Status for a χEFT P4 current |
|---|---|---|
| AFH response database | Seven isotopes × eight responses × 2,401 (q_A) points; native Wolfram/Python comparison passed | Valid AFH input; not a χEFT many-body calculation |
| P3-C one-body amplitude | Xe-129/Xe-131 AFH OBDMEs embedded in positive-energy free Dirac spinors | Model-specific impulse result only |
| Quark-source matching | Universal (u+d+s) vector source and effective-(g_{VN}) convention in the P4 benchmark | Need coefficients, scale and matching carried into χEFT operators |
| χEFT NN+3N Hamiltonian | No pinned interaction or matrix elements for the seven isotopes in this repository | Missing order, regulator/cutoff, LECs and valence-space Hamiltonian |
| χEFT nuclear states | No consistently generated states for Xe-128,129,130,131,132,134,136 | Missing; AFH one-body densities do not establish force/current consistency |
| Two-body vector operators/densities | No matching `VV` higher-order operator or Xe pair-density archive | Needed if the chosen order includes `VV` exchange currents |
| Current conservation and relativistic terms | P4 longitudinal reduction is explicitly conditional; P3-C is a free-spinor impulse embedding | Must be derived/checked against the selected Hamiltonian and truncation |

Chiral-force VS-IMSRG calculations exist for Xe-129 and Xe-131 spin-dependent
responses, including axial two-body-current effects, but those results do not
provide the required seven-isotope vector-current inputs. See
[Hu et al. (2022)](https://arxiv.org/abs/2109.00193). The requirement that
currents and nuclear dynamics use a consistent scheme/regulator is also
discussed in [Pastore et al. (2009)](https://arxiv.org/abs/0906.1800) and
[Krebs (2020)](https://arxiv.org/abs/2008.00974).

## What would close the calculation

1. Pin the full flavor-source matching, renormalization scale, and the chosen
   `VV` chiral order; state whether the calculation is truncated before its
   higher-order two-body vector current or includes that operator.
2. Supply or generate the matching NN+3N interaction, order, regulator,
   fitted LECs and effective valence-space Hamiltonian.
3. Supply consistently generated Xe ground states and one-body densities for
   all seven AFH isotopes; if the selected `VV` order contains a two-body
   current, also supply consistently evolved two-body operators/densities.
4. Recompute complex Xe-129/Xe-131 P3-C amplitudes at the pinned kinematics,
   compare direct, multipole and response contractions before squaring, and
   check the continuity relation with the selected force/current pair.
5. Regenerate P4 using those same inputs, then report chiral truncation,
   regulator, basis and many-body uncertainties separately from halo and
   quadrature uncertainty.

Until those inputs and checks exist, retain
`native_afh_validated=true` for the new native AFH fold and
`full_current_complete=false` for its leading-static one-body candidate.


## 2026-11-24: current T3 replay status and source-specific χEFT gate

The exact-mass Borah2020+Sufian2017 one-body impulse model now has an independently replayed full-complex T3 contraction: 24 benchmark/rank rows and 74 common numeric scalar fields are exactly identical between primary and independent contraction outputs. The maximum signed (M)-(Phi''), signed (Delta)-(Sigma'), kernel decomposition and imaginary residuals are recorded in `results/2026-11-23-p3c-independent-contraction-replay/output/independent_fullcomplex_contraction_o22_summary.json`.

The independent replay consumes the primary P3-C direct-current multipoles; it independently reconstructs the magnetic-substate current matrices and full complex T3 Lorentz contraction. Therefore the T3 contraction machinery is validated for this impulse model, while the direct-current multipole generation is not independently regenerated by that replay. The four-point spherical quadrature replay separately validates the current integration rule at Xe-129/Xe-131, 200/267 keV.

The next χEFT calculation is source-specific: the universal (u+d+s) vector source is (3J_B^mu) in the three-flavor theory. The published (O(p^3)) `AV` two-body vector current is not the `VV` operator required here. The new [VV source-matching audit](CHI_EFT_VV_SOURCE_MATCHING_AUDIT_2026-11-24.md) identifies the exact operator derivation and the reusable VS-IMSRG xenon method.

Once the χEFT (VV) current and its consistently evolved Xe nuclear matrix elements are available, reuse the validated T3 contraction implementation, but generate a new run manifest and independent full-complex replay using the χEFT multipoles. Do not transfer the impulse-model PASS to the χEFT calculation.
