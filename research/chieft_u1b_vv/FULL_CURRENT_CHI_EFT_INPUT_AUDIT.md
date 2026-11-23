# Full nuclear vector current: chiEFT model-input audit

**Updated:** 2026-11-22
**Selected direction:** a current and nuclear Hamiltonian derived consistently in chiral EFT (chiEFT).
**Promotion state:** `full_current_complete=false`.

## Decision captured

The model direction is now chiEFT, following the user's selection. This is a
framework choice, not yet a fully specified numerical Hamiltonian/current
model. The present P3-C result remains the completed direct covariant
comparison for its stated positive-energy free-spinor impulse embedding. It
does not become a chiEFT result by relabeling the one-body current.

The P4 nuclear-source convention to be matched is the flavor-universal vector
source used by the current benchmark, proportional to
\(\bar u\gamma^\mu u+\bar d\gamma^\mu d+\bar s\gamma^\mu s\). The actual
matching coefficients and their renormalization-scale convention must be
carried into the chiral source before a two-nucleon operator is evaluated.
Electromagnetic exchange-current formulae are not a drop-in replacement for
this source: their flavor/isospin weights and the dark-matter-side Lorentz
channel have to match the present interaction.

## What the repository currently supplies

| Needed object | Present input | Gap for the selected chiEFT calculation |
|---|---|---|
| AFH one-body nuclear responses | AFH v6.0 default-density responses; seven isotopes; eight native-validated response columns | These are one-body AFH responses, not a consistent chiEFT force/current calculation |
| P3-C one-body amplitude | Xe-129/Xe-131 AFH OBDMEs embedded in positive-energy free Dirac spinors at selected recoil points | No chiEFT eigenstates or model-matched bound-nucleon wavefunctions |
| Nuclear Hamiltonian | AFH shell-model density tables and existing response code | No pinned chiEFT NN+3N matrix elements, order, regulator/cutoff, fitted LEC set, or effective valence-space Hamiltonian for the seven isotopes |
| Two-body vector current | No chiEFT two-body operator implementation in the P4/P3-C path | Need operators matched to the universal vector source and the chosen dark-matter Lorentz channel, at a declared order and with the same regularization as the force |
| Nuclear pair densities | No Xe two-body transition-density archive for the selected states | Need two-body densities (or consistently evolved two-body operators) to contract exchange currents; one-body OBDMEs cannot reconstruct them |
| Current conservation | Current P4 Coulomb reduction is explicitly conditional on nuclear-current continuity | Need to verify \(q_\mu J^\mu=[H,\rho]\) (up to the selected convention/order) with the chosen Hamiltonian and current |
| Relativistic effects | Exact free-spinor one-body impulse diagnostic | In a nonrelativistic chiEFT calculation, relativistic effects enter as ordered corrections to currents/forces; a relativistic mean-field Dirac orbital cannot be mixed in without a separately consistent covariant model |

Repository inventory found the AFH one-body response/OBDME and P3-C impulse
paths, but no matching chiEFT NN+3N interaction, Xe two-body density, or
two-body vector-current implementation. Thus the native eight-response AFH
archive remains valid as an AFH response validation; it is not sufficient
input for the selected chiEFT many-body current.

## Why an explicit chiEFT specification is necessary

1. Chiral power counting for dark-matter currents is channel-dependent. For
   the repository's vector–vector branch, the Hoferichter–Klos–Schwenk
   calculation states that the only two-body vector current through
   \(\mathcal O(p^3)\) occurs in the axial-DM/vector-nucleon (`AV`) channel.
   It therefore cannot be inserted into this vector–vector (`VV`) fold as an
   exchange-current correction. A `VV` two-body term would require a declared
   higher-order treatment and consistent matching. The paper derives the
   channel-dependent one- and two-body currents through third order; see
   [Hoferichter, Klos & Schwenk (2015)](https://arxiv.org/abs/1503.04811).
2. Nuclear forces and currents must share the chosen scheme and regulator for
   a consistent continuity relation. This relation is central to replacing
   the longitudinal nuclear current by a charge response. See
   [Pastore et al. (2009)](https://arxiv.org/abs/0906.1800) for a potential
   and consistent electromagnetic-current construction, and
   [Krebs (2020)](https://arxiv.org/abs/2008.00974) for the regularization
   requirement in chiEFT currents. These electromagnetic operators still
   require rematching before use for the present flavor-universal dark-matter
   source.
3. Heavy xenon calculations using chiEFT interactions and IMSRG demonstrate
   that a consistent nuclear-structure route is feasible for selected
   isotopes, but the cited Xe study is a spin-dependent calculation with
   axial-vector two-body-current effects. It is not the required vector
   response archive for all seven AFH isotopes. See
   [Hu et al. (2022)](https://arxiv.org/abs/2109.00193).

## Work required before computing a chiEFT P3-C/P4 result

1. Pin the quark-level Wilson coefficients and source convention, including
   the strange-vector contribution, and derive the corresponding chiral
   source operators.
2. Choose and pin an NN+3N interaction, chiral order, regulator/cutoff, LEC
   set, many-body solver, and valence space/basis. In particular, specify
   whether the `VV` calculation is truncated where its two-body vector current
   is absent, or whether a higher-order `VV` operator is included. Use the
   same regulator and compatible order for the current operators.
3. Derive the one- and two-body vector charge/current operators for this
   source, including the declared relativistic corrections and exchange
   currents. Record the continuity identity and any transverse unconstrained
   pieces.
4. Generate or obtain consistently evolved nuclear states and one-/two-body
   densities for Xe-128,129,130,131,132,134,136. Xe-129 and Xe-131 must also
   support complex amplitude-level P3-C comparisons at the existing selected
   kinematic points.
5. Contract the one- and two-body operators before squaring, preserving
   signed interference. Build any response functions needed beyond the
   current eight one-body AFH channels, and compare direct, multipole, and
   response-level amplitudes.
6. Re-run P4 with the same halo, Pole-width, mediator-transfer, isotope
   abundance, and recoil-grid conventions. Report chiral truncation, regulator,
   basis, and many-body uncertainties separately from halo and quadrature
   uncertainties.

## Current status

The selected chiEFT framework is recorded, but the numerical model is not yet
defined: no interaction/order/regulator/LECs or matching Xe many-body states
and densities are present. The current branch is `VV`; the published
\(\mathcal O(p^3)\) `AV` exchange current is not its missing term. The present
work can therefore close native AFH P4 production while preserving
`full_current_complete=false`. A chiEFT full-current fold cannot be generated
from the archived AFH response CSVs alone.


## 2026-11-24 source-specific VV matching update

The universal three-flavor source used by the production benchmark is the (U(1)_B) vector source, (J^mu_{uds}=3J^mu_B), in the three-flavor low-energy theory. The source-specific chiral matching is now recorded in [`CHI_EFT_VV_SOURCE_MATCHING_AUDIT_2026-11-24.md`](CHI_EFT_VV_SOURCE_MATCHING_AUDIT_2026-11-24.md).

The channel distinction is explicit: the two-body vector current derived through (O(p^3)) in Hoferichter–Klos–Schwenk belongs to `AV`, not this project's `VV` branch. It must not be inserted into the VV fold. Higher-order electromagnetic charge operators are construction references only; their source, isospin projection, regulator and unitary convention must be matched to (U(1)_B) before reuse.

The Hu et al. VS-IMSRG/NAT Xe-129/Xe-131 workflow is a reusable nuclear-structure method reference. Its public notebook/repository supplies spin-dependent response fits, not the vector-charge/current operators or VV matrix elements required here.

The independent full-complex T3 contraction replay is already PASS for the current Borah2020+Sufian2017 exact-mass one-body impulse multipoles: 24 rows and 74 common numeric fields agree exactly with the primary output. That validates the contraction replay, not an independent regeneration of the current multipoles. Once the (U(1)_B) χEFT current is generated, the same T3 contraction must be rerun with those new complex multipoles and independently replayed.

**Current boundary:** retain `full_current_complete=false`. The next physics gate is the explicit (U(1)_B) χEFT operator derivation for a pinned NN+3N interaction/regulator, followed by Xe-129/Xe-131 states and consistently transformed operators/densities.

The public software audit identifies IMSRG++ as a reusable VS-IMSRG/NAT engine and NuHamil as a complementary NN+3N matrix-element and multipole-generation tool. One source-level correction is important: NuHamil's `Tmag_2B_J1_Tz0_Q100-NLO` is an isovector operator (its code uses the rank-one isospin tensor); the `Tz0` label denotes the projection, not isospin rank. It is not the (U(1)_B) current. NuHamil also does not implement the required higher-order vector current or two-body charge operator. See the detailed code audit in `CHI_EFT_VV_SOURCE_MATCHING_AUDIT_2026-11-24.md`.
