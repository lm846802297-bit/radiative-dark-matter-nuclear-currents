"""Three-body recoil kinematics and halo-speed helpers for Fig. 2.

This module implements only the tree-level kinematic support condition
printed in the manuscript and the standard untruncated Maxwell-Boltzmann
lab-frame speed distribution used as a documented DMFormFactor default.

Units:
  masses, recoil energy: GeV
  speeds passed to the halo function: km/s
  vmin_three_body_beta returns v/c (dimensionless)
"""
from __future__ import annotations

import math

C_KM_S = 299_792.458


def nuclear_recoil_momentum(mA_GeV: float, ER_GeV: float) -> float:
    """Return Q_A=|bold(q_A)|=sqrt(2 m_A E_R+E_R^2), in GeV."""
    if mA_GeV <= 0.0 or ER_GeV < 0.0:
        raise ValueError("Require mA>0 and ER>=0")
    return math.sqrt(2.0 * mA_GeV * ER_GeV + ER_GeV**2)


def momentum_transfer(mA_GeV: float, ER_GeV: float) -> float:
    """Backward-compatible alias for the spatial nuclear transfer Q_A."""
    return nuclear_recoil_momentum(mA_GeV, ER_GeV)


FORMAL_NUCLEON_TRANSFER_SCHEME = (
    "one-body impulse-current external-transfer assignment q_N^mu=q_A^mu; t_N=t_A=-2*m_A*E_R"
)
STATIC_NUCLEON_TRANSFER_SCHEME = (
    "static diagnostic q_N^0=0; t_N=-Q_A^2; not formal B"
)


def transfer_kinematics_from_q_omega(
    qA_GeV: float,
    omega_GeV: float,
    *,
    nucleon_transfer_scheme: str = FORMAL_NUCLEON_TRANSFER_SCHEME,
) -> dict:
    """Build the external nuclear transfer and one-body vertex transfer.

    In the q-aligned frame with metric (+---), q_A^mu=(omega,0,0,Q_A).
    Formal B assigns the same external four-transfer to the one-body
    nucleon-current vertex, q_N^mu=q_A^mu. For a free on-shell impulse
    benchmark this can be realized as k_N'^mu-k_N^mu=q_A^mu. In a bound
    nucleus this is the external transfer entering the one-body current;
    it does not reconstruct each bound nucleon's microscopic off-shell
    initial and final four-momenta. The q_N^0=0 static choice is diagnostic.
    """
    qA = float(qA_GeV)
    omega = float(omega_GeV)
    if not math.isfinite(qA) or not math.isfinite(omega):
        raise ValueError("q_A and omega must be finite")
    if qA < 0.0 or omega < 0.0:
        raise ValueError("q_A and omega must be non-negative")
    if nucleon_transfer_scheme == FORMAL_NUCLEON_TRANSFER_SCHEME:
        qN0 = omega
        scheme_id = "one_body_impulse_external_transfer_qN_equals_qA"
        qN_definition = (
            "external four-momentum at the one-body nucleon-current vertex; "
            "set q_N=q_A by impulse matching, not a microscopic bound-nucleon "
            "off-shell momentum reconstruction"
        )
    elif nucleon_transfer_scheme == STATIC_NUCLEON_TRANSFER_SCHEME:
        qN0 = 0.0
        scheme_id = "static_qN0_zero_diagnostic"
        qN_definition = (
            "static diagnostic vertex assignment q_N^0=0 and "
            "boldsymbol q_N=boldsymbol q_A; not the physical nuclear energy transfer"
        )
    else:
        raise ValueError(f"unknown nucleon transfer scheme: {nucleon_transfer_scheme}")

    Q_A2 = qA * qA
    tA = omega * omega - Q_A2
    tN = qN0 * qN0 - Q_A2
    if tA > 1.0e-14 or tN > 1.0e-14:
        raise ValueError("nuclear and nucleon transfers must be spacelike or zero")

    qA4 = [omega, 0.0, 0.0, qA]
    qN4 = [qN0, 0.0, 0.0, qA]
    return {
        "qA_four_vector_GeV": qA4,
        "qN_four_vector_GeV": qN4,
        "qA_energy_GeV": omega,
        "qN_energy_GeV": qN0,
        "Q_A_GeV": qA,
        "Q_A2_GeV2": Q_A2,
        "tA_GeV2": tA,
        "tN_GeV2": tN,
        "metric_signature": "+---",
        "four_vector_frame": "q-aligned frame with boldsymbol q_A along +z",
        "qA_definition": "P_A' - P_A",
        "qN_definition": qN_definition,
        "nucleon_transfer_scheme": nucleon_transfer_scheme,
        "nucleon_transfer_scheme_id": scheme_id,
    }


def recoil_transfer_kinematics(
    mA_GeV: float,
    ER_keV: float,
    *,
    nucleon_transfer_scheme: str = FORMAL_NUCLEON_TRANSFER_SCHEME,
) -> dict:
    """Return exact recoil q_A^mu and the selected one-body q_N^mu prescription.

    mA_GeV is in GeV and ER_keV is in keV. The formal q_N=q_A assignment
    is the default; the static q_N^0=0 prescription is diagnostic-only.
    """
    mA = float(mA_GeV)
    ER_keV = float(ER_keV)
    ER = ER_keV * 1.0e-6
    if not math.isfinite(mA) or not math.isfinite(ER_keV) or mA <= 0.0 or ER < 0.0:
        raise ValueError("invalid nuclear mass or recoil energy")
    Q_A = math.sqrt(2.0 * mA * ER + ER * ER)
    result = transfer_kinematics_from_q_omega(
        Q_A, ER, nucleon_transfer_scheme=nucleon_transfer_scheme
    )
    expected_tA = -2.0 * mA * ER
    if not math.isclose(
        result["tA_GeV2"], expected_tA, rel_tol=2.0e-12, abs_tol=2.0e-18
    ):
        raise ArithmeticError("q_A four-vector violates the exact nuclear mass-shell relation")
    result.update({
        "mA_GeV": mA,
        "ER_keV": ER_keV,
        "mass_shell_tA_GeV2": expected_tA,
    })
    return result



def wmax_squared(
    mchi_GeV: float,
    mA_GeV: float,
    ER_GeV: float,
    speed_km_s: float,
) -> float:
    """W_max^2 from the manuscript, in GeV^2."""
    if mchi_GeV <= 0.0 or mA_GeV <= 0.0 or ER_GeV < 0.0:
        raise ValueError("Require positive masses and ER>=0")
    beta = speed_km_s / C_KM_S
    if not 0.0 <= beta < 1.0:
        raise ValueError("Incident speed must satisfy 0<=v<c")
    gamma = 1.0 / math.sqrt(1.0 - beta * beta)
    Echi = mchi_GeV * gamma
    pchi = Echi * beta
    Q_A = nuclear_recoil_momentum(mA_GeV, ER_GeV)
    return (
        mchi_GeV**2
        + 2.0 * pchi * Q_A
        - 2.0 * (Echi + mA_GeV) * ER_GeV
    )


def vmin_three_body_beta(
    mchi_GeV: float,
    ms_GeV: float,
    mA_GeV: float,
    ER_GeV: float,
) -> float:
    r"""Minimum incident speed for chi1+A -> chi1+A+s.

    It solves W_max^2(E_R,v)=(m_chi+m_s)^2 using the W_max expression
    printed in the manuscript. The result is dimensionless beta=v/c.

    Define q=sqrt(2 m_A E_R+E_R^2) and
      a=m_s+(m_A E_R+m_s^2/2)/m_chi.
    Then
      beta_min=(q E_R+a sqrt(q^2-E_R^2+a^2))/(q^2+a^2).

    This is not the usual elastic-recoil v_min formula.
    """
    if mchi_GeV <= 0.0 or mA_GeV <= 0.0 or ms_GeV < 0.0 or ER_GeV < 0.0:
        raise ValueError("Require mchi,mA>0 and ms,ER>=0")
    if ER_GeV == 0.0 and ms_GeV == 0.0:
        return 0.0
    Q_A = nuclear_recoil_momentum(mA_GeV, ER_GeV)
    a = ms_GeV + (mA_GeV * ER_GeV + 0.5 * ms_GeV**2) / mchi_GeV
    denom = Q_A * Q_A + a * a
    if denom == 0.0:
        return 0.0
    radicand = Q_A * Q_A - ER_GeV * ER_GeV + a * a
    beta = (Q_A * ER_GeV + a * math.sqrt(radicand)) / denom
    return beta


def vmin_three_body_km_s(
    mchi_GeV: float,
    ms_GeV: float,
    mA_GeV: float,
    ER_GeV: float,
) -> float:
    """Three-body minimum speed in km/s."""
    return C_KM_S * vmin_three_body_beta(mchi_GeV, ms_GeV, mA_GeV, ER_GeV)


def mb_lab_speed_pdf(
    speed_km_s: float,
    v0_km_s: float = 220.0,
    vE_km_s: float = 232.0,
) -> float:
    """Normalized lab-frame speed PDF for an untruncated MB halo.

    This is the analytic speed distribution corresponding to
    f_gal(u) proportional to exp(-u^2/v0^2), boosted by vE.
    It integrates to one over speed in [0,infinity). No escape-speed
    cutoff is applied here.
    """
    if speed_km_s < 0.0 or v0_km_s <= 0.0 or vE_km_s <= 0.0:
        raise ValueError("Require speed>=0 and positive v0,vE")
    v = speed_km_s
    return (
        v
        / (math.sqrt(math.pi) * v0_km_s * vE_km_s)
        * (
            math.exp(-((v - vE_km_s) / v0_km_s) ** 2)
            - math.exp(-((v + vE_km_s) / v0_km_s) ** 2)
        )
    )
