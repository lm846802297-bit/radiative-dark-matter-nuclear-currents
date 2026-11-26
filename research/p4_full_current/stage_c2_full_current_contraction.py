#!/usr/bin/env python3
"""P4-FC Stage C2: contract the frozen C1 Lorentz grid with Stage B currents.

C1 and Stage B are immutable inputs. This stage performs no phase-space
integration and no nuclear-current regeneration. It constructs the rank-
resolved four-current matrices, contracts L^{mu nu} J_mu J_nu^*, applies the
AFH physical factor four and C0(A) exactly once, and folds the stored speed
nodes into the two halo spectra.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from sympy import Integer, Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j

ROOT = Path(__file__).resolve().parents[3]
TOPIC = ROOT / "research" / "three_body_full_current_nucleon_match"
HERE = TOPIC / "p4_full_current"
C1_DIR = HERE / "output" / "stage-c1-full-lorentz-grid"
STAGE_B_CSV = HERE / "output" / "stage-b-direct-current" / "direct_bound_current_o14.csv"
STAGE_A_RESPONSE_CSV = HERE / "output" / "native-seven-isotope-81-point" / "afh_v6_p4_fc_7isotope_81point_native_responses.csv"
STAGE_B_VALIDATION_REPORT = HERE / "output" / "stage-b-direct-current" / "stage_b_validation_report.json"
B_FREEZE = ROOT / "research" / "nucleon_form_factor_B_baseline" / "B_DATA_VERSION_FREEZE_2026-11-24.json"
ISOTOPE_INPUT = ROOT / "research" / "afh_default_density_natural_xenon_halo" / "input" / "natural_xenon_isotope_weights.csv"
STAGE_B_SCRIPT = HERE / "direct_bound_current_seven_isotope.py"
STAGE_B_TEST_SCRIPT = HERE / "test_stage_b_direct_current.py"
C1_SCRIPT = HERE / "stage_c1_full_lorentz_tensor_grid.py"
RATE_REFERENCE = TOPIC / "native_charge_sector_halo_fold.py"
C2_TEST_SCRIPT = HERE / "test_stage_c2_full_current_contraction.py"

ISOTOPES = (128, 129, 130, 131, 132, 134, 136)
HALOS = ("shm_cutoff", "mb_untruncated")
RECOIL_COUNT = 81
SPEED_ORDER = 24
M_N_GEV = 0.938918
MCHI_GEV = 1000.0
MS_GEV = 200.0e-6
DELTA_M_GEV = 1000.0e-6
MV_GEV = 10.0
G2 = 1.0
COUPLING_PRODUCT_SQUARED = 1.0e-8
RHO_GEV_CM3 = 0.3
WIDTH_SCHEME = "pole"
MAPPING = "threshold_u2"

# Immutable fingerprints from the C1 production manifest and B freeze.
EXPECTED_SHA256 = {
    "c1_tensor_npz": "98898fbdac7a08a81453972d013ad027e46679a02e67d41a3d59a529f1671c76",
    "c1_tensor_index": "3647a0b4261be85b146ffe040f4e0241fd1508d85cc0d5f1f9b502dc1b33c06c",
    "c1_run_manifest": "18d36a643205ad0f8cde581cb04f1ff6dbd3f2aa6e86e1877837b6b8678fc8e8",
    "c1_convergence_report": "1cb9df212451374652e0713a477f86e8954e3f9d8aee9aa71123e2910b3a0878",
    "stage_b_csv": "0b335b39cd288beda1367fd977e311d595c2a88c2c5e7dfcaa9456cf8e26d355",
    "stage_b_validation_report": "5f0272f2c976ab894d8e69fc6930c609c13a0aba93ee445c5b83710c0ff4f313",
    "stage_a_responses": "39a5b8c2401eb5f4e5acabb0a4f38b9a9d1fa7103d12d208462b3fcefef57c48",
    "isotope_input": "fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991",
    "b_freeze": "2b0b8c58baf559c176fe9e64d14f75973d7987b072fa6af2275f09430d5fbe88",
    "stage_b_script": "78058780862acb06d5c3fd86957dead6d16d9a8e3041c9b4bdc3d5450189b75e",
    "c1_script": "c82b5a47c12fb59967dd41cf0de54c01a3c2a42c0e16152e2e0484babaa91bbc",
    "rate_reference": "71ed1c7f4b68d243f7a6b7f65727066e13bb958852630a0ab61f9d690d95f5c1",
}

METRIC = np.asarray([1.0, -1.0, -1.0, -1.0], dtype=float)
CHANNELS = ("M", "PhiPP", "Delta", "SigmaPrime")
C_KM_S = 299792.458
C_CM_S = 2.99792458e10
GEV2_TO_CM2 = 0.3893793721e-27
AVOGADRO = 6.02214076e23
U_GEV = 0.93149410242
M_E_U = 0.000548579909065
Z_XE = 54


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_blob_sha1(data: bytes) -> str:
    payload = b"blob " + str(len(data)).encode("ascii") + b"\0" + data
    return hashlib.sha1(payload).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def cval(row: dict, name: str) -> complex:
    return complex(float(row[name + "_re"]), float(row[name + "_im"]))


def pair_array(value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape[-1] != 2 or not np.all(np.isfinite(arr)):
        raise ValueError("complex-pair array must be finite and end in [real,imag]")
    return arr[..., 0] + 1j * arr[..., 1]


@lru_cache(maxsize=None)
def w3(j2: int, K: int, mf2: int, mu: int, mi2: int) -> float:
    return float(wigner_3j(
        Rational(j2, 2), Integer(K), Rational(j2, 2),
        Rational(-mf2, 2), Integer(mu), Rational(mi2, 2)))


@lru_cache(maxsize=None)
def cg_vector(L: int, mu: int, K: int) -> float:
    return float(clebsch_gordan(L, 1, K, 0, mu, mu))


def zero_matrix(j2: int) -> np.ndarray:
    n = j2 + 1
    return np.zeros((4, n, n), dtype=complex)


def build_rank_components(row: dict, representation: str) -> dict[str, np.ndarray]:
    """Rebuild rank-K current matrices, including spin-zero K=0,L=1 rows.

    The only generalization relative to the Xe-129/131 P3C helper is that
    2*J_A is read from the immutable Stage B row instead of inferred from A.
    """
    if representation not in {"exact", "leading_nr"}:
        raise ValueError(representation)
    K = int(row["J"])
    j2_float = 2.0 * float(row["J_A"])
    j2 = int(round(j2_float))
    if abs(j2_float - j2) > 1.0e-12 or j2 < 0:
        raise ValueError(f"invalid target spin J_A={row['J_A']}")
    m2_values = tuple(range(-j2, j2 + 1, 2))
    m_index = {m2: i for i, m2 in enumerate(m2_values)}
    out = {
        "charge_identity": zero_matrix(j2),
        "charge_spin": zero_matrix(j2),
        "spatial_identity_LJm1": zero_matrix(j2),
        "spatial_identity_LJ": zero_matrix(j2),
        "spatial_identity_LJp1": zero_matrix(j2),
        "spatial_spin_LJm1": zero_matrix(j2),
        "spatial_spin_LJ": zero_matrix(j2),
        "spatial_spin_LJp1": zero_matrix(j2),
    }
    for mi2 in m2_values:
        ii = m_index[mi2]
        mf2 = mi2
        if mf2 in m_index:
            ff = m_index[mf2]
            phase = (-1.0) ** ((j2 - mf2) // 2)
            weight = w3(j2, K, mf2, 0, mi2)
            plane = (1j ** K) * math.sqrt(4.0 * math.pi * (2 * K + 1))
            for piece, column in (
                ("charge_identity", f"{representation}_charge_identity"),
                ("charge_spin", f"{representation}_charge_spin"),
            ):
                out[piece][0, ff, ii] += phase * weight * plane * cval(row, column)

    for mi2 in m2_values:
        ii = m_index[mi2]
        sph_parts = {
            name: np.zeros((3, len(m2_values)), dtype=complex)
            for name in out if name.startswith("spatial_")
        }
        for mu in (-1, 0, 1):
            mf2 = mi2 + 2 * mu
            if mf2 not in m_index:
                continue
            ff = m_index[mf2]
            phase = (-1.0) ** ((j2 - mf2) // 2)
            weight = w3(j2, K, mf2, mu, mi2)
            if weight == 0.0:
                continue
            for L, suffix in ((K - 1, "LJm1"), (K, "LJ"), (K + 1, "LJp1")):
                if L < 0:
                    continue
                plane = 4.0 * math.pi * (1j ** L) * math.sqrt((2 * L + 1) / (4.0 * math.pi))
                basis = phase * weight * plane * cg_vector(L, mu, K)
                for kind in ("identity", "spin"):
                    column = f"{representation}_vector_{kind}_{suffix}"
                    name = f"spatial_{kind}_{suffix}"
                    sph_parts[name][mu + 1, ff] += basis * cval(row, column)
        for name, sph in sph_parts.items():
            for ff in range(len(m2_values)):
                # Condon-Shortley inverse, spherical order [-1,0,+1].
                out[name][1, ff, ii] = (sph[0, ff] - sph[2, ff]) / math.sqrt(2.0)
                out[name][2, ff, ii] = 1j * (sph[0, ff] + sph[2, ff]) / math.sqrt(2.0)
                out[name][3, ff, ii] = sph[1, ff]

    out["charge"] = out["charge_identity"] + out["charge_spin"]
    out["spatial_identity"] = sum(
        (out[f"spatial_identity_{suffix}"] for suffix in ("LJm1", "LJ", "LJp1")),
        zero_matrix(j2),
    )
    out["spatial_spin"] = sum(
        (out[f"spatial_spin_{suffix}"] for suffix in ("LJm1", "LJ", "LJp1")),
        zero_matrix(j2),
    )
    out["spatial"] = out["spatial_identity"] + out["spatial_spin"]
    out["full"] = out["charge"] + out["spatial"]
    out["delta_LJ"] = out["spatial_identity_LJ"]
    out["sigma_LJ"] = out["spatial_spin_LJ"]
    return out


def bilinear(L: np.ndarray, A: np.ndarray, B: np.ndarray, JA: float) -> complex:
    lower_A = METRIC[:, None, None] * A
    lower_B = METRIC[:, None, None] * B
    return np.einsum("mn,mfi,nfi->", L, lower_A, lower_B.conjugate()) / (2.0 * JA + 1.0)


def bilinear_explicit(L: np.ndarray, A: np.ndarray, B: np.ndarray, JA: float) -> complex:
    """Independent loop implementation used for contraction cross-checks."""
    if A.shape != B.shape or A.shape[0] != 4 or L.shape != (4, 4):
        raise ValueError("invalid Lorentz/current dimensions")
    total = 0.0j
    n_final, n_initial = A.shape[1], A.shape[2]
    for mu in range(4):
        for nu in range(4):
            for mf in range(n_final):
                for mi in range(n_initial):
                    total += (
                        L[mu, nu] * METRIC[mu] * A[mu, mf, mi]
                        * METRIC[nu] * B[nu, mf, mi].conjugate()
                    )
    return total / (2.0 * JA + 1.0)


def proca_numerator_mixed_matrix(ER_GeV: float, Q_A_GeV: float, mV_GeV: float) -> np.ndarray:
    """Return M^mu{}_nu for the Proca numerator g_{mu nu}-q_mu q_nu/mV^2.

    The fixed convention is q_A^mu=(E_R,0,0,+Q_A), metric (+---).
    This is the massive-vector numerator, not a transverse projector. The
    scalar propagator denominator is applied separately in the rate kernel.
    """
    if not all(math.isfinite(x) for x in (ER_GeV, Q_A_GeV, mV_GeV)):
        raise ValueError("Proca kinematics must be finite")
    if Q_A_GeV < 0.0 or mV_GeV <= 0.0:
        raise ValueError("Proca Q_A must be nonnegative and mV positive")
    q_upper = np.asarray([ER_GeV, 0.0, 0.0, Q_A_GeV], dtype=float)
    q_lower = METRIC * q_upper
    return np.eye(4, dtype=float) - np.outer(q_lower, q_upper) / (mV_GeV * mV_GeV)


def apply_proca_numerator(
    L: np.ndarray, ER_GeV: float, Q_A_GeV: float, mV_GeV: float
) -> np.ndarray:
    """Contract the Proca numerator into both indices of L^{mu nu}.

    From M=T^mu(g_{mu nu}-q_mu q_nu/mV^2)J^nu, define
    T_eff^rho=T^mu M^mu{}_rho. The tensor entering the nuclear-current
    bilinear is therefore L_eff=M^T L M. This applies the numerator exactly
    once on the amplitude and once on its conjugate.
    """
    tensor = np.asarray(L, dtype=complex)
    if tensor.shape != (4, 4) or not np.isfinite(tensor).all():
        raise ValueError("Proca input L must be a finite 4x4 complex tensor")
    mixed = proca_numerator_mixed_matrix(ER_GeV, Q_A_GeV, mV_GeV)
    return mixed.T @ tensor @ mixed


def contract_components(L: np.ndarray, current: dict[str, np.ndarray], JA: float) -> dict[str, float]:
    full = bilinear(L, current["full"], current["full"], JA)
    charge = bilinear(L, current["charge"], current["charge"], JA)
    spatial = bilinear(L, current["spatial"], current["spatial"], JA)
    charge_spatial = 2.0 * bilinear(L, current["charge"], current["spatial"], JA).real
    m_phi = 2.0 * bilinear(L, current["charge_identity"], current["charge_spin"], JA).real
    delta_sigma = 2.0 * bilinear(L, current["delta_LJ"], current["sigma_LJ"], JA).real
    spatial_identity_spin = 2.0 * bilinear(
        L, current["spatial_identity"], current["spatial_spin"], JA
    ).real
    longitudinal = current["full"].copy()
    longitudinal[1:3] = 0.0
    transverse = current["full"].copy()
    transverse[[0, 3]] = 0.0
    long_kernel = bilinear(L, longitudinal, longitudinal, JA)
    transverse_kernel = bilinear(L, transverse, transverse, JA)
    longitudinal_transverse_cross = 2.0 * bilinear(L, longitudinal, transverse, JA).real
    decomposition = charge.real + spatial.real + charge_spatial
    lt_decomposition = long_kernel.real + transverse_kernel.real + longitudinal_transverse_cross
    return {
        "K_full_re": float(full.real),
        "K_full_im": float(full.imag),
        "K_charge_re": float(charge.real),
        "K_charge_im": float(charge.imag),
        "K_spatial_re": float(spatial.real),
        "K_spatial_im": float(spatial.imag),
        "K_charge_spatial_cross": float(charge_spatial),
        "K_M_Phi_signed_cross": float(m_phi),
        "K_Delta_Sigma_LJ_signed_cross": float(delta_sigma),
        "K_spatial_identity_spin_signed_cross": float(spatial_identity_spin),
        "K_longitudinal_0z_re": float(long_kernel.real),
        "K_transverse_xy_re": float(transverse_kernel.real),
        "K_longitudinal_transverse_signed_cross": float(longitudinal_transverse_cross),
        "decomposition_LT_absolute_residual": float(abs(full.real - lt_decomposition)),
        "decomposition_LT_relative_residual": float(abs(full.real - lt_decomposition) / max(abs(full.real), abs(lt_decomposition), 1.0)),
        "decomposition_absolute_residual": float(abs(full.real - decomposition)),
        "decomposition_relative_residual": float(abs(full.real - decomposition) / max(abs(full.real), abs(decomposition), 1.0)),
        "imaginary_residual": float(abs(full.imag)),
        "imaginary_relative_residual": float(abs(full.imag) / max(abs(full.real), 1.0)),
    }



def contract_afh_continuity_projection(
    L: np.ndarray, current: dict[str, np.ndarray], JA: float, alpha: float
) -> dict[str, float]:
    """Separate the AFH six-response projection from the retained full J^3 current."""
    full = np.asarray(current["full"], dtype=complex)
    if L.shape != (4, 4) or full.ndim != 3 or full.shape[0] != 4:
        raise ValueError("invalid Lorentz tensor/current dimensions")
    rho, jz = full[0], full[3]
    d = jz - alpha * rho
    long_actual = np.zeros_like(full)
    long_actual[0], long_actual[3] = rho, jz
    long_projected = np.zeros_like(full)
    long_projected[0], long_projected[3] = rho, alpha * rho
    transverse = np.zeros_like(full)
    transverse[1:3] = full[1:3]
    k0z_full = bilinear(L, long_actual, long_actual, JA)
    k0z_lc = bilinear(L, long_projected, long_projected, JA)
    k_transverse = bilinear(L, transverse, transverse, JA)
    k_lt_projected = 2.0 * bilinear(L, long_projected, transverse, JA).real
    k_lt_actual = 2.0 * bilinear(L, long_actual, transverse, JA).real
    k_full = bilinear(L, full, full, JA)
    residual_terms = (
        -(L[0, 3] - alpha * L[3, 3]) * rho * d.conjugate()
        -(L[3, 0] - alpha * L[3, 3]) * d * rho.conjugate()
        + L[3, 3] * np.abs(d) ** 2    )
    residual_formula = float(np.real(np.sum(residual_terms)) / (2.0 * JA + 1.0))
    residual_direct = float(np.real(k0z_full - k0z_lc))
    formula_relative_residual = abs(residual_direct - residual_formula) / max(
        abs(residual_direct), abs(residual_formula), 1.0
    )
    # At large K_0z the direct difference K_0z-full - K_0z-LC loses digits to
    # cancellation even when the analytic residual formula is exact.  Check the
    # formula against the scale of the two large kernels as well, so that the
    # conditioning is explicit rather than hidden in a small-denominator ratio.
    formula_absolute_residual = abs(residual_direct - residual_formula)
    formula_conditioned_residual = formula_absolute_residual / max(
        abs(float(k0z_full.real)), abs(float(k0z_lc.real)), 1.0
    )
    k_afh6 = float(np.real(k0z_lc + k_transverse))
    k_projected_lt = float(k_lt_projected)
    k_nonconserved_lt = float(k_lt_actual - k_lt_projected)
    k_nonafh_residual = float(k_full.real - k_afh6)
    reconstructed = k_afh6 + k_projected_lt + residual_direct + k_nonconserved_lt
    decomposition_relative_residual = abs(float(k_full.real) - reconstructed) / max(
        abs(float(k_full.real)), abs(reconstructed), 1.0
    )
    return {
        "K_full_raw": float(k_full.real),
        "K_LC_longitudinal_projection_raw": float(k0z_lc.real),
        "K_transverse_xy_projection_raw": float(k_transverse.real),
        "K_AFH6_projection_raw": k_afh6,
        "K_projected_LT_cross_raw": k_projected_lt,
        "K_longitudinal_0z_nonconserved_residual_raw": residual_direct,
        "K_longitudinal_transverse_nonconserved_residual_raw": k_nonconserved_lt,
        "K_nonAFH_projection_residual_raw": k_nonafh_residual,
        "nonconserved_formula_relative_residual": float(formula_relative_residual),
        "nonconserved_formula_conditioned_residual": float(formula_conditioned_residual),
        "nonconserved_formula_absolute_residual": float(formula_absolute_residual),
        "projection_decomposition_relative_residual": float(decomposition_relative_residual),
        "imaginary_relative_residual": float(abs(k_full.imag) / max(abs(k_full.real), 1.0)),
    }


def sum_parts(parts: list[dict[str, float]]) -> dict[str, float]:
    if not parts:
        raise ValueError("cannot sum an empty rank list")
    keys = tuple(parts[0])
    return {
        key: float(max(part[key] for part in parts) if key.endswith("_relative_residual") else sum(part[key] for part in parts))
        for key in keys
    }


def build_afh_reference(row: dict, L: np.ndarray, JA: float, C0: float, Q: float, ER: float) -> dict[str, float]:
    M = cval(row, "AFH_M")
    Phi = cval(row, "AFH_PhiPP")
    Delta = cval(row, "AFH_Delta")
    Sigma = cval(row, "AFH_SigmaPrime")
    F1 = float(row["F1_effective"])
    F2 = float(row["F2_effective"])
    tN = float(row["tN_GeV2"])
    cM = F1 + (F1 + 2.0 * F2) * tN / (8.0 * M_N_GEV * M_N_GEV)
    cPhi = Q * Q * (F1 + 2.0 * F2) / (4.0 * M_N_GEV * M_N_GEV)
    LT = float(np.real(L[1, 1] + L[2, 2]))
    L00 = float(np.real(L[0, 0]))
    LC = float(np.real(L[0, 0] - (ER / Q) * (L[0, 3] + L[3, 0]) + (ER / Q) ** 2 * L[3, 3]))
    CD = F1 * F1 * Q * Q * LT / (2.0 * M_N_GEV * M_N_GEV)
    CS = (F1 + F2) ** 2 * Q * Q * LT / (8.0 * M_N_GEV * M_N_GEV)
    CDS = -F1 * (F1 + F2) * Q * Q * LT / (2.0 * M_N_GEV * M_N_GEV)
    prefactor = C0 * 4.0 * math.pi / (2.0 * JA + 1.0)
    charge_amp = cM * M + cPhi * Phi
    charge_L00 = prefactor * L00 * abs(charge_amp) ** 2
    charge_LC = prefactor * LC * abs(charge_amp) ** 2
    spatial = prefactor * (
        CD * abs(Delta) ** 2 + CS * abs(Sigma) ** 2
        + CDS * (Sigma * Delta.conjugate()).real
    )
    return {
        "K_charge_AFH_L00_only": float(charge_L00),
        "K_charge_AFH_LC": float(charge_LC),
        "K_spatial_AFH": float(spatial),
        "K_total_AFH_L00_plus_spatial": float(charge_L00 + spatial),
        "K_total_AFH_LC_plus_spatial": float(charge_LC + spatial),
        "signed_M_Phi_kernel_AFH_L00_only": float(prefactor * L00 * 2.0 * (cM * M * (cPhi * Phi).conjugate()).real),
        "signed_Delta_Sigma_kernel_AFH": float(prefactor * CDS * (Sigma * Delta.conjugate()).real),
    }


def relative_difference(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def physical_kernel(raw_kernel: float, C0: float) -> float:
    """Physical equal-p/n normalization: AFH factor four and C0 exactly once."""
    return 4.0 * C0 * raw_kernel


def rate_factor(mchi: float, rho: float, mean_atomic_weight: float) -> float:
    target_nuclei_per_kg = AVOGADRO * 1000.0 / mean_atomic_weight
    return (
        (rho / mchi) * C_CM_S * target_nuclei_per_kg * GEV2_TO_CM2
        * 1.0e-6 * 86400.0
    )


def verify_hashes() -> dict[str, str]:
    paths = {
        "c1_tensor_npz": C1_DIR / "full_lorentz_tensor_grid.npz",
        "c1_tensor_index": C1_DIR / "tensor_index.csv",
        "c1_run_manifest": C1_DIR / "run_manifest.json",
        "c1_convergence_report": C1_DIR / "convergence_report.json",
        "stage_b_csv": STAGE_B_CSV,
        "stage_b_validation_report": STAGE_B_VALIDATION_REPORT,
        "stage_a_responses": STAGE_A_RESPONSE_CSV,
        "isotope_input": ISOTOPE_INPUT,
        "b_freeze": B_FREEZE,
        "stage_b_script": STAGE_B_SCRIPT,
        "c1_script": C1_SCRIPT,
        "rate_reference": RATE_REFERENCE,
    }
    actual = {}
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"required C2 input/source missing: {name}: {path}")
        actual[name] = sha256_file(path)
        expected = EXPECTED_SHA256[name]
        if actual[name] != expected:
            raise ValueError(f"{name} SHA-256 mismatch: {actual[name]} != {expected}")
    return actual


def load_c0_map() -> tuple[dict[int, float], int]:
    rows = read_csv(STAGE_A_RESPONSE_CSV)
    if len(rows) != len(ISOTOPES) * RECOIL_COUNT:
        raise ValueError(f"Stage A response row count {len(rows)} != 567")
    values: dict[int, list[float]] = defaultdict(list)
    for row in rows:
        values[int(row["A"])].append(float(row["C0"]))
    if set(values) != set(ISOTOPES):
        raise ValueError(f"Stage A C0 isotope coverage mismatch: {sorted(values)}")
    c0 = {}
    for A, items in values.items():
        if max(items) - min(items) > 2.0e-14:
            raise ValueError(f"C0(A) is not constant across recoil points for Xe-{A}")
        c0[A] = items[0]
    return c0, len(rows)


def load_isotope_inputs() -> tuple[dict[int, dict[str, float]], float, float]:
    rows = read_csv(ISOTOPE_INPUT)
    if not rows:
        raise ValueError("natural-xenon isotope input is empty")
    total_natural_abundance = sum(float(r["natural_atom_fraction"]) for r in rows)
    if abs(total_natural_abundance - 1.0) > 2.0e-6:
        raise ValueError(f"full natural-xenon abundance sum is not unity: {total_natural_abundance}")
    mean_atomic_weight = sum(float(r["atomic_mass_u"]) * float(r["natural_atom_fraction"]) for r in rows) / total_natural_abundance
    data = {}
    for row in rows:
        if row["AFH_default_density_available"].strip().lower() != "true":
            continue
        A = int(row["A"])
        atomic_mass_u = float(row["atomic_mass_u"])
        mA = (atomic_mass_u - Z_XE * M_E_U) * U_GEV
        data[A] = {
            "abundance": float(row["natural_atom_fraction"]),
            "atomic_mass_u": atomic_mass_u,
            "mA_GeV": mA,
        }
    if tuple(sorted(data)) != ISOTOPES:
        raise ValueError(f"natural xenon AFH isotope coverage mismatch: {sorted(data)}")
    coverage = sum(x["abundance"] for x in data.values())
    if abs(coverage - 0.998158) > 5.0e-7:
        raise ValueError(f"natural abundance coverage changed: {coverage}")
    return data, mean_atomic_weight, coverage

def load_c1_archive() -> tuple[dict[str, np.ndarray], list[dict[str, str]], dict, dict]:
    manifest = json.loads((C1_DIR / "run_manifest.json").read_text(encoding="utf-8"))
    convergence = json.loads((C1_DIR / "convergence_report.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS" or manifest.get("production") is not True:
        raise ValueError("C2 requires a production C1 PASS archive")
    if manifest.get("scope_boundary", {}).get("stage_c1_complete") is not True:
        raise ValueError("C1 manifest does not mark stage_c1_complete")
    if manifest.get("scope_boundary", {}).get("full_current_complete") is not False:
        raise ValueError("C1 must not claim full_current_complete")
    if manifest.get("grid", {}).get("isotopes") != list(ISOTOPES):
        raise ValueError("C1 isotope grid mismatch")
    if manifest.get("grid", {}).get("recoil_points") != RECOIL_COUNT:
        raise ValueError("C1 recoil grid mismatch")
    if manifest.get("grid", {}).get("rows") != 27216 or manifest.get("grid", {}).get("open_rows") != 27216:
        raise ValueError("C1 tensor row/open coverage mismatch")
    if manifest.get("grid", {}).get("speed_order") != SPEED_ORDER:
        raise ValueError("C1 speed order mismatch")
    if manifest.get("grid", {}).get("phase_space_order") != [14, 6, 4]:
        raise ValueError("C1 phase-space order mismatch")
    benchmark = manifest.get("benchmark", {})
    for key, expected in (("mchi_GeV", MCHI_GEV), ("ms_GeV", MS_GEV), ("delta_m_GeV", DELTA_M_GEV), ("g2", G2)):
        if not math.isclose(float(benchmark.get(key, math.nan)), expected, rel_tol=0.0, abs_tol=1.0e-14):
            raise ValueError(f"C1 benchmark mismatch for {key}")
    if benchmark.get("width_scheme") != WIDTH_SCHEME or benchmark.get("mapping") != MAPPING:
        raise ValueError("C1 width/mapping convention mismatch")
    if manifest.get("mass_scheme") != "AME2020 nuclear mass=(atomic_mass_u-54*m_e_u)*u_GeV":
        raise ValueError("C1 mass scheme mismatch")
    halo_names = [x.get("name") for x in manifest.get("grid", {}).get("halo_models", [])]
    if halo_names != list(HALOS):
        raise ValueError(f"C1 halo models mismatch: {halo_names}")
    if manifest.get("normalization", {}).get("nucleon_form_factors_included") is not False:
        raise ValueError("C1 must not have Stage B nucleon form factors pre-applied")
    if manifest.get("normalization", {}).get("AFH_C0_applied") is not False:
        raise ValueError("C1 must not have C0 pre-applied")
    if manifest.get("normalization", {}).get("AFH_factor_four_applied") is not False:
        raise ValueError("C1 must not have AFH factor four pre-applied")

    with np.load(C1_DIR / "full_lorentz_tensor_grid.npz", allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    required = {
        "isotope_A", "ER_keV", "halo_model", "speed_index", "speed_kms",
        "speed_weight_kms", "halo_pdf_kms_inv", "beta", "pchi_GeV",
        "phase_space_open", "L_total", "L_before", "L_after", "L_interference",
    }
    missing = required - set(arrays)
    if missing:
        raise ValueError(f"C1 NPZ missing arrays: {sorted(missing)}")
    n = len(arrays["isotope_A"])
    if n != 27216:
        raise ValueError(f"C1 NPZ row count {n} != 27216")
    for key in ("L_total", "L_before", "L_after", "L_interference"):
        if arrays[key].shape != (n, 4, 4) or not np.isfinite(arrays[key]).all():
            raise ValueError(f"C1 {key} must be finite with shape (27216,4,4)")
    if not bool(np.all(arrays["phase_space_open"])):
        raise ValueError("C1 production archive contains closed speed nodes")

    index = read_csv(C1_DIR / "tensor_index.csv")
    if len(index) != n:
        raise ValueError("C1 tensor_index row count does not match NPZ")
    for i, row in enumerate(index):
        if int(row["row_id"]) != i:
            raise ValueError(f"C1 tensor_index row_id discontinuity at {i}")
        comparisons = (
            ("A", int(row["A"]), int(arrays["isotope_A"][i])),
            ("ER_keV", float(row["ER_keV"]), float(arrays["ER_keV"][i])),
            ("speed_index", int(row["speed_index"]), int(arrays["speed_index"][i])),
            ("speed_kms", float(row["speed_kms"]), float(arrays["speed_kms"][i])),
            ("speed_weight_kms", float(row["speed_weight_kms"]), float(arrays["speed_weight_kms"][i])),
            ("halo_pdf_kms_inv", float(row["halo_pdf_kms_inv"]), float(arrays["halo_pdf_kms_inv"][i])),
            ("beta", float(row["beta"]), float(arrays["beta"][i])),
            ("pchi_GeV", float(row["pchi_GeV"]), float(arrays["pchi_GeV"][i])),
        )
        for label, left, right in comparisons:
            if isinstance(left, int):
                if left != right:
                    raise ValueError(f"C1 NPZ/index mismatch {label} at row {i}")
            elif not math.isclose(left, right, rel_tol=2.0e-14, abs_tol=2.0e-14):
                raise ValueError(f"C1 NPZ/index mismatch {label} at row {i}")
        if row["halo_model"] != str(arrays["halo_model"][i]):
            raise ValueError(f"C1 NPZ/index mismatch halo_model at row {i}")
    return arrays, index, manifest, convergence


def run_unit_tests() -> dict:
    if not C2_TEST_SCRIPT.is_file():
        raise FileNotFoundError(f"C2 unit-test script missing: {C2_TEST_SCRIPT}")
    completed = subprocess.run(
        [sys.executable, str(C2_TEST_SCRIPT)], cwd=HERE,
        text=True, capture_output=True, check=False,
    )
    result = {
        "status": "PASS" if completed.returncode == 0 else "FAIL",
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "test_script_sha256": sha256_file(C2_TEST_SCRIPT),
    }
    if completed.returncode != 0:
        raise RuntimeError(f"C2 unit tests failed: {result}")
    return result


def group_stage_b_rows(rows: list[dict[str, str]]) -> dict[tuple[int, float], list[dict[str, str]]]:
    grouped: dict[tuple[int, float], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (int(row["A"]), round(float(row["ER_GeV"]), 9))
        grouped[key].append(row)
    expected = {
        (A, round(float(er), 9))
        for A in ISOTOPES
        for er in np.linspace(40.0e-6, 800.0e-6, RECOIL_COUNT)
    }
    if set(grouped) != expected:
        missing = sorted(expected - set(grouped))[:8]
        extra = sorted(set(grouped) - expected)[:8]
        raise ValueError(f"Stage B point coverage mismatch; missing={missing}, extra={extra}")
    rank_sets = defaultdict(set)
    expected_spin = {128: 0.0, 129: 0.5, 130: 0.0, 131: 1.5, 132: 0.0, 134: 0.0, 136: 0.0}
    for key, point_rows in grouped.items():
        point_rows.sort(key=lambda r: int(r["J"]))
        A = key[0]
        ranks = [int(r["J"]) for r in point_rows]
        rank_sets[A].add(tuple(ranks))
        if 0 not in ranks:
            raise ValueError(f"Stage B K=0,L=1 rank row missing at {key}")
        if len(ranks) != len(set(ranks)):
            raise ValueError(f"duplicate Stage B rank rows at {key}")
        ja_values = {float(r["J_A"]) for r in point_rows}
        if len(ja_values) != 1:
            raise ValueError(f"mixed target spin at {key}")
        if abs(next(iter(ja_values)) - expected_spin[A]) > 1.0e-12:
            raise ValueError(f"Stage B target spin mismatch for Xe-{A}")
        mass_values = {float(r["mA_GeV"]) for r in point_rows}
        q_values = {float(r["Q_A_GeV"]) for r in point_rows}
        tn_values = {float(r["tN_GeV2"]) for r in point_rows}
        f1_values = {float(r["F1_effective"]) for r in point_rows}
        f2_values = {float(r["F2_effective"]) for r in point_rows}
        if len(mass_values) != 1 or len(q_values) != 1 or len(tn_values) != 1:
            raise ValueError(f"mixed Stage B kinematics within rank set at {key}")
        if len(f1_values) != 1 or len(f2_values) != 1:
            raise ValueError(f"mixed form factors within Stage B rank set at {key}")
    for A, sets in rank_sets.items():
        if len(sets) != 1:
            raise ValueError(f"Stage B multipole rank set changes with recoil for Xe-{A}: {sets}")
    return grouped


def run(out_dir: Path) -> dict:
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    timer = time.perf_counter()
    unit_test_result = run_unit_tests()
    input_hashes = verify_hashes()
    freeze = json.loads(B_FREEZE.read_text(encoding="utf-8"))
    if freeze.get("version") != "B-FORMAL-BORAH2020-SUFIAN2017-AME2020-20261124":
        raise ValueError("B freeze version mismatch")
    if not str(freeze.get("status", "")).startswith("FROZEN"):
        raise ValueError("B freeze is not frozen")
    if freeze.get("kinematics", {}).get("t_N_GeV2") != "-2*m_A*E_R":
        raise ValueError("B freeze t_N prescription mismatch")
    if freeze.get("form_factor_inputs", {}).get("raw_to_effective") != "Fhat_i=F_i^(uds)/3 exactly once":
        raise ValueError("B freeze raw-to-effective form-factor convention mismatch")
    if freeze.get("nuclear_response_inputs", {}).get("physical_normalization") != "AFH factor 4 and isotope C0(A) correction applied once; W_M^phys(0)=A^2":
        raise ValueError("B freeze AFH physical normalization mismatch")
    benchmark = freeze.get("benchmark", {})
    expected_benchmark = {
        "m_chi1_GeV": MCHI_GEV, "m_s_keV": MS_GEV * 1.0e6,
        "delta_m_keV": DELTA_M_GEV * 1.0e6, "m_V_GeV": MV_GEV,
        "g2": G2, "coupling_product_squared": COUPLING_PRODUCT_SQUARED,
        "rho_chi_GeV_cm3": RHO_GEV_CM3,
    }
    for key, value in expected_benchmark.items():
        if not math.isclose(float(benchmark.get(key, math.nan)), value, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"B benchmark mismatch for {key}")

    c0_by_A, stage_a_rows = load_c0_map()
    isotope_data, mean_atomic_weight, abundance_coverage = load_isotope_inputs()
    arrays, tensor_index, c1_manifest, c1_convergence = load_c1_archive()
    stage_b_validation = json.loads(STAGE_B_VALIDATION_REPORT.read_text(encoding="utf-8"))
    if stage_b_validation.get("status") != "PASS_P4_FC_STAGE_B_DIRECT_CURRENT_AUDIT":
        raise ValueError("Stage B independent direct-current audit is not PASS")
    stage_b_coverage = stage_b_validation.get("coverage", {})
    if stage_b_coverage.get("isotope_recoil_points") != 567 or stage_b_coverage.get("rank_resolved_rows") != 891:
        raise ValueError("Stage B independent audit coverage mismatch")
    if stage_b_coverage.get("rank0_longitudinal_rows_all_targets") != 567 or stage_b_coverage.get("spin_zero_rank0_longitudinal_rows") != 405:
        raise ValueError("Stage B K=0,L=1 longitudinal-current coverage mismatch")
    stage_b_checks = stage_b_validation.get("checks", {})
    if stage_b_checks.get("max_direct_multipole_fit_residual", math.inf) > 1.0e-10:
        raise ValueError("Stage B direct multipole-fit audit failed")
    if max(stage_b_checks.get("max_leading_NR_vs_Native_AFH_amplitude_residuals", {}).values(), default=math.inf) > 1.0e-10:
        raise ValueError("Stage B leading-NR/native-AFH amplitude audit failed")
    if stage_b_checks.get("max_rank0_longitudinal_alias_residual", math.inf) > 1.0e-12:
        raise ValueError("Stage B K=0,L=1 longitudinal alias audit failed")
    stage_b_rows = read_csv(STAGE_B_CSV)
    stage_b_by_point = group_stage_b_rows(stage_b_rows)
    if len(stage_b_by_point) != 567:
        raise ValueError("Stage B must contain exactly 567 isotope/recoil points")

    # Precompute immutable current matrices once per isotope/recoil/rank.
    current_cache = {}
    for point_key, point_rows in stage_b_by_point.items():
        A, ER_round = point_key
        JA = float(point_rows[0]["J_A"])
        for row in point_rows:
            if abs(float(row["mA_GeV"]) - isotope_data[A]["mA_GeV"]) > 2.0e-12:
                raise ValueError(f"Stage B/AME2020 mass mismatch at {point_key}")
            Q = float(row["Q_A_GeV"])
            ER = float(row["ER_GeV"])
            expected_Q = math.sqrt(2.0 * float(row["mA_GeV"]) * ER + ER * ER)
            if abs(Q - expected_Q) > 2.0e-13:
                raise ValueError(f"Stage B exact Q_A mismatch at {point_key}")
            expected_tA = -2.0 * float(row["mA_GeV"]) * ER
            if abs(float(row["tN_GeV2"]) - expected_tA) > 2.0e-13:
                raise ValueError(f"Stage B t_N != t_A at {point_key}")
            K = int(row["J"])
            current_cache[(A, ER_round, K, "exact")] = build_rank_components(row, "exact")
            current_cache[(A, ER_round, K, "leading_nr")] = build_rank_components(row, "leading_nr")

    # Validate C1/Stage-B kinematic alignment and speed-grid completeness.
    speed_groups = defaultdict(list)
    for i, idx in enumerate(tensor_index):
        A = int(idx["A"])
        ER_keV = float(idx["ER_keV"])
        ER_GeV = ER_keV * 1.0e-6
        point_key = (A, round(ER_GeV, 9))
        if point_key not in stage_b_by_point:
            raise ValueError(f"C1 point has no Stage B current: {point_key}")
        first = stage_b_by_point[point_key][0]
        for label, a, b, tol in (
            ("mA", float(idx["mA_GeV"]), float(first["mA_GeV"]), 2.0e-12),
            ("Q_A", float(idx["Q_A_GeV"]), float(first["Q_A_GeV"]), 2.0e-13),
            ("t_A", float(idx["t_A_GeV2"]), float(first["tN_GeV2"]), 2.0e-13),
            ("t_N", float(idx["t_N_GeV2"]), float(first["tN_GeV2"]), 2.0e-13),
        ):
            if abs(a - b) > tol:
                raise ValueError(f"C1/Stage-B {label} mismatch at {point_key}: {a} vs {b}")
        if idx["phase_space_open"].lower() != "true":
            raise ValueError(f"closed C1 speed row in production index: {i}")
        speed_groups[(A, round(ER_keV, 9), idx["halo_model"])].append(i)
    if len(speed_groups) != 1134 or any(len(v) != SPEED_ORDER for v in speed_groups.values()):
        raise ValueError("C1 speed-node coverage must be 1134 points x 24 nodes")
    for key, indices in speed_groups.items():
        speed_indices = sorted(int(tensor_index[i]["speed_index"]) for i in indices)
        if speed_indices != list(range(SPEED_ORDER)):
            raise ValueError(f"C1 speed indices are incomplete for {key}: {speed_indices}")

    out_dir = out_dir.resolve()
    try:
        out_dir.relative_to(ROOT.resolve())
    except ValueError:
        pass
    else:
        raise ValueError("C2 run output must be outside the repository; pass an external --out-dir to preserve a clean source tree")
    out_dir.mkdir(parents=True, exist_ok=True)
    if any(out_dir.iterdir()):
        raise FileExistsError(f"C2 output directory is not empty: {out_dir}")

    per_speed_rows = []
    isotope_halo_rows = []
    natural_rate_rows = []
    validation_samples = []
    max_herm = 0.0
    max_psd_violation = 0.0
    max_decomp = 0.0
    max_proca_herm = 0.0
    max_proca_psd_violation = 0.0
    max_proca_decomp = 0.0
    max_imag = 0.0
    max_diagram_residual = 0.0
    max_lt_decomposition = 0.0
    max_afh_match_spinful = 0.0
    max_afh_match_spinzero_charge = 0.0
    max_explicit_residual = 0.0
    min_exact_kernel = math.inf
    max_exact_vs_nr = 0.0
    max_rank0_L1_transverse = 0.0
    max_longitudinal_formula_residual = 0.0
    max_longitudinal_formula_relative = 0.0
    max_projection_decomposition_residual = 0.0
    max_projected_lt_cross_relative = 0.0
    max_nonafh_projection_relative = 0.0
    max_longitudinal_nonconserved_relative = 0.0

    # Cache AFH rank data by point and contract the C1 tensors at each speed.
    afh_cache = {}
    for i, idx in enumerate(tensor_index):
        A = int(idx["A"])
        ER_keV = float(idx["ER_keV"])
        ER_GeV = ER_keV * 1.0e-6
        ER_round = round(ER_GeV, 9)
        point_key = (A, ER_round)
        halo = idx["halo_model"]
        point_rows = stage_b_by_point[point_key]
        JA = float(point_rows[0]["J_A"])
        Q = float(idx["Q_A_GeV"])
        mA = float(idx["mA_GeV"])
        tA = float(idx["t_A_GeV2"])
        C0 = c0_by_A[A]

        L = np.asarray(arrays["L_total"][i], dtype=complex)
        L_before = np.asarray(arrays["L_before"][i], dtype=complex)
        L_after = np.asarray(arrays["L_after"][i], dtype=complex)
        L_interference = np.asarray(arrays["L_interference"][i], dtype=complex)
        norm_L = max(float(np.linalg.norm(L)), 1.0e-300)
        herm = float(np.linalg.norm(L - L.conjugate().T) / norm_L)
        eig_min = float(np.linalg.eigvalsh(0.5 * (L + L.conjugate().T)).min())
        psd_violation = max(0.0, -eig_min / norm_L)
        max_herm = max(max_herm, herm)
        max_psd_violation = max(max_psd_violation, psd_violation)
        if herm > 1.0e-10 or psd_violation > 1.0e-10:
            raise ValueError(f"C1 Lorentz tensor fails Hermiticity/PSD at row {i}")
        L_decomp = L_before + L_after + L_interference
        ldecomp = float(np.linalg.norm(L - L_decomp) / norm_L)
        max_decomp = max(max_decomp, ldecomp)
        if ldecomp > 1.0e-12:
            raise ValueError(f"C1 diagram tensor decomposition fails at row {i}")

        # Apply the full Proca numerator to the C1 particle tensor before
        # every exact, leading-NR, AFH-projection, and diagram contraction.
        L_proca = apply_proca_numerator(L, ER_GeV, Q, MV_GEV)
        L_before_proca = apply_proca_numerator(L_before, ER_GeV, Q, MV_GEV)
        L_after_proca = apply_proca_numerator(L_after, ER_GeV, Q, MV_GEV)
        L_interference_proca = apply_proca_numerator(L_interference, ER_GeV, Q, MV_GEV)
        norm_L_proca = max(float(np.linalg.norm(L_proca)), 1.0e-300)
        herm_proca = float(np.linalg.norm(L_proca - L_proca.conjugate().T) / norm_L_proca)
        eig_min_proca = float(np.linalg.eigvalsh(0.5 * (L_proca + L_proca.conjugate().T)).min())
        psd_violation_proca = max(0.0, -eig_min_proca / norm_L_proca)
        max_proca_herm = max(max_proca_herm, herm_proca)
        max_proca_psd_violation = max(max_proca_psd_violation, psd_violation_proca)
        if herm_proca > 1.0e-10 or psd_violation_proca > 1.0e-10:
            raise ValueError(f"Proca-dressed Lorentz tensor fails Hermiticity/PSD at row {i}")
        L_proca_decomp = L_before_proca + L_after_proca + L_interference_proca
        ldecomp_proca = float(np.linalg.norm(L_proca - L_proca_decomp) / norm_L_proca)
        max_proca_decomp = max(max_proca_decomp, ldecomp_proca)
        if ldecomp_proca > 1.0e-12:
            raise ValueError(f"Proca-dressed diagram tensor decomposition fails at row {i}")

        exact_currents = [current_cache[(A, ER_round, int(r["J"]), "exact")] for r in point_rows]
        nr_currents = [current_cache[(A, ER_round, int(r["J"]), "leading_nr")] for r in point_rows]
        exact_total = sum_parts([contract_components(L_proca, cur, JA) for cur in exact_currents])
        nr_total = sum_parts([contract_components(L_proca, cur, JA) for cur in nr_currents])
        alpha = ER_GeV / Q
        exact_afh_projection = sum_parts([
            contract_afh_continuity_projection(L_proca, cur, JA, alpha) for cur in exact_currents
        ])
        nr_afh_projection = sum_parts([
            contract_afh_continuity_projection(L_proca, cur, JA, alpha) for cur in nr_currents
        ])
        max_longitudinal_formula_residual = max(
            max_longitudinal_formula_residual,
            exact_afh_projection["nonconserved_formula_conditioned_residual"],
            nr_afh_projection["nonconserved_formula_conditioned_residual"],
        )
        max_longitudinal_formula_relative = max(
            max_longitudinal_formula_relative,
            exact_afh_projection["nonconserved_formula_relative_residual"],
            nr_afh_projection["nonconserved_formula_relative_residual"],
        )
        max_projection_decomposition_residual = max(
            max_projection_decomposition_residual,
            exact_afh_projection["projection_decomposition_relative_residual"],
            nr_afh_projection["projection_decomposition_relative_residual"],
        )
        projected_lt_fraction = abs(exact_afh_projection["K_projected_LT_cross_raw"]) / max(
            abs(exact_afh_projection["K_AFH6_projection_raw"]), 1.0
        )
        nonafh_fraction = abs(exact_afh_projection["K_nonAFH_projection_residual_raw"]) / max(
            abs(exact_afh_projection["K_full_raw"]), 1.0
        )
        longitudinal_residual_fraction = abs(
            exact_afh_projection["K_longitudinal_0z_nonconserved_residual_raw"]
        ) / max(abs(exact_afh_projection["K_full_raw"]), 1.0)
        max_projected_lt_cross_relative = max(max_projected_lt_cross_relative, projected_lt_fraction)
        max_nonafh_projection_relative = max(max_nonafh_projection_relative, nonafh_fraction)
        max_longitudinal_nonconserved_relative = max(
            max_longitudinal_nonconserved_relative, longitudinal_residual_fraction
        )
        if max_longitudinal_formula_residual > 1.0e-12:
            raise ValueError(f"nonconserved longitudinal residual formula fails at row {i}")
        if max_projection_decomposition_residual > 1.0e-12:
            raise ValueError(f"AFH projection plus retained residuals do not reconstruct full kernel at row {i}")
        rank0_current = current_cache[(A, ER_round, 0, "exact")]
        rank0_parts = contract_components(L_proca, rank0_current, JA)
        rank0_transverse_fraction = abs(rank0_parts["K_transverse_xy_re"]) / max(abs(rank0_parts["K_spatial_re"]), 1.0e-300)
        max_rank0_L1_transverse = max(max_rank0_L1_transverse, rank0_transverse_fraction)
        if rank0_transverse_fraction > 1.0e-10:
            raise ValueError(f"K=0,L=1 current is not purely longitudinal at row {i}")
        max_lt_decomposition = max(max_lt_decomposition, exact_total["decomposition_LT_relative_residual"])
        if exact_total["decomposition_LT_relative_residual"] > 1.0e-12:
            raise ValueError(f"longitudinal/transverse current decomposition fails at row {i}")
        exact_div_num = sum(float(np.linalg.norm(ER_GeV * cur["full"][0] - Q * cur["full"][3]) ** 2) for cur in exact_currents)
        exact_div_den = sum(float(np.linalg.norm(cur["full"]) ** 2) for cur in exact_currents)
        nr_div_num = sum(float(np.linalg.norm(ER_GeV * cur["full"][0] - Q * cur["full"][3]) ** 2) for cur in nr_currents)
        nr_div_den = sum(float(np.linalg.norm(cur["full"]) ** 2) for cur in nr_currents)
        exact_current_divergence = math.sqrt(exact_div_num) / max((ER_GeV + Q) * math.sqrt(exact_div_den), 1.0e-300)
        nr_current_divergence = math.sqrt(nr_div_num) / max((ER_GeV + Q) * math.sqrt(nr_div_den), 1.0e-300)
        before_kernel = sum(
            bilinear(L_before_proca, cur["full"], cur["full"], JA).real for cur in exact_currents
        )
        after_kernel = sum(
            bilinear(L_after_proca, cur["full"], cur["full"], JA).real for cur in exact_currents
        )
        interference_kernel = sum(
            bilinear(L_interference_proca, cur["full"], cur["full"], JA).real for cur in exact_currents
        )
        diagram_residual = relative_difference(
            before_kernel + after_kernel + interference_kernel,
            exact_total["K_full_re"],
        )
        max_diagram_residual = max(max_diagram_residual, diagram_residual)
        if diagram_residual > 1.0e-11:
            raise ValueError(f"contracted before/after/interference do not reconstruct total at row {i}")

        # The AFH six-response comparison uses the continuity-projected current.
        # The actual J^3 remains in the full contraction and its nonconserved
        # residual is retained separately.
        afh_parts = [
            build_afh_reference(r, L_proca, JA, C0, Q, ER_GeV)
            for r in point_rows
        ]
        afh_total = {
            key: sum(part[key] for part in afh_parts)
            for key in afh_parts[0]
        }
        nr_reference = C0 * nr_afh_projection["K_AFH6_projection_raw"]
        afh_reference = afh_total["K_total_AFH_LC_plus_spatial"]
        match = relative_difference(nr_reference, afh_reference)
        if JA == 0.0:
            max_afh_match_spinzero_charge = max(max_afh_match_spinzero_charge, match)
            if match > 1.0e-8:
                raise ValueError(f"spin-zero continuity-projected leading-NR charge fails AFH LC check at row {i}: {match}")
        else:
            max_afh_match_spinful = max(max_afh_match_spinful, match)
            if match > 1.0e-8:
                raise ValueError(f"spinful continuity-projected leading-NR six-response kernel fails AFH check at row {i}: {match}")

        # Explicit nested-loop contraction at deterministic coverage points.
        if int(idx["speed_index"]) == 12 and ER_keV in (201.5, 268.0):
            explicit = sum(
                bilinear_explicit(L_proca, cur["full"], cur["full"], JA)
                for cur in exact_currents
            )
            residual = relative_difference(explicit.real, exact_total["K_full_re"])
            max_explicit_residual = max(max_explicit_residual, residual)
            validation_samples.append({
                "A": A, "ER_keV": ER_keV, "halo_model": halo,
                "speed_index": int(idx["speed_index"]),
                "einsum_kernel": exact_total["K_full_re"],
                "explicit_loop_kernel": float(explicit.real),
                "relative_residual": residual,
            })
            if residual > 1.0e-12 or abs(explicit.imag) > 1.0e-12 * max(abs(explicit.real), 1.0):
                raise ValueError(f"explicit loop contraction mismatch at row {i}")

        if exact_total["imaginary_relative_residual"] > 1.0e-12:
            raise ValueError(f"full exact contraction has an imaginary residual at row {i}")
        max_imag = max(max_imag, exact_total["imaginary_relative_residual"])
        if exact_total["decomposition_relative_residual"] > 1.0e-12:
            raise ValueError(f"full current component decomposition fails at row {i}")
        if exact_total["K_full_re"] < -1.0e-10 * max(abs(exact_total["K_charge_re"]), abs(exact_total["K_spatial_re"]), 1.0):
            raise ValueError(f"negative full-current kernel at row {i}: {exact_total['K_full_re']}")
        min_exact_kernel = min(min_exact_kernel, exact_total["K_full_re"])
        max_exact_vs_nr = max(
            max_exact_vs_nr,
            relative_difference(exact_total["K_full_re"], nr_total["K_full_re"]),
        )

        # AFH physical normalization: raw multipole contraction -> C0 -> x4
        # exactly once. The C0*NR comparison above deliberately remains in
        # the raw AFH normalization; this is the factor-four regression gate.
        physical_factor = 4.0 * C0
        K_phys = physical_kernel(exact_total["K_full_re"], C0)
        K_before_phys = physical_factor * before_kernel
        K_after_phys = physical_factor * after_kernel
        K_interference_phys = physical_factor * interference_kernel
        if relative_difference(K_before_phys + K_after_phys + K_interference_phys, K_phys) > 1.0e-11:
            raise ValueError(f"physical factor-four diagram sum mismatch at row {i}")

        speed_kms = float(idx["speed_kms"])
        speed_weight = float(idx["speed_weight_kms"])
        pdf = float(idx["halo_pdf_kms_inv"])
        beta = float(idx["beta"])
        pchi = float(idx["pchi_GeV"])
        denom = (MV_GEV * MV_GEV - tA) ** 2
        common = mA * COUPLING_PRODUCT_SQUARED / (256.0 * math.pi ** 3 * pchi * pchi * denom)
        halo_weight = speed_weight * pdf * beta * common
        node_rate = halo_weight * K_phys
        node_before_rate = halo_weight * K_before_phys
        node_after_rate = halo_weight * K_after_phys
        node_interference_rate = halo_weight * K_interference_phys
        node_nr_rate = halo_weight * physical_factor * nr_total["K_full_re"]

        per_speed_rows.append({
            "row_id": i, "A": A, "ER_keV": ER_keV, "halo_model": halo,
            "speed_index": int(idx["speed_index"]), "speed_kms": speed_kms,
            "speed_weight_kms": speed_weight, "halo_pdf_kms_inv": pdf,
            "beta": beta, "pchi_GeV": pchi, "mA_GeV": mA,
            "Q_A_GeV": Q, "t_A_GeV2": tA, "t_N_GeV2": float(idx["t_N_GeV2"]),
            "J_A": JA, "C0": C0, "AFH_factor_four": 4.0,
            "K_exact_raw": exact_total["K_full_re"],
            "K_exact_C0_raw_normalization": C0 * exact_total["K_full_re"],
            "K_exact_physical": K_phys,
            "K_charge_physical": physical_factor * exact_total["K_charge_re"],
            "K_spatial_physical": physical_factor * exact_total["K_spatial_re"],
            "K_charge_spatial_cross_physical": physical_factor * exact_total["K_charge_spatial_cross"],
            "K_M_Phi_signed_cross_physical": physical_factor * exact_total["K_M_Phi_signed_cross"],
            "K_Delta_Sigma_LJ_signed_cross_physical": physical_factor * exact_total["K_Delta_Sigma_LJ_signed_cross"],
            "K_spatial_identity_spin_signed_cross_physical": physical_factor * exact_total["K_spatial_identity_spin_signed_cross"],
            "K_longitudinal_0z_physical": physical_factor * exact_total["K_longitudinal_0z_re"],
            "K_transverse_xy_physical": physical_factor * exact_total["K_transverse_xy_re"],
            "K_longitudinal_transverse_signed_cross_physical": physical_factor * exact_total["K_longitudinal_transverse_signed_cross"],
            "K_rank0_L1_spatial_physical": physical_factor * rank0_parts["K_spatial_re"],
            "K_rank0_charge_L1_interference_physical": physical_factor * rank0_parts["K_charge_spatial_cross"],
            "K_rank0_full_physical": physical_factor * rank0_parts["K_full_re"],
            "K_rank0_L1_transverse_fraction": rank0_transverse_fraction,
            "K_before_physical": K_before_phys,
            "K_after_physical": K_after_phys,
            "K_interference_physical": K_interference_phys,
            "K_leading_nr_physical": physical_factor * nr_total["K_full_re"],
            "K_leading_nr_charge_C0_raw": C0 * nr_total["K_charge_re"],
            "K_exact_AFH6_projection_physical": physical_factor * exact_afh_projection["K_AFH6_projection_raw"],
            "K_exact_LC_longitudinal_projection_physical": physical_factor * exact_afh_projection["K_LC_longitudinal_projection_raw"],
            "K_exact_transverse_xy_projection_physical": physical_factor * exact_afh_projection["K_transverse_xy_projection_raw"],
            "K_exact_projected_LT_cross_physical": physical_factor * exact_afh_projection["K_projected_LT_cross_raw"],
            "K_exact_longitudinal_0z_nonconserved_residual_physical": physical_factor * exact_afh_projection["K_longitudinal_0z_nonconserved_residual_raw"],
            "K_exact_longitudinal_transverse_nonconserved_residual_physical": physical_factor * exact_afh_projection["K_longitudinal_transverse_nonconserved_residual_raw"],
            "K_exact_nonAFH_projection_residual_physical": physical_factor * exact_afh_projection["K_nonAFH_projection_residual_raw"],
            "K_leading_nr_AFH6_projection_physical": physical_factor * nr_afh_projection["K_AFH6_projection_raw"],
            "K_leading_nr_nonAFH_projection_residual_physical": physical_factor * nr_afh_projection["K_nonAFH_projection_residual_raw"],
            "AFH_reference_raw": afh_reference,
            "leading_nr_AFH_relative_residual": match,
            "K_decomposition_relative_residual": exact_total["decomposition_relative_residual"],
            "diagram_decomposition_relative_residual": diagram_residual,
            "imaginary_relative_residual": exact_total["imaginary_relative_residual"],
            "longitudinal_transverse_decomposition_relative_residual": exact_total["decomposition_LT_relative_residual"],
            "AFH_projection_decomposition_relative_residual": exact_afh_projection["projection_decomposition_relative_residual"],
            "longitudinal_residual_formula_relative_residual": exact_afh_projection["nonconserved_formula_relative_residual"],
            "projected_LT_cross_over_AFH6": projected_lt_fraction,
            "nonAFH_projection_residual_over_full": nonafh_fraction,
            "longitudinal_0z_nonconserved_residual_over_full": longitudinal_residual_fraction,
            "exact_current_divergence_over_full_current_scale": exact_current_divergence,
            "leading_nr_current_divergence_over_full_current_scale": nr_current_divergence,
            "L_hermiticity_relative_residual": herm,
            "L_psd_negative_eigenvalue_relative": psd_violation,
            "L_diagram_decomposition_relative_residual": ldecomp,
            "L_proca_hermiticity_relative_residual": herm_proca,
            "L_proca_psd_negative_eigenvalue_relative": psd_violation_proca,
            "L_proca_diagram_decomposition_relative_residual": ldecomp_proca,
            "Proca_numerator_applied_once": True,
            "Proca_mV_GeV": MV_GEV,
            "rate_node_integrand": node_rate,
            "rate_node_before_integrand": node_before_rate,
            "rate_node_after_integrand": node_after_rate,
            "rate_node_interference_integrand": node_interference_rate,
            "rate_node_leading_nr_integrand": node_nr_rate,
        })

    if not validation_samples:
        raise ValueError("no explicit-loop validation samples were selected")

    # Integrate the stored speed nodes only; no phase-space or speed nodes are regenerated.
    node_groups = defaultdict(list)
    for row in per_speed_rows:
        key = (int(row["A"]), round(float(row["ER_keV"]), 9), row["halo_model"])
        node_groups[key].append(row)
    isotope_integrals = {}
    for key, rows in sorted(node_groups.items()):
        A, ER_keV, halo = key
        if len(rows) != SPEED_ORDER:
            raise ValueError(f"speed-node count mismatch for {key}: {len(rows)}")
        sums = {
            name: sum(float(r[name]) for r in rows)
            for name in (
                "rate_node_integrand", "rate_node_before_integrand",
                "rate_node_after_integrand", "rate_node_interference_integrand",
                "rate_node_leading_nr_integrand",
            )
        }
        if relative_difference(
            sums["rate_node_before_integrand"] + sums["rate_node_after_integrand"] + sums["rate_node_interference_integrand"],
            sums["rate_node_integrand"],
        ) > 1.0e-11:
            raise ValueError(f"speed-integrated diagram contributions fail for {key}")
        abundance = isotope_data[A]["abundance"]
        isotope_integrals[key] = sums
        isotope_halo_rows.append({
            "A": A, "ER_keV": ER_keV, "halo_model": halo,
            "natural_atom_fraction": abundance,
            "mA_GeV": isotope_data[A]["mA_GeV"],
            "halo_integral_exact_physical": sums["rate_node_integrand"],
            "halo_integral_before_physical": sums["rate_node_before_integrand"],
            "halo_integral_after_physical": sums["rate_node_after_integrand"],
            "halo_integral_interference_physical": sums["rate_node_interference_integrand"],
            "halo_integral_leading_nr_physical": sums["rate_node_leading_nr_integrand"],
            "abundance_weighted_exact_integral": abundance * sums["rate_node_integrand"],
        })

    mean_atomic_weight = float(mean_atomic_weight)
    rate_prefactor = rate_factor(MCHI_GEV, RHO_GEV_CM3, mean_atomic_weight)
    for halo in HALOS:
        for ER_keV in sorted({float(r["ER_keV"]) for r in per_speed_rows}):
            selected = [r for r in isotope_halo_rows if r["halo_model"] == halo and abs(float(r["ER_keV"]) - ER_keV) < 1.0e-9]
            if len(selected) != len(ISOTOPES):
                raise ValueError(f"natural xenon isotope coverage mismatch at {halo}, ER={ER_keV}")
            weighted = {
                name: sum(float(r["natural_atom_fraction"]) * float(r[key]) for r in selected)
                for name, key in (
                    ("exact", "halo_integral_exact_physical"),
                    ("before", "halo_integral_before_physical"),
                    ("after", "halo_integral_after_physical"),
                    ("interference", "halo_integral_interference_physical"),
                    ("leading_nr", "halo_integral_leading_nr_physical"),
                )
            }
            rates = {name: rate_prefactor * value for name, value in weighted.items()}
            if relative_difference(rates["before"] + rates["after"] + rates["interference"], rates["exact"]) > 1.0e-11:
                raise ValueError(f"natural rate diagram reconstruction fails at {halo}, ER={ER_keV}")
            natural_rate_rows.append({
                "ER_keV": ER_keV, "halo_model": halo,
                "rate_full_current_events_kg_day_keV": rates["exact"],
                "rate_before_events_kg_day_keV": rates["before"],
                "rate_after_events_kg_day_keV": rates["after"],
                "rate_signed_interference_events_kg_day_keV": rates["interference"],
                "rate_leading_nr_control_events_kg_day_keV": rates["leading_nr"],
                "natural_atom_fraction_coverage": abundance_coverage,
                "mean_atomic_weight_u": mean_atomic_weight,
                "rate_prefactor": rate_prefactor,
            })

    # Acceptance summary.
    validation = {
        "status": "PASS",
        "C2_unit_tests": unit_test_result["status"],
        "C1_production_archive_reused": True,
        "Stage_A_native_responses_reused_read_only": True,
        "Stage_B_direct_current_reused_read_only": True,
        "Stage_B_independent_nuclear_amplitude_audit_reused": True,
        "Stage_B_independent_audit_status": stage_b_validation["status"],
        "Stage_B_independent_audit_max_direct_multipole_fit_residual": stage_b_checks["max_direct_multipole_fit_residual"],
        "Stage_B_independent_audit_max_leading_NR_vs_Native_AFH_amplitude_residual": max(stage_b_checks["max_leading_NR_vs_Native_AFH_amplitude_residuals"].values()),
        "phase_space_integrals_recomputed": False,
        "speed_nodes_regenerated": False,
        "C1_tensor_rows": len(per_speed_rows),
        "C1_point_groups": len(speed_groups),
        "Stage_B_recoil_points": len(stage_b_by_point),
        "Stage_A_response_rows": stage_a_rows,
        "C0_isotope_count": len(c0_by_A),
        "max_L_hermiticity_relative": max_herm,
        "max_L_psd_negative_eigenvalue_relative": max_psd_violation,
        "max_L_diagram_decomposition_relative": max_decomp,
        "max_Proca_L_hermiticity_relative": max_proca_herm,
        "max_Proca_L_psd_negative_eigenvalue_relative": max_proca_psd_violation,
        "max_Proca_L_diagram_decomposition_relative": max_proca_decomp,
        "Proca_numerator_applied_once": True,
        "Proca_numerator_convention": "q_A^mu=(E_R,0,0,+Q_A), metric (+---), L_Proca=M^T L M",
        "max_contracted_component_decomposition_relative": max(float(r["K_decomposition_relative_residual"]) for r in per_speed_rows),
        "max_contracted_diagram_decomposition_relative": max_diagram_residual,
        "max_longitudinal_transverse_decomposition_relative": max_lt_decomposition,
        "max_longitudinal_residual_formula_relative": max_longitudinal_formula_relative,
        "max_longitudinal_residual_formula_conditioned": max_longitudinal_formula_residual,
        "max_AFH_projection_decomposition_relative": max_projection_decomposition_residual,
        "max_projected_longitudinal_transverse_cross_over_AFH6": max_projected_lt_cross_relative,
        "max_nonAFH_projection_residual_over_full": max_nonafh_projection_relative,
        "max_longitudinal_0z_nonconserved_residual_over_full": max_longitudinal_nonconserved_relative,
        "AFH_six_response_comparison_uses_continuity_projected_current": True,
        "actual_longitudinal_current_retained_in_full_contraction": True,
        "nonconserved_longitudinal_residual_retained_separately": True,
        "max_imaginary_relative_residual": max_imag,
        "max_spinful_leading_NR_to_AFH_relative_difference": max_afh_match_spinful,
        "max_spinzero_leading_NR_charge_to_AFH_relative_difference": max_afh_match_spinzero_charge,
        "max_explicit_loop_vs_einsum_relative_difference": max_explicit_residual,
        "explicit_loop_sample_count": len(validation_samples),
        "max_rank0_L1_transverse_fraction_all_isotopes": max_rank0_L1_transverse,
        "minimum_exact_raw_kernel": min_exact_kernel,
        "maximum_exact_vs_leading_NR_relative_difference": max_exact_vs_nr,
        "AFH_factor_four_applied_once": True,
        "C0_applied_once": True,
        "C0_times_four_applied_to_physical_kernel": True,
        "C0_times_four_not_applied_to_raw_AFH_matching_check": True,
        "independent_nuclear_amplitude_validation": "Stage B independent one-body impulse audit plus C2 explicit Lorentz contraction and all-grid continuity-projected leading-NR/AFH six-response match; actual longitudinal residual retained separately",
        "t_N_equals_t_A": True,
        "mediator_uses_t_A_not_Q_A_squared": True,
        "natural_rate_rows": len(natural_rate_rows),
        "isotope_halo_rows": len(isotope_halo_rows),
        "tolerances": {
            "hermiticity": 1.0e-10, "PSD_negative_eigenvalue_relative": 1.0e-10,
            "Proca_tensor_hermiticity": 1.0e-10, "Proca_tensor_PSD_negative_eigenvalue_relative": 1.0e-10,
            "tensor_diagram_decomposition": 1.0e-12, "Proca_tensor_diagram_decomposition": 1.0e-12,
            "contracted_component_decomposition": 1.0e-12,
            "contracted_diagram_decomposition": 1.0e-11,
            "longitudinal_transverse_decomposition": 1.0e-12,
            "longitudinal_residual_formula_conditioned": 1.0e-12,
            "AFH_projection_decomposition": 1.0e-12,
            "imaginary_relative": 1.0e-12,
            "leading_NR_to_AFH": 1.0e-8,
            "explicit_loop_vs_einsum": 1.0e-12,
            "rank0_L1_transverse_fraction_all_isotopes": 1.0e-10,
        },
    }
    if max_afh_match_spinful > 1.0e-8 or max_afh_match_spinzero_charge > 1.0e-8:
        validation["status"] = "FAIL_AFH_MATCH"
    if (
        max_proca_herm > 1.0e-10
        or max_proca_psd_violation > 1.0e-10
        or max_proca_decomp > 1.0e-12
        or max_explicit_residual > 1.0e-12
        or max_rank0_L1_transverse > 1.0e-10
        or max_longitudinal_formula_residual > 1.0e-12
        or max_projection_decomposition_residual > 1.0e-12
    ):
        validation["status"] = "FAIL_INDEPENDENT_CONTRACTION"
    if validation["status"] != "PASS":
        raise RuntimeError(f"C2 acceptance failed: {validation}")

    write_csv(out_dir / "per_speed_full_current_contractions.csv", per_speed_rows)
    write_csv(out_dir / "isotope_halo_contributions.csv", isotope_halo_rows)
    write_csv(out_dir / "natural_xenon_full_current_rates.csv", natural_rate_rows)
    (out_dir / "validation_samples.json").write_text(json.dumps(validation_samples, indent=2) + "\n", encoding="utf-8")
    (out_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    (out_dir / "unit_test_output.txt").write_text(unit_test_result["stdout"] + unit_test_result["stderr"], encoding="utf-8")

    runtime = {
        "started_utc": started,
        "wall_time_seconds": time.perf_counter() - timer,
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "numpy": np.__version__,
        "sympy": __import__("sympy").__version__,
        "command_line": sys.argv,
        "github_actions_used": False,
        "phase_space_or_speed_integration_recomputed": False,
    }
    try:
        import scipy
        runtime["scipy"] = scipy.__version__
    except Exception:
        runtime["scipy"] = None
    try:
        runtime["gpu_info"] = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            text=True, stderr=subprocess.STDOUT,
        ).strip()
    except Exception:
        runtime["gpu_info"] = None
    try:
        import psutil
        runtime["memory_total_bytes"] = int(psutil.virtual_memory().total)
    except Exception:
        runtime["memory_total_bytes"] = None
    try:
        runtime["git_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        runtime["git_branch"] = subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip()
        runtime["git_status_porcelain"] = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()
    except Exception:
        runtime["git_commit"] = None
        runtime["git_branch"] = None
        runtime["git_status_porcelain"] = None
    if not runtime["git_commit"] or runtime["git_status_porcelain"] is None:
        raise RuntimeError("C2 production run requires a Git checkout")
    if runtime["git_status_porcelain"]:
        raise RuntimeError("C2 production run requires a clean Git working tree")

    output_files = [
        out_dir / "per_speed_full_current_contractions.csv",
        out_dir / "isotope_halo_contributions.csv",
        out_dir / "natural_xenon_full_current_rates.csv",
        out_dir / "validation_samples.json",
        out_dir / "validation.json",
        out_dir / "unit_test_output.txt",
    ]
    output_hashes = {p.name: sha256_file(p) for p in output_files}
    manifest = {
        "schema": "p4-fc-stage-c2-full-current-contraction-v2-proca",
        "status": validation["status"],
        "production": True,
        "repository": "lm846802297-bit/LZ-3-body-for-hdm",
        "parent_stage": "C1",
        "scope": "Full one-body four-current contraction with the Proca numerator g_{mu nu}-q_mu q_nu/mV^2 over the complete seven-isotope, 81-recoil, two-halo C1 speed grid; includes the physical AFH factor four and C0(A) exactly once, then folds stored speed nodes into natural-xenon rates.",
        "benchmark": {
            "mchi_GeV": MCHI_GEV, "ms_GeV": MS_GEV,
            "delta_m_GeV": DELTA_M_GEV, "mV_GeV": MV_GEV,
            "g2": G2, "coupling_product_squared": COUPLING_PRODUCT_SQUARED,
            "rho_chi_GeV_cm3": RHO_GEV_CM3,
        },
        "kinematics": {
            "mass_scheme": "AME2020 exact nuclear mass=(atomic_mass_u-54*m_e_u)*u_GeV",
            "Q_A_squared_GeV2": "2*m_A*E_R+E_R^2",
            "t_A_GeV2": "-2*m_A*E_R",
            "t_N_GeV2": "t_A; formal external-transfer one-body impulse prescription",
            "AFH_argument": "Q_A",
            "form_factor_argument": "t_N",
            "mediator_argument": "t_A",
        },
        "propagator": {
            "metric": "(+---)",
            "q_A_upper_GeV": "(E_R,0,0,+Q_A)",
            "numerator": "g_{mu nu}-q_{A,mu}q_{A,nu}/mV^2",
            "numerator_contraction": "L_Proca=M^T L M with M^mu{}_nu=delta^mu_nu-q_A,nu*q_A^mu/mV^2",
            "numerator_applied_once": True,
            "scalar_denominator": "(mV^2-t_A)^2, applied once in the rate kernel",
        },
        "normalization": {
            "Stage_B_form_factors": "effective Fhat_i=F_i^(uds)/3 already included in Stage B current matrices",
            "DM_tensor_coupling_convention": "unit scalar-emission vertex amplitude; g2 enters the pole width, while the external coupling product carries (g1*g2*gVN)^2",
            "external_coupling_basis": "(g1*g2*gVN)^2 with gVN=3*gVq",
            "C0_applied_once": True,
            "AFH_factor_four_applied_once": True,
            "physical_kernel": "4*C0(A)*K_raw",
            "leading_NR_AFH_comparison": "C0(A)*K_AFH6 is compared to the raw AFH multipole contraction after J^3 is continuity-projected to (E_R/Q_A)J^0; actual J^3 and nonconserved residuals remain in the full exact-current kernel; factor four is not included in this raw matching check",
            "no_additional_nucleon_form_factor_factor": True,
        },
        "grid": {
            "isotopes": list(ISOTOPES), "recoil_points": RECOIL_COUNT,
            "halo_models": list(HALOS), "speed_order": SPEED_ORDER,
            "phase_space_order_reused": [14, 6, 4],
            "tensor_rows": len(per_speed_rows), "isotope_halo_points": len(isotope_halo_rows),
            "natural_rate_rows": len(natural_rate_rows),
        },
        "source_sha256": {
            "C2_script": sha256_file(Path(__file__).resolve()),
            "C1_script": input_hashes["c1_script"],
            "Stage_B_direct_current_script": input_hashes["stage_b_script"],
            "rate_formula_reference": input_hashes["rate_reference"],
        },
        "input_sha256": input_hashes,
        "upstream_C1_source_sha256": c1_manifest.get("source_sha256", {}),
        "upstream_C1_input_sha256": c1_manifest.get("input_sha256", {}),
        "upstream_C1_convergence_report": c1_convergence,
        "Stage_B_independent_validation_report": stage_b_validation,
        "unit_tests": unit_test_result,
        "validation": validation,
        "runtime": runtime,
        "outputs": output_hashes,
        "scope_boundary": {
            "Proca_numerator_applied": True,
            "Proca_numerator_applied_once": True,
            "C1_recomputed": False,
            "Stage_A_recomputed": False,
            "Stage_B_recomputed": False,
            "phase_space_integrals_recomputed": False,
            "speed_nodes_regenerated": False,
            "stage_c2_complete": True,
            "one_body_impulse_current_complete": True,
            "independent_one_body_impulse_amplitude_validation_complete": True,
            "independent_nuclear_amplitude_validation_complete": False,
            "P3_tensor_validation_complete": True,
            "full_current_complete": False,
            "full_current_completion_boundary": "The archived one-body impulse current is contracted completely; generic full-current completion is not claimed without a consistent nuclear Hamiltonian and exchange-current closure.",
            "detector_efficiency_or_exposure_applied": False,
        },    }
    manifest_path = out_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    sums_path = out_dir / "SHA256SUMS.txt"
    all_files = output_files + [manifest_path]
    sums_path.write_text("".join(f"{sha256_file(p)}  {p.name}\n" for p in all_files), encoding="ascii")
    print(json.dumps({
        "status": manifest["status"],
        "Proca_numerator_applied_once": True,
        "tensor_rows": len(per_speed_rows),
        "natural_rate_rows": len(natural_rate_rows),
        "max_spinful_AFH_match": max_afh_match_spinful,
        "max_spinzero_charge_AFH_match": max_afh_match_spinzero_charge,
        "max_explicit_loop_residual": max_explicit_residual,
        "max_diagram_residual": max_diagram_residual,
        "wall_time_seconds": runtime["wall_time_seconds"],
        "output_dir": str(out_dir),
    }, indent=2))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True, help="external empty directory for the native C2 run outputs")
    args = parser.parse_args()
    run(args.out_dir)


if __name__ == "__main__":
    main()
