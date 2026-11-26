#!/usr/bin/env python3
"""Direct seven-isotope HO impulse-current multipoles for P4-FC Stage B.

The bound basis is the AFH v6.0 harmonic-oscillator single-particle basis.
Positive-energy free Dirac spinors are folded with those HO momentum-space
states; the one-body vertex uses the fixed external q_A four-vector.  This is
an explicit impulse embedding, not a unique off-shell nuclear current.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
import platform
import time
import json
import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TOPIC_DIR = ROOT / "research" / "three_body_full_current_nucleon_match"
if str(TOPIC_DIR) not in sys.path:
    sys.path.insert(0, str(TOPIC_DIR))

import numpy as np
from scipy.special import eval_genlaguerre, roots_hermite, sph_harm_y
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j

from nucleon_vector_form_factors import (
    borah_sufian2017, effective_nucleon_form_factors,
)

BASE_RESULT_DIR = TOPIC_DIR / "results" / "2026-11-22" / "p3c_bound_current_amplitude"
OUT_DIR = TOPIC_DIR / "p4_full_current" / "output" / "stage-b-direct-current"
AFH_NATIVE_CSV_AU = BASE_RESULT_DIR / "native" / "afh_v6_p3c_xe129_xe131_native_multipoles.csv"
NATIVE_CSV_P4 = TOPIC_DIR / "p4_full_current" / "output" / "native-seven-isotope-81-point" / "afh_v6_p4_fc_7isotope_81point_native_multipoles.csv"
DENSITY_DIR = TOPIC_DIR / "results" / "2026-11-22" / "native_afh_eight_response" / "inputs" / "density_matrices"
NUCLEAR_MASS_INPUT = ROOT / "research" / "afh_default_density_natural_xenon_halo" / "input" / "natural_xenon_isotope_weights.csv"
AFH_SHA = "dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372"
NUCLEAR_MASS_SHA256 = "fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991"
NATIVE_CSV_P4_SHA256 = "26113db76a969dba98b0fb7b4d6b5389189d2dec33ec680d699741fb0b43a2dc"
DENSITY_SHA256 = {
    128: "fa2017a7ea1951cad6407078aeebede34b9182d0554f38968d5033c6c16bbeb4",
    129: "a0bd848d443a073b80ad731cc87fcfb29d7c899da932d5a230e2e245295f158e",
    130: "38eedadeecc45acc637264b386f8e0d2b63c8dc9e811395e348fbb3e4da8beea",
    131: "4ca110f6828ef6cebe4373fb0c9704220d8227c4a131bc2b5801dc2884e00a00",
    132: "5a367bb5a4d278bf9c3c9f7e2aae8cd7252efe8c3371f1da2d06516c766f07f3",
    134: "843bd0d85d77f904c41b6a8ecb349369cadfec9d3b79f477406208e8dc410a7c",
    136: "9069106f11e2efd537380c62d7e4e63114ed84ec6bbfde687bd5da493d988df0",
}
U_GEV = 0.93149410242
M_E_U = 0.000548579909065
Z_XE = 54
AU_NUCLEAR_MASS_SCHEME = "historical mA=A*u"
AME2020_NUCLEAR_MASS_SCHEME = "P4 exact nuclear mass=(atomic_mass_u-54*m_e_u)*u_GeV"
M_N = 0.938918
F1 = 1.0
F2 = -0.184
RECOIL_KEV = tuple(40.0 + 9.5 * i for i in range(81))
ER_GEV = tuple(x * 1.0e-6 for x in RECOIL_KEV)
ISOTOPES = (128, 129, 130, 131, 132, 134, 136)
TARGET_SPIN = {128: 0.0, 129: 0.5, 130: 0.0, 131: 1.5, 132: 0.0, 134: 0.0, 136: 0.0}
METRIC = np.diag([1.0, -1.0, -1.0, -1.0])
PAULI = np.array([
    [[0, 1], [1, 0]],
    [[0, -1j], [1j, 0]],
    [[1, 0], [0, -1]],
], dtype=complex)
I2 = np.eye(2, dtype=complex)


def _gamma_matrices():
    z = np.zeros((2, 2), dtype=complex)
    g0 = np.block([[I2, z], [z, -I2]])
    spatial = [np.block([[z, s], [-s, z]]) for s in PAULI]
    return np.asarray([g0, *spatial])


GAMMA = _gamma_matrices()
SIGMA = np.empty((4, 4, 4, 4), dtype=complex)
for _mu in range(4):
    for _nu in range(4):
        SIGMA[_mu, _nu] = 0.5j * (GAMMA[_mu] @ GAMMA[_nu] - GAMMA[_nu] @ GAMMA[_mu])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def b_parameter_GeV_inv(A: int) -> float:
    b_fm = math.sqrt(41.467 / (45.0 * A ** (-1.0 / 3.0) - 25.0 * A ** (-2.0 / 3.0)))
    return b_fm * 5.0677


@lru_cache(maxsize=1)
def _nuclear_mass_table() -> dict[int, float]:
    """Use the exact P4 production formula, not the rounded display column."""
    masses: dict[int, float] = {}
    with NUCLEAR_MASS_INPUT.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if row.get("AFH_default_density_available", "").strip().lower() != "true":
                continue
            A = int(row["A"])
            atomic_mass_u = float(row["atomic_mass_u"])
            masses[A] = (atomic_mass_u - Z_XE * M_E_U) * U_GEV
    expected_isotopes = (128, 129, 130, 131, 132, 134, 136)
    if tuple(sorted(masses)) != expected_isotopes:
        raise ValueError(
            f"expected P4 AME2020 isotope set {expected_isotopes}, got {tuple(sorted(masses))}"
        )
    return masses


def nuclear_mass_GeV(A: int, mass_scheme: str = "ame2020") -> float:
    """Return the nuclear mass under an explicitly selected comparison scheme."""
    if mass_scheme == "au":
        return int(A) * U_GEV
    if mass_scheme != "ame2020":
        raise ValueError(f"unknown nuclear mass scheme: {mass_scheme}")
    try:
        return _nuclear_mass_table()[int(A)]
    except KeyError as exc:
        raise ValueError(f"no pinned AME2020 nuclear mass for Xe-{A}") from exc


def recoil_Q_GeV(A: int, ER: float, mass_scheme: str = "ame2020") -> float:
    m_A = nuclear_mass_GeV(A, mass_scheme)
    return math.sqrt(2.0 * m_A * ER + ER * ER)


def _rational_twice(x2: int) -> Rational:
    return Rational(x2, 2)


@lru_cache(maxsize=None)
def _cg(l: int, ml2: int, s2: int, j2: int, m2: int) -> float:
    ml = Rational(ml2, 2)
    ss = Rational(s2, 2)
    jj = Rational(j2, 2)
    mm = Rational(m2, 2)
    return float(clebsch_gordan(l, Rational(1, 2), jj, ml, ss, mm))


@lru_cache(maxsize=None)
def _w3j(jp2: int, K: int, j2: int, mf2: int, M: int, mi2: int) -> float:
    return float(wigner_3j(
        Rational(jp2, 2), K, Rational(j2, 2),
        Rational(-mf2, 2), M, Rational(mi2, 2)))


@lru_cache(maxsize=None)
def _cg_vector(L: int, mu: int, K: int) -> float:
    return float(clebsch_gordan(L, 1, K, 0, mu, mu))


@dataclass(frozen=True)
class Orbit:
    N: int
    j2: int

    @property
    def l(self) -> int:
        # AFH v6.0 Lnumber[N,j]: choose l=j±1/2 so N-l is even.
        return (self.j2 + 1) // 2 if (self.N - (self.j2 + 1) // 2) % 2 == 0 else (self.j2 - 1) // 2

    @property
    def n(self) -> int:
        return (self.N - self.l) // 2

    @property
    def j(self) -> float:
        return self.j2 / 2.0

    def ms2(self) -> tuple[int, ...]:
        return tuple(range(-self.j2, self.j2 + 1, 2))


def _angles(p: np.ndarray):
    r = np.linalg.norm(p, axis=1)
    safe = np.maximum(r, 1e-300)
    theta = np.arccos(np.clip(p[:, 2] / safe, -1.0, 1.0))
    phi = np.arctan2(p[:, 1], p[:, 0])
    return r, theta, phi


def _ho_spin_angular(orbit: Orbit, p_dimless: np.ndarray) -> np.ndarray:
    """Return Fourier HO spin-angular amplitudes, without Gaussian or HO phase."""
    r, theta, phi = _angles(p_dimless)
    out = np.zeros((len(r), len(orbit.ms2()), 2), dtype=complex)
    for im, m2 in enumerate(orbit.ms2()):
        # Dirac/Pauli matrices use the standard spin basis (+1/2, -1/2).
        for spin_index, s2 in enumerate((1, -1)):
            ml2 = m2 - s2
            if ml2 % 2:
                continue
            ml = ml2 // 2
            if abs(ml) > orbit.l:
                continue
            coeff = _cg(orbit.l, ml2, s2, orbit.j2, m2)
            if coeff == 0.0:
                continue
            y = sph_harm_y(orbit.l, ml, theta, phi)
            lag = eval_genlaguerre(orbit.n, orbit.l + 0.5, r * r)
            radial_poly = (
                math.sqrt(2.0 * math.factorial(orbit.n) / math.gamma(orbit.n + orbit.l + 1.5))
                * np.power(r, orbit.l) * lag
            )
            out[:, im, spin_index] += coeff * radial_poly * y
    phase = ((-1) ** orbit.n) * ((-1j) ** orbit.l)
    return phase * out


def _positive_energy_spinor(p: np.ndarray) -> np.ndarray:
    """NR-normalized u(p)/sqrt(2E), shape (node, Dirac, Pauli-spin)."""
    p2 = np.einsum("ni,ni->n", p, p)
    e = np.sqrt(M_N * M_N + p2)
    scale = np.sqrt((e + M_N) / (2.0 * e))
    lower = np.zeros((len(p), 2, 2), dtype=complex)
    for i in range(3):
        lower += (p[:, i] / (e + M_N))[:, None, None] * PAULI[i][None, :, :]
    u = np.zeros((len(p), 4, 2), dtype=complex)
    u[:, :2] = scale[:, None, None] * I2[None]
    u[:, 2:] = scale[:, None, None] * lower
    return u


def exact_covariant_current(p_i: np.ndarray, p_f: np.ndarray, Q: float, q0: float, f1=F1, f2=F2):
    """Return normalized covariant Dirac-Pauli current matrices per quadrature node."""
    u_i = _positive_energy_spinor(p_i)
    u_f = _positive_energy_spinor(p_f)
    q4 = np.asarray([q0, 0.0, 0.0, Q])
    qlow = METRIC @ q4
    out = np.empty((len(p_i), 4, 2, 2), dtype=complex)
    for mu in range(4):
        vertex = np.broadcast_to(f1 * GAMMA[mu], (len(p_i), 4, 4)).copy()
        for nu in range(4):
            vertex += (1j * f2 / (2.0 * M_N) * qlow[nu]) * SIGMA[mu, nu][None]
        ubar = np.swapaxes(u_f.conj(), 1, 2) @ GAMMA[0]
        out[:, mu] = np.einsum("nar,nrt,ntb->nab", ubar, vertex, u_i)
    return out


def leading_nr_current(k: np.ndarray, Q: float, q0: float, f1=F1, f2=F2):
    """Project-matched Pauli reduction for external q_A=(q0,0,0,Q)."""
    out = np.zeros((len(k), 4, 2, 2), dtype=complex)
    t = q0 * q0 - Q * Q
    out[:, 0] = (f1 + (f1 + 2.0 * f2) * t / (8.0 * M_N * M_N)) * I2[None]
    q = np.asarray([0.0, 0.0, Q])
    spin_orbit = np.cross(q[None, :], k)
    for i in range(3):
        out[:, 0] += (1j * (f1 + 2.0 * f2) / (4.0 * M_N * M_N) * spin_orbit[:, i])[:, None, None] * PAULI[i][None]
        out[:, i + 1] = (f1 * k[:, i] / M_N)[:, None, None] * I2[None]
    # i(F1+F2)/(2m) (sigma x q)
    for i in range(3):
        axis = np.eye(3)[i]
        for a in range(3):
            coefficient = np.cross(q, axis)[a]
            out[:, i + 1] += (1j * (f1 + f2) * coefficient / (2.0 * M_N)) * PAULI[a][None]
    return out


def momentum_grid(Q: float, b: float, order: int):
    nodes, weights = roots_hermite(order)
    x, y, z = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    wx, wy, wz = np.meshgrid(weights, weights, weights, indexing="ij")
    k_dimless = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    w = (wx * wy * wz).ravel() * math.exp(-((Q * b / 2.0) ** 2))
    a = Q / 2.0
    p_i = k_dimless / b - np.asarray([0.0, 0.0, a])
    p_f = k_dimless / b + np.asarray([0.0, 0.0, a])
    k = k_dimless / b
    return k_dimless, p_i, p_f, k, w


def _allowed_J(jp2: int, j2: int, lp: int, l: int, *, spatial=False):
    lo = int(abs(jp2 - j2) // 2)
    hi = (jp2 + j2) // 2
    out = []
    for J in range(lo, hi + 1):
        parity = (lp + l) % 2
        if spatial or (J - parity) % 2 == 0:
            out.append(J)
    return out


def _extract_charge_reduced(mat0: np.ndarray, orbit_f: Orbit, orbit_i: Orbit):
    """Fit q-z charge matrix to AFH j_JY_J0 reduced matrix elements."""
    Js = _allowed_J(orbit_f.j2, orbit_i.j2, orbit_f.l, orbit_i.l, spatial=False)
    rows, values = [], []
    fmap = {m2: i for i, m2 in enumerate(orbit_f.ms2())}
    imap = {m2: i for i, m2 in enumerate(orbit_i.ms2())}
    for mi2, ii in imap.items():
        if mi2 not in fmap:
            continue
        mf2 = mi2; ff = fmap[mf2]
        if abs(mf2) > orbit_f.j2:
            continue
        row = []
        for J in Js:
            w3 = _w3j(orbit_f.j2, J, orbit_i.j2, mf2, 0, mi2)
            phase = (-1.0) ** ((orbit_f.j2 - mf2) // 2)
            plane = (1j ** J) * math.sqrt(4.0 * math.pi * (2 * J + 1))
            row.append(phase * w3 * plane)
        if any(abs(x) > 1e-14 for x in row):
            rows.append(row); values.append(mat0[ff, ii])
    if not Js:
        return {}
    design = np.asarray(rows, complex); values = np.asarray(values, complex)
    solution, *_ = np.linalg.lstsq(design, values, rcond=1e-12)
    error_norm = np.linalg.norm(design @ solution - values)
    value_norm = np.linalg.norm(values)
    # Several forbidden spin-charge transitions are exact zeros.  Their
    # quadrature residue is below 1e-12 in the dimensionless current units;
    # treating that residue as a relative error makes a null channel report 1.
    residual = (float(error_norm) if value_norm <= 1e-12 and error_norm <= 1e-12
                else float(error_norm / max(value_norm, 1e-300)))
    return {J: complex(solution[i]) for i, J in enumerate(Js)}, float(residual)


def _extract_vector_reduced(mat_spherical: dict[int, np.ndarray], orbit_f: Orbit, orbit_i: Orbit):
    """Fit q-z vector-current matrix to coupled [j_L Y_L x J_1]_K tensors."""
    Js = _allowed_J(orbit_f.j2, orbit_i.j2, orbit_f.l, orbit_i.l, spatial=True)
    unknowns = []
    for K in Js:
        for L in (K - 1, K, K + 1):
            if L < 0:
                continue
            # Do not impose a parity rule on L alone: the external q vector
            # in the covariant current also carries parity and couples to the
            # spin/current vector.  The complete L=K-1,K,K+1 basis is needed
            # to separate transverse and longitudinal pieces at fixed q.
            unknowns.append((K, L))
    rows, values = [], []
    for mu in (-1, 0, 1):
        mat = mat_spherical[mu]
        fmap = {m2: i for i, m2 in enumerate(orbit_f.ms2())}
        imap = {m2: i for i, m2 in enumerate(orbit_i.ms2())}
        for mi2, ii in imap.items():
            mf2 = mi2 + 2 * mu
            if mf2 not in fmap:
                continue
            ff = fmap[mf2]
            row = []
            for K, L in unknowns:
                w3 = _w3j(orbit_f.j2, K, orbit_i.j2, mf2, mu, mi2)
                phase = (-1.0) ** ((orbit_f.j2 - mf2) // 2)
                plane = 4.0 * math.pi * (1j ** L) * math.sqrt((2 * L + 1) / (4.0 * math.pi))
                row.append(phase * w3 * plane * _cg_vector(L, mu, K))
            if any(abs(x) > 1e-14 for x in row):
                rows.append(row); values.append(mat[ff, ii])
    if not unknowns:
        return {}, 0.0
    design = np.asarray(rows, complex); values = np.asarray(values, complex)
    solution, *_ = np.linalg.lstsq(design, values, rcond=1e-12)
    error_norm = np.linalg.norm(design @ solution - values)
    value_norm = np.linalg.norm(values)
    residual = (float(error_norm) if value_norm <= 1e-12 and error_norm <= 1e-12
                else float(error_norm / max(value_norm, 1e-300)))
    return {key: complex(solution[i]) for i, key in enumerate(unknowns)}, float(residual)


def _read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def _read_density(A: int):
    rows = _read_csv(DENSITY_DIR / f"afh_default_xe{A}_density_matrix.csv")
    return [{
        "coefficient": float(r["coefficient"]), "rank": int(r["rank"]),
        "isospin": int(r["isospin"]), "Nprime": int(r["Nprime"]),
        "two_jprime": int(r["two_jprime"]), "N": int(r["N"]),
        "two_j": int(r["two_j"]),
    } for r in rows]


def prepare_current(A, ER, order, *, current_kind="exact", f1=F1, f2=F2, vertex_energy="external", mass_scheme="ame2020"):
    Q = recoil_Q_GeV(A, ER, mass_scheme); b = b_parameter_GeV_inv(A)
    nodes, weights_1d = roots_hermite(order)
    x, y, z = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    wx, wy, wz = np.meshgrid(weights_1d, weights_1d, weights_1d, indexing="ij")
    kdim = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    weights = (wx * wy * wz).ravel() * math.exp(-((Q * b / 2.0) ** 2))
    shift_dim = np.asarray([0.0, 0.0, Q * b / 2.0])
    pi_dim, pf_dim = kdim - shift_dim, kdim + shift_dim
    p_i, p_f = pi_dim / b, pf_dim / b
    k = kdim / b
    q0_internal = np.sqrt(M_N * M_N + np.einsum("ni,ni->n", p_f, p_f)) - np.sqrt(M_N * M_N + np.einsum("ni,ni->n", p_i, p_i))
    q0_internal_rms = math.sqrt(float(np.average(q0_internal * q0_internal, weights=weights)))
    q0 = ER if vertex_energy == "external" else q0_internal
    if current_kind == "exact":
        # exact_covariant_current accepts one q0; external is constant. The
        # internal q0 alternative is evaluated node by node below.
        if vertex_energy == "external":
            cur = exact_covariant_current(p_i, p_f, Q, ER, f1, f2)
        else:
            # Vectorized node-dependent q0 vertex, preserving q_N=p_f-p_i.
            cur = np.empty((len(p_i), 4, 2, 2), complex)
            u_i, u_f = _positive_energy_spinor(p_i), _positive_energy_spinor(p_f)
            ubar = np.swapaxes(u_f.conj(), 1, 2) @ GAMMA[0]
            for mu in range(4):
                vertex = np.broadcast_to(f1 * GAMMA[mu], (len(p_i), 4, 4)).copy()
                for nu in range(4):
                    qnu = q0_internal if nu == 0 else np.full(len(p_i), -Q if nu == 3 else 0.0)
                    vertex += (1j * f2 / (2.0 * M_N) * qnu)[:, None, None] * SIGMA[mu, nu][None]
                cur[:, mu] = np.einsum("nar,nrt,ntb->nab", ubar, vertex, u_i)
    elif current_kind == "leading_nr":
        cur = leading_nr_current(k, Q, ER, f1, f2)
    elif current_kind == "convection_only":
        cur = leading_nr_current(k, Q, ER, 1.0, -1.0)
    elif current_kind == "magnetization_only":
        cur = leading_nr_current(k, Q, ER, 0.0, 1.0)
    else:
        raise ValueError(current_kind)
    identity_cur = np.zeros_like(cur)
    spin_cur = np.zeros_like(cur)
    for mu in range(4):
        scalar = 0.5 * np.trace(cur[:, mu], axis1=1, axis2=2)
        identity_cur[:, mu] = scalar[:, None, None] * I2[None]
        for a in range(3):
            coeff = 0.5 * np.einsum("ab,nba->n", PAULI[a], cur[:, mu])
            spin_cur[:, mu] += coeff[:, None, None] * PAULI[a][None]
    return {"A": A, "ER": ER, "Q": Q, "b": b, "kdim": kdim, "pi_dim": pi_dim, "pf_dim": pf_dim,
            "weights": weights, "p_i": p_i, "p_f": p_f, "k": k, "q0_internal": q0_internal,
            "q0_internal_rms": q0_internal_rms,
            "current": cur, "identity_current": identity_cur, "spin_current": spin_cur,
            "current_kind": current_kind, "vertex_energy": vertex_energy, "wf_cache": {}}


def direct_current_for_pair(row, prepared, component="full"):
    Q = prepared["Q"]
    ER = prepared["ER"]
    w = prepared["weights"]
    oi = Orbit(row["N"], row["two_j"]); of = Orbit(row["Nprime"], row["two_jprime"])
    cache = prepared["wf_cache"]
    if (oi, "i") not in cache:
        cache[(oi, "i")] = _ho_spin_angular(oi, prepared["pi_dim"])
    if (of, "f") not in cache:
        cache[(of, "f")] = _ho_spin_angular(of, prepared["pf_dim"])
    psi_i = cache[(oi, "i")]; psi_f = cache[(of, "f")]
    cur = prepared[{"full": "current", "identity": "identity_current", "spin": "spin_current"}[component]]
    # Include the Fourier HO phases and the Gaussian already removed by GH.
    phase_i = ((-1) ** oi.n) * ((-1j) ** oi.l)
    phase_f = ((-1) ** of.n) * ((-1j) ** of.l)
    phase = np.conjugate(phase_f) * phase_i
    current_m = np.einsum("n,nfs,nmst,nit->mfi", w, psi_f.conj(), cur, psi_i, optimize=True)
    current_m = phase * current_m
    # spherical spatial current, with Condon-Shortley basis.
    cart = {mu: current_m[mu] for mu in range(1, 4)}
    spherical = {
        0: cart[3],
        1: -(cart[1] + 1j * cart[2]) / math.sqrt(2.0),
        -1: (cart[1] - 1j * cart[2]) / math.sqrt(2.0),
    }
    c_reduced, c_resid = _extract_charge_reduced(current_m[0], of, oi)
    v_reduced, v_resid = _extract_vector_reduced(spherical, of, oi)
    q0_internal = prepared["q0_internal"]
    return c_reduced, v_reduced, {"charge_fit_relative_residual": c_resid, "vector_fit_relative_residual": v_resid,
                                  "qN0_min_GeV": float(np.min(q0_internal)), "qN0_max_GeV": float(np.max(q0_internal)),
                                  "qN0_weighted_rms_GeV": prepared["q0_internal_rms"]}


def run(order: int = 14, vertex_energy: str = "external", mass_scheme: str = "ame2020",
        native_csv: Path | None = None, *, isotopes=None, recoil_energies=None):
    if vertex_energy not in {"external", "internal"}:
        raise ValueError("vertex_energy must be 'external' or 'internal'")
    if mass_scheme not in {"au", "ame2020"}:
        raise ValueError("mass_scheme must be 'au' or 'ame2020'")
    if native_csv is None:
        native_csv = AFH_NATIVE_CSV_AU if mass_scheme == "au" else NATIVE_CSV_P4
    native_csv = Path(native_csv)
    if not native_csv.is_file():
        raise FileNotFoundError(native_csv)
    if sha256(ROOT / "v6" / "dmformfactor-V6.m") != AFH_SHA:
        raise ValueError("AFH source SHA-256 mismatch")
    if mass_scheme == "ame2020" and sha256(NUCLEAR_MASS_INPUT) != NUCLEAR_MASS_SHA256:
        raise ValueError("P4 exact nuclear-mass input SHA-256 mismatch")
    if mass_scheme == "ame2020" and native_csv.resolve() == NATIVE_CSV_P4.resolve():
        if sha256(native_csv) != NATIVE_CSV_P4_SHA256:
            raise ValueError("Stage A Native multipole CSV SHA-256 mismatch")
    for isotope in ISOTOPES:
        density_path = DENSITY_DIR / f"afh_default_xe{isotope}_density_matrix.csv"
        if sha256(density_path) != DENSITY_SHA256[isotope]:
            raise ValueError(f"AFH Xe-{isotope} density-matrix SHA-256 mismatch")
    native_rows = _read_csv(native_csv)
    native = {}
    native_kinematics = {}
    for r in native_rows:
        point = (int(r["A"]), round(float(r["ER_GeV"]), 9))
        native[(point[0], point[1], int(r["J"]), r["channel"])] = complex(float(r["amplitude_real"]), float(r["amplitude_imag"]))
        native_kinematics[point] = {
            "mA_GeV": float(r["mA_GeV"]),
            "Q_A_GeV": float(r["Q_A_GeV"]),
            "t_A_GeV2": float(r["t_A_GeV2"]),
        }
    isotope_set = tuple(ISOTOPES if isotopes is None else isotopes)
    recoil_set = tuple(ER_GEV if recoil_energies is None else recoil_energies)
    if not isotope_set or any(A not in ISOTOPES for A in isotope_set):
        raise ValueError(f"isotopes must be a nonempty subset of {ISOTOPES}")
    if not recoil_set or any(not (40.0e-6 - 1e-15 <= x <= 800.0e-6 + 1e-15) for x in recoil_set):
        raise ValueError("recoil energies must lie in the Stage A 40--800 keV grid range")
    run_started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    run_timer = time.perf_counter()
    results = []
    for A in isotope_set:
        density = [r for r in _read_density(A) if r["isospin"] == 0]
        ranks = sorted({r["rank"] for r in density})
        # AFH SetIsotope sets T_iso=|A/2-Z|; only the density-matrix row's
        # isospin index is zero, so preserve this separate target-isospin factor.
        T_iso = abs(A / 2.0 - 54.0)
        iso_factor = math.sqrt(2.0 / (2.0 * T_iso + 1.0))
        for ER in recoil_set:
            point = (A, round(ER, 9))
            if point not in native_kinematics:
                raise ValueError(f"Stage A Native AFH archive has no point Xe-{A} ER={ER}")
            # Formal B impulse prescription: q_N^mu=q_A^mu, hence
            # t_N=t_A=-2*m_A*E_R. Do not substitute -Q_A^2 here.
            mA_GeV = nuclear_mass_GeV(A, mass_scheme)
            tN_GeV2 = -2.0 * mA_GeV * ER
            QN2_GeV2 = -tN_GeV2
            native_kin = native_kinematics[point]
            if abs(native_kin["mA_GeV"] - mA_GeV) > 2e-12:
                raise ValueError(f"Stage A m_A mismatch at Xe-{A} ER={ER}")
            if abs(native_kin["Q_A_GeV"] - recoil_Q_GeV(A, ER, mass_scheme)) > 2e-13:
                raise ValueError(f"Stage A Q_A mismatch at Xe-{A} ER={ER}")
            if abs(native_kin["t_A_GeV2"] - tN_GeV2) > 2e-13:
                raise ValueError(f"Stage A t_A mismatch at Xe-{A} ER={ER}")
            F1_raw, F2_raw = borah_sufian2017(QN2_GeV2)
            f1, f2 = effective_nucleon_form_factors(F1_raw, F2_raw)
            by_kind = {}
            for kind in ("exact", "leading_nr"):
                prepared = prepare_current(A, ER, order, current_kind=kind, f1=f1, f2=f2,
                                           vertex_energy=vertex_energy if kind == "exact" else "external",
                                           mass_scheme=mass_scheme)
                by_component = {}
                qn_min, qn_max, qn_rms = math.inf, -math.inf, 0.0
                for component in ("identity", "spin"):
                    charge_amp = {J: 0j for J in ranks}
                    vector_amp = {}
                    fit_residuals = []
                    for row in density:
                        c, v, diag = direct_current_for_pair(row, prepared, component=component)
                        rank = row["rank"]
                        charge_amp[rank] = charge_amp.get(rank, 0j) + row["coefficient"] * iso_factor * c.get(rank, 0j)
                        for (K, L), amp in v.items():
                            if K == rank:
                                vector_amp[(rank, L)] = vector_amp.get((rank, L), 0j) + row["coefficient"] * iso_factor * amp
                        fit_residuals.extend([diag["charge_fit_relative_residual"], diag["vector_fit_relative_residual"]])
                        qn_min = min(qn_min, diag["qN0_min_GeV"]); qn_max = max(qn_max, diag["qN0_max_GeV"])
                        qn_rms = max(qn_rms, diag["qN0_weighted_rms_GeV"])
                    by_component[component] = {"charge": charge_amp, "vector": vector_amp, "fit_residual_max": max(fit_residuals)}
                by_kind[kind] = {**by_component, "qN0_min": qn_min, "qN0_max": qn_max, "qN0_rms": qn_rms}
            for J in ranks:
                row = {"A": A, "J_A": TARGET_SPIN[A], "ER_GeV": ER, "ER_keV": ER * 1.0e6,
                       "Q_A_GeV": recoil_Q_GeV(A, ER, mass_scheme), "J": J,
                       "mass_scheme": AME2020_NUCLEAR_MASS_SCHEME if mass_scheme == "ame2020" else AU_NUCLEAR_MASS_SCHEME,
                       "mA_GeV": mA_GeV, "tN_GeV2": tN_GeV2, "QN2_GeV2": QN2_GeV2,
                       "F1_uds_raw": F1_raw, "F2_uds_raw": F2_raw,
                       "F1_effective": f1, "F2_effective": f2,
                       "form_factor_prescription": "Borah2020_EM_plus_Sufian2017_unit_strange_current"}
                for channel in ("M", "PhiPP", "Delta", "SigmaPrime"):
                    row[f"AFH_{channel}"] = native.get((A, round(ER, 9), J, channel), 0j)
                for channel in ("M", "PhiPP"):
                    if (A, round(ER, 9), J, channel) not in native:
                        raise ValueError(f"Missing Native AFH {channel} multipole at Xe-{A}, ER={ER}, J={J}")
                if J > 0:
                    for channel in ("Delta", "SigmaPrime"):
                        if (A, round(ER, 9), J, channel) not in native:
                            raise ValueError(f"Missing Native AFH {channel} multipole at Xe-{A}, ER={ER}, J={J}")
                Q = row["Q_A_GeV"]
                cM = f1 + (f1 + 2 * f2) * tN_GeV2 / (8 * M_N * M_N)
                cPhi = Q * Q * (f1 + 2 * f2) / (4 * M_N * M_N)
                expected = {"M": cM * row["AFH_M"], "PhiPP": cPhi * row["AFH_PhiPP"],
                            "Delta": f1 * row["AFH_Delta"], "SigmaPrime": (f1 + f2) * row["AFH_SigmaPrime"]}
                row["expected_AFH_current_M"] = expected["M"]
                row["expected_AFH_current_PhiPP"] = expected["PhiPP"]
                row["expected_AFH_current_Delta"] = expected["Delta"]
                row["expected_AFH_current_SigmaPrime"] = expected["SigmaPrime"]
                for name in ("exact", "leading_nr"):
                    row[f"{name}_charge_identity"] = by_kind[name]["identity"]["charge"].get(J, 0j)
                    row[f"{name}_charge_spin"] = by_kind[name]["spin"]["charge"].get(J, 0j)
                    row[f"{name}_vector_identity_LJ"] = by_kind[name]["identity"]["vector"].get((J, J), 0j)
                    row[f"{name}_vector_spin_LJ"] = by_kind[name]["spin"]["vector"].get((J, J), 0j)
                    row[f"{name}_vector_identity_LJm1"] = by_kind[name]["identity"]["vector"].get((J, J - 1), 0j)
                    row[f"{name}_vector_identity_LJp1"] = by_kind[name]["identity"]["vector"].get((J, J + 1), 0j)
                    row[f"{name}_vector_spin_LJm1"] = by_kind[name]["spin"]["vector"].get((J, J - 1), 0j)
                    row[f"{name}_vector_spin_LJp1"] = by_kind[name]["spin"]["vector"].get((J, J + 1), 0j)
                # AFH current-channel conversion fixed by the project q_A/q_AFH
                # conventions and confirmed independently with pure convection
                # and pure magnetization controls.
                row["exact_Delta_effective"] = 1j * M_N / Q * row["exact_vector_identity_LJ"]
                row["leading_NR_Delta_effective"] = 1j * M_N / Q * row["leading_nr_vector_identity_LJ"]
                row["exact_SigmaPrime_effective"] = -2j * M_N / Q * row["exact_vector_spin_LJ"]
                row["leading_NR_SigmaPrime_effective"] = -2j * M_N / Q * row["leading_nr_vector_spin_LJ"]
                row["exact_vs_leading_vector_identity_LJm1"] = row["exact_vector_identity_LJm1"] - row["leading_nr_vector_identity_LJm1"]
                row["exact_vs_leading_vector_identity_LJ"] = row["exact_vector_identity_LJ"] - row["leading_nr_vector_identity_LJ"]
                row["exact_vs_leading_vector_identity_LJp1"] = row["exact_vector_identity_LJp1"] - row["leading_nr_vector_identity_LJp1"]
                row["exact_vs_leading_vector_spin_LJm1"] = row["exact_vector_spin_LJm1"] - row["leading_nr_vector_spin_LJm1"]
                row["exact_vs_leading_vector_spin_LJ"] = row["exact_vector_spin_LJ"] - row["leading_nr_vector_spin_LJ"]
                row["exact_vs_leading_vector_spin_LJp1"] = row["exact_vector_spin_LJp1"] - row["leading_nr_vector_spin_LJp1"]
                row["exact_rank0_longitudinal_identity"] = row["exact_vector_identity_LJp1"] if J == 0 else 0j
                row["exact_rank0_longitudinal_spin"] = row["exact_vector_spin_LJp1"] if J == 0 else 0j
                row["leading_nr_rank0_longitudinal_identity"] = row["leading_nr_vector_identity_LJp1"] if J == 0 else 0j
                row["leading_nr_rank0_longitudinal_spin"] = row["leading_nr_vector_spin_LJp1"] if J == 0 else 0j
                row["leading_NR_charge_M_residual"] = row["leading_nr_charge_identity"] - expected["M"]
                row["leading_NR_charge_PhiPP_residual"] = row["leading_nr_charge_spin"] - expected["PhiPP"]
                row["leading_NR_Delta_residual"] = row["leading_NR_Delta_effective"] - expected["Delta"]
                row["leading_NR_SigmaPrime_residual"] = row["leading_NR_SigmaPrime_effective"] - expected["SigmaPrime"]
                row["exact_vs_leading_charge_identity"] = row["exact_charge_identity"] - row["leading_nr_charge_identity"]
                row["exact_vs_leading_charge_spin"] = row["exact_charge_spin"] - row["leading_nr_charge_spin"]
                row["exact_vs_leading_Delta"] = row["exact_Delta_effective"] - row["leading_NR_Delta_effective"]
                row["exact_vs_leading_SigmaPrime"] = row["exact_SigmaPrime_effective"] - row["leading_NR_SigmaPrime_effective"]
                row["signed_charge_interference_AFH"] = (row["AFH_M"] * row["AFH_PhiPP"].conjugate()).real
                row["signed_spatial_interference_AFH"] = (row["AFH_Delta"] * row["AFH_SigmaPrime"].conjugate()).real
                row["signed_spatial_interference_exact"] = (row["exact_Delta_effective"] * row["exact_SigmaPrime_effective"].conjugate()).real
                row["signed_charge_interference_expected_NR"] = 2.0 * (expected["M"] * expected["PhiPP"].conjugate()).real
                row["signed_charge_interference_exact"] = 2.0 * (row["exact_charge_identity"] * row["exact_charge_spin"].conjugate()).real
                row["signed_spatial_interference_expected_NR"] = 2.0 * (expected["Delta"] * expected["SigmaPrime"].conjugate()).real
                row["signed_spatial_interference_exact_weighted"] = 2.0 * (row["exact_Delta_effective"] * row["exact_SigmaPrime_effective"].conjugate()).real
                row["exact_charge_identity_relative_to_expected"] = abs(row["exact_charge_identity"] - expected["M"]) / max(abs(expected["M"]), 1e-300)
                row["exact_charge_spin_relative_to_expected"] = abs(row["exact_charge_spin"] - expected["PhiPP"]) / max(abs(expected["PhiPP"]), 1e-300)
                row["exact_Delta_relative_to_expected"] = abs(row["exact_Delta_effective"] - expected["Delta"]) / max(abs(expected["Delta"]), 1e-300)
                row["exact_SigmaPrime_relative_to_expected"] = abs(row["exact_SigmaPrime_effective"] - expected["SigmaPrime"]) / max(abs(expected["SigmaPrime"]), 1e-300)
                row["max_fit_residual"] = max(by_kind["exact"]["identity"]["fit_residual_max"], by_kind["exact"]["spin"]["fit_residual_max"],
                                               by_kind["leading_nr"]["identity"]["fit_residual_max"], by_kind["leading_nr"]["spin"]["fit_residual_max"])
                row["qN0_internal_quadrature_node_range_GeV"] = [by_kind["exact"]["qN0_min"], by_kind["exact"]["qN0_max"]]
                row["qN0_internal_weighted_rms_GeV"] = by_kind["exact"]["qN0_rms"]
                results.append(row)
    f1_zero_raw, f2_zero_raw = borah_sufian2017(0.0)
    return {"status": "COMPLETED_P4_FC_STAGE_B_DIRECT_CURRENT_RUN", "quadrature_order_per_axis": order,
            "mass_scheme_key": mass_scheme,
            "nuclear_mass_scheme": AME2020_NUCLEAR_MASS_SCHEME if mass_scheme == "ame2020" else AU_NUCLEAR_MASS_SCHEME,
            "nuclear_mass_input_sha256": sha256(NUCLEAR_MASS_INPUT) if mass_scheme == "ame2020" else None,
            "nuclear_mass_constants": {"u_GeV": U_GEV, "m_e_u": M_E_U if mass_scheme == "ame2020" else None, "Z_Xe": Z_XE if mass_scheme == "ame2020" else None},
            "native_AFH_multipoles_sha256": sha256(native_csv),
            "native_AFH_points_used": len({(int(r["A"]), round(float(r["ER_GeV"]), 9)) for r in native_rows}),
            "vertex_energy_prescription": vertex_energy,
            "form_factors": {"prescription": "Borah2020_EM_plus_Sufian2017_unit_strange_current", "strange_model": "sufian2017_unit_current", "raw_zero_transfer": {"F1_uds": f1_zero_raw, "F2_uds": f2_zero_raw}, "effective_zero_transfer": {"F1": f1_zero_raw / 3.0, "F2": f2_zero_raw / 3.0}, "normalization": "raw F1uds,F2uds divided by 3 once for external gVN=3*gVq", "argument": "QN2=-tN with tN=tA=-2*mA*ER", "mN_GeV": M_N},
            "quadrature_rule": "tensor-product Gauss-Hermite with exp(-x^2) weight; shifted HO momenta p_i=k-q_A/2, p_f=k+q_A/2",
            "model": f"AFH HO momentum-space orbitals embedded into positive-energy free Dirac spinors; finite-tN Borah2020 proton/neutron EM plus Sufian2017 unit-strange F1,F2; qN=qA, tN=tA; exact-current vertex q0 uses {vertex_energy} energy prescription; leading-NR control uses external ER",
            "null_fit_amplitude_threshold": 1e-12,
            "isotopes": list(isotope_set), "recoil_energies_GeV": list(recoil_set),
            "runtime": {
                "run_started_utc": run_started_utc,
                "wall_time_seconds": time.perf_counter() - run_timer,
                "python": sys.version,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "processor": platform.processor(),
                "logical_cpu_count": os.cpu_count(),
                "numpy": np.__version__,
                "scipy": __import__("scipy").__version__,
                "sympy": __import__("sympy").__version__,
            },
            "source_sha256": {"AFH": AFH_SHA,
                              "native_multipoles": sha256(native_csv),
                              "density_matrices": {str(A): sha256(DENSITY_DIR / f"afh_default_xe{A}_density_matrix.csv") for A in ISOTOPES},
                              "form_factor_module": sha256(TOPIC_DIR / "nucleon_vector_form_factors.py"),
                              "borah_source": sha256(ROOT / "research" / "nucleon_form_factor_B_baseline" / "borah2020_lowq2.py"),
                              "borah_parameters": sha256(ROOT / "research" / "nucleon_form_factor_B_baseline" / "data" / "borah2020" / "parameters.json"),
                              "borah_covariance": sha256(ROOT / "research" / "nucleon_form_factor_B_baseline" / "data" / "borah2020" / "ELZExpFF.dat"),
                              "nuclear_mass_input": sha256(NUCLEAR_MASS_INPUT) if mass_scheme == "ame2020" else None,
                              "rate_kinematics": sha256(ROOT / "research" / "three_body_calculation_and_validation" / "rate_kinematics.py"),
                              "direct_script": sha256(Path(__file__).resolve())},
            "results": results}


def _complex_fields(row):
    out = {}
    for k, v in row.items():
        if isinstance(v, complex):
            out[k + "_re"] = v.real; out[k + "_im"] = v.imag
        else:
            out[k] = v
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=14)
    parser.add_argument("--vertex-energy", choices=("external", "internal"), default="external")
    parser.add_argument("--mass-scheme", choices=("au", "ame2020"), default="ame2020")
    parser.add_argument("--native-csv", type=Path, default=None)
    parser.add_argument("--isotopes", type=int, nargs="+", default=None)
    parser.add_argument("--recoil-kev", type=float, nargs="+", default=None)
    parser.add_argument("--out", type=Path, default=OUT_DIR / "direct_bound_current_o14.json")
    args = parser.parse_args()
    recoil_energies = None if args.recoil_kev is None else tuple(x * 1.0e-6 for x in args.recoil_kev)
    payload = run(args.order, args.vertex_energy, args.mass_scheme, args.native_csv,
                  isotopes=args.isotopes, recoil_energies=recoil_energies)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    csv_out = args.out.with_suffix(".csv")
    rows = [_complex_fields(x) for x in payload["results"]]
    if not rows:
        raise RuntimeError("direct-current run produced no rows")
    json_payload = dict(payload)
    json_payload["results"] = rows
    args.out.write_text(json.dumps(json_payload, indent=2) + "\n", encoding="utf-8")
    with csv_out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"status": payload["status"], "order": args.order, "row_count": len(rows),
                      "max_fit_residual": max(r["max_fit_residual"] for r in rows),
                      "max_leading_nr_channel_residual": max(
                          abs(r[key]) for r in payload["results"] for key in (
                              "leading_NR_charge_M_residual", "leading_NR_charge_PhiPP_residual",
                              "leading_NR_Delta_residual", "leading_NR_SigmaPrime_residual")),
                      "max_exact_vs_leading_nr_channel_difference": max(
                          abs(r[key]) for r in payload["results"] for key in (
                              "exact_vs_leading_charge_identity", "exact_vs_leading_charge_spin",
                              "exact_vs_leading_Delta", "exact_vs_leading_SigmaPrime"))}, indent=2))


if __name__ == "__main__":
    main()
