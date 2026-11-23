#!/usr/bin/env python3
"""Direct HO impulse-current multipoles for the P3-C AFH comparison.

The bound basis is the AFH v6.0 harmonic-oscillator single-particle basis.
Positive-energy free Dirac spinors are folded with those HO momentum-space
states; the one-body vertex uses the fixed external q_A four-vector.  This is
an explicit impulse embedding, not a unique off-shell nuclear current.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.special import eval_genlaguerre, roots_hermite, sph_harm_y
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j

ROOT = Path(__file__).resolve().parents[6]
OUT_DIR = Path(__file__).resolve().parents[1]
AFH_NATIVE_CSV = OUT_DIR / "native" / "afh_v6_p3c_xe129_xe131_native_multipoles.csv"
DENSITY_DIR = ROOT / "research" / "three_body_full_current_nucleon_match" / "results" / "2026-11-22" / "native_afh_eight_response" / "inputs" / "density_matrices"
AFH_SHA = "dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372"
U_GEV = 0.93149410242
M_N = 0.938918
F1 = 1.0
F2 = -0.184
ER_GEV = (0.0002, 0.000267)
ISOTOPES = (129, 131)
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


def recoil_Q_GeV(A: int, ER: float) -> float:
    m_A = A * U_GEV
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


def prepare_current(A, ER, order, *, current_kind="exact", f1=F1, f2=F2, vertex_energy="external"):
    Q = recoil_Q_GeV(A, ER); b = b_parameter_GeV_inv(A)
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


def run(order: int = 20, vertex_energy: str = "external"):
    if vertex_energy not in {"external", "internal"}:
        raise ValueError("vertex_energy must be 'external' or 'internal'")
    if not AFH_NATIVE_CSV.is_file():
        raise FileNotFoundError(AFH_NATIVE_CSV)
    if sha256(ROOT / "v6" / "dmformfactor-V6.m") != AFH_SHA:
        raise ValueError("AFH source SHA-256 mismatch")
    native_rows = _read_csv(AFH_NATIVE_CSV)
    native = {}
    for r in native_rows:
        native[(int(r["A"]), round(float(r["ER_GeV"]), 9), int(r["J"]), r["channel"])] = complex(float(r["amplitude_real"]), float(r["amplitude_imag"]))
    results = []
    for A in ISOTOPES:
        density = [r for r in _read_density(A) if r["isospin"] == 0]
        ranks = sorted({r["rank"] for r in density})
        # AFH SetIsotope sets T_iso=|A/2-Z|; only the density-matrix row's
        # isospin index is zero, so preserve this separate target-isospin factor.
        T_iso = abs(A / 2.0 - 54.0)
        iso_factor = math.sqrt(2.0 / (2.0 * T_iso + 1.0))
        for ER in ER_GEV:
            by_kind = {}
            for kind in ("exact", "leading_nr", "convection_only", "magnetization_only"):
                prepared = prepare_current(A, ER, order, current_kind=kind,
                                           vertex_energy=vertex_energy if kind == "exact" else "external")
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
                row = {"A": A, "ER_GeV": ER, "Q_A_GeV": recoil_Q_GeV(A, ER), "J": J}
                for channel in ("M", "PhiPP", "Delta", "SigmaPrime"):
                    row[f"AFH_{channel}"] = native.get((A, round(ER, 9), J, channel), 0j)
                Q = row["Q_A_GeV"]
                cM = F1 + (F1 + 2 * F2) * (ER * ER - Q * Q) / (8 * M_N * M_N)
                cPhi = Q * Q * (F1 + 2 * F2) / (4 * M_N * M_N)
                expected = {"M": cM * row["AFH_M"], "PhiPP": cPhi * row["AFH_PhiPP"],
                            "Delta": F1 * row["AFH_Delta"], "SigmaPrime": (F1 + F2) * row["AFH_SigmaPrime"]}
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
    return {"status": "COMPLETED_MODEL_SPECIFIC_DIRECT_COVARIANT_VS_AFH_COMPARISON", "quadrature_order_per_axis": order,
            "vertex_energy_prescription": vertex_energy,
            "form_factors": {"F1": F1, "F2": F2, "mN_GeV": M_N},
            "quadrature_rule": "tensor-product Gauss-Hermite with exp(-x^2) weight; shifted HO momenta p_i=k-q_A/2, p_f=k+q_A/2",
            "model": f"AFH HO momentum-space orbitals embedded into positive-energy free Dirac spinors; spatial shift p_f-p_i=+q_A; exact-current vertex q0 uses {vertex_energy} energy prescription; leading-NR control uses external ER",
            "null_fit_amplitude_threshold": 1e-12,
            "source_sha256": {"AFH": AFH_SHA, "density_Xe129": sha256(DENSITY_DIR / "afh_default_xe129_density_matrix.csv"),
                              "density_Xe131": sha256(DENSITY_DIR / "afh_default_xe131_density_matrix.csv"),
                              "native_multipoles": sha256(AFH_NATIVE_CSV),
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
    parser.add_argument("--order", type=int, default=20)
    parser.add_argument("--vertex-energy", choices=("external", "internal"), default="external")
    parser.add_argument("--out", type=Path, default=OUT_DIR / "direct_bound_current_pilot.json")
    args = parser.parse_args()
    payload = run(args.order, args.vertex_energy)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    csv_out = args.out.with_suffix(".csv")
    rows = [_complex_fields(x) for x in payload["results"]]
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
