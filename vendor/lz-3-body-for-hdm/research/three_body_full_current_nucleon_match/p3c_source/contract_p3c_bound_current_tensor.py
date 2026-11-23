#!/usr/bin/env python3
"""Contract direct bound-current multipoles with Pole-width T3 tensors."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from sympy import Integer, Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j

ROOT = Path(__file__).resolve().parents[6]
OUT_DIR = Path(__file__).resolve().parents[1]
C0_SOURCE = ROOT / "research" / "three_body_full_current_nucleon_match" / "results" / "2026-11-22" / "native_afh_eight_response" / "p3_python_reference" / "metadata.json"
C0 = {}
JA = {129: 0.5, 131: 1.5}
M_N = 0.938918
F1 = 1.0
F2 = -0.184
METRIC = np.asarray([1.0, -1.0, -1.0, -1.0])
CHANNELS = ("M", "PhiPP", "Delta", "SigmaPrime")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def cval(row, name: str) -> complex:
    return complex(float(row[name + "_re"]), float(row[name + "_im"]))


def pair_array(value):
    a = np.asarray(value, dtype=float)
    if a.shape[-1] != 2:
        raise ValueError("expected [real,imag] pairs")
    return a[..., 0] + 1j * a[..., 1]


@lru_cache(maxsize=None)
def w3(j2: int, K: int, mf2: int, mu: int, mi2: int) -> float:
    return float(wigner_3j(
        Rational(j2, 2), Integer(K), Rational(j2, 2),
        Rational(-mf2, 2), Integer(mu), Rational(mi2, 2)))


@lru_cache(maxsize=None)
def cg_vector(L: int, mu: int, K: int) -> float:
    return float(clebsch_gordan(L, 1, K, 0, mu, mu))


def zero_matrix(j2: int):
    n = j2 + 1
    return np.zeros((4, n, n), dtype=complex)


def build_rank_components(row, A: int, representation: str):
    """Rebuild four-current transition matrices for one target tensor rank."""
    if representation not in {"exact", "leading_nr"}:
        raise ValueError(representation)
    K = int(row["J"])
    j2 = 1 if A == 129 else 3
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
            w = w3(j2, K, mf2, 0, mi2)
            plane = (1j ** K) * math.sqrt(4.0 * math.pi * (2 * K + 1))
            for piece, col in (("charge_identity", f"{representation}_charge_identity"),
                               ("charge_spin", f"{representation}_charge_spin")):
                out[piece][0, ff, ii] += phase * w * plane * cval(row, col)
    # Construct the vector current in spherical projections, then convert to
    # Cartesian components for contraction with the q-aligned Lorentz tensor.
    for mi2 in m2_values:
        ii = m_index[mi2]
        sph_parts = {name: np.zeros((3, len(m2_values)), dtype=complex) for name in out if name.startswith("spatial_")}
        for mu in (-1, 0, 1):
            mf2 = mi2 + 2 * mu
            if mf2 not in m_index:
                continue
            ff = m_index[mf2]
            phase = (-1.0) ** ((j2 - mf2) // 2)
            w = w3(j2, K, mf2, mu, mi2)
            if w == 0.0:
                continue
            for L, suffix in ((K - 1, "LJm1"), (K, "LJ"), (K + 1, "LJp1")):
                if L < 0:
                    continue
                plane = 4.0 * math.pi * (1j ** L) * math.sqrt((2 * L + 1) / (4.0 * math.pi))
                basis = phase * w * plane * cg_vector(L, mu, K)
                for kind, comp in (("identity", f"{representation}_vector_identity_{suffix}"),
                                   ("spin", f"{representation}_vector_spin_{suffix}")):
                    sph_name = f"spatial_{kind}_{suffix}"
                    # Per-component storage uses Cartesian slots only after
                    # all three spherical projections have been evaluated.
                    sph_parts.setdefault(sph_name, np.zeros((3, len(m2_values)), complex))
                    sph_parts[sph_name][mu + 1, ff] += basis * cval(row, comp)
        for name, sph in sph_parts.items():
            # Spherical slots are ordered [-1,0,+1]. Condon-Shortley inverse:
            # Jx=(J_-1-J_+1)/sqrt(2), Jy=i(J_-1+J_+1)/sqrt(2), Jz=J_0.
            for ff in range(len(m2_values)):
                out[name][1, ff, ii] = (sph[0, ff] - sph[2, ff]) / math.sqrt(2.0)
                out[name][2, ff, ii] = 1j * (sph[0, ff] + sph[2, ff]) / math.sqrt(2.0)
                out[name][3, ff, ii] = sph[1, ff]

    out["charge"] = out["charge_identity"] + out["charge_spin"]
    out["spatial_identity"] = sum((out[f"spatial_identity_{s}"] for s in ("LJm1", "LJ", "LJp1")), zero_matrix(j2))
    out["spatial_spin"] = sum((out[f"spatial_spin_{s}"] for s in ("LJm1", "LJ", "LJp1")), zero_matrix(j2))
    out["spatial"] = out["spatial_identity"] + out["spatial_spin"]
    out["full"] = out["charge"] + out["spatial"]
    out["delta_LJ"] = out["spatial_identity_LJ"]
    out["sigma_LJ"] = out["spatial_spin_LJ"]
    return out


def bilinear(L: np.ndarray, A: np.ndarray, B: np.ndarray, JA_value: float) -> complex:
    lower_A = METRIC[:, None, None] * A
    lower_B = METRIC[:, None, None] * B
    return np.einsum("mn,mfi,nfi->", L, lower_A, lower_B.conjugate()) / (2.0 * JA_value + 1.0)


def build_afh_reference(direct_row, native, L, JA_value, Q, ER, A):
    def get(channel):
        return native.get((A, round(ER, 9), int(direct_row["J"]), channel), 0j)
    M, Phi, Delta, Sigma = (get(x) for x in CHANNELS)
    t = ER * ER - Q * Q
    cM = F1 + (F1 + 2.0 * F2) * t / (8.0 * M_N * M_N)
    cPhi = Q * Q * (F1 + 2.0 * F2) / (4.0 * M_N * M_N)
    LT = float(np.real(L[1, 1] + L[2, 2]))
    L00 = float(np.real(L[0, 0]))
    LC = float(np.real(L[0, 0] - (ER / Q) * (L[0, 3] + L[3, 0]) + (ER / Q) ** 2 * L[3, 3]))
    CD = F1 * F1 * Q * Q * LT / (2.0 * M_N * M_N)
    CS = (F1 + F2) ** 2 * Q * Q * LT / (8.0 * M_N * M_N)
    CDS = -F1 * (F1 + F2) * Q * Q * LT / (2.0 * M_N * M_N)
    prefactor = C0[A] * 4.0 * math.pi / (2.0 * JA_value + 1.0)
    charge_amp = cM * M + cPhi * Phi
    charge = prefactor * LC * abs(charge_amp) ** 2
    charge_L00 = prefactor * L00 * abs(charge_amp) ** 2
    spatial = prefactor * (CD * abs(Delta) ** 2 + CS * abs(Sigma) ** 2 + CDS * (Sigma * Delta.conjugate()).real)
    return {
        "AFH_M_re": M.real, "AFH_M_im": M.imag,
        "AFH_PhiPP_re": Phi.real, "AFH_PhiPP_im": Phi.imag,
        "AFH_Delta_re": Delta.real, "AFH_Delta_im": Delta.imag,
        "AFH_SigmaPrime_re": Sigma.real, "AFH_SigmaPrime_im": Sigma.imag,
        "AFH_cM_M_re": (cM * M).real, "AFH_cM_M_im": (cM * M).imag,
        "AFH_cPhi_PhiPP_re": (cPhi * Phi).real, "AFH_cPhi_PhiPP_im": (cPhi * Phi).imag,
        "L_C": LC, "L_T": LT,
        "L00": L00, "K_charge_AFH": charge, "K_charge_AFH_L00_only": charge_L00,
        "K_spatial_AFH": spatial,
        "K_total_AFH_four_channel": charge + spatial,
        "K_total_AFH_L00_plus_spatial": charge_L00 + spatial,
        "signed_M_Phi_kernel_AFH": prefactor * LC * 2.0 * (cM * M * (cPhi * Phi).conjugate()).real,
        "signed_M_Phi_kernel_AFH_L00_only": prefactor * L00 * 2.0 * (cM * M * (cPhi * Phi).conjugate()).real,
        "signed_Delta_Sigma_kernel_AFH": prefactor * CDS * (Sigma * Delta.conjugate()).real,
        "normalization_prefactor_C0_4pi_over_2JA1": prefactor,
    }


def run(direct_csv: Path, native_csv: Path, tensor_json: Path, vertex_energy: str = "external"):
    if vertex_energy not in {"external", "internal"}:
        raise ValueError("vertex_energy must be external or internal")
    c0_meta = json.loads(C0_SOURCE.read_text(encoding="utf-8"))
    C0.clear()
    C0.update({int(x["A"]): float(x["common_matrix_normalization"]) for x in c0_meta["isotope_metadata"]})
    direct_rows = read_csv(direct_csv)
    native_rows = read_csv(native_csv)
    native = {
        (int(r["A"]), round(float(r["ER_GeV"]), 9), int(r["J"]), r["channel"]):
            complex(float(r["amplitude_real"]), float(r["amplitude_imag"]))
        for r in native_rows
    }
    tensor_archive = json.loads(tensor_json.read_text(encoding="utf-8"))
    if tensor_archive.get("status") != "PASS_SPIN_RESOLVED_DM_TENSOR" or tensor_archive.get("formal_width_scheme") != "pole":
        raise ValueError("requires accepted Pole-width spin-resolved tensors")
    tensors = {}
    for b in tensor_archive["benchmarks"]:
        if b.get("width_scheme") != "pole":
            raise ValueError("non-Pole tensor row in production archive")
        A = int(b["isotope_A"])
        ER = float(b["ER_GeV"])
        G = pair_array(b["G_q_aligned_real_imag"])
        L = np.einsum("mana->mn", G)
        tensors[(A, round(ER, 9), b["benchmark_base"])] = (b, L)
    grouped = {}
    for r in direct_rows:
        key = (int(r["A"]), round(float(r["ER_GeV"]), 9), int(r["J"]))
        grouped[key] = r

    results = []
    for (A, ERround, K), row in sorted(grouped.items()):
        ER = float(row["ER_GeV"])
        Q = float(row["Q_A_GeV"])
        JA_value = JA[A]
        rank_data = {rep: build_rank_components(row, A, rep) for rep in ("exact", "leading_nr")}
        matching_tensors = sorted((name, data) for (a, e, name), data in tensors.items() if a == A and e == ERround)
        if not matching_tensors:
            raise KeyError(f"no Pole-width tensors for Xe-{A}, ER={ER}")
        for bname, (tensor, L) in matching_tensors:
            afh = build_afh_reference(row, native, L, JA_value, Q, ER, A)
            out = {
                "A": A, "ER_GeV": ER, "ER_keV": ER * 1e6, "Q_A_GeV": Q,
                "J_A": JA_value, "J": K, "benchmark": tensor["benchmark"],
                "benchmark_base": bname, "width_scheme": "pole",
                "quadrature": {k: tensor[k] for k in ("nw", "nz", "nphi", "quadrature_points", "mapping")},
                "C0": C0[A], "DM_L_hermiticity_residual": float(np.linalg.norm(L - L.conjugate().T) / max(np.linalg.norm(L), 1e-300)),
                **afh,
            }
            for rep, currents in rank_data.items():
                parts = {name: bilinear(L, currents[name], currents[name], JA_value) for name in (
                    "full", "charge", "spatial", "charge_identity", "charge_spin",
                    "spatial_identity", "spatial_spin", "delta_LJ", "sigma_LJ")}
                longitudinal_current = currents["full"].copy()
                longitudinal_current[1:3] = 0.0
                transverse_current = currents["full"].copy()
                transverse_current[[0, 3]] = 0.0
                K_longitudinal = bilinear(L, longitudinal_current, longitudinal_current, JA_value)
                K_transverse = bilinear(L, transverse_current, transverse_current, JA_value)
                total = parts["full"]
                cross_charge_space = 2.0 * bilinear(L, currents["charge"], currents["spatial"], JA_value).real
                cross_M_Phi = 2.0 * bilinear(L, currents["charge_identity"], currents["charge_spin"], JA_value).real
                cross_Delta_Sigma = 2.0 * bilinear(L, currents["delta_LJ"], currents["sigma_LJ"], JA_value).real
                expected_decomposition = parts["charge"].real + parts["spatial"].real + cross_charge_space
                out[f"K_{rep}_matrix_raw"] = float(total.real)
                out[f"K_{rep}_longitudinal_0z_matrix_raw"] = float(K_longitudinal.real)
                out[f"K_{rep}_transverse_xy_matrix_raw"] = float(K_transverse.real)
                out[f"K_{rep}_charge_matrix_raw"] = float(parts["charge"].real)
                out[f"K_{rep}_spatial_matrix_raw"] = float(parts["spatial"].real)
                out[f"K_{rep}_charge_spatial_cross_matrix_raw"] = float(cross_charge_space)
                out[f"K_{rep}_M_Phi_cross_matrix_raw"] = float(cross_M_Phi)
                out[f"K_{rep}_Delta_Sigma_LJ_cross_matrix_raw"] = float(cross_Delta_Sigma)
                out[f"K_{rep}_matrix_decomposition_absolute_residual"] = float(abs(total.real - expected_decomposition))
                out[f"K_{rep}_matrix_imaginary_residual"] = float(abs(total.imag))
                out[f"K_{rep}_matrix_physical_candidate_4C0"] = float(4.0 * C0[A] * total.real)
                out[f"K_{rep}_matrix_physical_candidate_C0"] = float(C0[A] * total.real)
                out[f"K_{rep}_longitudinal_0z_physical_C0"] = float(C0[A] * K_longitudinal.real)
                out[f"K_{rep}_transverse_xy_physical_C0"] = float(C0[A] * K_transverse.real)
                out[f"{rep}_longitudinal_0z_vs_AFH_charge_relative_difference"] = float(
                    abs(C0[A] * K_longitudinal.real - afh["K_charge_AFH"])
                    / max(abs(afh["K_charge_AFH"]), 1e-300)
                ) if afh["K_charge_AFH"] != 0.0 else 0.0
                out[f"rank_kernel_contribution_{rep}"] = float(total.real)
            out["leading_nr_vs_AFH_four_channel_raw_ratio"] = (
                out["K_leading_nr_matrix_raw"] / out["K_total_AFH_four_channel"]
                if out["K_total_AFH_four_channel"] != 0.0 else None)
            out["exact_vs_leading_nr_kernel_relative_difference"] = abs(
                out["K_exact_matrix_raw"] - out["K_leading_nr_matrix_raw"]
            ) / max(abs(out["K_leading_nr_matrix_raw"]), 1e-300)
            out["leading_nr_vs_AFH_four_channel_raw_relative_difference"] = abs(
                out["K_leading_nr_matrix_raw"] - out["K_total_AFH_four_channel"]
            ) / max(abs(out["K_total_AFH_four_channel"]), 1e-300)
            out["leading_nr_vs_AFH_four_channel_relative_difference"] = abs(
                out["K_leading_nr_matrix_physical_candidate_C0"] - out["K_total_AFH_four_channel"]
            ) / max(abs(out["K_total_AFH_four_channel"]), 1e-300)
            out["leading_nr_vs_AFH_L00_plus_spatial_relative_difference"] = abs(
                out["K_leading_nr_matrix_physical_candidate_C0"] - out["K_total_AFH_L00_plus_spatial"]
            ) / max(abs(out["K_total_AFH_L00_plus_spatial"]), 1e-300)
            out["leading_nr_charge_vs_AFH_L00_relative_difference"] = abs(
                C0[A] * out["K_leading_nr_charge_matrix_raw"] - out["K_charge_AFH_L00_only"]
            ) / max(abs(out["K_charge_AFH_L00_only"]), 1e-300) if out["K_charge_AFH_L00_only"] != 0.0 else 0.0
            out["leading_nr_M_Phi_cross_vs_AFH_L00_relative_difference"] = abs(
                C0[A] * out["K_leading_nr_M_Phi_cross_matrix_raw"] - out["signed_M_Phi_kernel_AFH_L00_only"]
            ) / max(abs(out["signed_M_Phi_kernel_AFH_L00_only"]), 1e-300) if out["signed_M_Phi_kernel_AFH_L00_only"] != 0.0 else 0.0
            out["leading_nr_Delta_Sigma_cross_vs_AFH_relative_difference"] = abs(
                C0[A] * out["K_leading_nr_Delta_Sigma_LJ_cross_matrix_raw"] - out["signed_Delta_Sigma_kernel_AFH"]
            ) / max(abs(out["signed_Delta_Sigma_kernel_AFH"]), 1e-300) if out["signed_Delta_Sigma_kernel_AFH"] != 0.0 else 0.0
            # Exact-current nuclear-current divergence under external q_A.
            for rep, currents in rank_data.items():
                div = ER * currents["full"][0] - Q * currents["full"][3]
                den = max((ER + Q) * np.linalg.norm(currents["full"]), 1e-300)
                out[f"{rep}_nuclear_current_divergence_over_full_current_scale"] = float(np.linalg.norm(div) / den)
            results.append(out)

    max_leading_nr_amplitude_residual = 0.0
    for r in direct_rows:
        for direct_name, afh_name in (
            ("leading_nr_charge_identity", "expected_AFH_current_M"),
            ("leading_nr_charge_spin", "expected_AFH_current_PhiPP"),
            ("leading_NR_Delta_effective", "expected_AFH_current_Delta"),
            ("leading_NR_SigmaPrime_effective", "expected_AFH_current_SigmaPrime"),
        ):
            max_leading_nr_amplitude_residual = max(
                max_leading_nr_amplitude_residual,
                abs(cval(r, direct_name) - cval(r, afh_name)),
            )
    max_fit_residual = max(float(r["max_fit_residual"]) for r in direct_rows)
    max_NR_AFH_l00_kernel = max(x["leading_nr_vs_AFH_L00_plus_spatial_relative_difference"] for x in results)
    max_mp_cross_residual = max(x["leading_nr_M_Phi_cross_vs_AFH_L00_relative_difference"] for x in results)
    max_ds_cross_residual = max(x["leading_nr_Delta_Sigma_cross_vs_AFH_relative_difference"] for x in results)
    max_kernel_decomposition_residual = max(
        x[f"K_{rep}_matrix_decomposition_absolute_residual"]
        / max(abs(x[f"K_{rep}_matrix_raw"]), 1e-300)
        for x in results for rep in ("exact", "leading_nr")
    )
    max_kernel_imaginary_residual = max(
        x[f"K_{rep}_matrix_imaginary_residual"]
        / max(abs(x[f"K_{rep}_matrix_raw"]), 1e-300)
        for x in results for rep in ("exact", "leading_nr")
    )
    max_kernel_scale = max(abs(x["K_exact_matrix_raw"]) for x in results)
    min_exact_kernel = min(x["K_exact_matrix_raw"] for x in results)
    acceptance = {
        "leading_nr_vs_native_AFH_max_absolute_amplitude_residual": max_leading_nr_amplitude_residual,
        "leading_nr_amplitude_absolute_tolerance": 1e-10,
        "max_direct_multipole_fit_residual": max_fit_residual,
        "multipole_fit_tolerance": 1e-10,
        "leading_nr_vs_AFH_L00_plus_spatial_max_relative_kernel_difference": max_NR_AFH_l00_kernel,
        "leading_nr_kernel_tolerance": 1e-8,
        "signed_M_Phi_interference_max_relative_residual": max_mp_cross_residual,
        "signed_Delta_Sigma_interference_max_relative_residual": max_ds_cross_residual,
        "signed_interference_tolerance": 1e-8,
        "kernel_decomposition_max_relative_residual": max_kernel_decomposition_residual,
        "kernel_decomposition_tolerance": 1e-10,
        "kernel_imaginary_max_relative_residual": max_kernel_imaginary_residual,
        "kernel_imaginary_tolerance": 1e-10,
        "minimum_exact_kernel_over_max_kernel_scale": min_exact_kernel / max(max_kernel_scale, 1e-300),
        "kernel_nonnegativity_tolerance": -1e-10,
    }
    acceptance["pass"] = (
        max_leading_nr_amplitude_residual < acceptance["leading_nr_amplitude_absolute_tolerance"]
        and max_fit_residual < acceptance["multipole_fit_tolerance"]
        and max_NR_AFH_l00_kernel < acceptance["leading_nr_kernel_tolerance"]
        and max_mp_cross_residual < acceptance["signed_interference_tolerance"]
        and max_ds_cross_residual < acceptance["signed_interference_tolerance"]
        and max_kernel_decomposition_residual < acceptance["kernel_decomposition_tolerance"]
        and max_kernel_imaginary_residual < acceptance["kernel_imaginary_tolerance"]
        and acceptance["minimum_exact_kernel_over_max_kernel_scale"] > acceptance["kernel_nonnegativity_tolerance"]
    )
    payload = {
        "status": "PASS_DIRECT_BOUND_CURRENT_T3_CONTRACTION_WITH_MODEL_SCOPE" if acceptance["pass"] else "FAIL_DIRECT_BOUND_CURRENT_T3_CONTRACTION",
        "acceptance": acceptance,
        "direct_model": f"positive-energy free Dirac spinors folded with AFH v6 HO momentum orbitals; one-body impulse current; exact-current q_N^0 uses {vertex_energy} energy prescription; p_f-p_i=+q_A; F1=1, F2=-0.184",
        "vertex_energy_prescription": vertex_energy,
        "formal_width_scheme": "pole",
        "spin_average": "unpolarized DM Pauli-basis trace from G; target average 1/(2J_A+1)",
        "lorentz_contraction": "G-traced L^{mu nu} contracted with lower-index nuclear current using metric (+---)",
        "normalization_audit": "records both C0 and 4*C0 candidates; the AFH response formula includes a 1/4 in the package response which cancels the physical equal-p/n factor four",
        "benchmark_count": len(results),
        "input_sha256": {"direct_multipoles": sha256(direct_csv), "native_AFH_multipoles": sha256(native_csv), "pole_T3_tensors": sha256(tensor_json), "C0_calibration_metadata": sha256(C0_SOURCE), "contraction_script": sha256(Path(__file__).resolve())},
        "results": results,
    }
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct-csv", type=Path, default=OUT_DIR / "direct_bound_current_o22.csv")
    parser.add_argument("--native-csv", type=Path, default=OUT_DIR / "native" / "afh_v6_p3c_xe129_xe131_native_multipoles.csv")
    parser.add_argument("--tensor-json", type=Path, default=OUT_DIR / "pole_width_t3_tensor_xe129_xe131.json")
    parser.add_argument("--vertex-energy", choices=("external", "internal"), default="external")
    parser.add_argument("--out", type=Path, default=OUT_DIR / "direct_bound_current_t3_contraction.json")
    args = parser.parse_args()
    result = run(args.direct_csv, args.native_csv, args.tensor_json, args.vertex_energy)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"], "benchmark_count": result["benchmark_count"],
        "acceptance": result["acceptance"],
        "max_NR_vs_AFH_four_channel_relative_difference": max(x["leading_nr_vs_AFH_four_channel_relative_difference"] for x in result["results"]),
        "max_NR_vs_AFH_L00_plus_spatial_relative_difference": max(x["leading_nr_vs_AFH_L00_plus_spatial_relative_difference"] for x in result["results"]),
        "max_signed_M_Phi_interference_relative_residual": result["acceptance"]["signed_M_Phi_interference_max_relative_residual"],
        "max_signed_Delta_Sigma_interference_relative_residual": result["acceptance"]["signed_Delta_Sigma_interference_max_relative_residual"],
        "max_exact_vs_NR_kernel_relative_difference": max(x["exact_vs_leading_nr_kernel_relative_difference"] for x in result["results"]),
        "output": str(args.out),
    }, indent=2))


if __name__ == "__main__":
    main()
