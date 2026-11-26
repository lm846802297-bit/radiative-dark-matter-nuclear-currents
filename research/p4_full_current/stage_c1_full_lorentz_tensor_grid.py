#!/usr/bin/env python3
"""P4-FC Stage C1: archive the full q-aligned L^{mu nu}(v,E_R) grid.

This stage computes only the missing three-body Lorentz tensors.  It reuses
the accepted Stage A Native AFH and Stage B direct-current archives, and does
not regenerate nuclear multipoles, Stage B currents, scalar halo folds, or
historical phase-space scans.  The output is an uncontracted tensor archive;
Stage C2 performs the nuclear-current contraction and physical rate fold.
"""
from __future__ import annotations

import argparse
import concurrent.futures
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
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[3]
TOPIC = ROOT / "research" / "three_body_full_current_nucleon_match"
HERE = Path(__file__).resolve().parent
if str(TOPIC) not in sys.path:
    sys.path.insert(0, str(TOPIC))

import native_charge_sector_halo_fold as halo  # noqa: E402
from three_body_phase_space_tensor import integrate_spin_resolved_tensor  # noqa: E402

ISOTOPES = (128, 129, 130, 131, 132, 134, 136)
RECOIL_KEV = tuple(40.0 + 9.5 * i for i in range(81))
HALOS = (
    ("shm_cutoff", halo.shm_truncated_pdf, 776.0),
    ("mb_untruncated", halo.mb_untruncated_pdf, 1600.0),
)
MCHI_GEV = 1000.0
MS_GEV = 0.000200
DELTA_M_GEV = 0.001000
G2 = 1.0
WIDTH_SCHEME = "pole"
MAPPING = "threshold_u2"
NV = 24
# Revised C1 production quadrature.  The original frozen (8,6,4) trial is
# archived separately as a FAIL_CONVERGENCE preflight; the u-direction order
# was the limiting quadrature for the untruncated MB halo anchors.
NW, NZ, NPHI = 14, 6, 4
PHASE_COMPARISON = (16, 14, 12)
SPEED_CONVERGENCE_TOL = 1.0e-5
PHASE_SPACE_TENSOR_TOL = 1.0e-6
DIRAC_TRACE_TOL = 1.0e-8
HERMITICITY_TOL = 1.0e-10
PSD_RELATIVE_TOL = 1.0e-10
PARTIAL_TRACE_TOL = 1.0e-12
DECOMPOSITION_TOL = 1.0e-12

MASS_INPUT = ROOT / "research" / "afh_default_density_natural_xenon_halo" / "input" / "natural_xenon_isotope_weights.csv"
STAGE_A_DIR = HERE / "output" / "native-seven-isotope-81-point"
STAGE_A_MULTIPOLES = STAGE_A_DIR / "afh_v6_p4_fc_7isotope_81point_native_multipoles.csv"
STAGE_A_RESPONSES = STAGE_A_DIR / "afh_v6_p4_fc_7isotope_81point_native_responses.csv"
STAGE_B_DIR = HERE / "output" / "stage-b-direct-current"
STAGE_B_CSV = STAGE_B_DIR / "direct_bound_current_o14.csv"
STAGE_B_SCRIPT = HERE / "direct_bound_current_seven_isotope.py"
B_FREEZE = ROOT / "research" / "nucleon_form_factor_B_baseline" / "B_DATA_VERSION_FREEZE_2026-11-24.json"

EXPECTED_INPUT_SHA256 = {
    "mass_input": "fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991",
    "stage_a_multipoles": "26113db76a969dba98b0fb7b4d6b5389189d2dec33ec680d699741fb0b43a2dc",
    "stage_a_responses": "39a5b8c2401eb5f4e5acabb0a4f38b9a9d1fa7103d12d208462b3fcefef57c48",
    "stage_b_csv": "0b335b39cd288beda1367fd977e311d595c2a88c2c5e7dfcaa9456cf8e26d355",
    "stage_b_script": "78058780862acb06d5c3fd86957dead6d16d9a8e3041c9b4bdc3d5450189b75e",
}

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def complex_from_pairs(value: Any, shape: tuple[int, ...]) -> np.ndarray:
    raw = np.asarray(value, dtype=float)
    if raw.shape != shape + (2,) or not np.isfinite(raw).all():
        raise ValueError(f"expected finite complex-pair array {shape + (2,)}, got {raw.shape}")
    return raw[..., 0] + 1j * raw[..., 1]


def tensor_diagnostics(result: dict[str, Any]) -> dict[str, Any]:
    L = complex_from_pairs(result["L_q_aligned_real_imag"], (4, 4))
    Lb = complex_from_pairs(result["L_before_real_imag"], (4, 4))
    La = complex_from_pairs(result["L_after_real_imag"], (4, 4))
    Li = complex_from_pairs(result["L_interference_real_imag"], (4, 4))
    G = complex_from_pairs(result["G_q_aligned_real_imag"], (4, 4, 4, 4))
    norm = max(float(np.linalg.norm(L)), 1.0e-300)
    herm = float(np.linalg.norm(L - L.conj().T) / norm)
    eig = np.linalg.eigvalsh(0.5 * (L + L.conj().T))
    psd_relative_min = float(np.min(eig) / norm)
    L_from_G = np.einsum("mana->mn", G)
    partial_trace = float(np.linalg.norm(L - L_from_G) / norm)
    decomposition = float(np.linalg.norm(L - Lb - La - Li) / norm)
    return {
        "L": L, "L_before": Lb, "L_after": La, "L_interference": Li,
        "hermiticity_relative": herm,
        "psd_min_eigenvalue_relative": psd_relative_min,
        "gram_partial_trace_relative": partial_trace,
        "decomposition_relative": decomposition,
        "rotation_orthogonality": float(result["q_frame_rotation_max_orthogonality_residual"]),
        "rotation_determinant": float(result["q_frame_rotation_max_determinant_residual"]),
        "quadrature_points": int(result["quadrature_points"]),
    }


def make_speed_nodes(vmin: float, vmax: float, nv: int, pdf) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not (math.isfinite(vmin) and math.isfinite(vmax)) or vmin < 0.0 or vmax <= vmin:
        raise ValueError(f"invalid speed interval [{vmin}, {vmax}]")
    x, w = leggauss(nv)
    speeds = 0.5 * (vmax + vmin) + 0.5 * (vmax - vmin) * x
    weights = 0.5 * (vmax - vmin) * w
    densities = np.asarray([pdf(float(v)) for v in speeds], dtype=float)
    if not np.isfinite(densities).all() or np.any(densities < 0.0):
        raise ValueError("halo PDF returned invalid speed density")
    return speeds, weights, densities


def run_one_task(task: tuple) -> dict[str, Any]:
    (A, ER_keV, mA, halo_name, speed_index, speed, speed_weight, pdf_value,
     nw, nz, nphi) = task
    beta = float(speed) / halo.C_KM_S
    pchi = MCHI_GEV * beta / math.sqrt(1.0 - beta * beta)
    transfer = halo.recoil_transfer_kinematics(float(mA), float(ER_keV))
    base = {
        "A": int(A), "ER_keV": float(ER_keV), "halo_model": str(halo_name),
        "speed_index": int(speed_index), "speed_kms": float(speed),
        "speed_weight_kms": float(speed_weight), "halo_pdf_kms_inv": float(pdf_value),
        "beta": beta, "pchi_GeV": pchi, "mA_GeV": float(mA),
        "Q_A_GeV": float(transfer["Q_A_GeV"]), "t_A_GeV2": float(transfer["tA_GeV2"]),
        "t_N_GeV2": float(transfer["tN_GeV2"]),
    }
    if pdf_value == 0.0:
        return {**base, "phase_space_open": False, "skip_reason": "zero_halo_pdf",
                "L": np.zeros((4, 4), complex), "L_before": np.zeros((4, 4), complex),
                "L_after": np.zeros((4, 4), complex), "L_interference": np.zeros((4, 4), complex),
                "hermiticity_relative": 0.0, "psd_min_eigenvalue_relative": 0.0,
                "gram_partial_trace_relative": 0.0, "decomposition_relative": 0.0,
                "rotation_orthogonality": 0.0, "rotation_determinant": 0.0,
                "quadrature_points": 0}
    result = integrate_spin_resolved_tensor(
        mchi=MCHI_GEV, ms=MS_GEV, delta_m=DELTA_M_GEV, mA=float(mA),
        ER=float(ER_keV) * 1.0e-6, speed_kms=float(speed), g2=G2,
        width_scheme=WIDTH_SCHEME, mapping=MAPPING,
        nw=int(nw), nz=int(nz), nphi=int(nphi), crosscheck_dirac_trace=False,
    )
    if not result.get("phase_space_open", False):
        return {**base, "phase_space_open": False, "skip_reason": "phase_space_closed",
                "L": np.zeros((4, 4), complex), "L_before": np.zeros((4, 4), complex),
                "L_after": np.zeros((4, 4), complex), "L_interference": np.zeros((4, 4), complex),
                "hermiticity_relative": 0.0, "psd_min_eigenvalue_relative": 0.0,
                "gram_partial_trace_relative": 0.0, "decomposition_relative": 0.0,
                "rotation_orthogonality": float(result.get("q_frame_rotation_max_orthogonality_residual", 0.0)),
                "rotation_determinant": float(result.get("q_frame_rotation_max_determinant_residual", 0.0)),
                "quadrature_points": 0}
    diag = tensor_diagnostics(result)
    return {**base, "phase_space_open": True, "skip_reason": "",
            **diag}


def build_tasks(isotopes, recoil_energies, isotope_data, nv=NV, nw=NW, nz=NZ, nphi=NPHI):
    tasks = []
    coverage = []
    for halo_name, pdf, vmax in HALOS:
        for A in isotopes:
            mA = float(isotope_data[int(A)]["mA_GeV"])
            for ER_keV in recoil_energies:
                vmin = halo.vmin_three_body_km_s(MCHI_GEV, MS_GEV, mA, float(ER_keV))
                point = {
                    "A": int(A), "ER_keV": float(ER_keV), "halo_model": halo_name,
                    "mA_GeV": mA, "vmin_kms": float(vmin), "vmax_kms": float(vmax),
                    "open_speed_interval": bool(vmin < vmax), "speed_node_count": 0,
                }
                if vmin < vmax:
                    speeds, weights, densities = make_speed_nodes(vmin, vmax, nv, pdf)
                    point["speed_node_count"] = int(nv)
                    for j, (v, w, f) in enumerate(zip(speeds, weights, densities)):
                        tasks.append((int(A), float(ER_keV), mA, halo_name, j,
                                      float(v), float(w), float(f), nw, nz, nphi))
                coverage.append(point)
    return tasks, coverage


def speed_integral(rows: list[dict[str, Any]]) -> np.ndarray:
    total = np.zeros((4, 4), dtype=complex)
    for row in rows:
        if not row["phase_space_open"]:
            continue
        vweight = (row["speed_weight_kms"] * row["halo_pdf_kms_inv"]
                   * row["beta"] / (row["pchi_GeV"] ** 2))
        total += vweight * row["L"]
    return total


def relative_matrix_difference(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    scale = max(float(np.linalg.norm(a)), float(np.linalg.norm(b)), 1.0e-300)
    return float(np.linalg.norm(a - b) / scale), float(np.max(np.abs(a - b)) / scale)


def run_convergence_checks(records, isotope_data, recoil_energies, enabled=True, nv_low=NV):
    if not enabled:
        return {"status": "NOT_RUN_PILOT_ONLY", "speed_checks": [], "phase_space_checks": []}
    anchor_energies = (201.5, 268.0)
    speed_checks = []
    phase_checks = []
    dirac_trace_checks = []
    lookup: dict[tuple, list[dict[str, Any]]] = {}
    for row in records:
        lookup.setdefault((row["A"], row["ER_keV"], row["halo_model"]), []).append(row)
    for A in (129, 131):
        for ER in anchor_energies:
            if ER not in recoil_energies:
                continue
            mA = float(isotope_data[A]["mA_GeV"])
            vmin = halo.vmin_three_body_km_s(MCHI_GEV, MS_GEV, mA, ER)
            for halo_name, pdf, vmax in HALOS:
                base_rows = sorted(lookup.get((A, ER, halo_name), []), key=lambda x: x["speed_index"])
                if len(base_rows) != nv_low:
                    raise ValueError(f"missing C1 production speed rows for convergence anchor {(A, ER, halo_name)}")
                base_integral = speed_integral(base_rows)
                speeds48, weights48, densities48 = make_speed_nodes(vmin, vmax, 48, pdf)
                high_rows = []
                for j, (v, w, f) in enumerate(zip(speeds48, weights48, densities48)):
                    high_rows.append(run_one_task((A, ER, mA, halo_name, j, float(v), float(w),
                                                   float(f), NW, NZ, NPHI)))
                high_integral = speed_integral(high_rows)
                rel, max_component = relative_matrix_difference(base_integral, high_integral)
                speed_checks.append({
                    "A": A, "ER_keV": ER, "halo_model": halo_name,
                    "nv_low": nv_low, "nv_high": 48, "relative_frobenius": rel,
                    "max_component_relative_to_tensor_norm": max_component,
                })
                median = base_rows[len(base_rows) // 2]
                high_phase = run_one_task((A, ER, mA, halo_name, median["speed_index"],
                                           median["speed_kms"], median["speed_weight_kms"],
                                           median["halo_pdf_kms_inv"],
                                           *PHASE_COMPARISON))
                if not median["phase_space_open"] or not high_phase["phase_space_open"]:
                    phase_rel = 0.0 if (not median["phase_space_open"] and not high_phase["phase_space_open"]) else math.inf
                    phase_component = phase_rel
                else:
                    phase_rel, phase_component = relative_matrix_difference(median["L"], high_phase["L"])
                phase_checks.append({
                    "A": A, "ER_keV": ER, "halo_model": halo_name,
                    "speed_kms": median["speed_kms"], "production_order": [NW, NZ, NPHI],
                    "comparison_order": list(PHASE_COMPARISON), "relative_frobenius": phase_rel,
                    "max_component_relative_to_tensor_norm": phase_component,
                })
                trace_result = integrate_spin_resolved_tensor(
                    mchi=MCHI_GEV, ms=MS_GEV, delta_m=DELTA_M_GEV, mA=mA,
                    ER=ER * 1.0e-6, speed_kms=median["speed_kms"], g2=G2,
                    width_scheme=WIDTH_SCHEME, mapping=MAPPING,
                    nw=NW, nz=NZ, nphi=NPHI, crosscheck_dirac_trace=True,
                )
                trace_crosscheck = trace_result.get("dirac_trace_crosscheck") or {}
                trace_residual = float(trace_crosscheck.get("max_relative_residual", math.inf))
                dirac_trace_checks.append({
                    "A": A, "ER_keV": ER, "halo_model": halo_name,
                    "speed_kms": median["speed_kms"],
                    "max_relative_residual": trace_residual,
                    "residuals": trace_crosscheck.get("relative_residuals", {}),
                })
    max_speed = max((x["relative_frobenius"] for x in speed_checks), default=0.0)
    max_phase = max((x["relative_frobenius"] for x in phase_checks), default=0.0)
    max_trace = max((x["max_relative_residual"] for x in dirac_trace_checks), default=0.0)
    return {
        "status": "PASS" if (max_speed <= SPEED_CONVERGENCE_TOL
                             and max_phase <= PHASE_SPACE_TENSOR_TOL
                             and max_trace <= DIRAC_TRACE_TOL) else "FAIL",
        "speed_checks": speed_checks, "phase_space_checks": phase_checks,
        "dirac_trace_checks": dirac_trace_checks,
        "max_speed_relative_frobenius": max_speed,
        "speed_tolerance": SPEED_CONVERGENCE_TOL,
        "max_phase_space_relative_frobenius": max_phase,
        "phase_space_tolerance": PHASE_SPACE_TENSOR_TOL,
        "max_independent_dirac_trace_residual": max_trace,
        "dirac_trace_tolerance": DIRAC_TRACE_TOL,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def total_memory_bytes() -> int | None:
    try:
        import psutil
        return int(psutil.virtual_memory().total)
    except Exception:
        if hasattr(os, "sysconf"):
            try:
                return int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"))
            except (ValueError, OSError, AttributeError):
                return None
        return None


def git_state() -> dict[str, Any]:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        branch = subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip()
        dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()
        return {"commit": commit, "branch": branch, "working_tree_dirty": bool(dirty), "status_porcelain": dirty}
    except Exception as exc:
        return {"commit": None, "branch": None, "working_tree_dirty": None, "git_error": str(exc)}


def parse_float_list(values, default):
    if not values:
        return tuple(default)
    return tuple(float(x) for x in values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=HERE / "output" / "stage-c1-full-lorentz-grid")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--isotopes", nargs="*", type=int, default=list(ISOTOPES))
    parser.add_argument("--recoil-kev", nargs="*", type=float, default=list(RECOIL_KEV))
    parser.add_argument("--nv", type=int, default=NV)
    parser.add_argument("--nw", type=int, default=NW)
    parser.add_argument("--nz", type=int, default=NZ)
    parser.add_argument("--nphi", type=int, default=NPHI)
    parser.add_argument("--pilot", action="store_true", help="mark reduced isotope/recoil runs as non-production")
    parser.add_argument("--skip-convergence", action="store_true", help="allowed only with --pilot")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or min(args.nv, args.nw, args.nz, args.nphi) < 2:
        parser.error("workers and all quadrature orders must be positive (orders >= 2)")
    if args.skip_convergence and not args.pilot:
        parser.error("--skip-convergence is permitted only for a non-production --pilot run")
    isotopes = tuple(args.isotopes)
    recoil_energies = tuple(args.recoil_kev)
    if not set(isotopes).issubset(ISOTOPES) or not isotopes or not recoil_energies:
        parser.error("isotopes must be a non-empty subset of the seven AFH isotopes and recoil grid must be non-empty")
    if not args.pilot and (isotopes != ISOTOPES or recoil_energies != RECOIL_KEV
                           or (args.nv, args.nw, args.nz, args.nphi) != (NV, NW, NZ, NPHI)):
        parser.error("production mode requires the complete frozen 7-isotope x 81-point grid and (24,14,6,4) quadrature")
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        parser.error(f"output directory is not empty: {args.out_dir}; use a new path or --overwrite")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    started = dt.datetime.now(dt.timezone.utc).isoformat()
    timer = time.perf_counter()
    required = {
        "mass_input": MASS_INPUT, "stage_a_multipoles": STAGE_A_MULTIPOLES,
        "stage_a_responses": STAGE_A_RESPONSES, "stage_b_csv": STAGE_B_CSV,
        "stage_b_script": STAGE_B_SCRIPT, "B_freeze": B_FREEZE,
    }
    for name, path in required.items():
        if not path.is_file():
            raise FileNotFoundError(f"required frozen input is missing: {name}: {path}")
    actual_hashes = {name: sha256_file(path) for name, path in required.items()}
    for name, expected in EXPECTED_INPUT_SHA256.items():
        if actual_hashes[name] != expected:
            raise ValueError(f"{name} SHA-256 mismatch: {actual_hashes[name]} != {expected}")
    isotope_data = halo.load_isotope_inputs(MASS_INPUT, mass_scheme="ame2020")
    if tuple(sorted(isotope_data)) != ISOTOPES:
        raise ValueError("AME2020 isotope input does not contain the frozen seven-isotope set")
    freeze = json.loads(B_FREEZE.read_text(encoding="utf-8"))
    if freeze.get("version") != "B-FORMAL-BORAH2020-SUFIAN2017-AME2020-20261124":
        raise ValueError("B data freeze version does not match the adopted production baseline")
    if not str(freeze.get("status", "")).startswith("FROZEN"):
        raise ValueError("B data freeze manifest is not in a frozen state")
    if "pole" not in str(freeze.get("benchmark", {}).get("width", "")).lower():
        raise ValueError("C1 requires the frozen asymmetric Pole width prescription")
    frozen_benchmark = freeze.get("benchmark", {})
    expected_benchmark = {
        "m_chi1_GeV": MCHI_GEV,
        "m_s_keV": MS_GEV * 1.0e6,
        "delta_m_keV": DELTA_M_GEV * 1.0e6,
        "g2": G2,
    }
    for key, expected in expected_benchmark.items():
        if not math.isclose(float(frozen_benchmark.get(key, math.nan)), expected,
                            rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"B freeze benchmark mismatch for {key}")

    tasks, coverage = build_tasks(isotopes, recoil_energies, isotope_data,
                                  nv=args.nv, nw=args.nw, nz=args.nz, nphi=args.nphi)
    if not tasks:
        raise ValueError("C1 grid contains no open speed nodes")
    print(f"C1: {len(tasks)} tensor nodes; workers={args.workers}", flush=True)
    if args.workers == 1:
        records = [run_one_task(task) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
            records = list(pool.map(run_one_task, tasks, chunksize=4))

    max_herm = max((x["hermiticity_relative"] for x in records), default=0.0)
    max_psd_violation = max((max(0.0, -x["psd_min_eigenvalue_relative"]) for x in records), default=0.0)
    max_trace = max((x["gram_partial_trace_relative"] for x in records), default=0.0)
    max_decomp = max((x["decomposition_relative"] for x in records), default=0.0)
    max_rotation = max((max(x["rotation_orthogonality"], x["rotation_determinant"]) for x in records), default=0.0)
    open_count = sum(bool(x["phase_space_open"]) for x in records)
    expected_coverage_keys = {
        (int(A), float(ER), str(halo_name))
        for halo_name, _pdf, _vmax in HALOS
        for A in isotopes
        for ER in recoil_energies
    }
    observed_coverage_keys = {
        (int(row["A"]), float(row["ER_keV"]), str(row["halo_model"]))
        for row in coverage
    }
    coverage_complete = (
        len(coverage) == len(expected_coverage_keys)
        and observed_coverage_keys == expected_coverage_keys
        and all(
            int(row["speed_node_count"]) == (args.nv if row["open_speed_interval"] else 0)
            for row in coverage
        )
        and len(tasks) == sum(int(row["speed_node_count"]) for row in coverage)
    )
    validation = {
        "status": "PASS" if (coverage_complete
                             and max_herm <= HERMITICITY_TOL and max_psd_violation <= PSD_RELATIVE_TOL
                             and max_trace <= PARTIAL_TRACE_TOL and max_decomp <= DECOMPOSITION_TOL
                             and max_rotation <= 1.0e-12) else "FAIL",
        "coverage_complete": coverage_complete,
        "point_coverage_rows": len(coverage),
        "expected_point_coverage_rows": len(expected_coverage_keys),
        "tensor_rows": len(records), "expected_tensor_rows_from_speed_grid": len(tasks),
        "phase_space_open_rows": open_count,
        "max_hermiticity_relative": max_herm, "hermiticity_tolerance": HERMITICITY_TOL,
        "max_psd_negative_eigenvalue_relative": max_psd_violation, "psd_relative_tolerance": PSD_RELATIVE_TOL,
        "max_gram_partial_trace_relative": max_trace, "partial_trace_tolerance": PARTIAL_TRACE_TOL,
        "max_current_diagram_decomposition_relative": max_decomp, "decomposition_tolerance": DECOMPOSITION_TOL,
        "max_q_frame_rotation_residual": max_rotation, "rotation_tolerance": 1.0e-12,
    }
    convergence = run_convergence_checks(
        records, isotope_data, recoil_energies,
        enabled=not args.skip_convergence, nv_low=args.nv,
    )
    if convergence["status"] not in {"PASS", "NOT_RUN_PILOT_ONLY"}:
        validation["status"] = "FAIL"
    if args.pilot:
        validation["status"] = "PILOT_NOT_PRODUCTION"

    arrays = {
        "isotope_A": np.asarray([r["A"] for r in records], dtype=np.int16),
        "ER_keV": np.asarray([r["ER_keV"] for r in records], dtype=np.float64),
        "halo_model": np.asarray([r["halo_model"] for r in records], dtype="U16"),
        "speed_index": np.asarray([r["speed_index"] for r in records], dtype=np.int16),
        "speed_kms": np.asarray([r["speed_kms"] for r in records], dtype=np.float64),
        "speed_weight_kms": np.asarray([r["speed_weight_kms"] for r in records], dtype=np.float64),
        "halo_pdf_kms_inv": np.asarray([r["halo_pdf_kms_inv"] for r in records], dtype=np.float64),
        "beta": np.asarray([r["beta"] for r in records], dtype=np.float64),
        "pchi_GeV": np.asarray([r["pchi_GeV"] for r in records], dtype=np.float64),
        "phase_space_open": np.asarray([r["phase_space_open"] for r in records], dtype=bool),
        "L_total": np.stack([r["L"] for r in records]),
        "L_before": np.stack([r["L_before"] for r in records]),
        "L_after": np.stack([r["L_after"] for r in records]),
        "L_interference": np.stack([r["L_interference"] for r in records]),
    }
    tensor_path = args.out_dir / "full_lorentz_tensor_grid.npz"
    np.savez_compressed(tensor_path, **arrays)
    index_rows = []
    for i, r in enumerate(records):
        index_rows.append({
            "row_id": i, "A": r["A"], "ER_keV": r["ER_keV"],
            "halo_model": r["halo_model"], "speed_index": r["speed_index"],
            "speed_kms": r["speed_kms"], "speed_weight_kms": r["speed_weight_kms"],
            "halo_pdf_kms_inv": r["halo_pdf_kms_inv"], "beta": r["beta"],
            "pchi_GeV": r["pchi_GeV"], "mA_GeV": r["mA_GeV"],
            "Q_A_GeV": r["Q_A_GeV"], "t_A_GeV2": r["t_A_GeV2"],
            "t_N_GeV2": r["t_N_GeV2"], "phase_space_open": r["phase_space_open"],
            "skip_reason": r["skip_reason"], "quadrature_points": r["quadrature_points"],
            "hermiticity_relative": r["hermiticity_relative"],
            "psd_min_eigenvalue_relative": r["psd_min_eigenvalue_relative"],
            "gram_partial_trace_relative": r["gram_partial_trace_relative"],
            "decomposition_relative": r["decomposition_relative"],
        })
    index_path = args.out_dir / "tensor_index.csv"
    write_csv(index_path, index_rows)
    coverage_path = args.out_dir / "point_coverage.csv"
    write_csv(coverage_path, coverage)
    convergence_path = args.out_dir / "convergence_report.json"
    convergence_path.write_text(json.dumps(convergence, indent=2) + "\n", encoding="utf-8")

    source_paths = {
        "C1_script": Path(__file__).resolve(),
        "phase_space_tensor": TOPIC / "three_body_phase_space_tensor.py",
        "halo_fold": TOPIC / "native_charge_sector_halo_fold.py",
        "stage_b_direct_current": STAGE_B_SCRIPT,
    }
    source_hashes = {name: sha256_file(path) for name, path in source_paths.items()}
    input_hashes = {name: actual_hashes[name] for name in actual_hashes}
    input_hashes["B_freeze_manifest"] = actual_hashes["B_freeze"]
    run_manifest = {
        "schema": "p4-fc-stage-c1-full-lorentz-tensor-grid-v1",
        "status": validation["status"] if convergence["status"] in {"PASS", "NOT_RUN_PILOT_ONLY"} else "FAIL_CONVERGENCE",
        "production": not args.pilot,
        "scope": "Uncontracted q-aligned three-body L^{mu nu}(v,E_R) tensors; no nuclear-current contraction or rate fold.",
        "repository": "lm846802297-bit/LZ-3-body-for-hdm",
        "git": git_state(),
        "benchmark": {
            "mchi_GeV": MCHI_GEV, "ms_GeV": MS_GEV, "delta_m_GeV": DELTA_M_GEV,
            "g2": G2, "width_scheme": WIDTH_SCHEME, "mapping": MAPPING,
        },
        "mass_scheme": "AME2020 nuclear mass=(atomic_mass_u-54*m_e_u)*u_GeV",
        "grid": {
            "isotopes": list(isotopes), "recoil_points": len(recoil_energies),
            "recoil_keV": list(recoil_energies),
            "halo_models": [{"name": h[0], "vmax_kms": h[2]} for h in HALOS],
            "speed_order": args.nv, "phase_space_order": [args.nw, args.nz, args.nphi],
            "tensor_shape": [4, 4], "tensor_frame": "right-handed q_A-aligned frame",
            "rows": len(records), "open_rows": open_count,
        },
        "normalization": {
            "nucleon_form_factors_included": False,
            "AFH_C0_applied": False,
            "AFH_factor_four_applied": False,
            "current_vertex_form_factor_convention": "reserved for Stage C2; Stage B archive is pinned unchanged",
        },
        "source_sha256": source_hashes,
        "input_sha256": input_hashes,
        "validation": validation,
        "convergence": convergence,
        "runtime": {
            "started_utc": started, "wall_time_seconds": time.perf_counter() - timer,
            "python": sys.version, "platform": platform.platform(),
            "machine": platform.machine(), "processor": platform.processor(),
            "logical_cpu_count": os.cpu_count(), "memory_total_bytes": total_memory_bytes(),
            "numpy": np.__version__, "scipy": __import__("scipy").__version__,
            "workers": args.workers, "command_line": sys.argv,
            "gpu_used": False, "github_actions_used": False,
        },
        "outputs": {
            "full_lorentz_tensor_grid.npz": sha256_file(tensor_path),
            "tensor_index.csv": sha256_file(index_path),
            "point_coverage.csv": sha256_file(coverage_path),
            "convergence_report.json": sha256_file(convergence_path),
        },
        "scope_boundary": {
            "stage_c0_anchor_validation_reused": True,
            "stage_a_native_afh_regenerated": False,
            "stage_b_direct_current_regenerated": False,
            "scalar_charge_halo_fold_repeated": False,
            "full_current_complete": False,
            "stage_c1_complete": validation["status"] == "PASS" and convergence["status"] == "PASS" and not args.pilot,
        },
    }
    manifest_path = args.out_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")
    output_paths = (tensor_path, index_path, coverage_path, convergence_path, manifest_path)
    sums_path = args.out_dir / "SHA256SUMS.txt"
    sums_path.write_text("".join(f"{sha256_file(p)}  {p.name}\n" for p in output_paths), encoding="utf-8")
    print(json.dumps({
        "status": run_manifest["status"], "production": run_manifest["production"],
        "tensor_rows": len(records), "open_rows": open_count,
        "validation": validation, "convergence": convergence,
        "output_dir": str(args.out_dir), "wall_time_seconds": time.perf_counter() - timer,
    }, indent=2), flush=True)
    return 0 if run_manifest["status"] == "PASS" or args.pilot else 1


if __name__ == "__main__":
    raise SystemExit(main())
