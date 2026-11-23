#!/usr/bin/env python3
"""Recompute kinematic, Ward, charge-block, q-frame, and T3 acceptance checks."""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[6]
OUT_DIR = Path(__file__).resolve().parents[1]
SOURCE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE_DIR))
import p3c_bound_current_amplitude as direct  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    contraction = json.loads((OUT_DIR / "direct_bound_current_t3_contraction.json").read_text())
    tensor = json.loads((OUT_DIR / "pole_width_t3_tensor_xe129_xe131.json").read_text())

    kinematics = []
    for A in (129, 131):
        for ER in (0.0002, 0.000267):
            mA = A * direct.U_GEV
            Q2 = 2.0 * mA * ER + ER * ER
            tA = -2.0 * mA * ER
            residual = (ER * ER - Q2) - tA
            kinematics.append({
                "A": A,
                "ER_GeV": ER,
                "t_A_GeV2": tA,
                "Q_A2_GeV2": Q2,
                "identity_residual_GeV2": residual,
            })
    max_kinematic_residual = max(abs(x["identity_residual_GeV2"]) for x in kinematics)

    charge_max = {"relative_residual": -1.0}
    for row in contraction["results"]:
        afh = row["K_charge_AFH_L00_only"] + row["signed_M_Phi_kernel_AFH_L00_only"]
        nr = (row["K_leading_nr_charge_matrix_raw"]
              + row["K_leading_nr_M_Phi_cross_matrix_raw"]) * row["C0"]
        relative = abs(nr - afh) / max(abs(nr), abs(afh), 1e-300)
        if relative > charge_max["relative_residual"]:
            charge_max = {
                "relative_residual": relative,
                "A": row["A"],
                "ER_keV": row["ER_keV"],
                "J": row["J"],
                "AFH_L00_charge_plus_signed_M_Phi": afh,
                "leading_NR_charge_plus_signed_M_Phi_times_C0": nr,
            }

    ward_by_case = []
    for A in (129, 131):
        for ER in (0.0002, 0.000267):
            prepared = direct.prepare_current(A, ER, 22, vertex_energy="internal")
            q0 = prepared["q0_internal"]
            Q = prepared["Q"]
            current = prepared["current"]
            q_lower = np.stack(
                [q0, np.zeros_like(q0), np.zeros_like(q0), -np.full_like(q0, Q)],
                axis=1,
            )
            divergence = np.einsum("nm,nmab->nab", q_lower, current)
            scale = (np.abs(q0) + Q) * np.sqrt(
                np.sum(np.abs(current) ** 2, axis=(1, 2, 3))
            )
            residual = np.sqrt(np.sum(np.abs(divergence) ** 2, axis=(1, 2))) / np.maximum(
                scale, 1e-300
            )
            ward_by_case.append({
                "A": A,
                "ER_GeV": ER,
                "quadrature_order_per_axis": 22,
                "quadrature_node_count": int(len(q0)),
                "max_relative_q_mu_J_mu_residual": float(np.max(residual)),
            })
    max_ward_residual = max(x["max_relative_q_mu_J_mu_residual"] for x in ward_by_case)

    q_frame_max = max(
        entry["q_frame_rotation_max_orthogonality_residual"]
        for entry in tensor["benchmarks"]
    )
    external_divergence_rows = contraction["results"]
    max_external_divergence_row = max(
        external_divergence_rows,
        key=lambda x: x["exact_nuclear_current_divergence_over_full_current_scale"],
    )
    max_external_divergence = max_external_divergence_row[
        "exact_nuclear_current_divergence_over_full_current_scale"
    ]

    checks = {
        "kinematic_identities": {
            "max_absolute_residual_GeV2": max_kinematic_residual,
            "tolerance_absolute_GeV2": 1e-13,
            "pass": max_kinematic_residual < 1e-13,
            "cases": kinematics,
        },
        "charge_only_leading_NR_vs_AFH": {
            **charge_max,
            "tolerance_relative": 5e-7,
            "pass": charge_max["relative_residual"] < 5e-7,
        },
        "on_shell_nucleon_Ward_identity_internal_q0": {
            "max_relative_residual": max_ward_residual,
            "tolerance_relative": 1e-10,
            "pass": max_ward_residual < 1e-10,
            "cases": ward_by_case,
        },
        "q_frame_tensor_basis": {
            "max_orthogonality_residual": q_frame_max,
            "tolerance": 1e-12,
            "pass": q_frame_max < 1e-12,
        },
        "external_q0_bound_nuclear_current_continuity": {
            "max_normalized_divergence": max_external_divergence,
            "A": max_external_divergence_row["A"],
            "ER_keV": max_external_divergence_row["ER_keV"],
            "J": max_external_divergence_row["J"],
            "interpretation": (
                "Reported model-specific nonconservation for bound HO orbitals with "
                "external qN0=ER; this is a diagnostic, not a pass criterion."
            ),
        },
        "T3_tensor_acceptance": tensor["acceptance"],
        "T3_direct_contraction_acceptance": contraction["acceptance"],
    }
    scoped_pass = (
        checks["kinematic_identities"]["pass"]
        and checks["charge_only_leading_NR_vs_AFH"]["pass"]
        and checks["on_shell_nucleon_Ward_identity_internal_q0"]["pass"]
        and checks["q_frame_tensor_basis"]["pass"]
        and checks["T3_tensor_acceptance"]["pass"]
        and checks["T3_direct_contraction_acceptance"]["pass"]
    )
    payload = {
        "status": "PASS_SCOPED_P3C_CHECKS" if scoped_pass else "FAIL_SCOPED_P3C_CHECKS",
        "scope": (
            "Selected Xe-129/Xe-131 P3-C direct covariant comparison; "
            "external-q0 nuclear-current divergence is reported separately."
        ),
        "source_sha256": {
            "direct_amplitude_script": sha256(SOURCE_DIR / "p3c_bound_current_amplitude.py"),
            "audit_script": sha256(Path(__file__).resolve()),
        },
        "checks": checks,
    }
    output = OUT_DIR / "p3c_current_identity_checks.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "max_kinematic_residual_GeV2": max_kinematic_residual,
        "max_charge_only_relative_residual": charge_max["relative_residual"],
        "max_on_shell_Ward_residual": max_ward_residual,
        "max_q_frame_orthogonality_residual": q_frame_max,
        "max_external_bound_current_divergence": max_external_divergence,
        "output": str(output),
    }, indent=2))
    if not scoped_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
