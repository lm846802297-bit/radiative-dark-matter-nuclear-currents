#!/usr/bin/env python3
"""Independent audit of the P4-FC Stage B direct finite-tN Dirac-Pauli current."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TOPIC_DIR = ROOT / "research" / "three_body_full_current_nucleon_match"
if str(TOPIC_DIR) not in sys.path:
    sys.path.insert(0, str(TOPIC_DIR))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import direct_bound_current_seven_isotope as direct
from nucleon_vector_form_factors import borah_sufian2017, effective_nucleon_form_factors

ISOTOPES = (128, 129, 130, 131, 132, 134, 136)
ENERGIES_KEV = tuple(40.0 + 9.5 * i for i in range(81))
CHARGE_RANKS = {
    128: (0,), 129: (0, 1), 130: (0,), 131: (0, 1, 2, 3),
    132: (0,), 134: (0,), 136: (0,),
}
EXPECTED_NATIVE_SHA256 = "26113db76a969dba98b0fb7b4d6b5389189d2dec33ec680d699741fb0b43a2dc"
EXPECTED_MASS_SHA256 = "fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991"
EXPECTED_DENSITY_SHA256 = dict(direct.DENSITY_SHA256)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def close(a: float, b: float, *, rel: float = 5e-7, abs_: float = 1e-12) -> bool:
    return abs(a - b) <= abs_ + rel * max(abs(a), abs(b))


def complex_field(row: dict[str, str], name: str) -> complex:
    return complex(float(row[name + "_re"]), float(row[name + "_im"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, default=direct.OUT_DIR / "direct_bound_current_o14.json")
    parser.add_argument("--csv", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    json_path = args.json.resolve()
    csv_path = args.csv.resolve() if args.csv else json_path.with_suffix(".csv")
    out_path = args.out.resolve() if args.out else json_path.parent / "stage_b_validation_report.json"
    if not json_path.is_file() or not csv_path.is_file():
        raise FileNotFoundError("Stage B JSON/CSV output is missing")

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    rows = read_csv(csv_path)
    if payload.get("status") != "COMPLETED_P4_FC_STAGE_B_DIRECT_CURRENT_RUN":
        raise ValueError(f"unexpected Stage B status: {payload.get('status')}")
    if payload.get("mass_scheme_key") != "ame2020":
        raise ValueError("production Stage B requires the P4 AME2020 nuclear-mass scheme")
    if payload.get("vertex_energy_prescription") != "external":
        raise ValueError("production Stage B requires the formal external q_N^0=E_R prescription")
    if tuple(payload.get("isotopes", ())) != ISOTOPES:
        raise ValueError("production Stage B does not cover all seven AFH isotopes")
    if tuple(round(float(x) * 1e6, 9) for x in payload.get("recoil_energies_GeV", ())) != tuple(round(x, 9) for x in ENERGIES_KEV):
        raise ValueError("production Stage B does not use the canonical 81-point recoil grid")

    native_path = direct.NATIVE_CSV_P4
    mass_path = direct.NUCLEAR_MASS_INPUT
    script_path = HERE / "direct_bound_current_seven_isotope.py"
    native_sha = sha256(native_path)
    mass_sha = sha256(mass_path)
    script_sha = sha256(script_path)
    assert native_sha == EXPECTED_NATIVE_SHA256, native_sha
    assert mass_sha == EXPECTED_MASS_SHA256, mass_sha
    source = payload["source_sha256"]
    assert source["native_multipoles"] == native_sha
    assert source["nuclear_mass_input"] == mass_sha
    assert source["direct_script"] == script_sha
    density_hashes = {
        str(A): sha256(direct.DENSITY_DIR / f"afh_default_xe{A}_density_matrix.csv")
        for A in ISOTOPES
    }
    assert density_hashes == {str(A): EXPECTED_DENSITY_SHA256[A] for A in ISOTOPES}
    assert source["density_matrices"] == density_hashes
    source_files = {
        "form_factor_module": TOPIC_DIR / "nucleon_vector_form_factors.py",
        "borah_source": ROOT / "research/nucleon_form_factor_B_baseline/borah2020_lowq2.py",
        "borah_parameters": ROOT / "research/nucleon_form_factor_B_baseline/data/borah2020/parameters.json",
        "borah_covariance": ROOT / "research/nucleon_form_factor_B_baseline/data/borah2020/ELZExpFF.dat",
        "rate_kinematics": ROOT / "research/three_body_calculation_and_validation/rate_kinematics.py",
    }
    for key, source_path in source_files.items():
        if key in source:
            assert source[key] == sha256(source_path), (key, source[key], sha256(source_path))

    expected_rows = sum(len(CHARGE_RANKS[A]) for A in ISOTOPES) * len(ENERGIES_KEV)
    assert len(rows) == expected_rows, (len(rows), expected_rows)
    assert len(payload["results"]) == expected_rows

    keys = set()
    by_point = {}
    max_fit = 0.0
    max_tn = 0.0
    max_q = 0.0
    max_q2_identity = 0.0
    max_ff_residual = {k: 0.0 for k in ("F1_uds_raw", "F2_uds_raw", "F1_effective", "F2_effective")}
    max_nr_match = {k: 0.0 for k in ("M", "PhiPP", "Delta", "SigmaPrime")}
    max_rank0_alias = 0.0
    rank0_rows = 0
    spin_zero_rank0_rows = 0
    rank0_longitudinal_max_by_A = {A: 0.0 for A in ISOTOPES}
    max_exact_vs_nr = {k: 0.0 for k in (
        "vector_identity_LJm1", "vector_identity_LJ", "vector_identity_LJp1",
        "vector_spin_LJm1", "vector_spin_LJ", "vector_spin_LJp1",
    )}

    for row in rows:
        A = int(row["A"])
        ER = float(row["ER_GeV"])
        ER_keV = float(row["ER_keV"])
        J = int(row["J"])
        JA = float(row["J_A"])
        key = (A, round(ER_keV, 9), J)
        assert key not in keys, f"duplicate Stage B rank row {key}"
        keys.add(key)
        assert A in ISOTOPES
        assert any(abs(ER_keV - e) < 1e-8 for e in ENERGIES_KEV)
        assert J in CHARGE_RANKS[A]
        assert JA == direct.TARGET_SPIN[A]
        by_point.setdefault((A, round(ER_keV, 9)), set()).add(J)

        mA = float(row["mA_GeV"])
        tN = float(row["tN_GeV2"])
        Q = float(row["Q_A_GeV"])
        QN2 = float(row["QN2_GeV2"])
        expected_tN = -2.0 * mA * ER
        expected_Q = math.sqrt(2.0 * mA * ER + ER * ER)
        max_tn = max(max_tn, abs(tN - expected_tN))
        max_q = max(max_q, abs(Q - expected_Q))
        max_q2_identity = max(max_q2_identity, abs(Q * Q - (QN2 + ER * ER)))
        assert close(tN, expected_tN, rel=0.0, abs_=2e-13)
        assert close(Q, expected_Q, rel=0.0, abs_=2e-13)
        assert close(QN2, -tN, rel=0.0, abs_=2e-13)
        assert close(Q * Q, QN2 + ER * ER, rel=0.0, abs_=2e-13)

        raw_f1, raw_f2 = borah_sufian2017(-tN)
        eff_f1, eff_f2 = effective_nucleon_form_factors(raw_f1, raw_f2)
        for name, expected in (
            ("F1_uds_raw", raw_f1), ("F2_uds_raw", raw_f2),
            ("F1_effective", eff_f1), ("F2_effective", eff_f2),
        ):
            residual = abs(float(row[name]) - expected)
            max_ff_residual[name] = max(max_ff_residual[name], residual)
            assert close(float(row[name]), expected, rel=0.0, abs_=2e-13), (A, ER_keV, name)

        max_fit = max(max_fit, float(row["max_fit_residual"]))
        assert float(row["max_fit_residual"]) < 1e-10
        for direct_name, expected_name, channel in (
            ("leading_nr_charge_identity", "expected_AFH_current_M", "M"),
            ("leading_nr_charge_spin", "expected_AFH_current_PhiPP", "PhiPP"),
            ("leading_NR_Delta_effective", "expected_AFH_current_Delta", "Delta"),
            ("leading_NR_SigmaPrime_effective", "expected_AFH_current_SigmaPrime", "SigmaPrime"),
        ):
            residual = abs(complex_field(row, direct_name) - complex_field(row, expected_name))
            max_nr_match[channel] = max(max_nr_match[channel], residual)
            assert residual < 1e-10, (A, ER_keV, J, channel, residual)

        for component in max_exact_vs_nr:
            name = "exact_vs_leading_" + component
            value = abs(complex_field(row, name))
            max_exact_vs_nr[component] = max(max_exact_vs_nr[component], value)

        if J == 0:
            rank0_rows += 1
            if JA == 0.0:
                spin_zero_rank0_rows += 1
            for rep in ("exact", "leading_nr"):
                for kind in ("identity", "spin"):
                    alias = complex_field(row, f"{rep}_rank0_longitudinal_{kind}")
                    expected = complex_field(row, f"{rep}_vector_{kind}_LJp1")
                    max_rank0_alias = max(max_rank0_alias, abs(alias - expected))
                    assert close(abs(alias - expected), 0.0, rel=0.0, abs_=1e-14)
                    if rep == "exact":
                        rank0_longitudinal_max_by_A[A] = max(
                            rank0_longitudinal_max_by_A[A], abs(alias)
                        )

    assert rank0_rows == len(ISOTOPES) * len(ENERGIES_KEV), rank0_rows
    assert spin_zero_rank0_rows == 5 * len(ENERGIES_KEV), spin_zero_rank0_rows
    assert max_rank0_alias == 0.0
    expected_points = {(A, round(E, 9)) for A in ISOTOPES for E in ENERGIES_KEV}
    assert set(by_point) == expected_points
    for A, E in expected_points:
        assert by_point[(A, E)] == set(CHARGE_RANKS[A]), (A, E, by_point[(A, E)])

    report = {
        "status": "PASS_P4_FC_STAGE_B_DIRECT_CURRENT_AUDIT",
        "scope": "Seven-isotope finite-tN direct Dirac-Pauli impulse-current multipoles; no T3 contraction or halo fold.",
        "source_sha256": {
            "direct_script": script_sha,
            "native_stage_a_multipoles": native_sha,
            "nuclear_mass_input": mass_sha,
            "density_matrices": density_hashes,
            "form_factor_module": sha256(TOPIC_DIR / "nucleon_vector_form_factors.py"),
            "afh_v6": sha256(ROOT / "v6/dmformfactor-V6.m"),
        },
        "coverage": {
            "isotopes": list(ISOTOPES),
            "recoil_points_per_isotope": len(ENERGIES_KEV),
            "isotope_recoil_points": len(by_point),
            "rank_resolved_rows": len(rows),
            "expected_rank_rows": expected_rows,
            "rank0_longitudinal_rows_all_targets": rank0_rows,
            "spin_zero_rank0_longitudinal_rows": spin_zero_rank0_rows,
            "max_rank0_longitudinal_amplitude_by_isotope": rank0_longitudinal_max_by_A,
            "rank0_L1_presence_check": "all five spin-zero isotopes retain explicit L=1 fields, including any exact zero values; no channel is dropped",
        },
        "checks": {
            "max_direct_multipole_fit_residual": max_fit,
            "max_tN_identity_residual_GeV2": max_tn,
            "max_QA_identity_residual_GeV": max_q,
            "max_QA2_minus_QN2_minus_ER2_residual_GeV2": max_q2_identity,
            "max_form_factor_recalculation_absolute_residuals": max_ff_residual,
            "max_leading_NR_vs_Native_AFH_amplitude_residuals": max_nr_match,
            "max_rank0_longitudinal_alias_residual": max_rank0_alias,
            "max_exact_vs_leading_NR_multipole_differences": max_exact_vs_nr,
            "fit_tolerance": 1e-10,
            "leading_NR_AFH_match_tolerance_absolute": 1e-10,
            "kinematic_tolerance_absolute": 2e-13,
        },
        "model": {
            "nucleon_transfer": "q_N^mu=q_A^mu; t_N=t_A=-2*m_A*E_R; Q_A^2=Q_N^2+E_R^2",
            "form_factors": "Borah2020 proton/neutron electromagnetic form factors plus Sufian2017 unit strange current; raw flavor sums divided by 3 once",
            "vertex_energy": payload["vertex_energy_prescription"],
            "bound_basis": "AFH v6.0 HO momentum orbitals embedded in positive-energy free Dirac spinors",
            "current": "exact covariant Dirac-Pauli and leading-NR control",
            "normalization": "Native AFH raw amplitudes; no factor-four or C0 rescaling",
            "even_even_targets": "retain direct K=0,L=1 longitudinal current; do not map it onto AFH Delta or SigmaPrime",
        },
        "native_run_runtime": payload.get("runtime", {}),
        "audit_runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "processor": platform.processor(),
            "logical_cpu_count": __import__("os").cpu_count(),
            "quadrature_order_per_axis": payload["quadrature_order_per_axis"],
        },
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    sums = [
        (json_path, sha256(json_path)),
        (csv_path, sha256(csv_path)),
        (out_path, sha256(out_path)),
    ]
    sums_path = out_path.parent / "SHA256SUMS_stage_b.txt"
    sums_path.write_text("".join(f"{digest}  {path.name}\n" for path, digest in sums), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
