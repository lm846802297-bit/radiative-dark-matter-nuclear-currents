#!/usr/bin/env python3
"""Independent integrity and response-reconstruction audit for the P4-FC Native AFH export."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

ISOTOPES = (128, 129, 130, 131, 132, 134, 136)
ENERGIES_KEV = tuple(40.0 + 9.5 * i for i in range(81))
CHARGE_RANKS = {
    128: (0,), 129: (0, 1), 130: (0,), 131: (0, 1, 2, 3),
    132: (0,), 134: (0,), 136: (0,),
}
SPATIAL_RANKS = {128: (), 129: (1,), 130: (), 131: (1, 2, 3), 132: (), 134: (), 136: ()}
DENSITY_ISOSCALAR_ROWS = {128: 16, 129: 25, 130: 16, 131: 39, 132: 16, 134: 16, 136: 16}
EXPECTED_AFH_SHA256 = "dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372"
EXPECTED_MASS_SHA256 = "fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991"
U_GEV = 0.93149410242
M_E_U = 0.000548579909065
Z_XE = 54
HBARC_GEV_FM = 0.1973269804


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def parse_metadata(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def f(row: dict[str, str], key: str) -> float:
    text = row[key].strip().replace("*^", "e")
    # Wolfram CSV can serialize arbitrary-precision reals as 0.5`20.
    if "`" in text:
        text = text.split("`", 1)[0]
    return float(text)


def close(a: float, b: float, *, rel: float = 5e-7, abs_: float = 2e-12) -> bool:
    return abs(a - b) <= abs_ + rel * max(abs(a), abs(b))


def b_parameter_GeV_inv(A: int) -> float:
    b_fm = math.sqrt(41.467 / (45.0 * A ** (-1.0 / 3.0) - 25.0 * A ** (-2.0 / 3.0)))
    return b_fm * 5.0677


def reconstruct_response(
    multipoles: dict[tuple[int, float, int, str], complex],
    A: int,
    ER_keV: float,
    y: float,
    channel: str,
) -> float:
    ranks = CHARGE_RANKS[A] if channel in ("M", "PhiPP") else SPATIAL_RANKS[A]
    total = 0.0
    for J in ranks:
        amp = multipoles.get((A, ER_keV, J, channel), 0j)
        total += abs(amp) ** 2 if channel not in ("MPhiPP", "SigmaPrimeDelta") else 0.0
    if channel == "MPhiPP":
        total = sum(
            (multipoles.get((A, ER_keV, J, "M"), 0j)
             * multipoles.get((A, ER_keV, J, "PhiPP"), 0j).conjugate()).real
            for J in CHARGE_RANKS[A]
        )
    elif channel == "SigmaPrimeDelta":
        total = sum(
            (multipoles.get((A, ER_keV, J, "SigmaPrime"), 0j)
             * multipoles.get((A, ER_keV, J, "Delta"), 0j).conjugate()).real
            for J in SPATIAL_RANKS[A]
        )
    # ResponseNuclear already carries E^(-2y); the E^y inside AFH Response
    # cancels it after the modulus square. Do not apply a second exponential.
    return total / 4.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent / "output/native-seven-isotope-81-point")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    root = args.root.resolve()

    multipole_path = outdir / "afh_v6_p4_fc_7isotope_81point_native_multipoles.csv"
    response_path = outdir / "afh_v6_p4_fc_7isotope_81point_native_responses.csv"
    term_path = outdir / "afh_v6_p4_fc_7isotope_81point_native_obdme_terms.csv"
    metadata_path = outdir / "native_export_metadata.txt"
    required = (multipole_path, response_path, term_path, metadata_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing Native export files: " + ", ".join(missing))

    afh_path = root / "v6" / "dmformfactor-V6.m"
    mass_path = root / "research/afh_default_density_natural_xenon_halo/input/natural_xenon_isotope_weights.csv"
    script_path = Path(__file__).resolve().parent / "export_p4_fc_native_multipoles_7isotope_81point.wls"
    metadata = parse_metadata(metadata_path)

    afh_sha = sha256(afh_path)
    mass_sha = sha256(mass_path)
    script_sha = sha256(script_path)
    audit_sha = sha256(Path(__file__).resolve())
    test_sha = sha256(Path(__file__).resolve().parent / "test_audit_p4_fc_native_export.py")
    assert afh_sha == EXPECTED_AFH_SHA256, (afh_sha, EXPECTED_AFH_SHA256)
    assert mass_sha == EXPECTED_MASS_SHA256, (mass_sha, EXPECTED_MASS_SHA256)
    assert metadata["afh_source_sha256"] == afh_sha
    assert metadata["mass_input_sha256"] == mass_sha
    assert metadata["script_sha256"] == script_sha
    assert metadata["multipole_csv_sha256"] == sha256(multipole_path)
    assert metadata["response_csv_sha256"] == sha256(response_path)
    assert metadata["obdme_terms_csv_sha256"] == sha256(term_path)

    multipole_rows = read_csv(multipole_path)
    response_rows = read_csv(response_path)
    term_rows = read_csv(term_path)
    expected_points = len(ISOTOPES) * len(ENERGIES_KEV)
    expected_charge_rows = 2 * sum(len(CHARGE_RANKS[A]) for A in ISOTOPES) * len(ENERGIES_KEV)
    expected_spatial_rows = 2 * sum(len(SPATIAL_RANKS[A]) for A in ISOTOPES) * len(ENERGIES_KEV)
    expected_term_rows = sum(DENSITY_ISOSCALAR_ROWS.values()) * len(ENERGIES_KEV)
    assert len(response_rows) == expected_points, (len(response_rows), expected_points)
    assert len(multipole_rows) == expected_charge_rows + expected_spatial_rows
    assert len(term_rows) == expected_term_rows
    assert int(metadata["recoil_isotope_points"]) == expected_points
    assert int(metadata["multipole_rows"]) == len(multipole_rows)
    assert int(metadata["native_response_rows"]) == len(response_rows)
    assert int(metadata["obdme_term_rows"]) == len(term_rows)

    expected_energy = set(ENERGIES_KEV)
    response_by_point = {}
    for row in response_rows:
        A = int(row["A"])
        ER = f(row, "ER_keV")
        assert A in ISOTOPES and any(abs(ER - e) < 1e-9 for e in expected_energy)
        key = (A, round(ER, 9))
        assert key not in response_by_point, f"duplicate native response point {key}"
        response_by_point[key] = row

    multipoles = {}
    max_imag = 0.0
    channel_counts = defaultdict(int)
    for row in multipole_rows:
        A = int(row["A"])
        ER = round(f(row, "ER_keV"), 9)
        J = int(row["J"])
        channel = row["channel"]
        amp = complex(f(row, "amplitude_real"), f(row, "amplitude_imag"))
        key = (A, ER, J, channel)
        response_reference = response_by_point[(A, ER)]
        assert close(f(row, "t_A_GeV2"), f(response_reference, "t_A_GeV2"), rel=0.0, abs_=2e-13)
        assert close(f(row, "Q_A_GeV"), f(response_reference, "Q_A_GeV"), rel=0.0, abs_=2e-13)
        assert key not in multipoles, f"duplicate multipole {key}"
        multipoles[key] = amp
        max_imag = max(max_imag, abs(amp.imag))
        channel_counts[channel] += 1
    assert set(channel_counts) == {"M", "PhiPP", "Delta", "SigmaPrime"}
    assert channel_counts["M"] == channel_counts["PhiPP"] == expected_charge_rows // 2
    assert channel_counts["Delta"] == channel_counts["SigmaPrime"] == expected_spatial_rows // 2
    assert max_imag <= 2e-12, max_imag

    max_response_imag = 0.0
    max_response_residual = {name: 0.0 for name in ("M", "PhiPP", "MPhiPP", "Delta", "SigmaPrime", "SigmaPrimeDelta")}
    max_cauchy_excess = 0.0
    max_kinematic_residual = 0.0
    max_tA_residual = 0.0
    max_y_residual = 0.0
    c0_by_A = {}
    for key, row in response_by_point.items():
        A, ER_keV = key
        ER = ER_keV * 1e-6
        mA = f(row, "mA_GeV")
        Q = f(row, "Q_A_GeV")
        y = f(row, "y")
        JA = f(row, "J_A")
        c0 = f(row, "C0")
        expected_Q = math.sqrt(2.0 * mA * ER + ER * ER)
        expected_tA = -2.0 * mA * ER
        max_tA_residual = max(max_tA_residual, abs(f(row, "t_A_GeV2") - expected_tA))
        max_kinematic_residual = max(max_kinematic_residual, abs(Q - expected_Q))
        expected_y = (expected_Q * b_parameter_GeV_inv(A) / 2.0) ** 2
        max_y_residual = max(max_y_residual, abs(y - expected_y))
        assert close(f(row, "ER_GeV"), ER, rel=0.0, abs_=2e-15)
        assert close(Q, expected_Q, rel=0.0, abs_=2e-13)
        assert close(f(row, "t_A_GeV2"), expected_tA, rel=0.0, abs_=2e-13)
        assert close(y, expected_y, rel=0.0, abs_=2e-12)
        assert JA == (0.0 if A not in (129, 131) else (0.5 if A == 129 else 1.5))
        if A in c0_by_A:
            assert c0_by_A[A] == c0
        assert math.isfinite(c0) and c0 > 0.0, (A, c0)
        c0_by_A[A] = c0

        for name in ("M", "SigmaDoublePrime", "SigmaPrime", "PhiPP", "PhiTildePrime", "Delta", "MPhiPP", "SigmaPrimeDelta"):
            max_response_imag = max(max_response_imag, abs(f(row, "Im_W_" + name + "_raw")))
        for channel in max_response_residual:
            reconstructed = reconstruct_response(multipoles, A, ER_keV, y, channel)
            stored = f(row, "W_" + channel + "_raw")
            residual = abs(reconstructed - stored)
            max_response_residual[channel] = max(max_response_residual[channel], residual)
            assert close(reconstructed, stored), (A, ER_keV, channel, reconstructed, stored, residual)

        wm = f(row, "W_M_raw")
        wphi = f(row, "W_PhiPP_raw")
        wdelta = f(row, "W_Delta_raw")
        wsigma = f(row, "W_SigmaPrime_raw")
        wmp = f(row, "W_MPhiPP_raw")
        wsd = f(row, "W_SigmaPrimeDelta_raw")
        max_cauchy_excess = max(
            max_cauchy_excess,
            max(0.0, abs(wmp) - math.sqrt(max(0.0, wm * wphi))),
            max(0.0, abs(wsd) - math.sqrt(max(0.0, wdelta * wsigma))),
        )
        assert abs(wmp) <= math.sqrt(max(0.0, wm * wphi)) + 5e-12
        assert abs(wsd) <= math.sqrt(max(0.0, wdelta * wsigma)) + 5e-12

    assert max_response_imag <= 2e-12, max_response_imag
    assert set(response_by_point) == {
        (A, round(ER, 9)) for A in ISOTOPES for ER in ENERGIES_KEV
    }
    assert metadata["status"] == "NATIVE_AFH_P4_FC_SEVEN_ISOTOPE_EXPORT_COMPLETE"

    report = {
        "status": "PASS_P4_FC_NATIVE_SEVEN_ISOTOPE_EXPORT_AUDIT",
        "scope": "Native AFH response/multipole export integrity and rank reconstruction only; no direct bound-current contraction or halo fold.",
        "source_sha256": {
            "audit_script": audit_sha,
            "audit_unit_tests": test_sha,
            "afh_v6": afh_sha,
            "isotope_mass_input": mass_sha,
            "native_export_wls": script_sha,
            "native_multipoles_csv": sha256(multipole_path),
            "native_responses_csv": sha256(response_path),
            "obdme_terms_csv": sha256(term_path),
        },
        "coverage": {
            "isotopes": list(ISOTOPES),
            "recoil_points_per_isotope": len(ENERGIES_KEV),
            "isotope_recoil_points": len(response_rows),
            "charge_multipole_rows": expected_charge_rows,
            "spatial_multipole_rows": expected_spatial_rows,
            "obdme_term_rows": expected_term_rows,
            "multipole_channel_rows": dict(channel_counts),
        },
        "checks": {
            "max_multipole_imaginary_residual": max_imag,
            "max_native_response_imaginary_residual": max_response_imag,
            "max_Q_A_reconstruction_residual_GeV": max_kinematic_residual,
            "max_t_A_reconstruction_residual_GeV2": max_tA_residual,
            "max_y_reconstruction_residual": max_y_residual,
            "max_response_absolute_residuals": max_response_residual,
            "max_Cauchy_excess": max_cauchy_excess,
            "response_reconstruction_tolerance": "abs <= 2e-12 + 5e-7*max(abs(reconstructed),abs(native))",
            "Cauchy_tolerance_absolute": 5e-12,
        },
        "normalization": {
            "native_response_convention": "raw AFH package tau1=tau2=0",
            "physical_factor_four_and_C0_applied": False,
            "C0_by_isotope": {str(A): c0_by_A[A] for A in ISOTOPES},
        },
        "native_runtime": {
            "wolfram_version": metadata.get("wolfram_version"),
            "system": metadata.get("system"),
            "system_id": metadata.get("system_id"),
            "processor_type": metadata.get("processor_type"),
            "processor_count": metadata.get("processor_count"),
            "machine_name": metadata.get("machine_name"),
        },
        "audit_runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "processor": platform.processor(),
            "logical_cpu_count": __import__("os").cpu_count(),
        },
    }
    report_path = outdir / "validation_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    sums_path = outdir / "SHA256SUMS.txt"
    outputs = (multipole_path, response_path, term_path, metadata_path, report_path)
    sums_path.write_text("".join(f"{sha256(path)}  {path.name}\n" for path in outputs), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
