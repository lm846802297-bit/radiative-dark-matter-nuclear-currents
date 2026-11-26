#!/usr/bin/env python3
"""Contract the Stage B finite-tN current with the archived exact-mass T3 tensors.

No three-body phase-space integral is evaluated here.  The T3 Gram/Lorentz
tensors and exact-q Native AFH multipoles are read from the immutable
P3C/P4 exact-mass closeout archive.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
import platform
import sys
import tarfile
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TOPIC = ROOT / "research" / "three_body_full_current_nucleon_match"
HERE = Path(__file__).resolve().parent
if str(TOPIC) not in sys.path:
    sys.path.insert(0, str(TOPIC))

import direct_bound_current_seven_isotope as stage_b
import contract_p3c_finite_tN_bound_current as legacy_contract

T3_ARCHIVE = (
    TOPIC / "results" / "2026-11-23" / "p3c_p4_exact_mass_rebase"
    / "native_full_current_closeout"
    / "P3C_P4_Native_ExactMass_FullCurrent_Closeout_2026-11-23.tar.gz"
)
STAGE_A_RESPONSE_CSV = (
    TOPIC / "p4_full_current" / "output" / "native-seven-isotope-81-point"
    / "afh_v6_p4_fc_7isotope_81point_native_responses.csv"
)
STAGE_B_SCRIPT = HERE / "direct_bound_current_seven_isotope.py"
STAGE_B_PRODUCTION_CSV = (
    TOPIC / "p4_full_current" / "output" / "stage-b-direct-current"
    / "direct_bound_current_o14.csv"
)
MASS_INPUT = (
    ROOT / "research" / "afh_default_density_natural_xenon_halo"
    / "input" / "natural_xenon_isotope_weights.csv"
)
OUT_DIR = HERE / "output" / "stage-b-existing-t3-anchor-integration"

EXPECTED_T3_ARCHIVE_GIT_BLOB_SHA1 = "d09de743947b16e62f1a2b183c1785aec83d3847"
EXPECTED_NATIVE_EXACT_MASS_SHA256 = "f23f84d90050044560602d74e75d385b46f851b932202f42a508be4061a22a91"
EXPECTED_STAGE_A_RESPONSE_SHA256 = "39a5b8c2401eb5f4e5acabb0a4f38b9a9d1fa7103d12d208462b3fcefef57c48"
EXPECTED_STAGE_B_SCRIPT_SHA256 = "78058780862acb06d5c3fd86957dead6d16d9a8e3041c9b4bdc3d5450189b75e"
EXPECTED_STAGE_B_PRODUCTION_CSV_SHA256 = "0b335b39cd288beda1367fd977e311d595c2a88c2c5e7dfcaa9456cf8e26d355"
EXPECTED_AFH_SHA256 = "dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372"
EXPECTED_FORM_FACTOR_MODULE_SHA256 = "d409248b3a6d99a8fa9119bb6e5ca906e67360c9442bc37562ff385b8ee85b5b"
EXPECTED_BORAH_SOURCE_SHA256 = "c1acaca81c65e981612958e793804bde3489e5307640e0e9a44a9dad199a83cb"
EXPECTED_BORAH_PARAMETERS_SHA256 = "6dfba4a914259e4dd48bb670b26ab892d53dde75a6e013a79914dde1e40f7281"
EXPECTED_BORAH_COVARIANCE_SHA256 = "044097e0190d815a2f37c6f3192c3ce6c6fb6fce014da6428a862dc92985dbbd"
EXPECTED_RATE_KINEMATICS_SHA256 = "daf409e8246c04c290511dfb78f60992131057f75044cacaeebb159d3f9239a2"
EXPECTED_LEGACY_CONTRACTION_GIT_BLOB_SHA1 = "3a3bc3fc31cb5b4ef5fea3f88a3f22d6652aa46d"
EXPECTED_MASS_SHA256 = "fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991"

EXPECTED_MASSES_GEV = {
    129: 120.04644920255713,
    131: 121.90971983538435,
}
ANCHOR_ER_GEV = (0.0002, 0.000267)
ANCHOR_ISOTOPES = (129, 131)
M_N = 0.938918
METRIC = np.asarray([1.0, -1.0, -1.0, -1.0])


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv_bytes(data: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))


def read_csv_file(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def complex_pairs(value: Any) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape[-1] != 2:
        raise ValueError("complex tensor data must use a final [real,imag] axis")
    return array[..., 0] + 1j * array[..., 1]


def gram_partial_trace_residual(L: np.ndarray, G: np.ndarray) -> float:
    """Check L^{mu nu}=sum_a G^{mu a nu a} for the archived DM-spin Gram tensor."""
    reconstructed = np.einsum("mana->mn", G)
    return float(np.linalg.norm(reconstructed - L) / max(float(np.linalg.norm(L)), 1.0))


def jsonable(value: Any) -> Any:
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def csvable(row: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in row.items():
        if isinstance(value, complex):
            out[key + "_re"] = value.real
            out[key + "_im"] = value.imag
        elif isinstance(value, (dict, list, tuple, np.ndarray)):
            out[key] = json.dumps(jsonable(value), separators=(",", ":"))
        elif isinstance(value, np.generic):
            out[key] = value.item()
        else:
            out[key] = value
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty CSV: {path}")
    flat = [csvable(row) for row in rows]
    fields = list(dict.fromkeys(key for row in flat for key in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(flat)


def read_archive_inputs(archive_path: Path) -> dict[str, Any]:
    archive_bytes = archive_path.read_bytes()
    archive_sha = sha256_bytes(archive_bytes)
    archive_git_blob_sha = git_blob_sha1(archive_bytes)
    if archive_git_blob_sha != EXPECTED_T3_ARCHIVE_GIT_BLOB_SHA1:
        raise ValueError(f"exact-mass closeout Git blob SHA-1 mismatch: {archive_git_blob_sha}")

    with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:gz") as archive:
        members = archive.getmembers()

        native_member = None
        native_bytes = None
        for member in members:
            if not member.isfile() or not member.name.lower().endswith(".csv"):
                continue
            stream = archive.extractfile(member)
            if stream is None:
                continue
            data = stream.read()
            if sha256_bytes(data) == EXPECTED_NATIVE_EXACT_MASS_SHA256:
                native_member, native_bytes = member, data
                break
        if native_member is None or native_bytes is None:
            raise FileNotFoundError(
                "the pinned exact-mass Native AFH multipole CSV was not found inside the closeout archive"
            )

        candidates = []
        for member in members:
            if not member.isfile() or not member.name.lower().endswith(".json"):
                continue
            stream = archive.extractfile(member)
            if stream is None:
                continue
            data = stream.read()
            try:
                payload = json.loads(data.decode("utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError):
                continue
            if not isinstance(payload, dict) or not isinstance(payload.get("benchmarks"), list):
                continue
            if not str(payload.get("status", "")).startswith("PASS_P3C_P4_EXACT_MASS_T3_TENSOR"):
                continue
            if not bool(payload.get("acceptance", {}).get("pass", False)):
                continue
            tensor_rows = [
                row for row in payload["benchmarks"]
                if isinstance(row, dict)
                and "L_q_aligned_real_imag" in row
                and "G_q_aligned_real_imag" in row
            ]
            exact_rows = []
            for row in tensor_rows:
                isotope = int(row.get("isotope_A", row.get("A", -1)))
                mass = row.get("mA_P4_GeV", row.get("mA_GeV"))
                if isotope in EXPECTED_MASSES_GEV and mass is not None:
                    if abs(float(mass) - EXPECTED_MASSES_GEV[isotope]) <= 2e-12:
                        exact_rows.append(row)
            if exact_rows:
                status = str(payload.get("status", ""))
                score = (
                    len(exact_rows),
                    int("P4" in status or "p4_exact_mass" in member.name.lower()),
                    max(int(row.get("quadrature_points", 0)) for row in exact_rows),
                )
                candidates.append((score, member, data, payload, exact_rows))

        if not candidates:
            raise FileNotFoundError(
                "no P4 exact-mass JSON containing both L_q_aligned and G_q_aligned tensors "
                "was found inside the closeout archive"
            )
        candidates.sort(key=lambda item: item[0], reverse=True)
        _, tensor_member, tensor_bytes, tensor_payload, tensor_rows = candidates[0]

        # If the archive carries an internal SHA256SUMS file, retain its exact
        # record for the selected tensor and Native multipole members.
        internal_sums = {}
        for member in members:
            base = Path(member.name).name.lower()
            if member.isfile() and base in {"sha256sums.txt", "sha256sums"}:
                stream = archive.extractfile(member)
                if stream is not None:
                    for line in stream.read().decode("utf-8", errors="replace").splitlines():
                        parts = line.strip().split(None, 1)
                        if len(parts) == 2:
                            internal_sums[Path(parts[1].strip().lstrip("*")).name] = parts[0]
                break

    native_sha = sha256_bytes(native_bytes)
    tensor_sha = sha256_bytes(tensor_bytes)
    if native_sha != EXPECTED_NATIVE_EXACT_MASS_SHA256:
        raise AssertionError("Native exact-mass CSV hash changed after archive extraction")
    expected_internal_native = internal_sums.get(Path(native_member.name).name)
    if expected_internal_native and expected_internal_native != native_sha:
        raise ValueError("archive-internal Native multipole SHA-256 does not match extracted bytes")
    expected_internal_tensor = internal_sums.get(Path(tensor_member.name).name)
    if expected_internal_tensor and expected_internal_tensor != tensor_sha:
        raise ValueError("archive-internal T3 tensor SHA-256 does not match extracted bytes")

    return {
        "archive_sha256": archive_sha,
        "archive_git_blob_sha1": archive_git_blob_sha,
        "native_member": native_member.name,
        "native_bytes": native_bytes,
        "native_sha256": native_sha,
        "tensor_member": tensor_member.name,
        "tensor_payload": tensor_payload,
        "tensor_rows": tensor_rows,
        "tensor_sha256": tensor_sha,
        "internal_sha256": internal_sums,
    }


def enrich_native_csv(native_bytes: bytes) -> bytes:
    rows = read_csv_bytes(native_bytes)
    if not rows:
        raise ValueError("exact-mass Native AFH multipole CSV is empty")
    fields = list(rows[0])
    for required in ("A", "ER_GeV", "mA_GeV", "Q_A_GeV", "J", "channel",
                     "amplitude_real", "amplitude_imag"):
        if required not in fields:
            raise ValueError(f"Native exact-mass CSV lacks required field {required}")
    if "t_A_GeV2" not in fields:
        fields.append("t_A_GeV2")
    if "ER_keV" not in fields:
        fields.append("ER_keV")
    for row in rows:
        er = float(row["ER_GeV"])
        mass = float(row["mA_GeV"])
        row["t_A_GeV2"] = f"{-2.0 * mass * er:.17g}"
        row["ER_keV"] = f"{er * 1.0e6:.17g}"
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def load_stage_a_c0() -> dict[int, float]:
    actual = sha256_file(STAGE_A_RESPONSE_CSV)
    if actual != EXPECTED_STAGE_A_RESPONSE_SHA256:
        raise ValueError(f"Stage A Native response CSV SHA-256 mismatch: {actual}")
    rows = read_csv_file(STAGE_A_RESPONSE_CSV)
    c0: dict[int, float] = {}
    for row in rows:
        isotope = int(row["A"])
        value = float(row["C0"])
        if isotope in c0 and abs(c0[isotope] - value) > 2e-14:
            raise ValueError(f"Stage A C0 is not constant across recoil points for Xe-{isotope}")
        c0[isotope] = value
    if not set(EXPECTED_MASSES_GEV).issubset(c0):
        raise ValueError("Stage A response archive does not contain C0 for Xe-129 and Xe-131")
    return c0


def select_anchor_tensor_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected = []
    for row in rows:
        if row.get("phase_space_open", True) is False:
            continue
        isotope = int(row.get("isotope_A", row.get("A", -1)))
        if isotope not in ANCHOR_ISOTOPES:
            continue
        er_value = row.get("ER_GeV")
        er = float(er_value) if er_value is not None else float(row["ER_keV"]) * 1e-6
        if not any(abs(er - anchor) <= 5e-13 for anchor in ANCHOR_ER_GEV):
            continue
        if "L_q_aligned_real_imag" not in row or "G_q_aligned_real_imag" not in row:
            continue
        mass = row.get("mA_P4_GeV", row.get("mA_GeV"))
        if mass is None or abs(float(mass) - EXPECTED_MASSES_GEV[isotope]) > 2e-12:
            continue
        width = str(row.get("width_scheme", "pole")).lower()
        if width != "pole":
            continue
        selected.append(row)
    if not selected:
        raise ValueError("no open P4 exact-mass Pole-width T3 anchor rows matched Xe-129/131 at 200/267 keV")
    return selected


def contract_components(L: np.ndarray, current: dict[str, np.ndarray], JA: float) -> dict[str, float]:
    bilinear = legacy_contract.bilinear
    full = bilinear(L, current["full"], current["full"], JA)
    charge = bilinear(L, current["charge"], current["charge"], JA)
    spatial = bilinear(L, current["spatial"], current["spatial"], JA)
    charge_spatial = 2.0 * bilinear(L, current["charge"], current["spatial"], JA).real
    m_phi = 2.0 * bilinear(
        L, current["charge_identity"], current["charge_spin"], JA
    ).real
    delta_sigma_lj = 2.0 * bilinear(
        L, current["delta_LJ"], current["sigma_LJ"], JA
    ).real
    spatial_identity_spin = 2.0 * bilinear(
        L, current["spatial_identity"], current["spatial_spin"], JA
    ).real

    longitudinal = current["full"].copy()
    longitudinal[1:3] = 0.0
    transverse = current["full"].copy()
    transverse[[0, 3]] = 0.0
    long_kernel = bilinear(L, longitudinal, longitudinal, JA)
    transverse_kernel = bilinear(L, transverse, transverse, JA)

    decomposition = charge.real + spatial.real + charge_spatial
    decomposition_absolute = float(abs(full.real - decomposition))
    decomposition_relative = decomposition_absolute / max(abs(full.real), abs(decomposition), 1.0)
    imaginary_absolute = float(abs(full.imag))
    imaginary_relative = imaginary_absolute / max(abs(full.real), 1.0)
    return {
        "K_full_re": float(full.real),
        "K_full_im": float(full.imag),
        "K_charge_re": float(charge.real),
        "K_charge_im": float(charge.imag),
        "K_spatial_re": float(spatial.real),
        "K_spatial_im": float(spatial.imag),
        "K_charge_spatial_cross": float(charge_spatial),
        "K_M_Phi_signed_cross": float(m_phi),
        "K_Delta_Sigma_LJ_signed_cross": float(delta_sigma_lj),
        "K_spatial_identity_spin_signed_cross": float(spatial_identity_spin),
        "K_longitudinal_0z_re": float(long_kernel.real),
        "K_transverse_xy_re": float(transverse_kernel.real),
        "decomposition_absolute_residual": decomposition_absolute,
        "decomposition_relative_residual": float(decomposition_relative),
        "imaginary_residual": imaginary_absolute,
        "imaginary_relative_residual": float(imaginary_relative),
    }


def sum_parts(parts: list[dict[str, float]]) -> dict[str, float]:
    if not parts:
        raise ValueError("cannot sum an empty rank list")
    keys = tuple(parts[0])
    return {
        key: float(max(part[key] for part in parts) if key.endswith("_relative_residual") else sum(part[key] for part in parts))
        for key in keys
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, default=T3_ARCHIVE)
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    parser.add_argument("--order", type=int, default=14)
    args = parser.parse_args()

    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    timer = time.perf_counter()
    archive_inputs = read_archive_inputs(args.archive.resolve())
    c0_by_A = load_stage_a_c0()

    stage_b_sha = sha256_file(STAGE_B_SCRIPT)
    if stage_b_sha != EXPECTED_STAGE_B_SCRIPT_SHA256:
        raise ValueError(f"Stage B direct-current source SHA-256 mismatch: {stage_b_sha}")
    if sha256_file(MASS_INPUT) != EXPECTED_MASS_SHA256:
        raise ValueError("P4 exact nuclear-mass input SHA-256 mismatch")
    if sha256_file(ROOT / "v6" / "dmformfactor-V6.m") != EXPECTED_AFH_SHA256:
        raise ValueError("AFH v6.0 source SHA-256 mismatch")
    ff_module = TOPIC / "nucleon_vector_form_factors.py"
    borah_source = ROOT / "research" / "nucleon_form_factor_B_baseline" / "borah2020_lowq2.py"
    borah_parameters = ROOT / "research" / "nucleon_form_factor_B_baseline" / "data" / "borah2020" / "parameters.json"
    borah_covariance = ROOT / "research" / "nucleon_form_factor_B_baseline" / "data" / "borah2020" / "ELZExpFF.dat"
    rate_kinematics = ROOT / "research" / "three_body_calculation_and_validation" / "rate_kinematics.py"
    for source_path, expected, label in (
        (ff_module, EXPECTED_FORM_FACTOR_MODULE_SHA256, "nucleon form-factor module"),
        (borah_source, EXPECTED_BORAH_SOURCE_SHA256, "Borah2020 source"),
        (borah_parameters, EXPECTED_BORAH_PARAMETERS_SHA256, "Borah2020 parameters"),
        (borah_covariance, EXPECTED_BORAH_COVARIANCE_SHA256, "Borah2020 covariance"),
        (rate_kinematics, EXPECTED_RATE_KINEMATICS_SHA256, "canonical rate kinematics module"),
    ):
        actual = sha256_file(source_path)
        if actual != expected:
            raise ValueError(f"{label} SHA-256 mismatch: {actual}")
    legacy_contract_path = TOPIC / "contract_p3c_finite_tN_bound_current.py"
    legacy_contract_bytes = legacy_contract_path.read_bytes()
    legacy_contract_git_sha1 = git_blob_sha1(legacy_contract_bytes)
    if legacy_contract_git_sha1 != EXPECTED_LEGACY_CONTRACTION_GIT_BLOB_SHA1:
        raise ValueError("validated P3C contraction module Git blob SHA-1 mismatch")
    legacy_contract_sha256 = sha256_bytes(legacy_contract_bytes)
    production_csv_sha = sha256_file(STAGE_B_PRODUCTION_CSV)
    if production_csv_sha != EXPECTED_STAGE_B_PRODUCTION_CSV_SHA256:
        raise ValueError("immutable Stage B 81-point production CSV fingerprint mismatch")

    enriched_native = enrich_native_csv(archive_inputs["native_bytes"])
    with tempfile.TemporaryDirectory(prefix="p4fc-existing-t3-") as tmp:
        native_path = Path(tmp) / "native_exact_mass_multipoles_enriched.csv"
        native_path.write_bytes(enriched_native)
        stage_b_payload = stage_b.run(
            order=args.order,
            vertex_energy="external",
            mass_scheme="ame2020",
            native_csv=native_path,
            isotopes=ANCHOR_ISOTOPES,
            recoil_energies=ANCHOR_ER_GEV,
        )

    stage_b_rows = stage_b_payload["results"]
    stage_b_by_point: dict[tuple[int, float], list[dict[str, Any]]] = {}
    for row in stage_b_rows:
        key = (int(row["A"]), round(float(row["ER_GeV"]), 9))
        stage_b_by_point.setdefault(key, []).append(row)
    expected_points = {
        (A, round(ER, 9)) for A in ANCHOR_ISOTOPES for ER in ANCHOR_ER_GEV
    }
    if set(stage_b_by_point) != expected_points:
        raise ValueError(f"Stage B exact-anchor current coverage mismatch: {sorted(stage_b_by_point)}")

    native_rows = read_csv_bytes(archive_inputs["native_bytes"])
    native = {
        (int(row["A"]), round(float(row["ER_GeV"]), 9), int(row["J"]), row["channel"]):
            complex(float(row["amplitude_real"]), float(row["amplitude_imag"]))
        for row in native_rows
    }

    max_fit_residual = max(float(row["max_fit_residual"]) for row in stage_b_rows)
    max_nr_native_residual = max(
        abs(row[key])
        for row in stage_b_rows
        for key in (
            "leading_NR_charge_M_residual",
            "leading_NR_charge_PhiPP_residual",
            "leading_NR_Delta_residual",
            "leading_NR_SigmaPrime_residual",
        )
    )

    tensor_rows = select_anchor_tensor_rows(archive_inputs["tensor_rows"])
    tensor_points = {
        (int(row.get("isotope_A", row.get("A"))), round(float(row.get("ER_GeV", float(row.get("ER_keV", 0.0)) * 1e-6)), 9))
        for row in tensor_rows
    }
    if tensor_points != expected_points:
        raise ValueError(f"archived T3 tensor anchor coverage mismatch: {sorted(tensor_points)}")
    contraction_results = []
    rank_results = []
    max_leading_nr_afh_relative = 0.0
    max_decomposition_absolute_residual = 0.0
    max_decomposition_relative_residual = 0.0
    max_imaginary_absolute_residual = 0.0
    max_imaginary_relative_residual = 0.0
    min_exact_kernel = float("inf")
    max_exact_vs_nr_relative = 0.0
    max_gram_partial_trace_residual = 0.0

    for tensor_index, tensor_row in enumerate(tensor_rows):
        isotope = int(tensor_row.get("isotope_A", tensor_row.get("A")))
        er_value = tensor_row.get("ER_GeV")
        er = float(er_value) if er_value is not None else float(tensor_row["ER_keV"]) * 1e-6
        point_key = (isotope, round(er, 9))
        rows = sorted(stage_b_by_point[point_key], key=lambda row: int(row["J"]))
        JA = float(rows[0]["J_A"])
        mass = float(tensor_row.get("mA_P4_GeV", tensor_row.get("mA_GeV")))
        q = float(
            tensor_row.get(
                "Q_A_P4_GeV",
                tensor_row.get("qA_magnitude_GeV", tensor_row.get("Q_A_GeV", rows[0]["Q_A_GeV"])),
            )
        )
        tA = float(tensor_row.get("tA_P4_GeV2", -2.0 * mass * er))
        if abs(mass - float(rows[0]["mA_GeV"])) > 2e-12:
            raise ValueError(f"P4 mass mismatch in archived T3 tensor for Xe-{isotope}, ER={er}")
        if abs(q - float(rows[0]["Q_A_GeV"])) > 2e-13:
            raise ValueError(f"Q_A mismatch in archived T3 tensor for Xe-{isotope}, ER={er}")
        if abs(tA - float(rows[0]["tN_GeV2"])) > 2e-13:
            raise ValueError(f"t_A/t_N mismatch in archived T3 tensor for Xe-{isotope}, ER={er}")

        L = complex_pairs(tensor_row["L_q_aligned_real_imag"])
        G = complex_pairs(tensor_row["G_q_aligned_real_imag"])
        if L.shape != (4, 4) or G.shape != (4, 4, 4, 4):
            raise ValueError(
                f"unexpected archived tensor dimensions for Xe-{isotope}, ER={er}: "
                f"L={L.shape}, G={G.shape}"
            )
        partial_trace_residual = gram_partial_trace_residual(L, G)
        max_gram_partial_trace_residual = max(max_gram_partial_trace_residual, partial_trace_residual)
        if partial_trace_residual > 1e-12:
            raise ValueError(
                f"archived L tensor does not equal the DM-spin partial trace of G: "
                f"Xe-{isotope}, ER={er}, residual={partial_trace_residual}"
            )

        legacy_contract.C0[isotope] = c0_by_A[isotope]
        exact_rank_parts = []
        nr_rank_parts = []
        afh_rank_parts = []
        for row in rows:
            exact = legacy_contract.build_rank_components(row, isotope, "exact")
            leading_nr = legacy_contract.build_rank_components(row, isotope, "leading_nr")
            exact_parts = contract_components(L, exact, JA)
            nr_parts = contract_components(L, leading_nr, JA)
            afh = legacy_contract.build_afh_reference(row, native, L, JA, q, er, isotope)
            exact_rank_parts.append(exact_parts)
            nr_rank_parts.append(nr_parts)
            afh_rank_parts.append(afh)
            rank_results.append({
                "tensor_index": tensor_index,
                "benchmark": tensor_row.get("benchmark", f"tensor_{tensor_index}"),
                "A": isotope,
                "ER_keV": er * 1e6,
                "J_A": JA,
                "K": int(row["J"]),
                "F1_effective": float(row["F1_effective"]),
                "F2_effective": float(row["F2_effective"]),
                "exact": exact_parts,
                "leading_nr": nr_parts,
                "AFH_reference": {
                    "K_total_AFH_L00_plus_spatial": float(afh["K_total_AFH_L00_plus_spatial"]),
                    "K_total_AFH_four_channel": float(afh["K_total_AFH_four_channel"]),
                    "signed_M_Phi_kernel_AFH": float(afh["signed_M_Phi_kernel_AFH"]),
                    "signed_Delta_Sigma_kernel_AFH": float(afh["signed_Delta_Sigma_kernel_AFH"]),
                },
            })

        exact_total = sum_parts(exact_rank_parts)
        nr_total = sum_parts(nr_rank_parts)
        afh_total = {
            key: float(sum(float(part[key]) for part in afh_rank_parts))
            for key in (
                "K_total_AFH_L00_plus_spatial",
                "K_total_AFH_four_channel",
                "K_charge_AFH",
                "K_spatial_AFH",
                "signed_M_Phi_kernel_AFH",
                "signed_Delta_Sigma_kernel_AFH",
            )
        }
        nr_c0 = c0_by_A[isotope] * nr_total["K_full_re"]
        exact_c0 = c0_by_A[isotope] * exact_total["K_full_re"]
        nr_afh_relative = abs(nr_c0 - afh_total["K_total_AFH_L00_plus_spatial"]) / max(
            abs(afh_total["K_total_AFH_L00_plus_spatial"]), 1e-300
        )
        exact_vs_nr_relative = abs(exact_c0 - nr_c0) / max(abs(nr_c0), 1e-300)
        max_leading_nr_afh_relative = max(max_leading_nr_afh_relative, nr_afh_relative)
        max_decomposition_absolute_residual = max(
            max_decomposition_absolute_residual,
            exact_total["decomposition_absolute_residual"],
            nr_total["decomposition_absolute_residual"],
        )
        max_decomposition_relative_residual = max(
            max_decomposition_relative_residual,
            exact_total["decomposition_relative_residual"],
            nr_total["decomposition_relative_residual"],
        )
        max_imaginary_absolute_residual = max(
            max_imaginary_absolute_residual,
            exact_total["imaginary_residual"],
            nr_total["imaginary_residual"],
        )
        max_imaginary_relative_residual = max(
            max_imaginary_relative_residual,
            exact_total["imaginary_relative_residual"],
            nr_total["imaginary_relative_residual"],
        )
        min_exact_kernel = min(min_exact_kernel, exact_total["K_full_re"])
        max_exact_vs_nr_relative = max(max_exact_vs_nr_relative, exact_vs_nr_relative)

        contraction_results.append({
            "tensor_index": tensor_index,
            "tensor_member": archive_inputs["tensor_member"],
            "tensor_member_sha256": archive_inputs["tensor_sha256"],
            "tensor_payload_status": archive_inputs["tensor_payload"].get("status"),
            "benchmark": tensor_row.get("benchmark", f"tensor_{tensor_index}"),
            "benchmark_base": tensor_row.get("benchmark_base"),
            "A": isotope,
            "J_A": JA,
            "ER_GeV": er,
            "ER_keV": er * 1e6,
            "mA_GeV": mass,
            "Q_A_GeV": q,
            "t_A_GeV2": tA,
            "mchi_GeV": tensor_row.get("mchi_GeV"),
            "ms_GeV": tensor_row.get("ms_GeV"),
            "delta_m_GeV": tensor_row.get("delta_m_GeV"),
            "speed_kms": tensor_row.get("speed_kms"),
            "g2": tensor_row.get("g2"),
            "width_scheme": tensor_row.get("width_scheme", "pole"),
            "mapping": tensor_row.get("mapping"),
            "quadrature": {
                "nw": tensor_row.get("nw"),
                "nz": tensor_row.get("nz"),
                "nphi": tensor_row.get("nphi"),
                "quadrature_points": tensor_row.get("quadrature_points"),
            },
            "archived_tensor_checks": {
                "gram_hermiticity_relative_residual": tensor_row.get("gram_hermiticity_relative_residual"),
                "tensor_hermiticity_relative_residual": tensor_row.get("tensor_hermiticity_relative_residual"),
                "gram_min_eigenvalue": tensor_row.get("gram_min_eigenvalue"),
                "tensor_min_eigenvalue": tensor_row.get("tensor_min_eigenvalue"),
                "dirac_trace_max_relative_residual": (
                    tensor_row.get("dirac_trace_crosscheck", {}).get("max_relative_residual")
                    if isinstance(tensor_row.get("dirac_trace_crosscheck"), dict) else None
                ),
                "L_from_G_partial_trace_relative_residual": partial_trace_residual,
            },
            "C0": c0_by_A[isotope],
            "exact_raw": exact_total,
            "leading_nr_raw": nr_total,
            "AFH_reference": afh_total,
            "normalization": {
                "exact_C0": exact_c0,
                "exact_4C0_diagnostic_only": 4.0 * exact_c0,
                "leading_nr_C0": nr_c0,
                "leading_nr_4C0_diagnostic_only": 4.0 * nr_c0,
                "factor_four_applied_to_production": False,
            },
            "leading_nr_vs_AFH_L00_plus_spatial_relative_difference": nr_afh_relative,
            "exact_vs_leading_nr_relative_difference": exact_vs_nr_relative,
        })

    if not contraction_results:
        raise ValueError("no Stage B current/T3 tensor contractions were produced")
    if max_fit_residual > 1e-10:
        raise ValueError(f"Stage B exact-anchor multipole fit residual failed: {max_fit_residual}")
    if max_nr_native_residual > 1e-10:
        raise ValueError(f"leading-NR exact-anchor Native AFH amplitude match failed: {max_nr_native_residual}")
    if max_leading_nr_afh_relative > 1e-8:
        raise ValueError(f"leading-NR full T3 contraction does not reproduce AFH: {max_leading_nr_afh_relative}")
    if max_gram_partial_trace_residual > 1e-12:
        raise ValueError(f"archived L/G partial-trace identity failed: {max_gram_partial_trace_residual}")
    if max_decomposition_relative_residual > 1e-12:
        raise ValueError(f"full-current kernel decomposition failed: {max_decomposition_relative_residual}")
    if max_imaginary_relative_residual > 1e-12:
        raise ValueError(f"full-current contraction imaginary residual failed: {max_imaginary_relative_residual}")
    if min_exact_kernel < -1e-10:
        raise ValueError(f"negative exact full-current kernel: {min_exact_kernel}")

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    current_json = outdir / "stage_b_exact_anchor_current_o14.json"
    current_csv = outdir / "stage_b_exact_anchor_current_o14.csv"
    contraction_json = outdir / "stage_b_existing_t3_anchor_contractions.json"
    contraction_csv = outdir / "stage_b_existing_t3_anchor_contractions.csv"
    rank_csv = outdir / "stage_b_existing_t3_anchor_rank_contributions.csv"
    validation_path = outdir / "stage_b_existing_t3_anchor_validation_report.json"
    manifest_path = outdir / "stage_b_existing_t3_anchor_run_manifest.json"

    current_payload = {
        "status": "COMPLETED_STAGE_B_EXACT_MASS_ANCHOR_CURRENT",
        "scope": "Stage B finite-tN direct current evaluated only at the four exact-q Native/T3 anchor points; no 567-point Stage B rerun.",
        "quadrature_order_per_axis": args.order,
        "runtime": stage_b_payload.get("runtime", {}),
        "source_sha256": stage_b_payload.get("source_sha256", {}),
        "native_exact_mass_multipoles_member": archive_inputs["native_member"],
        "native_exact_mass_multipoles_sha256": archive_inputs["native_sha256"],
        "results": stage_b_rows,
    }
    current_json.write_text(json.dumps(jsonable(current_payload), indent=2) + "\n", encoding="utf-8")
    write_csv(current_csv, stage_b_rows)

    contraction_payload = {
        "status": "COMPLETED_STAGE_B_CURRENT_CONTRACTED_WITH_EXISTING_P4_EXACT_MASS_T3",
        "scope": "Fixed-speed exact-mass Xe-129/Xe-131 anchor contractions using the archived Pole-width T3 tensors; no T3 phase-space integration or halo fold.",
        "normalization": "Direct current kernels use C0 once. The factor-four columns are diagnostics only and are not used in the production result.",
        "archive_sha256": archive_inputs["archive_sha256"],
        "archive_git_blob_sha1": archive_inputs["archive_git_blob_sha1"],
        "tensor_member": archive_inputs["tensor_member"],
        "tensor_member_sha256": archive_inputs["tensor_sha256"],
        "results": contraction_results,
        "rank_results": rank_results,
    }
    contraction_json.write_text(json.dumps(jsonable(contraction_payload), indent=2) + "\n", encoding="utf-8")
    write_csv(contraction_csv, contraction_results)
    write_csv(rank_csv, rank_results)

    validation = {
        "status": "PASS_STAGE_B_CURRENT_X_EXISTING_T3_ANCHOR_INTEGRATION",
        "scope": "Stage B Borah2020+Sufian2017 finite-tN direct current contracted with already archived P4 exact-mass Pole-width T3 tensors.",
        "coverage": {
            "isotopes": sorted({int(row["A"]) for row in contraction_results}),
            "recoil_energies_keV": sorted({round(float(row["ER_keV"]), 9) for row in contraction_results}),
            "exact_anchor_current_points": len(stage_b_by_point),
            "rank_resolved_current_rows": len(stage_b_rows),
            "existing_T3_tensor_rows_used": len(contraction_results),
            "rank_contraction_rows": len(rank_results),
        },
        "checks": {
            "max_direct_current_multipole_fit_residual": max_fit_residual,
            "max_leading_NR_to_exact_Native_AFH_amplitude_residual": max_nr_native_residual,
            "max_leading_NR_full_T3_to_AFH_relative_difference": max_leading_nr_afh_relative,
            "max_full_current_decomposition_absolute_residual": max_decomposition_absolute_residual,
            "max_full_current_decomposition_relative_residual": max_decomposition_relative_residual,
            "max_full_current_imaginary_absolute_residual": max_imaginary_absolute_residual,
            "max_full_current_imaginary_relative_residual": max_imaginary_relative_residual,
            "max_archived_L_from_G_partial_trace_relative_residual": max_gram_partial_trace_residual,
            "minimum_exact_full_current_raw_kernel": min_exact_kernel,
            "maximum_exact_vs_leading_NR_relative_kernel_difference": max_exact_vs_nr_relative,
            "leading_NR_AFH_tolerance": 1e-8,
            "multipole_fit_tolerance": 1e-10,
            "normalization": "C0 applied once; factor four not applied to direct current kernels",
        },
        "reuse_policy": {
            "Stage_A_native_7isotope_81point": "reused unchanged",
            "Stage_B_7isotope_81point_production": "reused unchanged as immutable reference; only four exact-q anchor points were evaluated for T3 alignment",
            "P3C_P4_exact_mass_T3_archive": "reused unchanged; no integrate_spin_resolved_tensor call is made",
            "historical_quadrature_scans": "not rerun",
        },
    }
    validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")

    manifest = {
        "record_type": "Stage B current x existing exact-mass T3 anchor integration",
        "status": validation["status"],
        "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "repository": "lm846802297-bit/LZ-3-body-for-hdm",
        "pull_request": 143,
        "scope": contraction_payload["scope"],
        "commands": [
            "python research/three_body_full_current_nucleon_match/p4_full_current/stage_b_existing_t3_anchor_integration.py",
        ],
        "inputs": {
            "exact_mass_t3_archive": str(T3_ARCHIVE.relative_to(ROOT)),
            "exact_mass_t3_archive_sha256": archive_inputs["archive_sha256"],
            "exact_mass_t3_archive_git_blob_sha1": archive_inputs["archive_git_blob_sha1"],
            "t3_tensor_member": archive_inputs["tensor_member"],
            "t3_tensor_member_sha256": archive_inputs["tensor_sha256"],
            "native_exact_mass_multipoles_member": archive_inputs["native_member"],
            "native_exact_mass_multipoles_sha256": archive_inputs["native_sha256"],
            "stage_a_native_responses_sha256": sha256_file(STAGE_A_RESPONSE_CSV),
            "stage_b_81point_production_csv_sha256": production_csv_sha,
            "p4_mass_input_sha256": sha256_file(MASS_INPUT),
            "afh_v6_sha256": sha256_file(ROOT / "v6" / "dmformfactor-V6.m"),
            "nucleon_form_factor_module_sha256": sha256_file(ff_module),
            "borah2020_source_sha256": sha256_file(borah_source),
            "borah2020_parameters_sha256": sha256_file(borah_parameters),
            "borah2020_covariance_sha256": sha256_file(borah_covariance),
            "rate_kinematics_sha256": sha256_file(rate_kinematics),
            "validated_p3c_contraction_module_git_blob_sha1": legacy_contract_git_sha1,
            "validated_p3c_contraction_module_sha256": legacy_contract_sha256,
        },
        "source_sha256": {
            "integration_script": sha256_file(Path(__file__).resolve()),
            "stage_b_direct_current_script": stage_b_sha,
            "legacy_validated_contraction_module": sha256_file(TOPIC / "contract_p3c_finite_tN_bound_current.py"),
        },
        "stage_b_exact_anchor_runtime": stage_b_payload.get("runtime", {}),
        "integration_runtime": {
            "started_utc": started_utc,
            "wall_time_seconds": time.perf_counter() - timer,
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "logical_cpu_count": os.cpu_count(),
            "numpy": np.__version__,
            "github_actions_used": False,
        },
        "validation_report": validation,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    output_files = (current_json, current_csv, contraction_json, contraction_csv, rank_csv, validation_path, manifest_path)
    sums_path = outdir / "SHA256SUMS_stage_b_existing_t3_integration.txt"
    sums_path.write_text(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in output_files),
        encoding="utf-8",
    )
    print(json.dumps({
        "status": validation["status"],
        "tensor_rows": len(contraction_results),
        "rank_rows": len(rank_results),
        "max_leading_nr_afh_relative": max_leading_nr_afh_relative,
        "max_exact_vs_leading_nr_relative": max_exact_vs_nr_relative,
        "output_dir": str(outdir),
        "wall_time_seconds": time.perf_counter() - timer,
    }, indent=2))


if __name__ == "__main__":
    main()
