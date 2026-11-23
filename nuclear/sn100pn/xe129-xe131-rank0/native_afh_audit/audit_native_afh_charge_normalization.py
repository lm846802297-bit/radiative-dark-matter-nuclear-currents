"""Audit exported Native AFH W_M00 at q_A=0 without altering the data."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ISOTOPES = {128: 0, 129: 1, 130: 0, 131: 3, 132: 0, 134: 0, 136: 0}
REL_TOL = 2.5e-4  # archival regression window, not a physical error estimate
BASE_COMMIT = "5204f1a78b32b6ff5b423cc5d278d781e88add5e"

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def collect(folder):
    manifest = folder / "SHA256SUMS.txt"
    pins = {line.split(maxsplit=1)[1].lstrip(" *"): line.split()[0]
            for line in manifest.read_text().splitlines() if line.strip()}
    rows = []
    for A, twice_J in ISOTOPES.items():
        relpath = f"native/afh_v6_xe{A}_eight_response_native.csv"
        csv_path = folder / relpath
        file_sha = digest(csv_path)
        if pins[relpath] != file_sha:
            raise ValueError(f"Input SHA-256 mismatch: Xe-{A}")
        with csv_path.open(newline="", encoding="utf-8-sig") as f:
            zero = [r for r in csv.DictReader(f) if float(r["q_GeV"]) == 0]
        if len(zero) != 1 or int(zero[0]["A"]) != A or float(zero[0]["y"]) != 0:
            raise ValueError(f"Zero-recoil row lock failed: Xe-{A}")
        W = float(zero[0]["W_M00_AFH"])
        if not math.isfinite(W) or W <= 0 or float(zero[0]["Im_W_M00_AFH"]) != 0:
            raise ValueError(f"Unexpected zero-recoil coherent response: Xe-{A}")
        factor = 4 * math.pi / (twice_J + 1)
        R_raw = factor * W
        R_equal = 4 * R_raw
        delta = R_equal / A**2 - 1
        if abs(delta) >= REL_TOL:
            raise ValueError(f"Xe-{A} coherent response outside archived-count window: {delta}")
        rows.append({"A": A, "twice_J": twice_J, "sha256": file_sha,
                     "W_native_at_zero": W, "spin_angle_factor": factor,
                     "R_raw_at_zero": R_raw, "ideal_R_raw": A**2 / 4,
                     "R_equal_pn_at_zero": R_equal, "ideal_R_equal_pn": A**2,
                     "fractional_count_residual": delta})
    return {"baseline_commit": BASE_COMMIT, "input_manifest_sha256": digest(manifest),
            "scope": "zero-momentum coherent charge normalization of published Native AFH CSVs only",
            "relative_regression_tolerance": REL_TOL, "isotopes": rows,
            "max_abs_fractional_count_residual": max(abs(r["fractional_count_residual"]) for r in rows)}

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("native_directory", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    a.output.write_text(json.dumps(collect(a.native_directory), indent=2, allow_nan=False) + "\n")
