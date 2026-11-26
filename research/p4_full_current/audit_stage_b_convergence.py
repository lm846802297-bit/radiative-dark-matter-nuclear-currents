#!/usr/bin/env python3
"""Compare Stage B direct-current multipoles across Gauss-Hermite orders."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

ORDERS = (10, 14, 18, 22)
ISOTOPES = (128, 129, 130, 131, 132, 134, 136)
ANCHORS_KEV = (40.0, 420.0, 800.0)
CHARGE_FIELDS = ("charge_identity", "charge_spin")
VECTOR_FIELDS = tuple(
    f"vector_{kind}_{suffix}"
    for kind in ("identity", "spin")
    for suffix in ("LJm1", "LJ", "LJp1")
)
COMPONENTS = CHARGE_FIELDS + VECTOR_FIELDS
REPRESENTATIONS = ("exact", "leading_nr")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cval(row: dict, name: str) -> complex:
    return complex(float(row[name + "_re"]), float(row[name + "_im"]))


def load(path: Path) -> tuple[dict, dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = {
        (int(row["A"]), round(float(row["ER_keV"]), 9), int(row["J"])): row
        for row in payload["results"]
    }
    if len(rows) != len(payload["results"]):
        raise ValueError(f"duplicate rank rows in {path}")
    return payload, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--o10", type=Path, required=True)
    parser.add_argument("--o14", type=Path, required=True)
    parser.add_argument("--o18", type=Path, required=True)
    parser.add_argument("--o22", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    files = {10: args.o10.resolve(), 14: args.o14.resolve(), 18: args.o18.resolve(), 22: args.o22.resolve()}
    loaded = {order: load(path) for order, path in files.items()}
    anchor_keys = set()
    for A in ISOTOPES:
        for ER in ANCHORS_KEV:
            # Stage A grid is 40 + 9.5*n keV; 420 and 800 are exact grid points.
            for J in (0, 1, 2, 3):
                anchor_keys.add((A, round(ER, 9), J))
    # Retain only physically present ranks, taken from the actual order-14 result.
    ref_rows = loaded[14][1]
    anchor_keys = {key for key in anchor_keys if key in ref_rows}
    expected_points = {(A, ER) for A in ISOTOPES for ER in ANCHORS_KEV}
    for order, (payload, rows) in loaded.items():
        if tuple(payload.get("isotopes", ())) != ISOTOPES:
            raise ValueError(f"order {order} run does not contain all seven isotopes")
        found_points = {(A, E) for A, E, _ in rows}
        if not expected_points.issubset(found_points):
            raise ValueError(f"order {order} run is missing one or more convergence anchors")
        if int(payload["quadrature_order_per_axis"]) != order:
            raise ValueError(f"order metadata mismatch for order {order}")

    comparisons = {}
    maximum = {}
    for lower, upper in ((10, 14), (14, 18), (18, 22)):
        label = f"order_{lower}_to_{upper}"
        max_abs = 0.0
        max_rel = 0.0
        worst = None
        compared = 0
        pair_pass = True
        for key in sorted(anchor_keys):
            row_a = loaded[lower][1][key]
            row_b = loaded[upper][1][key]
            for rep in REPRESENTATIONS:
                for component in COMPONENTS:
                    field = f"{rep}_{component}"
                    a = cval(row_a, field)
                    b = cval(row_b, field)
                    absolute = abs(a - b)
                    relative = absolute / max(abs(a), abs(b), 1e-300)
                    compared += 1
                    if absolute > 1e-10 + 1e-5 * max(abs(a), abs(b)):
                        pair_pass = False
                    if absolute > max_abs:
                        max_abs = absolute
                        worst = {
                            "A": key[0], "ER_keV": key[1], "J": key[2],
                            "representation": rep, "component": component,
                            "lower": [a.real, a.imag], "upper": [b.real, b.imag],
                            "absolute_difference": absolute,
                            "relative_difference": relative,
                        }
                    max_rel = max(max_rel, relative)
        comparisons[label] = {
            "rank_rows": len(anchor_keys),
            "complex_components_compared": compared,
            "max_absolute_difference": max_abs,
            "max_relative_difference": max_rel,
            "pass": pair_pass,
            "worst_component": worst,
        }
        maximum[label] = max_abs

    final = comparisons["order_18_to_22"]
    accepted = final["pass"]
    report = {
        "status": "PASS_STAGE_B_GAUSS_HERMITE_CONVERGENCE" if accepted else "FAIL_STAGE_B_GAUSS_HERMITE_CONVERGENCE",
        "scope": "All seven isotopes at 40, 420, and 800 keV; direct and leading-NR charge plus all vector L=K-1,K,K+1 multipoles.",
        "quadrature_orders": list(ORDERS),
        "acceptance": {
            "final_order_18_to_22": "max absolute difference <= 1e-10 + 1e-5*max amplitude scale",
            "pass": accepted,
        },
        "comparisons": comparisons,
        "source_sha256": {str(order): sha256(path) for order, path in files.items()},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
