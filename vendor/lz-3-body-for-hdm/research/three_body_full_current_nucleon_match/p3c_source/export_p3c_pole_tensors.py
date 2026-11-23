#!/usr/bin/env python3
"""Recompute the accepted Pole-width T3 Gram tensors for Xe-129 and Xe-131."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
OUT = Path(__file__).resolve().parents[1] / "pole_width_t3_tensor_xe129_xe131.json"
sys.path.insert(0, str(ROOT / "research" / "three_body_full_current_nucleon_match"))
sys.path.insert(0, str(ROOT / "research" / "vector_current_response_audit"))

import three_body_phase_space_tensor as phase  # noqa: E402
import vector_current_response_audit as vca  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    results = []
    for A in (129, 131):
        mA = A * phase.MASS_U_GEV
        for case in vca.DEFAULT_CASES:
            result = phase.integrate_spin_resolved_tensor(
                mchi=phase.MCHI,
                ms=case["ms"],
                delta_m=case["delta_m"],
                mA=mA,
                ER=case["ER"],
                speed_kms=case["speed_kms"],
                g2=case["g2"],
                width_scheme="pole",
                mapping=case["mapping"],
                nw=case["nw"],
                nz=case["nz"],
                nphi=case["nphi"],
                crosscheck_dirac_trace=True,
            )
            if not result.get("phase_space_open", False):
                raise RuntimeError(f"closed benchmark: Xe-{A} {case['benchmark']}")
            result["benchmark"] = f"{case['benchmark']}_Xe{A}"
            result["benchmark_base"] = case["benchmark"]
            result["isotope_A"] = A
            result["width_scheme"] = "pole"
            results.append(result)
    max_trace = max(x["dirac_trace_crosscheck"]["max_relative_residual"] for x in results)
    min_gram = min(x["gram_min_eigenvalue"] for x in results)
    min_tensor = min(x["tensor_min_eigenvalue"] for x in results)
    max_hermiticity = max(x["gram_hermiticity_relative_residual"] for x in results)
    acceptance = {
        "dirac_trace_max_relative_residual": max_trace,
        "dirac_trace_tolerance": 2e-9,
        "gram_minimum_eigenvalue": min_gram,
        "tensor_minimum_eigenvalue": min_tensor,
        "minimum_eigenvalue_tolerance": -1e-10,
        "gram_hermiticity_max_relative_residual": max_hermiticity,
        "hermiticity_tolerance": 1e-12,
    }
    acceptance["pass"] = (
        max_trace < acceptance["dirac_trace_tolerance"]
        and min_gram > acceptance["minimum_eigenvalue_tolerance"]
        and min_tensor > acceptance["minimum_eigenvalue_tolerance"]
        and max_hermiticity < acceptance["hermiticity_tolerance"]
    )
    output = {
        "status": "PASS_SPIN_RESOLVED_DM_TENSOR" if acceptance["pass"] else "FAIL_SPIN_RESOLVED_DM_TENSOR",
        "formal_width_scheme": "pole",
        "isotopes": [129, 131],
        "benchmarks": results,
        "acceptance": acceptance,
        "source_sha256": {
            "phase_space_tensor": sha256(ROOT / "research" / "three_body_full_current_nucleon_match" / "three_body_phase_space_tensor.py"),
            "vector_current_response_audit": sha256(ROOT / "research" / "vector_current_response_audit" / "vector_current_response_audit.py"),
            "isotope_expansion_script": sha256(Path(__file__).resolve()),
        },
        "scope": "The same four accepted benchmark parameter sets, re-integrated at the isotope-specific mA for Xe-129 and Xe-131.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": output["status"],
        "benchmark_count": len(results),
        "isotopes": output["isotopes"],
        "dirac_trace_max_relative_residual": max(
            x["dirac_trace_crosscheck"]["max_relative_residual"] for x in results
        ),
        "acceptance": acceptance,
        "output": str(OUT),
    }, indent=2))
    if not acceptance["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
