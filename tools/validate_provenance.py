#!/usr/bin/env python3
"""Validate the repository-level SHA256SUMS.txt manifest."""
from __future__ import annotations
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "provenance" / "SHA256SUMS.txt"

def main() -> int:
    if not MANIFEST.is_file():
        print("missing checksum manifest", file=sys.stderr)
        return 2
    seen: set[str] = set()
    failures: list[str] = []
    checked = 0
    for line_number, raw in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            failures.append(f"line {line_number}: malformed entry")
            continue
        expected, relative = parts
        if len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected.lower()):
            failures.append(f"line {line_number}: invalid SHA-256 for {relative}")
            continue
        if relative in seen:
            failures.append(f"line {line_number}: duplicate path {relative}")
            continue
        seen.add(relative)
        path = (ROOT / relative).resolve()
        if ROOT not in path.parents:
            failures.append(f"line {line_number}: path escapes repository: {relative}")
            continue
        if not path.is_file():
            failures.append(f"line {line_number}: missing file: {relative}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected.lower():
            failures.append(f"line {line_number}: hash mismatch: {relative}")
        else:
            checked += 1
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"PASS: verified {checked} files")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
