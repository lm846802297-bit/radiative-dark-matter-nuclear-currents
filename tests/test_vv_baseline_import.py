import csv
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "benchmarks" / "vv_recovery" / "paper_fig3_vector_natural_xenon.json"
AFH = ROOT / "nuclear" / "afh" / "data" / "native_v6_2026-11-22"
IMPORT_META = ROOT / "results" / "provenance" / "VV_BASELINE_IMPORT_2026-11-23.json"
ISOTOPES = (128, 129, 130, 131, 132, 134, 136)


def scalar_values(value):
    if isinstance(value, dict):
        for item in value.values():
            yield from scalar_values(item)
    elif isinstance(value, list):
        for item in value:
            yield from scalar_values(item)
    else:
        yield value


class VVBaselineImportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.benchmark = json.loads(BENCHMARK.read_text(encoding="utf-8"))
        cls.report = (AFH / "REPORT.md").read_text(encoding="utf-8")
        cls.import_meta = json.loads(IMPORT_META.read_text(encoding="utf-8"))

    def test_paper_fig3_benchmark_values_are_present(self):
        values = list(scalar_values(self.benchmark))
        for expected in (1000, 5, 300, 10, 0.1, 1.725, 0.04857, 2.84):
            self.assertTrue(
                any(isinstance(v, (int, float)) and abs(float(v) - expected) < 1e-12 for v in values),
                expected,
            )

    def test_native_afh_seven_isotope_coverage(self):
        total = 0
        for A in ISOTOPES:
            path = AFH / "native" / f"afh_v6_xe{A}_eight_response_native.csv"
            with path.open(newline="", encoding="utf-8-sig") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 2401, f"Xe-{A}")
            self.assertGreaterEqual(len(rows[0]), 10)
            total += len(rows)
        self.assertEqual(total, 16807)

    def test_native_spatial_rank_rows(self):
        path = AFH / "native" / "spatial_multipoles" / "afh_v6_xe129_xe131_spatial_multipoles_native.csv"
        with path.open(newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 9604)
        ranks = {(int(r["A"]), int(r["multipole_J"])) for r in rows}
        self.assertEqual({j for a, j in ranks if a == 129}, {1})
        self.assertEqual({j for a, j in ranks if a == 131}, {1, 2, 3})

    def test_report_preserves_current_validation_boundary(self):
        validation = self.import_meta["native_afh_validation"]
        self.assertIn("full_current_complete=false", self.report)
        self.assertIn("direct covariant bound-nucleus P3-C amplitude gate", self.report)
        self.assertEqual(validation["response_rows"], 16807)
        self.assertEqual(validation["max_native_imaginary_residual"], 0)
        self.assertFalse(validation["full_current_complete"])
        self.assertFalse(validation["physically_validated_full_current"])


if __name__ == "__main__":
    unittest.main()
