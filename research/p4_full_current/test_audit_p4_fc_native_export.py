import math
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_p4_fc_native_export as audit


class TestP4FCNativeExportAudit(unittest.TestCase):
    def test_expected_native_grid(self):
        self.assertEqual(len(audit.ENERGIES_KEV), 81)
        self.assertEqual(audit.ENERGIES_KEV[0], 40.0)
        self.assertEqual(audit.ENERGIES_KEV[-1], 800.0)

    def test_rank_maps_include_even_even_charge_only(self):
        self.assertEqual(audit.CHARGE_RANKS[128], (0,))
        self.assertEqual(audit.SPATIAL_RANKS[128], ())
        self.assertEqual(audit.SPATIAL_RANKS[129], (1,))
        self.assertEqual(audit.SPATIAL_RANKS[131], (1, 2, 3))

    def test_response_reconstruction_diagonal(self):
        multipoles = {(128, 40.0, 0, "M"): 2.0 + 0.0j}
        y = 0.1
        # Native ResponseNuclear's E^(-2y) cancels the E^y in AFH Response;
        # reconstruction must not apply a second exponential.
        expected = 4.0 / 4.0
        self.assertAlmostEqual(audit.reconstruct_response(multipoles, 128, 40.0, y, "M"), expected)

    def test_response_reconstruction_signed_interference(self):
        multipoles = {
            (129, 40.0, 0, "M"): 2.0 + 0.0j,
            (129, 40.0, 0, "PhiPP"): -3.0 + 0.0j,
        }
        y = 0.0
        self.assertEqual(audit.reconstruct_response(multipoles, 129, 40.0, y, "MPhiPP"), -1.5)


if __name__ == "__main__":
    unittest.main()
