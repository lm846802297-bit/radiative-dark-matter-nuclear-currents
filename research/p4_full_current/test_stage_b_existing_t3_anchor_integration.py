import csv
import io
import unittest

import numpy as np

import stage_b_existing_t3_anchor_integration as integration


class TestStageBExistingT3Integration(unittest.TestCase):
    def test_git_blob_sha1_matches_git_object_format(self):
        self.assertEqual(
            integration.git_blob_sha1(b"test"),
            "30d74d258442c7c65512eafab474568dd706c430",
        )

    def test_complex_pair_conversion(self):
        value = integration.complex_pairs([[[1.0, 2.0], [3.0, -4.0]]])
        self.assertEqual(value.shape, (1, 2))
        self.assertEqual(value[0, 0], 1.0 + 2.0j)
        self.assertEqual(value[0, 1], 3.0 - 4.0j)

    def test_native_csv_enrichment_adds_exact_tA(self):
        source = (
            "A,ER_GeV,mA_GeV,Q_A_GeV,y,J,channel,amplitude_real,amplitude_imag\n"
            "129,0.0002,120.04644920255713,0.21913151229575095,0.1,0,M,1,0\n"
        ).encode()
        enriched = integration.enrich_native_csv(source).decode()
        row = next(csv.DictReader(io.StringIO(enriched)))
        self.assertAlmostEqual(float(row["ER_keV"]), 200.0, places=12)
        self.assertAlmostEqual(
            float(row["t_A_GeV2"]),
            -2.0 * 120.04644920255713 * 0.0002,
            places=16,
        )

    def test_anchor_tensor_filter_uses_exact_p4_mass_and_pole_width(self):
        row = {
            "isotope_A": 129,
            "ER_GeV": 0.0002,
            "mA_P4_GeV": integration.EXPECTED_MASSES_GEV[129],
            "width_scheme": "pole",
            "L_q_aligned_real_imag": np.zeros((4, 4, 2)).tolist(),
            "G_q_aligned_real_imag": np.zeros((4, 4, 4, 4, 2)).tolist(),
        }
        self.assertEqual(integration.select_anchor_tensor_rows([row]), [row])

    def test_anchor_tensor_filter_rejects_historical_Au_mass(self):
        row = {
            "isotope_A": 129,
            "ER_GeV": 0.0002,
            "mA_GeV": 129.0 * 0.93149410242,
            "width_scheme": "pole",
            "L_q_aligned_real_imag": np.zeros((4, 4, 2)).tolist(),
            "G_q_aligned_real_imag": np.zeros((4, 4, 4, 4, 2)).tolist(),
        }
        self.assertEqual(integration.select_anchor_tensor_rows([row]), [])

    def test_gram_partial_trace_reconstructs_lorentz_tensor(self):
        G = np.zeros((4, 4, 4, 4), dtype=complex)
        for mu in range(4):
            for nu in range(4):
                for a in range(4):
                    G[mu, a, nu, a] = (mu + 1) * (nu + 1) / (a + 1)
        L = np.einsum("mana->mn", G)
        self.assertEqual(integration.gram_partial_trace_residual(L, G), 0.0)

    def test_gram_partial_trace_detects_mismatch(self):
        G = np.zeros((4, 4, 4, 4), dtype=complex)
        L = np.eye(4, dtype=complex)
        self.assertGreater(integration.gram_partial_trace_residual(L, G), 0.1)

    def test_full_kernel_decomposition(self):
        n = 2
        z = np.zeros((4, n, n), dtype=complex)
        charge_id = z.copy()
        charge_spin = z.copy()
        spatial_id = z.copy()
        spatial_spin = z.copy()
        charge_id[0, 0, 0] = 2.0
        charge_spin[0, 1, 1] = -0.5
        spatial_id[3, 0, 1] = 0.25
        spatial_spin[1, 1, 0] = 0.75
        current = {
            "charge_identity": charge_id,
            "charge_spin": charge_spin,
            "spatial_identity": spatial_id,
            "spatial_spin": spatial_spin,
            "charge": charge_id + charge_spin,
            "spatial": spatial_id + spatial_spin,
            "full": charge_id + charge_spin + spatial_id + spatial_spin,
            "delta_LJ": spatial_id,
            "sigma_LJ": spatial_spin,
        }
        parts = integration.contract_components(np.eye(4, dtype=complex), current, 0.5)
        self.assertLess(parts["decomposition_absolute_residual"], 1e-12)
        self.assertLess(parts["imaginary_residual"], 1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
