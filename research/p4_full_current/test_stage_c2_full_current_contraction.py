import importlib.util
import unittest
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("stage_c2_full_current_contraction.py")
spec = importlib.util.spec_from_file_location("stage_c2", MODULE_PATH)
stage_c2 = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = stage_c2
spec.loader.exec_module(stage_c2)


class StageC2ContractionTests(unittest.TestCase):
    def test_bilinear_matches_independent_loop_for_target_spins(self):
        rng = np.random.default_rng(20261124)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            L = x @ x.conjugate().T
            A = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            B = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            got = stage_c2.bilinear(L, A, B, JA)
            ref = stage_c2.bilinear_explicit(L, A, B, JA)
            self.assertLess(abs(got - ref) / max(abs(got), abs(ref), 1.0), 1e-13)

    def _row(self, JA, K):
        row = {"J_A": str(JA), "J": str(K)}
        names = [
            "charge_identity", "charge_spin",
            "vector_identity_LJm1", "vector_identity_LJ", "vector_identity_LJp1",
            "vector_spin_LJm1", "vector_spin_LJ", "vector_spin_LJp1",
        ]
        for rep in ("exact", "leading_nr"):
            for name in names:
                row[f"{rep}_{name}_re"] = 0.0
                row[f"{rep}_{name}_im"] = 0.0
        return row

    def test_spin_zero_rank0_L1_longitudinal_current_is_retained(self):
        row = self._row(0.0, 0)
        row["exact_vector_identity_LJp1_re"] = 1.0
        row["leading_nr_vector_identity_LJp1_re"] = 1.0
        current = stage_c2.build_rank_components(row, "exact")
        self.assertEqual(current["full"].shape, (4, 1, 1))
        self.assertGreater(abs(current["spatial_identity_LJp1"][3, 0, 0]), 0.0)
        self.assertAlmostEqual(abs(current["spatial_identity_LJp1"][1, 0, 0]), 0.0, places=14)
        self.assertAlmostEqual(abs(current["spatial_identity_LJp1"][2, 0, 0]), 0.0, places=14)
        L = np.diag([1.0, 1.0, 1.0, 1.0]).astype(complex)
        result = stage_c2.contract_components(L, current, 0.0)
        self.assertGreater(result["K_spatial_re"], 0.0)
        self.assertAlmostEqual(result["K_transverse_xy_re"], 0.0, places=14)
        self.assertAlmostEqual(result["K_longitudinal_transverse_signed_cross"], 0.0, places=14)

    def test_rank0_L1_is_longitudinal_for_all_target_spins(self):
        for JA in (0.0, 0.5, 1.5):
            row = self._row(JA, 0)
            row["exact_vector_identity_LJp1_re"] = 1.0
            current = stage_c2.build_rank_components(row, "exact")
            self.assertGreater(abs(current["spatial_identity_LJp1"][3]).max(), 0.0)
            self.assertAlmostEqual(abs(current["spatial_identity_LJp1"][1]).max(), 0.0, places=14)
            self.assertAlmostEqual(abs(current["spatial_identity_LJp1"][2]).max(), 0.0, places=14)

    def test_spinful_dimensions_follow_target_spin(self):
        row_half = self._row(0.5, 1)
        row_three_half = self._row(1.5, 1)
        self.assertEqual(stage_c2.build_rank_components(row_half, "exact")["full"].shape, (4, 2, 2))
        self.assertEqual(stage_c2.build_rank_components(row_three_half, "exact")["full"].shape, (4, 4, 4))

    def test_physical_afh_normalization_applies_factor_four_once(self):
        raw = 2.0
        c0 = 1.0001
        raw_matching = c0 * raw
        physical = stage_c2.physical_kernel(raw, c0)
        self.assertAlmostEqual(raw_matching, 2.0002, places=12)
        self.assertAlmostEqual(physical, 8.0008, places=12)
        self.assertAlmostEqual(physical / raw_matching, 4.0, places=14)

    def test_current_component_decomposition_is_exact(self):
        rng = np.random.default_rng(18)
        JA = 0.5
        n = 2
        current = {}
        for key in (
            "charge_identity", "charge_spin", "spatial_identity_LJm1", "spatial_identity_LJ",
            "spatial_identity_LJp1", "spatial_spin_LJm1", "spatial_spin_LJ", "spatial_spin_LJp1",
        ):
            current[key] = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
        current["charge"] = current["charge_identity"] + current["charge_spin"]
        current["spatial_identity"] = sum(current[f"spatial_identity_{s}"] for s in ("LJm1", "LJ", "LJp1"))
        current["spatial_spin"] = sum(current[f"spatial_spin_{s}"] for s in ("LJm1", "LJ", "LJp1"))
        current["spatial"] = current["spatial_identity"] + current["spatial_spin"]
        current["full"] = current["charge"] + current["spatial"]
        current["delta_LJ"] = current["spatial_identity_LJ"]
        current["sigma_LJ"] = current["spatial_spin_LJ"]
        x = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        L = x @ x.conjugate().T
        result = stage_c2.contract_components(L, current, JA)
        self.assertLess(result["decomposition_relative_residual"], 1e-13)
        self.assertLess(result["decomposition_LT_relative_residual"], 1e-13)
        self.assertLess(result["imaginary_relative_residual"], 1e-13)


    def test_conserved_current_projects_to_lc(self):
        rng = np.random.default_rng(20261125)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 7)) + 1j * rng.normal(size=(4, 7))
            L = x @ x.conjugate().T
            current = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            alpha = 0.0047
            current[3] = alpha * current[0]
            got = stage_c2.contract_afh_continuity_projection(
                L, {"full": current}, JA, alpha
            )
            self.assertLess(abs(got["K_longitudinal_0z_nonconserved_residual_raw"]), 1e-12)
            self.assertLess(got["nonconserved_formula_relative_residual"], 1e-12)
            self.assertLess(got["nonconserved_formula_conditioned_residual"], 1e-12)
            self.assertLess(got["projection_decomposition_relative_residual"], 1e-12)
            lc = np.real(L[0, 0] - alpha * (L[0, 3] + L[3, 0]) + alpha * alpha * L[3, 3])
            expected = lc * np.sum(np.abs(current[0]) ** 2) / (2.0 * JA + 1.0)
            self.assertLess(abs(got["K_LC_longitudinal_projection_raw"] - expected), 1e-11 * max(abs(expected), 1.0))

    def test_nonconserved_longitudinal_residual_matches_independent_formula(self):
        rng = np.random.default_rng(20261126)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 9)) + 1j * rng.normal(size=(4, 9))
            L = x @ x.conjugate().T
            current = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            alpha = 0.011
            rho, jz = current[0], current[3]
            d = jz - alpha * rho
            lower = np.stack((rho, -jz))
            block = L[np.ix_([0, 3], [0, 3])]
            k_full_0z = np.real(np.einsum("ab,afi,bfi->", block, lower, lower.conjugate())) / (2.0 * JA + 1.0)
            lc = np.real(L[0, 0] - alpha * (L[0, 3] + L[3, 0]) + alpha * alpha * L[3, 3])
            k_lc = lc * np.sum(np.abs(rho) ** 2) / (2.0 * JA + 1.0)
            expected_residual = (
                -(L[0, 3] - alpha * L[3, 3]) * rho * d.conjugate()
                -(L[3, 0] - alpha * L[3, 3]) * d * rho.conjugate()
                + L[3, 3] * np.abs(d) ** 2
            )
            expected_residual = np.real(np.sum(expected_residual)) / (2.0 * JA + 1.0)
            got = stage_c2.contract_afh_continuity_projection(
                L, {"full": current}, JA, alpha
            )
            self.assertLess(abs(got["K_longitudinal_0z_nonconserved_residual_raw"] - (k_full_0z-k_lc)), 1e-11 * max(abs(k_full_0z), abs(k_lc), 1.0))
            self.assertLess(abs(got["K_longitudinal_0z_nonconserved_residual_raw"] - expected_residual), 1e-11 * max(abs(expected_residual), 1.0))
            self.assertLess(got["nonconserved_formula_relative_residual"], 1e-12)
            self.assertLess(got["nonconserved_formula_conditioned_residual"], 1e-12)
            self.assertLess(got["projection_decomposition_relative_residual"], 1e-12)

    def test_afh_six_response_projection_and_all_residuals_reconstruct_full_kernel(self):
        rng = np.random.default_rng(20261127)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 8)) + 1j * rng.normal(size=(4, 8))
            L = x @ x.conjugate().T
            current = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            got = stage_c2.contract_afh_continuity_projection(
                L, {"full": current}, JA, 0.006
            )
            reconstructed = (
                got["K_AFH6_projection_raw"]
                + got["K_projected_LT_cross_raw"]
                + got["K_longitudinal_0z_nonconserved_residual_raw"]
                + got["K_longitudinal_transverse_nonconserved_residual_raw"]
            )
            self.assertLess(abs(got["K_full_raw"] - reconstructed), 1e-12 * max(abs(got["K_full_raw"]), 1.0))
            self.assertLess(got["projection_decomposition_relative_residual"], 1e-12)


    def test_proca_numerator_matches_direct_covariant_amplitude(self):
        rng = np.random.default_rng(20261125)
        metric = np.diag(stage_c2.METRIC)
        for _ in range(80):
            ER = rng.uniform(1.0e-5, 8.0e-4)
            QA = rng.uniform(0.02, 0.45)
            mV = rng.uniform(1.0, 15.0)
            T = rng.normal(size=4) + 1j * rng.normal(size=4)
            J = rng.normal(size=4) + 1j * rng.normal(size=4)
            L = np.outer(T, T.conjugate())
            L_proca = stage_c2.apply_proca_numerator(L, ER, QA, mV)
            J_matrix = J.reshape(4, 1, 1)
            got = stage_c2.bilinear(L_proca, J_matrix, J_matrix, 0.0).real

            q_upper = np.asarray([ER, 0.0, 0.0, QA])
            q_lower = stage_c2.METRIC * q_upper
            proca_lower = metric - np.outer(q_lower, q_lower) / (mV * mV)
            amplitude = T @ proca_lower @ J
            expected = abs(amplitude) ** 2
            self.assertLess(abs(got - expected) / max(abs(got), abs(expected), 1.0), 2e-13)

    def test_proca_dressed_tensor_preserves_hermiticity_psd_and_diagram_sum(self):
        rng = np.random.default_rng(20261126)
        for _ in range(30):
            before_amplitude = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            after_amplitude = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            L_before = before_amplitude @ before_amplitude.conjugate().T
            L_after = after_amplitude @ after_amplitude.conjugate().T
            L_interference = (
                before_amplitude @ after_amplitude.conjugate().T
                + after_amplitude @ before_amplitude.conjugate().T
            )
            L_total = (before_amplitude + after_amplitude) @ (
                before_amplitude + after_amplitude
            ).conjugate().T

            args = (2.0e-4, 0.23, 10.0)
            dressed_total = stage_c2.apply_proca_numerator(L_total, *args)
            dressed_before = stage_c2.apply_proca_numerator(L_before, *args)
            dressed_after = stage_c2.apply_proca_numerator(L_after, *args)
            dressed_interference = stage_c2.apply_proca_numerator(L_interference, *args)
            reconstructed = dressed_before + dressed_after + dressed_interference

            self.assertLess(
                np.linalg.norm(dressed_total - dressed_total.conjugate().T)
                / max(np.linalg.norm(dressed_total), 1.0),
                2e-13,
            )
            self.assertGreaterEqual(
                np.linalg.eigvalsh(0.5 * (dressed_total + dressed_total.conjugate().T)).min(),
                -2e-12,
            )
            self.assertLess(
                np.linalg.norm(dressed_total - reconstructed)
                / max(np.linalg.norm(dressed_total), 1.0),
                2e-13,
            )

    def test_zero_transfer_proca_numerator_reduces_to_metric_vertex(self):
        rng = np.random.default_rng(20261127)
        x = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        L = x @ x.conjugate().T
        got = stage_c2.apply_proca_numerator(L, 0.0, 0.0, 10.0)
        np.testing.assert_allclose(got, L, rtol=0.0, atol=1e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)
