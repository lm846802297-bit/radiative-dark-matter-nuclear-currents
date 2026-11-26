import unittest
from unittest.mock import patch

import numpy as np

import stage_c1_full_lorentz_tensor_grid as c1


def pairs(value):
    value = np.asarray(value, dtype=complex)
    return np.stack((value.real, value.imag), axis=-1).tolist()


class TestStageC1FullLorentzTensorGrid(unittest.TestCase):
    def test_speed_nodes_cover_interval_and_weights(self):
        speeds, weights, densities = c1.make_speed_nodes(
            250.0, 776.0, 24, c1.halo.shm_truncated_pdf
        )
        self.assertEqual(speeds.shape, (24,))
        self.assertEqual(weights.shape, (24,))
        self.assertEqual(densities.shape, (24,))
        self.assertGreater(float(speeds.min()), 250.0)
        self.assertLess(float(speeds.max()), 776.0)
        self.assertAlmostEqual(float(weights.sum()), 526.0, places=11)
        self.assertTrue(np.isfinite(densities).all())
        self.assertTrue((densities >= 0.0).all())

    def test_tensor_diagnostics_reconstruct_partial_trace_and_decomposition(self):
        rng = np.random.default_rng(741)
        X = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        L = X @ X.conj().T
        G = np.zeros((4, 4, 4, 4), dtype=complex)
        for a in range(4):
            G[:, a, :, a] = L / 4.0
        result = {
            "L_q_aligned_real_imag": pairs(L),
            "L_before_real_imag": pairs(0.25 * L),
            "L_after_real_imag": pairs(0.25 * L),
            "L_interference_real_imag": pairs(0.50 * L),
            "G_q_aligned_real_imag": pairs(G),
            "q_frame_rotation_max_orthogonality_residual": 0.0,
            "q_frame_rotation_max_determinant_residual": 0.0,
            "quadrature_points": 336,
        }
        diag = c1.tensor_diagnostics(result)
        self.assertLess(diag["hermiticity_relative"], 1e-13)
        self.assertGreaterEqual(diag["psd_min_eigenvalue_relative"], -1e-13)
        self.assertLess(diag["gram_partial_trace_relative"], 1e-13)
        self.assertLess(diag["decomposition_relative"], 1e-13)
        self.assertEqual(diag["L"].shape, (4, 4))

    def test_tensor_diagnostics_rejects_nonfinite_tensor(self):
        L = np.eye(4, dtype=complex)
        bad = L.copy()
        bad[0, 0] = np.nan
        result = {
            "L_q_aligned_real_imag": pairs(bad),
            "L_before_real_imag": pairs(L),
            "L_after_real_imag": pairs(L),
            "L_interference_real_imag": pairs(L),
            "G_q_aligned_real_imag": pairs(np.zeros((4, 4, 4, 4), complex)),
            "q_frame_rotation_max_orthogonality_residual": 0.0,
            "q_frame_rotation_max_determinant_residual": 0.0,
            "quadrature_points": 1,
        }
        with self.assertRaises(ValueError):
            c1.tensor_diagnostics(result)

    def test_build_tasks_uses_both_halo_grids(self):
        isotope_data = {129: {"mA_GeV": 120.04644920255713}}
        tasks, coverage = c1.build_tasks(
            (129,), (201.5,), isotope_data, nv=4, nw=8, nz=6, nphi=4
        )
        self.assertEqual(len(coverage), 2)
        self.assertEqual(len(tasks), 8)
        self.assertEqual({x[3] for x in tasks}, {"shm_cutoff", "mb_untruncated"})
        self.assertEqual({x[8:] for x in tasks}, {(8, 6, 4)})
        self.assertTrue(all(x["open_speed_interval"] for x in coverage))

    def test_speed_integral_keeps_full_complex_matrix(self):
        L = np.asarray([
            [2.0, 1.0j, 0.0, 0.0],
            [-1.0j, 3.0, 0.5, 0.0],
            [0.0, 0.5, 4.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ], dtype=complex)
        rows = [{
            "phase_space_open": True,
            "speed_weight_kms": 2.0,
            "halo_pdf_kms_inv": 0.5,
            "beta": 0.01,
            "pchi_GeV": 1.0,
            "L": L,
        }]
        got = c1.speed_integral(rows)
        np.testing.assert_allclose(got, 0.01 * L, rtol=0.0, atol=1e-14)
        self.assertEqual(got[0, 1], 0.01j)
        self.assertEqual(got[1, 0], -0.01j)

    def test_frozen_benchmark_keV_to_GeV_conversion(self):
        self.assertAlmostEqual(c1.MS_GEV, 200.0e-6, places=15)
        self.assertAlmostEqual(c1.DELTA_M_GEV, 1000.0e-6, places=15)
        self.assertEqual(c1.WIDTH_SCHEME, "pole")

    def test_c1_revised_quadrature_order(self):
        self.assertEqual((c1.NW, c1.NZ, c1.NPHI), (14, 6, 4))
        self.assertEqual(c1.PHASE_COMPARISON, (16, 14, 12))

    def test_run_one_task_returns_full_complex_lorentz_tensor(self):
        L = np.eye(4, dtype=complex)
        G = np.zeros((4, 4, 4, 4), dtype=complex)
        for a in range(4):
            G[:, a, :, a] = L / 4.0
        result = {
            "phase_space_open": True,
            "L_q_aligned_real_imag": pairs(L),
            "L_before_real_imag": pairs(0.25 * L),
            "L_after_real_imag": pairs(0.25 * L),
            "L_interference_real_imag": pairs(0.50 * L),
            "G_q_aligned_real_imag": pairs(G),
            "q_frame_rotation_max_orthogonality_residual": 0.0,
            "q_frame_rotation_max_determinant_residual": 0.0,
            "quadrature_points": 336,
        }
        task = (129, 201.5, 120.04644920255713, "shm_cutoff", 0,
                600.0, 10.0, 0.001, 14, 6, 4)
        with patch.object(c1, "integrate_spin_resolved_tensor", return_value=result) as mocked:
            row = c1.run_one_task(task)
        kwargs = mocked.call_args.kwargs
        self.assertAlmostEqual(kwargs["ms"], 200.0e-6, places=15)
        self.assertAlmostEqual(kwargs["delta_m"], 1000.0e-6, places=15)
        self.assertEqual(kwargs["width_scheme"], "pole")
        self.assertEqual(row["L"].shape, (4, 4))
        np.testing.assert_allclose(row["L"], L)

    def test_relative_matrix_difference_is_scale_invariant(self):
        A = np.eye(4, dtype=complex)
        B = (1.0 + 1.0e-7) * A
        rel, component = c1.relative_matrix_difference(A, B)
        self.assertAlmostEqual(rel, 1.0e-7, delta=1e-14)
        self.assertAlmostEqual(component, 5.0e-8, delta=1e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)
