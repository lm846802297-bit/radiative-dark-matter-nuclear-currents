import math
import unittest
from pathlib import Path

import numpy as np

import direct_bound_current_seven_isotope as direct


class TestP4FCStageBDirectCurrent(unittest.TestCase):
    def test_seven_isotope_recoil_grid(self):
        self.assertEqual(direct.ISOTOPES, (128, 129, 130, 131, 132, 134, 136))
        self.assertEqual(len(direct.RECOIL_KEV), 81)
        self.assertEqual(direct.RECOIL_KEV[0], 40.0)
        self.assertEqual(direct.RECOIL_KEV[-1], 800.0)

    def test_exact_nuclear_masses_and_target_spins(self):
        masses = direct._nuclear_mass_table()
        self.assertEqual(tuple(sorted(masses)), direct.ISOTOPES)
        self.assertEqual(direct.TARGET_SPIN, {
            128: 0.0, 129: 0.5, 130: 0.0, 131: 1.5,
            132: 0.0, 134: 0.0, 136: 0.0,
        })
        for A in direct.ISOTOPES:
            self.assertGreater(masses[A], 0.0)

    def test_spin_zero_rank_zero_longitudinal_basis_is_allowed(self):
        # A rank-K=0 vector multipole requires L=1; the q-aligned
        # spherical component is mu=0 and its Clebsch-Gordan coefficient
        # must be nonzero. It is not an AFH Delta/Sigma-prime rank.
        self.assertIn(0, direct._allowed_J(1, 1, 0, 0, spatial=True))
        self.assertNotEqual(direct._cg_vector(1, 0, 0), 0.0)
        self.assertEqual(direct._cg_vector(1, 1, 0), 0.0)

    def test_native_afh_is_stage_a_archive(self):
        path = direct.NATIVE_CSV_P4
        self.assertTrue(path.is_file())
        self.assertEqual(
            direct.sha256(path),
            "26113db76a969dba98b0fb7b4d6b5389189d2dec33ec680d699741fb0b43a2dc",
        )

    def test_direct_vector_fit_retains_K0_L1_longitudinal_component(self):
        initial = direct.Orbit(0, 1)
        final = direct.Orbit(0, 1)
        matrices = {
            -1: np.zeros((2, 2), dtype=complex),
             0: np.eye(2, dtype=complex),
             1: np.zeros((2, 2), dtype=complex),
        }
        reduced, residual = direct._extract_vector_reduced(matrices, final, initial)
        self.assertIn((0, 1), reduced)
        self.assertGreater(abs(reduced[(0, 1)]), 0.0)
        self.assertLess(residual, 1e-12)

    def test_one_point_spin_zero_pilot_preserves_rank_zero_L1(self):
        result = direct.run(
            order=2,
            vertex_energy="external",
            mass_scheme="ame2020",
            isotopes=(128,),
            recoil_energies=(40.0e-6,),
        )
        self.assertEqual(len(result["results"]), 1)
        row = result["results"][0]
        self.assertEqual(row["A"], 128)
        self.assertEqual(row["J_A"], 0.0)
        self.assertEqual(row["J"], 0)
        self.assertIn("exact_rank0_longitudinal_identity", row)
        self.assertIn("leading_nr_rank0_longitudinal_identity", row)
        self.assertTrue(math.isfinite(abs(row["exact_rank0_longitudinal_identity"])))
        self.assertAlmostEqual(row["tN_GeV2"], -2.0 * row["mA_GeV"] * row["ER_GeV"], places=15)
        self.assertAlmostEqual(
            row["QN2_GeV2"],
            -row["tN_GeV2"],
            places=15,
        )
        self.assertAlmostEqual(
            row["Q_A_GeV"] ** 2,
            row["QN2_GeV2"] + row["ER_GeV"] ** 2,
            places=13,
        )

    def test_finite_form_factor_current_shape(self):
        p_i = np.asarray([[0.01, -0.02, -0.10]])
        p_f = np.asarray([[0.01, -0.02, 0.12]])
        baseline = direct.exact_covariant_current(p_i, p_f, 0.22, 0.0002, 1.0, -0.184)
        finite = direct.exact_covariant_current(p_i, p_f, 0.22, 0.0002, 0.84, -0.184)
        self.assertEqual(finite.shape, (1, 4, 2, 2))
        self.assertGreater(float(np.linalg.norm(finite - baseline)), 0.0)


if __name__ == "__main__":
    unittest.main()
