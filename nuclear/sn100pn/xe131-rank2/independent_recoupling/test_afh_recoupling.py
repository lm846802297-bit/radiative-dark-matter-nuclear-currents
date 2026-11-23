import unittest
import math
import numpy as np
try:
    from afh_recoupling import operators
except ModuleNotFoundError as e:
    if e.name != 'afh_recoupling':
        raise
    operators = None


class TestAFHRecoupling(unittest.TestCase):
    def test_rank0_spin_orbit_and_charge_control(self):
        self.assertIsNotNone(operators, 'AFH recoupled operator implementation is absent')
        for orbit, sdl in (((4, 1), 0), ((4, 3), -3), ((4, 5), 2), ((4, 7), -5), ((5, 11), 5)):
            m, phi = operators(0.0, *orbit, *orbit, 0)
            self.assertAlmostEqual(m, math.sqrt((orbit[1] + 1)/(4*math.pi)), places=13)
            self.assertAlmostEqual(phi, -sdl*m/3, places=13)

    def test_rank2_zero_charge_and_hermiticity(self):
        self.assertIsNotNone(operators, 'AFH recoupled operator implementation is absent')
        self.assertEqual(operators(0.0, 4, 3, 4, 3, 2)[0], 0.0)
        self.assertGreater(abs(operators(0.0, 4, 3, 4, 3, 2)[1]), 1e-3)
        for y in (0.0, .2, 1.5, 5.0):
            np.testing.assert_allclose(operators(y, 4, 3, 4, 5, 2),
                                       -np.array(operators(y, 4, 5, 4, 3, 2)),
                                       atol=2e-13, rtol=2e-13)


if __name__ == '__main__':
    unittest.main()
