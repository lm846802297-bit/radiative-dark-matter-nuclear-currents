import json
from pathlib import Path
import unittest

import numpy as np

from audit import orbit, radial


class TestAudit(unittest.TestCase):
    def test_ho_0s_charge_exact(self):
        y = np.linspace(0, 12, 97)
        for method, nodes in [('analytic', 0), ('laguerre', 32), ('laguerre', 80)]:
            m, phi = radial(y, 0, 0, method, nodes)
            np.testing.assert_allclose(m, np.exp(-y), atol=3e-14, rtol=0)
            self.assertAlmostEqual(phi[0], 1/3, places=14)

    def test_spin_orbit_partner_signs(self):
        self.assertEqual(orbit(4, 7), (0, 4, -5))
        self.assertEqual(orbit(4, 9), (0, 4, 4))
        self.assertEqual(orbit(4, 3), (1, 2, -3))
        self.assertEqual(orbit(4, 5), (1, 2, 2))

    def test_release_gate(self):
        summary = json.loads((Path(__file__).parent / 'native/summary.json').read_text())
        for a in ('129', '131'):
            r = summary['isotopes'][a]
            self.assertLess(max(r['max_abs_amplitude_delta_analytic_vs_80']), 3e-11)
            self.assertGreater(r['min_reconstructed_gram_eigenvalue'], -5e-12)
            self.assertGreater(r['printed_negative_diagonal_nodes_at_1e-8'], 0)
            self.assertLess(abs(r['AFH_q0_physical_F_M_minus_A2'])/int(a)**2, 2e-4)
        mine = next(r for r in summary['isotopes']['129']['landmarks'] if r['q_GeV'] == .15)
        self.assertGreater(mine['SN100PN_rank0'][0], 0)
        self.assertLess(mine['SN100PN_printed_all_allowed_multipoles'][0], 0)


if __name__ == '__main__':
    unittest.main()
