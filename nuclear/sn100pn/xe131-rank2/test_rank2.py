import unittest
from pathlib import Path

import numpy as np

from rank2 import angular, radial


class RankTwoTests(unittest.TestCase):
    def test_angular_0s_monopole(self):
        got = angular((0, 1), (0, 1), 0, 24)
        np.testing.assert_allclose(got, [np.sqrt(2/(4*np.pi)), 0, 0, 0], atol=1e-14)

    def test_angular_selection_and_spin_orbit(self):
        np.testing.assert_allclose(angular((0, 1), (0, 1), 2, 24), 0)
        for orbit, sigma_l in [((4,7),-5),((4,9),4),((3,5),-4)]:
            m, a, b, c = angular(orbit, orbit, 0, 24)
            np.testing.assert_allclose([a,b,c], [sigma_l*m,0,0], atol=1e-12)

    def test_radial_zero_energy(self):
        for orth_method in ('analytic','laguerre'):
            terms = radial(np.array([0.]), (4, 7), (4, 7), 2, orth_method, 48)
            np.testing.assert_allclose(terms[:,0], [0., 2/15, -1/10, 1/15], atol=1e-12)

    def test_pinned_endpoint(self):
        summary = Path(__file__).with_name('native') / 'summary.json'
        if not summary.exists():
            self.skipTest('run full pinned-input replay first')
        import json
        brief = json.loads(summary.read_text())
        self.assertGreater(abs(brief['landmarks'][0]['J2_amp'][2]), 0.01)
        self.assertEqual(brief['landmarks'][0]['J2_amp'][:2], [0.,0.])
        self.assertLess(brief['amplitude_checks']['analytic_vs_GL80'], 3e-11)
