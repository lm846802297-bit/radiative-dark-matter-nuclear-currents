import json
from pathlib import Path
import unittest

import numpy as np

from compute import coefficient_vector, responses


class SpatialTests(unittest.TestCase):
    def test_interference_in_single_square(self):
        a=np.array([[2.,-1.],[3.,4.]])
        c,w,v=coefficient_vector(.2,15.,.81,-.184)
        r=responses(a)
        self.assertGreater(r[2],0)
        self.assertLess(c[2]*r[2],0)
        self.assertAlmostEqual(c@r,w*np.sum((a@v)**2))
        self.assertLessEqual(abs(r[2]),(r[0]*r[1])**.5)

    def test_zero_recoil_weight(self):
        c,w,v=coefficient_vector(0,15.,1.,-.184)
        np.testing.assert_array_equal(c,np.zeros(3))
        self.assertEqual(w,0)
        self.assertLess(v[1],0)

    def test_archived_odd_isotope_controls(self):
        s=json.loads((Path(__file__).parent/'native/summary.json').read_text())
        self.assertEqual([(r['A'],r['ER_keV']) for r in s['benchmarks']],[(129,200.),(129,267.),(131,200.),(131,267.)])
        for a in ('129','131'):
            self.assertTrue(s['checks']['published_q0_control'][a]['all_ten_q0_signs_agree'])
        for r in s['benchmarks']:
            if r['A']==131:
                np.testing.assert_array_equal(r['SN100PN_F_normalized_amplitudes_p_Delta_SigmaPrime'][1],[0.,0.])
                np.testing.assert_array_equal(r['SN100PN_F_normalized_amplitudes_n_Delta_SigmaPrime'][1],[0.,0.])


if __name__=='__main__':
    unittest.main()
