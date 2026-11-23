import unittest

import numpy as np

from finite_t import independent_ff, signed_terms, separate


class FiniteT(unittest.TestCase):
    def test_equal_species_and_anchors(self):
        for scheme in ('anchors','bishara_linear','isoscalar_zexp'):
            self.assertTrue(np.allclose(independent_ff(scheme,0),[1.,-.184],rtol=0,atol=1e-15))
        self.assertAlmostEqual(independent_ff('bishara_linear',-.06)[0],1-.06*8.392/3)
        self.assertAlmostEqual(independent_ff('bishara_linear',-.06)[1],-.184)

    def test_one_interference_product(self):
        np.testing.assert_allclose(signed_terms(np.array([3.,2.]),[[16.,4.,-8.]],2),[[288,32,-192]])

    def test_model_and_current_effect_close_on_both_paths(self):
        r=separate({'AFH_BD':10.,'SN100PN':20.},{'AFH_BD':8.,'SN100PN':12.})
        self.assertAlmostEqual(r['nuclear_only_anchors'],2.)
        self.assertAlmostEqual(r['nuclear_only_at_scheme'],1.5)
        self.assertAlmostEqual(r['model_by_ff_ratio_of_ratios'],.75)
        self.assertAlmostEqual(r['joint_vs_anchor_AFH_BD'],1.2)

if __name__=='__main__':
    unittest.main()
