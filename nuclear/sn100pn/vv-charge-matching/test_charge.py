import unittest

import numpy as np

from fixed_charge import charge_coefficients, contract, dirac_pauli_check, nuclear_transfer


class ChargeTests(unittest.TestCase):
    def test_exact_mass_shell(self):
        m, er = 122., .000267
        q, t = nuclear_transfer(m, er)
        self.assertAlmostEqual(q*q, 2*m*er + er*er)
        self.assertAlmostEqual(t, -2*m*er)
        self.assertGreater(q*q + t, 0)

    def test_spin_orbit_sign_and_fourth_order_remainder(self):
        result = dirac_pauli_check()
        self.assertGreater(min(result['consecutive_error_ratios']), 15)
        self.assertLess(max(result['consecutive_error_ratios']), 17)
        self.assertGreater(result['wrong_spin_orbit_sign_final_error_over_correct_error'], 100)

    def test_one_cross_term_and_equal_pn_normalization(self):
        c = charge_coefficients(.2,-.04,1.,-.184,.938918)
        self.assertGreater(c[1],0)
        ranks = np.array([[2.,-3.],[4.,1.]])
        response = np.array([(ranks[:,0]**2).sum(),(ranks[:,1]**2).sum(),np.prod(ranks,axis=1).sum()])
        terms = contract(c,response,3.)
        self.assertAlmostEqual(terms[2],6*c[0]*c[1]*response[2])
        self.assertAlmostEqual(sum(terms),3*sum((ranks@c)**2))


if __name__ == '__main__':
    unittest.main()
