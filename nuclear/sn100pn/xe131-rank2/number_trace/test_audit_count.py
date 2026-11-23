import unittest
from decimal import Decimal as D, localcontext
from audit_count import occupation, last_place

class Coefficients(unittest.TestCase):
    def test_independent_forward_inverse(self):
        with localcontext() as ctx:
            ctx.prec=50
            T=D('11.5')
            ji=4
            g=8
            # App. B Eq. (14): a(NN)=sqrt((2Ji+1)/g) n_N at J=0.
            ap=(D(ji)/g).sqrt()*D('2.8')
            an=(D(ji)/g).sqrt()*D('7.71')
            # App. B Eqs. (16)-(17) with C^{TM}_{00;TM}=1 and C^{TM}_{10;TM}=+sqrt(T/(T+1)).
            psi0=(2*T+1).sqrt()*(ap+an)/D(2).sqrt()
            psi1=(2*T+1).sqrt()*(an-ap)/(D(6).sqrt()*(T/(T+1)).sqrt())
            p,n,*_=occupation(psi0,psi1,g,T,ji*(2*T+1))
            self.assertLess(abs(p-D('2.8')),D('1e-47'))
            self.assertLess(abs(n-D('7.71')),D('1e-47'))

    def test_exact_closed_core(self):
        with localcontext() as ctx:
            ctx.prec=50
            T=D('10.5'); Ji=2; g=10
            psi0=(2*D(g)*Ji*(2*T+1)).sqrt()
            p,n,*_=occupation(psi0,None,g,T,Ji*(2*T+1))
            self.assertLess(abs(p-g),D('1e-47'))
            self.assertEqual(p,n)

    def test_decimal_precision_retains_trailing_zeros(self):
        self.assertEqual(last_place(D('43.81780')),D('1e-5'))
        self.assertEqual(last_place(D('2.121437835')),D('1e-9'))
        self.assertEqual(last_place(D('6.00')),D('1e-2'))
if __name__=='__main__': unittest.main()
