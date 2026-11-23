#!/usr/bin/env python3
"""Independent AFH v6 recoupling (3j/6j/9j plus HO 1F1) for M and Phi''.

Origin: DMFormFactor-v6.m, MJ and PhiPPJ, fixed source SHA-256 recorded
in the accompanying README. No numerical angular projection is used.
"""
from functools import lru_cache
import math
import numpy as np
from scipy.special import factorial2, hyp1f1
from sympy import Rational
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j

HALF = Rational(1, 2)


@lru_cache(None)
def w(fn, *a):
    try:
        return float(fn(*a))
    except ValueError:  # exact inadmissible triangles
        return 0.0


def qnorm(j):
    return math.sqrt(2*j+1)


def orbit(N, two_j):
    hi = (two_j+1)//2
    l = hi if (N-hi) % 2 == 0 else hi-1
    return (N-l)//2, l


def radial(y, np_, lp, n, l, L, grad=None):
    # n,n' are zero-based; grad is the l -> l±1 term. Standard HO phase.
    derivative = grad is not None
    if y == 0 and L > int(derivative):
        return 0.0
    pref = ((2*math.sqrt(y))**(L-int(derivative))/float(factorial2(2*L+1))
            * math.exp(-y)
            * math.sqrt(math.factorial(n)*math.factorial(np_) *
                        math.gamma(n+l+1.5)*math.gamma(np_+lp+1.5)))
    b = L+1.5
    pieces = []
    for m in range(n+1):
        for mp in range(np_+1):
            a = (L-l-lp-2*(m+mp)-int(derivative))/2
            z = (L+l+lp+2*(m+mp)+3-int(derivative))/2
            coef = ((-1)**(m+mp)*math.gamma(z) /
                    (math.factorial(m)*math.factorial(mp)*math.factorial(n-m)*
                     math.factorial(np_-mp)*math.gamma(l+m+1.5)*math.gamma(lp+mp+1.5)))
            f = hyp1f1(a,b,y)
            if derivative:
                f = -z*f + (2*m if grad==1 else 2*l+2*m+1)*hyp1f1(a+1,b,y)
            pieces.append(float(coef*f))
    return pref*math.fsum(pieces)


@lru_cache(None)
def structure(Np, tjp, N, tj, J):
    np_,lp = orbit(Np,tjp)
    n,l = orbit(N,tj)
    jp,j = Rational(tjp,2),Rational(tj,2)
    cfg = (np_,lp,n,l)
    if (lp+l+J) % 2 or not abs(jp-j)<=J<=jp+j or min(np_,n)<0:
        return cfg,0.,()
    overall = qnorm(lp)*qnorm(jp)*qnorm(j)/math.sqrt(4*math.pi)
    mw = ((-1)**int(HALF+j+J)*overall*qnorm(l)*qnorm(J)
          *w(wigner_3j,lp,J,l,0,0,0)*w(wigner_6j,lp,jp,HALF,j,l,J))
    terms = []
    for L in (J+1,J-1):
        if L<1:
            continue
        K = J+1 if L==J+1 else J
        for d in (1,-1):
            if l+d<0:
                continue
            recouple = math.fsum(
                (-1.)**(J-C+(d==1))*(2*C+1)
                *w(wigner_6j,L,1,C,1,J,1)*w(wigner_6j,L,1,C,l,lp,l+d)
                *w(wigner_9j,lp,l,C,HALF,HALF,1,jp,j,J)
                for C in sorted((J,L))
            )
            c = ((-1)**(lp+1)*6*overall*qnorm(L)*math.sqrt(K)
                 *qnorm(l+d)*math.sqrt(l+(d==1))*w(wigner_3j,lp,L,l+d,0,0,0)*recouple)
            if c:
                terms.append((L,d,c))
    return cfg,mw,tuple(terms)


def operators(y,Np,tjp,N,tj,J):
    cfg,mw,terms=structure(Np,tjp,N,tj,J)
    m = mw*radial(y,*cfg,J) if mw else 0.0
    phi = math.fsum(c*radial(y,*cfg,L,grad=d) for L,d,c in terms)
    return m,phi
