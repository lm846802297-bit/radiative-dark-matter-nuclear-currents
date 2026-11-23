#!/usr/bin/env python3
"""Xe131 signed SN100PN rank-2 charge and spin-orbit response replay."""
import argparse
import csv
from functools import lru_cache
import json
from pathlib import Path
import re
import sys
import tarfile

import numpy as np
from scipy.optimize import brentq
from scipy.special import gamma, genlaguerre, hyp1f1, roots_genlaguerre, sph_harm_y, spherical_jn
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j

# Works in staging and as a sibling of the pinned monopole audit in the project.
PAIR = ('xe129-xe131-rank0', 'radiative_sn100pn_stage')
for folder in PAIR:
    src = Path(__file__).resolve().parent.parent / folder
    if (src / 'audit.py').is_file():
        sys.path.insert(0, str(src))
        break
from audit import (AFH_DIR, ARXIV_SOURCE, AUTHORS_SHA, BASELINE_SHA,
                   CHANNELS, EXPECTED_SHA, RAW, TEX_SHA, amplitudes, iso_response,
                   load, occupations, orbit, printed_response, sha256)

A = 131
SIGMA_X = np.array([[0,1],[1,0]])
SIGMA_Y = np.array([[0,-1j],[1j,0]])
SIGMA_Z = np.array([[1,0],[0,-1]])


def read_densities(raw, J):
    result = {}
    active = None
    for line in raw.decode('ascii').splitlines():
        if line.startswith('THE'):
            rank, twice_tau = map(int, re.fullmatch(
                r'THE ONE-BODY DENSITY MATRIX FOR 2JO = (\d+) TO =, (\d+)', line).groups())
            active = twice_tau//2 if rank == 2*J else None
            if active is not None:
                if twice_tau not in (0,2) or active in result:
                    raise ValueError('bad OBDME channel')
                result[active] = {}
        elif active is not None:
            dot, na, ja, nb, jb, v = line.split()
            key = (int(na),int(ja)), (int(nb),int(jb))
            if dot != '.' or key in result[active]:
                raise ValueError(f'bad or repeated OBDME line: {line}')
            result[active][key] = float(v)
    if set(result) != {0,1}:
        raise ValueError('missing OBDME rank/isospin channel')
    return result


def spinor(ell, twice_j, twice_m, theta, sin_theta):
    # Positive spin-1/2 basis, Condon-Shortley Y_lm and C(l,m_l;1/2,m_s|j,m).
    result = np.zeros((3, len(theta), 2), dtype=complex)
    for s, twice_ms in enumerate((1,-1)):
        ml = (twice_m-twice_ms)//2
        if abs(ml) > ell:
            continue
        c = float(clebsch_gordan(ell, Rational(1,2), Rational(twice_j,2), ml,
                                 Rational(twice_ms,2), Rational(twice_m,2)))
        if c != 0:
            Y, d = sph_harm_y(ell, ml, theta, 0., diff_n=1)
            result[:, :, s] = c*np.array([Y,d[:,0],d[:,1]/sin_theta])
    return result


@lru_cache(None)
def angular(oa, ob, J, nodes):
    """Four reduced integrals for M, Phi_r, Phi_theta(radial), Phi_theta(azimuth)."""
    _,la,_ = orbit(*oa)
    _,lb,_ = orbit(*ob)
    ja,jb=oa[1],ob[1]
    if not abs(ja-jb) <= 2*J <= ja+jb or (la+lb+J)%2:
        return np.zeros(4)
    x,w = np.polynomial.legendre.leggauss(nodes)
    theta = np.arccos(x)
    sn = np.sqrt(1-x*x)
    yr, dy = sph_harm_y(J,0,theta,0.,diff_n=1)
    ytheta=dy[:,0]
    sr=sn[:,None,None]*SIGMA_X+x[:,None,None]*SIGMA_Z
    st=x[:,None,None]*SIGMA_X-sn[:,None,None]*SIGMA_Z
    def mat(op, a):
        return np.einsum('...ij,...j->...i',op,a)
    samples=[]
    we=[]
    for m in range(-min(ja,jb),min(ja,jb)+1,2):
        bra=spinor(la,ja,m,theta,sn)[0]
        ket,dth,dph=spinor(lb,jb,m,theta,sn)
        terms=np.array([yr[:,None]*ket,
                        1j*yr[:,None]*(mat(st,dph)-mat(SIGMA_Y,dth)),
                        1j*ytheta[:,None]*mat(SIGMA_Y,ket),
                        -1j*ytheta[:,None]*mat(sr,dph)])
        samples.append(2*np.pi*np.einsum('x,sxi,xi->s',w,terms,bra.conj()))
        we.append((-1)**((ja-m)//2)*float(wigner_3j(
            Rational(ja,2),J,Rational(jb,2),Rational(-m,2),0,Rational(m,2))))
    samples,we=np.array(samples),np.array(we)
    reduced=we@samples/np.dot(we,we)
    np.testing.assert_allclose(samples,we[:,None]*reduced,rtol=0,atol=2e-12,
                               err_msg=f'Wigner-Eckart check for {oa},{ob},J={J}')
    np.testing.assert_allclose(reduced.imag,0,rtol=0,atol=1e-14)
    return reduced.real


def shifted(poly):
    return np.r_[0.,poly]


def radial(y, oa, ob, J, method, nodes=0):
    n,l,_=orbit(*oa)
    m,k,_=orbit(*ob)
    norm=np.sqrt(gamma(n+1)*gamma(m+1)/(gamma(n+l+1.5)*gamma(m+k+1.5)))
    alpha=(l+k+1)/2
    ap=genlaguerre(n,l+0.5).c[::-1]
    bp=genlaguerre(m,k+0.5).c[::-1]
    P=np.polynomial.polynomial
    dp=P.polysub(P.polyadd(k*bp,2*shifted(P.polyder(bp))),shifted(bp))
    ab=P.polymul(ap,bp)
    abd=P.polymul(ap,dp)
    if method=='laguerre':
        x,w=roots_genlaguerre(nodes,alpha)
        z=2*np.sqrt(y[:,None]*x)
        normal=norm*w*P.polyval(x,ab)
        deriv=norm*w*P.polyval(x,abd)
        f=spherical_jn(J,z)
        if J==0:
            d=np.divide(-spherical_jn(1,z),z,out=np.full_like(z,-1/3),where=z!=0)
            c=np.zeros_like(d)  # corresponding angular factors vanish identically
        else:
            d=np.divide(spherical_jn(2,z,derivative=True),z,out=np.full_like(z,2/15),where=z!=0)
            c=np.divide(f,z*z,out=np.full_like(z,1/15),where=z!=0)
        return np.array([f@normal,d@normal,c@deriv,c@normal])
    if method!='analytic':
        raise ValueError(method)
    def integral(poly,b):
        a=alpha+1+np.arange(poly.size)
        return norm*sum(p*gamma(t)*hyp1f1(t,b,-y) for p,t in zip(poly,a))
    if J==0:
        m0=integral(ab,1.5)
        return np.array([m0,-integral(ab,2.5)/3,np.zeros_like(m0),np.zeros_like(m0)])
    return np.array([4*y*integral(shifted(ab),3.5)/15,
                     2*integral(ab,3.5)/15-4*y*integral(shifted(ab),4.5)/105,
                     integral(abd,3.5)/15,integral(ab,3.5)/15])


def nuclear(y, raw, J, method, r_nodes=0, a_nodes=32):
    densities=read_densities(raw,J)
    t=(A-108)/2
    weights=np.sqrt(4*np.pi/(4*(2*t+1)))*np.array([np.sqrt(2),np.sqrt(6*t/(t+1))])
    both=np.zeros((len(y),2,2))
    combos=sorted(set().union(*densities.values()))
    for a,b in combos:
        angulars=angular(a,b,J,a_nodes)
        terms=angulars[:,None]*radial(y,a,b,J,method,r_nodes)
        coefficients=weights*np.array([densities[tau].get((a,b),0.) for tau in (0,1)])
        both+=np.array([terms[0],terms[1:].sum(axis=0)]).T[:,:,None]*coefficients
    per_species=both@np.array([[.5,.5],[-.5,.5]])
    return per_species.reshape(len(y),4)


RANK0_COMMIT = '51da16ceb636413c0a0c5e0b98f64224fd4a216d'
RANK0_AUDIT_SHA = '7f83ac0ce395194d7b200710e5c7d8738ed8c16c1a44ba4d904dd94316d67e0a'


def outer(amp):
    return amp[:, :, None]*amp[:, None, :]


def validate_rmes(y, raw):
    two = read_densities(raw, 2)
    pairs = sorted(set(two[0]) | set(two[1]))
    residual = max(abs(two[tau][a,b] - (-1)**((a[1]-b[1])//2)*two[tau][b,a])
                   for tau in (0,1) for a,b in pairs)
    # Input values are used unchanged, including the tiny partner mismatches.
    assert residual < 5e-5
    rm, rp, ra = 0., 0., 0.
    for a,b in pairs:
        ang = angular(a,b,2,40)
        rev = angular(b,a,2,40)
        terms = radial(y,a,b,2,'analytic')
        rv = radial(y,b,a,2,'analytic')
        phase = (-1)**((a[1]-b[1])//2)
        rm = max(rm, np.max(abs(ang[0]*terms[0]-phase*rev[0]*rv[0])))
        rp = max(rp, np.max(abs(ang[1:]@terms[1:]-phase*rev[1:]@rv[1:])))
        ra = max(ra, np.max(abs(ang-angular(a,b,2,24))))
    assert max(rm,rp,ra) < 2e-12
    return {'ordered_entries_tau0':len(two[0]),'ordered_entries_tau1':len(two[1]),
            'published_density_hermiticity_max_abs_residual':residual,
            'M_RME_hermiticity_max_abs_residual':float(rm),
            'PhiPP_RME_hermiticity_max_abs_residual':float(rp),
            'angular_GL24_vs_GL40_max_abs_RME_factor':float(ra)}


def run(inputs, output):
    from audit import __file__ as audit_path
    if sha256(Path(audit_path).read_bytes()) != RANK0_AUDIT_SHA:
        raise ValueError('rank0 helper checkout SHA-256 mismatch')
    inputs.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    af = f'afh_v6_xe{A}_eight_response_native.csv'
    ob = f'{A}Xe_OBDME_SN100PN.dat'
    _, raw = load(inputs, ob, f'{RAW}/Raghda-AK/DM-Form-Factor/{AUTHORS_SHA}/{ob}', EXPECTED_SHA[A][0])
    _, base = load(inputs, af,
        f'{RAW}/lm846802297-bit/radiative-dark-matter-nuclear-currents/{BASELINE_SHA}/{AFH_DIR}/{af}',
        EXPECTED_SHA[A][1])
    archive, _ = load(inputs, '2311.15764v1', ARXIV_SOURCE)
    with tarfile.open(archive, 'r:gz') as tar:
        tex = tar.extractfile('tek/bodysupplementary.tex').read()
    if sha256(tex) != TEX_SHA:
        raise ValueError('Supplementary TeX member SHA-256 mismatch')
    base = np.genfromtxt(base.decode().splitlines(), names=True, delimiter=',')
    q,y = base['q_GeV'],base['y']
    assert q.size == 2401 and q[0] == 0 and q[-1] == 0.6
    np.testing.assert_allclose(q, np.arange(2401)/4000, rtol=0, atol=2e-15)
    prior_amp0 = amplitudes(y, occupations(raw,A), 'laguerre',80)
    amp0 = nuclear(y,raw,0,'analytic',a_nodes=40)
    amp2 = nuclear(y,raw,2,'laguerre',80,a_nodes=40)
    by_formula = nuclear(y,raw,2,'analytic',a_nodes=40)
    amp2_coarse = nuclear(y,raw,2,'laguerre',32,a_nodes=24)
    check = {'general_J0_vs_prior_J0':float(np.max(abs(amp0-prior_amp0))),
             'analytic_vs_GL80':float(np.max(abs(by_formula-amp2))),
             'GL32_angular24_vs_GL80_angular40':float(np.max(abs(amp2_coarse-amp2)))}
    assert max(check.values()) < 3e-11
    assert np.all(amp2[0,:2] == 0)
    assert abs(amp2[0,3]) > 0.1
    spectral = np.linalg.eigvalsh(outer(amp0)+outer(amp2))
    assert spectral.min() > -1e-11
    assert np.max(abs(spectral[:,:2])) < 1e-11  # rank <= 2
    for yk in (1e-8,2e-8):
        small = nuclear(np.array([yk,yk/2]),raw,2,'analytic',a_nodes=40)
        np.testing.assert_allclose(small[0,:2],2*small[1,:2],rtol=1e-7,atol=1e-15)
    p = range(4)
    ui = [(i,j) for i in p for j in p if i<=j]
    tags = ('Mp','Mn','PhiPPp','PhiPPn')
    g0,g2 = outer(amp0),outer(amp2)
    total = g0 + g2
    printed = printed_response(tex.decode(),A,y)
    reference = 4*np.pi*np.column_stack([base[f'W_{k}00_AFH'] for k in CHANNELS])
    e0,e2,et,ep = map(iso_response, (g0,g2,total,printed))
    assert np.max(abs((printed[0]-total[0])/np.maximum(abs(total[0]),1))) < 0.005
    normalized = np.max(abs(et-ep),axis=0)/abs(ep[0])
    assert max(normalized) < 0.002  # absolute scale relative to printed q=0, not relative at nodes
    sign = np.signbit(amp0[:,0]+amp0[:,1])
    ix = int(np.flatnonzero(sign[:-1]!=sign[1:])[0])
    def amp_at(k,J):
        yk=(k/q[-1])**2*y[-1]
        return nuclear(np.array([yk]),raw,J,'analytic',a_nodes=40)[0]
    physical_zero=brentq(lambda k: amp_at(k,0)[:2].sum(),q[ix],q[ix+1],xtol=1e-14)
    around=amp_at(physical_zero,2)
    def landmark(i):
        return {'q_GeV':q[i],'J0_amp':amp0[i].tolist(),'J2_amp':amp2[i].tolist(),
                'J0_equal_pn':e0[i].tolist(),'J2_equal_pn':e2[i].tolist(),
                'SN100PN_full_equal_pn':et[i].tolist(),'AFH_BD_full_equal_pn':reference[i].tolist(),
                'published_rounded_full_equal_pn':ep[i].tolist()}
    summary = {'scope':'Xe131 all allowed M/PhiPP elastic one-body multipoles J=0,2; no cross section',
               'pinned_inputs':{ob:EXPECTED_SHA[A][0],af:EXPECTED_SHA[A][1],'arxiv_tex_member':TEX_SHA,
                                'rank0_audit.py':RANK0_AUDIT_SHA},
               'source_commits':{'SN100PN_OBDME':AUTHORS_SHA,'AFH_BD_table':BASELINE_SHA,'rank0_audit':RANK0_COMMIT},
               'q_grid_nodes':len(q),'amplitude_checks':check,'single_particle_operator_checks':validate_rmes(y,raw),
               'minimum_full_gram_eigenvalue':float(spectral.min()),
               'max_abs_printed_vs_full_equal_pn_over_printed_q0':normalized.tolist(),
               'landmarks':[landmark(int(k*4000)) for k in [0,.05,.10,.15,.20,.30,.40,.60]],
               'first_J0_charge_zero_bracket':{'nodes':[landmark(ix),landmark(ix+1)],
                 'continuous_J0_charge_zero_GeV':physical_zero,
                 'J2_equal_pn_M_F_at_J0_zero':float(around[:2].sum()**2)}}
    with (output/'xe131_rank2.csv').open('w',newline='') as handle:
        writer=csv.writer(handle,lineterminator='\n')
        writer.writerow(['q_GeV','y']+[f'{rank}_amp_{s}' for rank in ('J0','J2') for s in tags]
                        +[f'{rank}_F_{tags[i]}_{tags[j]}' for rank in ('J2','full') for i,j in ui]
                        +[f'{src}_F_{ch}_eqpn' for src in ('SN100PN','AFH_BD','printed') for ch in CHANNELS])
        writer.writerows([q[k],y[k],*amp0[k],*amp2[k],*(g2[k,i,j] for i,j in ui),
                          *(total[k,i,j] for i,j in ui),*et[k],*reference[k],*ep[k]] for k in range(len(q)))
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.inputs,args.output)
