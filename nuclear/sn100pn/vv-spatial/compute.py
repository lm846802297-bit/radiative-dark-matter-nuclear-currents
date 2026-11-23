#!/usr/bin/env python3
"""Signed xenon spatial AFH responses with the frozen four-point VV T3 benchmark."""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import lzma
from pathlib import Path
import re
import subprocess
import sys
import tarfile

import mpmath as mp
import numpy as np

BASE='7f82851b777e21a8618a8bd5767448dcc58f4a70'
LZ_HELPER='ebe6e8cbd2e67d79e2a34305651d6b53ac565795'
AFH_REF='aac9b33d66f6a7fe1ca5d5d8cc648f5afe52a81b'
SN_REF='1395f0c161905341fc89a2978148b1ce7c9ee1dd'
SN_DIR=Path('nuclear/sn100pn')
SRC=Path('benchmarks/vv_recovery/source_snapshot')
P3=SRC/'research/three_body_full_current_nucleon_match'
NATIVE=Path('nuclear/afh/data/native_v6_2026-11-22/native/spatial_multipoles/afh_v6_xe129_xe131_spatial_multipoles_native.csv')
SUP='tek/bodysupplementary.tex'
EXPECTED={
 'afh':'dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372',
 'helper':'a9c0a68e354bd3b4504e6a65c24788662e95bdcc1466debeccdcbe5ef3a963e0',
 129:'26ef95e96f79bcd883cc3374e26aac62fdb80be23439e1604a96703468e8dd91',
 131:'2dd36c7aa74a88f8fc979b2423ac2970a605300d0eca973d94cb139bda438bd0',
 'sup':'00d6b232355a140d70b64fae5ef01c5a58203866b3137d9bae6a61e8ec21cc38',
 'native_rank':'cc1d89d0d3cce9cb138226a9dcc880e76bde8a3d2bc22698ac072900320726e1',
 'finite_t_summary':'8fa341408682ba3fd6ccf41054fe3ebd9cce82698a23014c1a6125cb26305e65',
}
A_VALUES=(129,131)
NORM_J={129:2,131:4}
MN=.938918
SP=('Delta','SigmaPrime','SigmaPrimeDelta')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def pinned(raw,expected,label):
    actual=sha(raw)
    if actual!=expected: raise ValueError(f'hash mismatch for {label}: {actual}')
    return raw


def load(path,expected):
    return pinned(path.read_bytes(),expected,path)


def read_sums(p):
    return {f:h for h,f in (l.split() for l in p.read_text().splitlines() if l and not l.startswith('#'))}


def af_rows(raw,a):
    arr=re.findall(rf'sdgXe{a}\[(\d+)\]=\{{[^=]+?=([-+0-9.]+),(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)\}};',
                   raw.decode())
    if len(arr)!={129:39,131:67}[a]: raise ValueError('wrong native density block size')
    return [(int(j),int(t),tuple(map(int,v)),rho) for _,rho,j,t,*v in arr if j!='0' and t=='0']


def sn_rows(raw,a):
    out=[]
    for line in raw.decode().splitlines():
        if line.startswith('THE'):
            j,t=[int(x)//2 for x in re.search(r'2JO = (\d+) TO =, (\d+)',line).groups()]
        elif j>0:
            dot,*data,rho=line.split()
            if dot!='.': raise ValueError('malformed density')
            out.append((j,t,tuple(map(int,data)),rho))
    if sorted({j for j,t,v,rho in out})!=list(range(1,NORM_J[a])):
        raise ValueError('wrong SN100PN density rank set')
    return out


def responses(amp):
    a=np.asarray(amp,dtype=float)
    d,s=a[...,0],a[...,1]
    return np.array([(d*d).sum(),(s*s).sum(),(s*d).sum()])


def coefficient_vector(q,Lt,f1,f2):
    weight=q*q*Lt/(2*MN*MN)
    v=np.array([f1,-(f1+f2)/2])
    return weight*np.array([v[0]**2,v[1]**2,2*v[0]*v[1]]),weight,v


def calculate(a,rows,q,obs,precision=50):
    with mp.workdps(precision):
        t=mp.mpf(a-108)/2
        b=mp.sqrt((mp.mpf(41467)/1000)/(45*a**(-mp.mpf(1)/3)-25*a**(-mp.mpf(2)/3)))*mp.mpf(50677)/10000
        y=(mp.mpf(str(q))*b/2)**2
        norm=(mp.sqrt(2/(2*t+1)),mp.sqrt(6*t/((2*t+1)*(t+1))))
        a0,a1=[],[]
        for j in range(1,NORM_J[a]):
            for tau,alist in ((0,a0),(1,a1)):
                terms=[tuple(mp.mpf(rho)*d for d in obs(y,*orb,j))
                       for rank,t,orb,rho in rows if (rank,t)==(j,tau)]
                alist.append([norm[tau]*mp.fsum(term[i] for term in terms) for i in range(2)])
        return y,a0,a1


def supplement(tex,a):
    section=tex.split(r'\label{nuclear responses}',1)[1].split(r'\subsection*{$^{'+str(a)+r'}$Xe}'+'\n')[1].split(r'\subsection*',1)[0]
    pat=r"\$F_\{([^{}]+)\}\^\{\(([pn]),([pn])\)\}\s*=\s*e\^\{-2\s*y\}\s*\(([+-]?\d+(?:\.\d*)?)(?:\\times10\^\{([+-]?\d+)\})?"
    symbols={r'\Delta':'Delta',r"\Sigma'":"SigmaPrime",r"\Sigma' \Delta":"SigmaPrimeDelta"}
    data={}
    for key,s,t,d,e in re.findall(pat,section):
        if key in symbols:
            data[symbols[key],s,t]=float(d)*10**int(e or 0)
    if len(data)!=10: raise ValueError(f'expected 10 published spatial formula entries, got {len(data)}')
    return data


def supplement_compare(pn,printed):
    by={key:[] for key in SP}
    sign=[]
    for (key,p,n),v in printed.items():
        ia,ib='pn'.index(p),'pn'.index(n)
        # Published interference order is Sigma' of first species, Delta of second.
        i,j={'Delta':(0,0),'SigmaPrime':(1,1),'SigmaPrimeDelta':(1,0)}[key]
        v2=float(pn[ia,:,i]@pn[ib,:,j])
        by[key].append(abs(v2-v)/max(abs(v),1e-30))
        sign.append(v*v2>0)
    if not all(sign): raise ValueError('literature mixed order/sign mismatch')
    return {k:max(v) for k,v in by.items()}


def main(args):
    root=args.repo.resolve()
    def git(*v):return subprocess.check_output(['git','-C',root,*v],text=True).strip()
    if git('rev-parse','HEAD')!=BASE or git('status','--short'):raise ValueError('wrong or dirty benchmark checkout')
    load(args.helper,EXPECTED['helper'])
    ps=importlib.util.spec_from_file_location('lz_spatial_operator_helper',args.helper)
    mod=importlib.util.module_from_spec(ps)
    ps.loader.exec_module(mod)
    af=load(args.afh,EXPECTED['afh'])
    sn={a:sn_rows(load(args.sn_inputs/f'{a}Xe_OBDME_SN100PN.dat',EXPECTED[a]),a) for a in A_VALUES}
    with tarfile.open(args.paper,'r:gz') as tar:
        tex=pinned(tar.extractfile(SUP).read(),EXPECTED['sup'],SUP).decode()
    arch=root/SN_DIR/'vv-charge-matching'
    sums=read_sums(arch/'SHA256SUMS.txt')
    ts=lzma.decompress(load(arch/'native/radiative_tensors.json.xz',sums['native/radiative_tensors.json.xz']))
    ht=read_sums(arch/'RAW_TENSOR_SHA256SUMS.txt')['native/radiative_tensors.json']
    t3=json.loads(pinned(ts,ht,'uncompressed archived T3'))
    finite=json.loads(load(root/SN_DIR/'vv-finite-t/native/finite_t_summary.json',EXPECTED['finite_t_summary']))
    names={(r['A'],r['ER_keV']):r for r in finite['benchmarks']}
    native=load(root/NATIVE,EXPECTED['native_rank'])
    known={}
    for row in csv.DictReader(io.StringIO(native.decode())):
        a=int(row['A']);i=round(float(row['q_GeV'])*4000);j=int(row['multipole_J'])
        assert float(row['Im_M_Delta_J_native'])==0 and float(row['Im_M_SigmaPrime_J_native'])==0
        known[a,i,j]=np.array([float(row['M_Delta_J_native']),float(row['M_SigmaPrime_J_native'])])
    if len(known)!=9604:raise ValueError('unexpected native signed spatial grid row count')
    manifest=read_sums(root/SRC/'SHA256SUMS.txt')
    hashes={}
    for x in ('three_body_multipole_leptonic_weights.py','nucleon_vector_form_factors.py'):
        p=P3/x
        h=manifest[str(p.relative_to(SRC))]
        load(root/p,h);hashes[str(p)]=h
    sys.path.insert(0,str(root/P3))
    import three_body_multipole_leptonic_weights as p3
    obs=mod.obs
    results=[]; checks={}
    active_err=source_err=precision_check=0.
    for a in A_VALUES:
        bd=af_rows(af,a)
        jnorm=NORM_J[a]
        norm=np.sqrt(4*np.pi/jnorm)
        def get(rows,q):
            y,z,o=calculate(a,rows,q,obs)
            return float(y),np.array(z,dtype=float),np.array(o,dtype=float)
        # Both species are independent reduced nuclear inputs even when the source couples equally.
        _y,zero0,zero1=get(sn[a],0)
        pns=np.array([norm*(zero0-zero1)/2,norm*(zero0+zero1)/2])
        pp=supplement_compare(pns,supplement(tex,a))
        checks[a]={'published_q0_species_response_max_relative_by_response':pp,'all_ten_q0_signs_agree':True}
        for er,case in ((200.,'green_450keV'),(267.,'near_pole_362keV_g2_1')):
            f=names[a,er]
            q=f['Q_A_GeV']
            # Native actual amplitudes are checked on adjacent native q nodes; target q is direct source evaluation.
            br=np.floor(q*4000).astype(int)
            for k in (0,br,br+1):
                _y,b0,_b1=get(bd,k/4000)
                objective=np.array([known[a,k,j] for j in range(1,jnorm)])
                e=float(np.max(abs(b0-objective)))
                source_err=max(source_err,e)
                if e>=3e-13: raise ValueError(f'native signed anchor mismatch Xe{a} {k}: {e}')
            sY,s0,s1=get(sn[a],q)
            _,b0,b1=get(bd,q)
            assert np.all(b1==0)
            with mp.workdps(100):
                q50=calculate(a,sn[a],q,obs,50)
                q80=calculate(a,sn[a],q,obs,80)
                e=max(abs(u-v) for block,block2 in zip(q50[1:],q80[1:])
                      for row,row2 in zip(block,block2) for u,v in zip(row,row2))
                precision_check=max(precision_check,float(e))
            p_amp=norm*(s0-s1)/2
            n_amp=norm*(s0+s1)/2
            amps={'AFH_BD':norm*b0,'SN100PN':p_amp+n_amp}
            L=p3.complex_from_pairs(t3[f'Xe{a}_{case}']['refined_grid']['L_q_aligned_real_imag'])
            LT=(L[1,1]+L[2,2]).real
            eig,vec=np.linalg.eigh(L[1:,1:])
            assert eig.min()>=-1e-14*np.linalg.norm(L[1:,1:])
            currents=(vec*np.sqrt(np.maximum(eig,0))).T
            row={'A':a,'ER_keV':er,'Q_A_GeV':q,'y':sY,'L_T':LT,
                 'allowed_density_ranks':list(range(1,jnorm)),
                 'raw_norm_AFH_iso0_signed':b0.tolist(),
                 'SN100PN_F_normalized_amplitudes_p_Delta_SigmaPrime':p_amp.tolist(),
                 'SN100PN_F_normalized_amplitudes_n_Delta_SigmaPrime':n_amp.tolist(),
                 'schemes':{}}
            for scheme in ('anchors','bishara_linear','isoscalar_zexp'):
                fc=f['cases'][scheme]
                F1,F2=fc['F1hat_p'],fc['F2hat_p']
                assert F1==fc['F1hat_n'] and F2==fc['F2hat_n']
                cs,w,v=coefficient_vector(q,LT,F1,F2)
                src=p3.extract_afh_vector_coefficients(L,q,omega_GeV=er*1e-6,f1=F1,f2=F2)
                src2=p3.afh_eq107_spatial_coefficients_from_amplitudes(currents,np.array([0.,0.,q]),
                           f1=F1,f2=F2,spin_average_denominator=1)
                for source in (src,src2):
                    expect=np.array([source[f'C_{p}'] for p in SP])
                    np.testing.assert_allclose(cs,expect,rtol=5e-14,atol=0)
                    active_err=max(active_err,float(np.max(abs(cs/expect-1))))
                scheme_out={}
                for model,ar in amps.items():
                    R=responses(ar)
                    assert abs(R[2])<=np.sqrt(R[0]*R[1])+1e-13
                    terms=cs*R
                    sq=w*np.sum((ar@v)**2)
                    np.testing.assert_allclose(terms.sum(),sq,rtol=3e-14,atol=1e-16)
                    c=fc['models'][model]['conditional_L_C']['kernel']
                    scheme_out[model]={'responses_equal_pn':R.tolist(),'spatial_terms_Delta_SigmaPrime_SigmaPrimeDelta':terms.tolist(),
                            'spatial_kernel':float(sq),'conditional_charge_plus_spatial':float(c+sq),
                            'spatial_over_conditional_charge':float(sq/c)}
                g,h=scheme_out['AFH_BD'],scheme_out['SN100PN']
                scheme_out['both_over_charge_only_model_ratio']=(h['conditional_charge_plus_spatial']/g['conditional_charge_plus_spatial'])/f['effects'][scheme]['nuclear_only_at_scheme']
                row['schemes'][scheme]=scheme_out
            results.append(row)
    if precision_check>=1e-45 or max(max(x['published_q0_species_response_max_relative_by_response'].values()) for x in checks.values())>0.04:
        raise ValueError(f'precision or published low-Q control failed: {precision_check}, {checks}')
    result={'scope':'leading-static transverse AFH diagnostic at four fixed radiative VV points, not full current or cross section',
            'parent':BASE,'source_pins':{'afh_package_ref':AFH_REF,'spatial_operator_helper_ref':LZ_HELPER,'SN100PN_OBDME_author_ref':SN_REF},
            'checks':{'signed_BD_native_anchor_max_abs':source_err,'particle_coefficient_crosscheck_max_relative':active_err,
                      'SN100PN_50_vs_80_digit_signed_input_max_abs':precision_check,'published_q0_control':checks},
            'inputs_sha256':{**EXPECTED,**hashes,'compressed_T3_tensor':sums['native/radiative_tensors.json.xz'],'raw_T3_tensor':ht},
            'benchmarks':results}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','afh','helper','sn-inputs','paper','out'):
        p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
