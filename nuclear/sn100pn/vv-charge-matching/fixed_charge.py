#!/usr/bin/env python3
"""Fixed-speed radiative VV charge-only response swap on pinned xenon inputs."""
import argparse
import hashlib
import io
import json
import lzma
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

PIN = 'e045901c684419ba70af9d39aacc434dbbbfbb33'
UPSTREAM_P3 = 'dd9a7387814cd157d9cde739fa80181e4cd7d428'
AFH_DATA = 'nuclear/afh/data/native_v6_2026-11-22/native'
SN_DIR = 'nuclear/sn100pn'
SOURCE = 'benchmarks/vv_recovery/source_snapshot'
P3 = 'research/three_body_full_current_nucleon_match'
F1, F2, MN = 1., -.184, .938918
MASS_U = .93149410242
TWOJPLUSONE = {129:2,131:4}
CHANNELS = ('M','PhiPP','MPhiPP')


def nuclear_transfer(mass, er):
    if not mass > 0 or not er >= 0:
        raise ValueError('nonphysical elastic mass or recoil')
    return math.sqrt(2*mass*er + er*er), -2*mass*er


def charge_coefficients(q,t,f1,f2,mass):
    if not q >= 0 or not t <= 0 or not mass > 0:
        raise ValueError('nonphysical transfer or nucleon mass')
    return np.array([f1+t*(f1+2*f2)/(8*mass*mass), q*q*(f1+2*f2)/(4*mass*mass)])


def contract(c,f,L):
    if not L >= 0:
        raise ValueError('charge tensor factor must be nonnegative')
    m,p = c
    return L*np.array([m*m*f[0],p*p*f[1],2*m*p*f[2]])


def dirac_pauli_check():
    # Exact normalized positive-energy Dirac spinors; this check uses no P3 helper.
    sig = np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
    def sigma(v):
        return np.einsum('a,aij->ij',v,sig)
    def energy(p):
        return math.sqrt(MN**2+np.dot(p,p))
    def spinor(p):
        e=energy(p)
        return math.sqrt((e+MN)/(2*e))*np.vstack((np.eye(2),sigma(p)/(e+MN)))
    rows=[]
    for lam in (1.,.5,.25,.125):
        q=lam*np.array([.12,-.07,.15])
        k=lam*np.array([.06,.11,-.03])
        pi,pf=k-q/2,k+q/2
        t=(energy(pf)-energy(pi))**2-np.dot(q,q)
        z=np.zeros((2,2))
        sq=sigma(q)
        # gamma^0 times the temporal Dirac-Pauli vertex, with q_N=p_f-p_i.
        operator=F1*np.eye(4)+(F2/(2*MN))*np.block([[z,sq],[-sq,z]])
        exact=spinor(pf).conj().T@operator@spinor(pi)
        cm,cphi=charge_coefficients(np.linalg.norm(q),t,F1,F2,MN)
        spin=1j*cphi*sigma(np.cross(q,k))/np.dot(q,q)
        correct=cm*np.eye(2)+spin
        wrong=cm*np.eye(2)-spin
        rows.append({'lambda':lam,'error':float(np.linalg.norm(exact-correct)),
                     'wrong_sign_error':float(np.linalg.norm(exact-wrong))})
    return {'rows':rows,
            'consecutive_error_ratios':[a['error']/b['error'] for a,b in zip(rows,rows[1:])],
            'wrong_spin_orbit_sign_final_error_over_correct_error':rows[-1]['wrong_sign_error']/rows[-1]['error']}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def checksum_map(path):
    rows=[l.split() for l in path.read_text().splitlines() if l and not l.startswith('#')]
    return {name:digest for digest,name in rows}


def load_tables(root,A,inputs):
    def get(rel,expected):
        raw=(root/rel).read_bytes()
        got=sha(raw)
        if got!=expected:
            raise ValueError(f'input SHA mismatch: {rel}')
        inputs[str(rel)]=got
        return raw
    af_dir=Path(AFH_DATA)
    af_name=f'afh_v6_xe{A}_eight_response_native.csv'
    # Already provenance-pinned in each rank-zero run and the frozen AFH archive.
    p0=Path(SN_DIR)/'xe129-xe131-rank0'
    meta=json.loads((root/p0/'native/summary.json').read_text())
    af_raw=get(af_dir/af_name,meta['isotopes'][str(A)]['input_sha256'][af_name])
    af=np.genfromtxt(io.BytesIO(af_raw),names=True,delimiter=',')
    if A==129:
        stage=p0
        parts=['native/xe129_rank0.csv.xz']
        raw_name='xe129_rank0.csv'
    else:
        stage=Path(SN_DIR)/'xe131-rank2'
        prefix='native/xe131_rank2.csv.xz.part-'
        parts=[prefix+'00',prefix+'01']
        raw_name='xe131_rank2.csv'
    compressed_hashes=checksum_map(root/stage/'SHA256SUMS.txt')
    raw=lzma.decompress(b''.join(get(stage/p,compressed_hashes[p]) for p in parts))
    raw_manifest=checksum_map(root/stage/'RAW_CSV_SHA256SUMS.txt')
    expected_raw=next(v for k,v in raw_manifest.items() if k.endswith(raw_name))
    if sha(raw)!=expected_raw:
        raise ValueError('decompressed SN100PN table mismatch')
    inputs[f'{stage}/{raw_name} (decompressed)']=expected_raw
    sn=np.genfromtxt(io.BytesIO(raw),names=True,delimiter=',')
    q=af['q_GeV']
    assert np.array_equal(q,sn['q_GeV'])
    assert q.size==2401 and q[0]==0 and q[-1]==.6
    ref=(16*np.pi/TWOJPLUSONE[A])*np.column_stack([af[f'W_{k}00_AFH'] for k in CHANNELS])
    if A==131:
        alt=np.column_stack([sn[f'SN100PN_F_{k}_eqpn'] for k in CHANNELS])
    else:
        m=sn['SN100PN_rank0_amp_Mp']+sn['SN100PN_rank0_amp_Mn']
        p=sn['SN100PN_rank0_amp_PhiPPp']+sn['SN100PN_rank0_amp_PhiPPn']
        alt=np.column_stack((m*m,p*p,m*p))
    return q,ref,alt


def sample(q,grid,rows):
    if not grid[0] <= q <= grid[-1]:
        raise ValueError('extrapolation not allowed')
    f=np.array([np.interp(q,grid,c) for c in rows.T])
    # A common linear interpolation of Gram entries is PSD; do not use separate fits.
    if min(f[0],f[1]) < 0 or f[0]*f[1]-f[2]**2 < -1e-10:
        raise ValueError('non-positive charge response')
    return f


def replay(root,out):
    def git(*cmd):
        return subprocess.check_output(['git','-C',root,*cmd],text=True).strip()
    if git('rev-parse','HEAD')!=PIN or git('status','--short'):
        raise ValueError(f'a clean exact checkout of {PIN} is required')
    out.mkdir(parents=True,exist_ok=True)
    inputs={}
    sroot=root/SOURCE
    s_hashes=checksum_map(sroot/'SHA256SUMS.txt')
    for file in (
        f'{P3}/three_body_phase_space_tensor.py',f'{P3}/three_body_multipole_leptonic_weights.py',
        f'{P3}/nucleon_vector_form_factors.py',
        'research/three_body_calculation_and_validation/eq6_kernel.py',
        'research/three_body_spin_resolved_current/spin_resolved_three_body_current.py',
        'research/vector_current_response_audit/vector_current_response_audit.py'):
        h=sha((sroot/file).read_bytes())
        assert h==s_hashes[file]
        inputs[f'{SOURCE}/{file}']=h
    sys.path.insert(0,str(sroot/P3))
    import three_body_multipole_leptonic_weights as p3weights
    import three_body_phase_space_tensor as p3t
    check=dirac_pauli_check()
    assert 15<min(check['consecutive_error_ratios'])<=max(check['consecutive_error_ratios'])<17
    results=[]
    tensors={}
    for A in (129,131):
        grid,af,sn=load_tables(root,A,inputs)
        for name in ('green_450keV','near_pole_362keV_g2_1'):
            cfg=dict(next(case for case in p3t.vca.DEFAULT_CASES if case['benchmark']==name))
            del cfg['benchmark']
            coarse=p3t.integrate_spin_resolved_tensor(
                mchi=p3t.MCHI,mA=A*MASS_U,width_scheme='pole',**cfg)
            fine_cfg=dict(cfg)
            fine_cfg.update(nw=cfg['nw']+16,nz=24,nphi=12)
            fine=p3t.integrate_spin_resolved_tensor(
                mchi=p3t.MCHI,mA=A*MASS_U,width_scheme='pole',**fine_cfg)
            tag=f'Xe{A}_{name}'
            tensors[tag]={'baseline_grid':coarse,'refined_grid':fine}
            er=cfg['ER']
            q,t=nuclear_transfer(A*MASS_U,er)
            assert math.isclose(fine['qA_magnitude_GeV'],q,abs_tol=1e-15,rel_tol=0)
            cc=charge_coefficients(q,t,F1,F2,MN)
            L=p3weights.complex_from_pairs(fine['L_q_aligned_real_imag'])
            L0=p3weights.complex_from_pairs(coarse['L_q_aligned_real_imag'])
            err=float(np.linalg.norm(L-L0)/np.linalg.norm(L))
            assert coarse['phase_space_open'] and fine['phase_space_open']
            p3=p3weights.extract_afh_vector_coefficients(L,q,omega_GeV=er,f1=F1,f2=F2,mN_GeV=MN)
            effective=p3['L_C_continuity_adjusted_charge_weight']
            coarse_LC,_=p3weights.continuity_adjusted_coulomb_weight(L0,q,er)
            active_errors={'L00':float(abs(L[0,0]-L0[0,0])/abs(L[0,0])),
                           'conditional_LC':abs(effective-coarse_LC)/effective}
            assert max(active_errors.values())<1e-8, (tag, active_errors)

            vv=np.array([p3['C_M'],p3['C_PhiPP'],p3['C_MPhiPP']])
            np.testing.assert_allclose(vv,contract(cc,np.ones(3),effective),rtol=2e-15,atol=0)
            np.testing.assert_allclose([t,t],[p3['tA_GeV2'],p3['tN_GeV2']],rtol=2e-15,atol=0)
            assert max(coarse['dirac_trace_crosscheck']['max_relative_residual'],
                       fine['dirac_trace_crosscheck']['max_relative_residual'])<5e-9, tag
            r={'isotope':A,'benchmark':name,'ER_keV':er*1e6,'Q_A_GeV':q,
               'Q_A2_GeV2':q*q,'t_A_equals_t_N_GeV2':t,'c_M':cc[0],'c_PhiPP':cc[1],
               'L_q_frame_refinement_relative_Frobenius':err,
               'charge_weights_refinement_relative':active_errors,
               'L00':float(L[0,0].real),'conditional_L_C':effective,'models':{}}
            for label,response in [('AFH_BD',af),('SN100PN',sn)]:
                f=sample(q,grid,response)
                terms=contract(cc,f,L[0,0].real)
                cond=contract(cc,f,effective)
                if sum(terms) < 0:
                    raise ValueError('negative VV charge kernel')
                r['models'][label]={'F_equal_pn':f.tolist(),'L00_terms':terms.tolist(),
                    'L00_charge_kernel':float(sum(terms)),
                    'conditional_LC_terms':cond.tolist(),
                    'conditional_LC_charge_kernel':float(sum(cond)),
                    'M_only_kernel':float(terms[0]),
                    'signed_charge_correction_over_M':float((terms[1]+terms[2])/terms[0])}
            a,b=r['models']['AFH_BD'],r['models']['SN100PN']
            for key in ('L00_charge_kernel','conditional_LC_charge_kernel','M_only_kernel'):
                r[f'SN100PN_over_AFH_BD_{key}']=b[key]/a[key]
            assert abs(r['SN100PN_over_AFH_BD_L00_charge_kernel']-
                       r['SN100PN_over_AFH_BD_conditional_LC_charge_kernel'])<1e-14
            results.append(r)
    header={'scope':'charge-only fixed radiative VV diagnostic; no cross section',
            'source_commit':PIN,'historical_P3_source':UPSTREAM_P3,'inputs_sha256':inputs,
            'q_transfer':'impulse q_N=q_A, omega=E_R, Q_A^2=2 m_A E_R+E_R^2, t_N=t_A=-2 m_A E_R',
            'form_factors':{'Fhat1_anchor':F1,'Fhat2_anchor':F2,'mN_GeV':MN},
            'response_normalization':'physical equal-p/n, native AFH raw spin factor removed once, no C0 re-normalization',
            'dirac_pauli_small_momentum':check,'benchmark_rows':results}
    (out/'summary.json').write_text(json.dumps(header,indent=2)+'\n')
    (out/'radiative_tensors.json').write_text(json.dumps(tensors,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    replay(args.repo.resolve(),args.output)
