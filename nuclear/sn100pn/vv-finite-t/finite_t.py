#!/usr/bin/env python3
"""Replay finite-invariant-t nucleon form factors on frozen radiative VV kernels."""
import argparse
import hashlib
import json
import lzma
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

PIN='112555411a63afb5e7a900f3748a03edfdbf3a10'
ARCHIVE=Path('nuclear/sn100pn/vv-charge-matching')
VV_SOURCE=Path('benchmarks/vv_recovery/source_snapshot')
MATCHER=VV_SOURCE/'research/three_body_full_current_nucleon_match'
MODELS=('AFH_BD','SN100PN')
SCHEMES=('anchors','bishara_linear','isoscalar_zexp')
M_N=.938918
MU=.880
KAPPA_S=-.064
SLOPES=(5.57,2.84,-.018)
T_CUT=9*.13957**2
T0=-.28
AE=(-1.809,.91,1.92,-.98)
BM=(-1.938,.78,3.71,-4.0)


def independent_ff(scheme,t):
    """Species-common effective universal charges in the external g_VN basis."""
    if t>0 or scheme not in SCHEMES:
        raise ValueError((scheme,t))
    if scheme=='anchors':
        return np.array([1.,-.184])
    if scheme=='bishara_linear':
        return np.array([1+t*math.fsum(SLOPES)/3, MU-1+KAPPA_S])
    # Independent polynomial with the q=0 constraint imposed as differences.
    def z(point):
        a,b=math.sqrt(T_CUT-point),math.sqrt(T_CUT-T0)
        return (a-b)/(a+b)
    zm,z0=z(t),z(0.)
    powers=[zm**k-z0**k for k in range(1,5)]
    ge=1+math.fsum(c*v for c,v in zip(AE,powers))
    gm=MU*(1+math.fsum(c*v for c,v in zip(BM,powers)))
    tau=-t/(4*M_N**2)
    return np.array([(ge+tau*gm)/(1+tau)+SLOPES[2]*t,
                     (gm-ge)/(1+tau)+KAPPA_S])


def coeff(q2,t,f):
    s=(f[0]+2*f[1])/(4*M_N**2)
    return np.array([f[0]+t*s/2,q2*s])


def signed_terms(c,f,L):
    f=np.asarray(f,dtype=float)
    cm,cp=c
    return f*np.array([cm**2,cp**2,2*cm*cp])*L


def separate(anchors,model):
    a,b=MODELS
    na=anchors[b]/anchors[a]
    nm=model[b]/model[a]
    fa=model[a]/anchors[a]
    fb=model[b]/anchors[b]
    combined=model[b]/anchors[a]
    np.testing.assert_allclose([na*fb,nm*fa,nm/na],[combined,combined,fb/fa],rtol=3e-15,atol=0)
    return {'nuclear_only_anchors':na,'nuclear_only_at_scheme':nm,
            'form_factor_only_AFH_BD':fa,'form_factor_only_SN100PN':fb,
            'model_by_ff_ratio_of_ratios':nm/na,'joint_vs_anchor_AFH_BD':combined}


def sha(raw): return hashlib.sha256(raw).hexdigest()


def read_sums(file):
    return {file:digest for digest,file in
            (line.split() for line in file.read_text().splitlines() if line and not line.startswith('#'))}


def checked(path,expected):
    raw=path.read_bytes()
    if sha(raw)!=expected:
        raise ValueError(f'SHA-256 mismatch {path}')
    return raw


def phase_factor(record,label):
    if label=='L00':
        return float(record['L00'])
    return float(record['conditional_L_C'])


def solve(root,output):
    rev=subprocess.check_output(['git','-C',root,'rev-parse','HEAD'],text=True).strip()
    dirty=subprocess.check_output(['git','-C',root,'status','--short'],text=True).strip()
    if rev!=PIN or dirty:
        raise ValueError("a clean exact checkout is required")
    archived=root/ARCHIVE
    frozen=read_sums(archived/'SHA256SUMS.txt')
    summed=checked(archived/'native/summary.json',frozen['native/summary.json'])
    points=checked(archived/'native/off_grid.json',frozen['native/off_grid.json'])
    compressed=checked(archived/'native/radiative_tensors.json.xz',frozen['native/radiative_tensors.json.xz'])
    tensors=lzma.decompress(compressed)
    expected_raw=read_sums(archived/'RAW_TENSOR_SHA256SUMS.txt')['native/radiative_tensors.json']
    if sha(tensors)!=expected_raw:
        raise ValueError('unpacked radiative tensor SHA mismatch')
    summary=json.loads(summed)
    offgrid=json.loads(points)
    vv=json.loads(tensors)
    sys.path.insert(0,str(root/MATCHER))
    pinmap=read_sums(root/VV_SOURCE/'SHA256SUMS.txt')
    input_hashes={str(ARCHIVE/k):frozen[k] for k in
                  ('native/summary.json','native/off_grid.json','native/radiative_tensors.json.xz')}
    input_hashes['radiative_tensors.json (uncompressed)']=expected_raw
    for name in ('nucleon_vector_form_factors.py','three_body_multipole_leptonic_weights.py'):
        path=MATCHER/name
        digest=pinmap[str(path.relative_to(VV_SOURCE))]
        checked(root/path,digest)
        input_hashes[str(path)]=digest
    import nucleon_vector_form_factors as source
    import three_body_multipole_leptonic_weights as p3
    off={(row['A'],row['E_R_keV']):row for row in offgrid['benchmark_rows']}
    cc_check=0.
    formf_check=0.
    kernel_check=0.
    anchor_check=0.
    t_check=0.
    out=[]
    for raw in summary['benchmark_rows']:
        a,er=raw['isotope'],raw['ER_keV']
        q,t=raw['Q_A_GeV'], raw['t_A_equals_t_N_GeV2']
        q2=raw['Q_A2_GeV2']
        arrays=off[a,er]
        ts=vv[f"Xe{a}_{raw['benchmark']}"]['refined_grid']
        l=p3.complex_from_pairs(ts['L_q_aligned_real_imag'])
        er_GeV=ts['ER_GeV']
        tr=source.recoil_transfer_kinematics(ts['mA_GeV'],er)
        t_check=max(t_check,abs(t-tr['tN_GeV2']),abs(-t+er_GeV**2-q2))
        np.testing.assert_allclose([q,q2,t],[tr['Q_A_GeV'],tr['Q_A2_GeV2'],tr['tN_GeV2']],rtol=2e-15,atol=0)
        results={'A':a,'ER_keV':er,'Q_A_GeV':q,'Q_A2_GeV2':q2,
                 't_A_equals_t_N_GeV2':t,'Q2_for_FF_positive':-t,
                 'L00':raw['L00'],'conditional_L_C':raw['conditional_L_C'],'cases':{}}
        for name in SCHEMES:
            fn=independent_ff(name,t)
            if name=='anchors':
                expected=np.array([source.F1N0,source.F2N0])
            else:
                ff=getattr(source,name)(-t)
                expected=np.asarray(source.effective_nucleon_form_factors(*ff))
            formf_check=max(formf_check,float(np.max(abs(fn-expected))))
            np.testing.assert_allclose(fn,expected,rtol=0,atol=8e-16)
            c=coeff(q2,t,fn)
            sourced=np.array(source.charge_coefficients(q,F1=fn[0],F2=fn[1],tN_GeV2=t))
            cc_check=max(cc_check,float(np.max(abs(c-sourced))))
            np.testing.assert_allclose(c,sourced,rtol=7e-16,atol=0)
            matched=p3.extract_afh_vector_coefficients(l,q,omega_GeV=er_GeV,f1=fn[0],f2=fn[1])
            np.testing.assert_allclose(
                [matched['C_M'],matched['C_PhiPP'],matched['C_MPhiPP']],
                signed_terms(c,np.ones(3),raw['conditional_L_C']),rtol=2e-15,atol=0)
            ff_info={'F1hat_p':fn[0],'F2hat_p':fn[1],'F1hat_n':fn[0],'F2hat_n':fn[1],
                     'c_M':c[0],'c_PhiPP':c[1],'models':{}}
            for model in MODELS:
                fs=np.array(arrays[f'exact_F_{model}'])
                gram=np.array([[fs[0],fs[2]],[fs[2],fs[1]]])
                assert np.linalg.eigvalsh(gram)[0]>-1e-12
                m={}
                for mode in ('L00','conditional_L_C'):
                    weight=phase_factor(raw,mode)
                    v=signed_terms(c,fs,weight)
                    total=float(v.sum())
                    other=float(weight*(c@gram@c))
                    err=abs(total-other)/abs(other)
                    kernel_check=max(kernel_check,err)
                    assert total>0 and err<3e-15
                    if name=='anchors':
                        key='L00_charge_kernel' if mode=='L00' else 'conditional_LC_charge_kernel'
                        before=arrays['exact_models_absolute_T3_kernels'][model][key]
                        e=abs(total/before-1)
                        anchor_check=max(anchor_check,e)
                        assert e<4e-15
                    m[mode]={'terms_M_PhiPP_MPhiPP':v.tolist(),'kernel':total,
                             'correction_PhiPP_plus_interference_over_M':float((v[1]+v[2])/v[0]),
                             'signed_interference_over_M':float(v[2]/v[0])}
                ff_info['models'][model]=m
            results['cases'][name]=ff_info
        def collapse(entry,mode):
            return {n:entry['models'][n][mode]['kernel'] for n in MODELS}
        cases=results['cases']
        anchor=collapse(cases['anchors'],'L00')
        results['effects']={name:separate(anchor,collapse(cases[name],'L00')) for name in SCHEMES}
        for name in SCHEMES:
            c0=collapse(cases['anchors'],'conditional_L_C')
            e0=separate(c0,collapse(cases[name],'conditional_L_C'))
            for key,v in results['effects'][name].items():
                assert math.isclose(v,e0[key],rel_tol=2e-15,abs_tol=0)
        out.append(results)
    if len(out)!=4:
        raise ValueError('expected four kinematic points')
    output.parent.mkdir(exist_ok=True,parents=True)
    result={'scope':'charge-only response swap with inherited invariant-t one-body inputs; no full cross section',
            'base':PIN,'inputs_sha256':input_hashes,
            'primary_sources':{'Bishara_Appendix_A_Eqs_A3_A6_A8':'https://arxiv.org/abs/1707.06998',
                               'Borah_isoscalar_Sachs_Eq16':'https://arxiv.org/abs/2003.13640'},
            'schemes':{'anchors':'constant inherited effective F1/F2',
                       'bishara_linear':'adopted first-order finite t_N Dirac, constant Pauli anchor',
                       'isoscalar_zexp':'alternate central sensitivity; finite EM isoscalar Dirac/Pauli, fixed strange Pauli, linear strange Dirac'},
            'verification':{'independent_vs_pinned_FF_max_abs':formf_check,
                            'independent_vs_pinned_coeff_max_abs':cc_check,
                            'kernel_gram_closure_max_relative':kernel_check,
                            'archived_anchors_replay_max_relative':anchor_check,
                            'exact_transfer_kinematics_max_abs_GeV2':t_check},
            'benchmarks':out}
    output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    cmd=argparse.ArgumentParser(description=__doc__)
    cmd.add_argument('--repo',type=Path,required=True)
    cmd.add_argument('--output',type=Path,required=True)
    options=cmd.parse_args()
    solve(options.repo.resolve(),options.output)
