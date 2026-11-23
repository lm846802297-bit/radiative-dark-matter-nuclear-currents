#!/usr/bin/env python3
"""Check fixed-speed charge ratios against source-level exact-momentum HO recoupling."""
import argparse
import json
import math
from pathlib import Path
import re
import sys
from urllib.request import urlopen

import numpy as np

from fixed_charge import MN, MASS_U, F1, F2, PIN, SN_DIR, TWOJPLUSONE, charge_coefficients, contract, load_tables, nuclear_transfer, sample, sha

RAW = 'https://raw.githubusercontent.com'
INPUTS = {
    '129Xe_OBDME_SN100PN.dat': ('Raghda-AK/DM-Form-Factor/1395f0c161905341fc89a2978148b1ce7c9ee1dd/129Xe_OBDME_SN100PN.dat',
        '26ef95e96f79bcd883cc3374e26aac62fdb80be23439e1604a96703468e8dd91'),
    '131Xe_OBDME_SN100PN.dat': ('Raghda-AK/DM-Form-Factor/1395f0c161905341fc89a2978148b1ce7c9ee1dd/131Xe_OBDME_SN100PN.dat',
        '2dd36c7aa74a88f8fc979b2423ac2970a605300d0eca973d94cb139bda438bd0'),
    'dmformfactor-V6.m': ('lm846802297-bit/LZ-3-body-for-hdm/aac9b33d66f6a7fe1ca5d5d8cc648f5afe52a81b/v6/dmformfactor-V6.m',
        'dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372'),
}
RECOUPLING=Path(SN_DIR)/'xe131-rank2/independent_recoupling/afh_recoupling.py'


def get(name,cache):
    p=cache/name
    upstream,digest=INPUTS[name]
    if p.exists():
        raw=p.read_bytes()
    else:
        with urlopen(f'{RAW}/{upstream}',timeout=40) as r:
            raw=r.read()
    if sha(raw)!=digest:
        raise ValueError(f'checksum mismatch: {name}')
    if not p.exists():
        p.write_bytes(raw)
    return raw.decode('ascii' if name.endswith('.dat') else 'utf-8')


def af_rows(raw,A):
    pattern=rf'sdgXe{A}\[(\d+)\]=\{{[^=]+?=([-+0-9.]+),(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)\}};'
    all_rows=re.findall(pattern,raw)
    assert len(all_rows)=={129:39,131:67}[A]
    return [(int(j),tuple(map(int,arr)),float(v))
            for _,v,j,t,*arr in all_rows if t=='0' and int(j)%2==0]


def sn_rows(raw,A):
    all_rows=[]
    for line in raw.splitlines():
        if line.startswith('THE'):
            j,t=[int(k)//2 for k in re.search(r'2JO = (\d+) TO =, (\d+)',line).groups()]
        elif t==0 and j%2==0:
            dot,*arr,v=line.split()
            assert dot=='.'
            all_rows.append((j,tuple(map(int,arr)),float(v)))
    assert len(all_rows)=={129:16,131:30}[A]
    return all_rows


def evaluate(root,cache,output,charge_summary=None):
    sys.path.insert(0,str(root/RECOUPLING.parent))
    from afh_recoupling import operators
    cache.mkdir(exist_ok=True,parents=True)
    af_source=get('dmformfactor-V6.m',cache)
    rows=[]
    for a in (129,131):
        norm=np.sqrt(8*np.pi/(TWOJPLUSONE[a]*(a-108+1)))
        b=5.0677*np.sqrt(41.467/(45*a**(-1/3)-25*a**(-2/3)))
        ranks=range(0,TWOJPLUSONE[a],2)
        def f_at(q,ob):
            y=(b*q/2)**2
            amps=norm*np.array([sum((v*np.array(operators(y,*s,j)) for j,s,v in ob if j==rank),
                                   start=np.zeros(2)) for rank in ranks])
            m,p=amps.T
            return np.array([m@m,p@p,m@p])
        astro={}
        grid,af_grid,sn_grid=load_tables(root,a,astro)
        af=af_rows(af_source,a)
        sn=sn_rows(get(f'{a}Xe_OBDME_SN100PN.dat',cache),a)
        for er in (200e-6,267e-6):
            q,t=nuclear_transfer(a*MASS_U,er)
            c=charge_coefficients(q,t,F1,F2,MN)
            lo=np.searchsorted(grid,q)-1
            for idx in [lo,lo+1]:
                for ob,table in ((af,af_grid),(sn,sn_grid)):
                    np.testing.assert_allclose(f_at(grid[idx],ob),table[idx],rtol=1e-8,atol=1e-9)
            exact=[f_at(q,ob) for ob in (af,sn)]
            interp=[sample(q,grid,table) for table in (af_grid,sn_grid)]
            x1,x2=[float(sum(contract(c,item,1))) for item in exact]
            l1,l2=[float(sum(contract(c,item,1))) for item in interp]
            rows.append({'A':a,'E_R_keV':er*1e6,'q_GeV':q,
                         'bracketing_q_GeV':grid[lo:lo+2].tolist(),
                         'exact_F_AFH_BD':exact[0].tolist(),'exact_F_SN100PN':exact[1].tolist(),
                         'exact_kernel_per_common_weight_AFH_BD':x1,
                         'exact_kernel_per_common_weight_SN100PN':x2,
                         'exact_SN100PN_over_AFH_BD':x2/x1,
                         'interpolated_SN100PN_over_AFH_BD':l2/l1,
                         'interpolation_ratio_error_relative':(l2/l1)/(x2/x1)-1,
                         'interpolation_AFH_kernel_error_relative':l1/x1-1,
                         'interpolation_SN100PN_kernel_error_relative':l2/x2-1})
    if charge_summary:
        summary=json.loads(charge_summary.read_text())
        assert summary['source_commit']==PIN
        weights={(r['isotope'],r['ER_keV']):r for r in summary['benchmark_rows']}
        for row in rows:
            r=weights[row['A'],row['E_R_keV']]
            assert math.isclose(r['Q_A_GeV'],row['q_GeV'],rel_tol=0,abs_tol=1e-15)
            c=np.array([r['c_M'],r['c_PhiPP']])
            terms={}
            for model in ('AFH_BD','SN100PN'):
                f=row[f'exact_F_{model}']
                p=contract(c,f,r['L00'])
                q=contract(c,f,r['conditional_L_C'])
                terms[model]={'L00_charge_terms':p.tolist(),'L00_charge_kernel':float(sum(p)),
                              'conditional_LC_terms':q.tolist(),'conditional_LC_charge_kernel':float(sum(q))}
            row['exact_models_absolute_T3_kernels']=terms
    info={'source_commit':PIN,'exact_recoupling_helper_sha256':sha((root/RECOUPLING).read_bytes()),
          'inputs_sha256':{name:p[1] for name,p in INPUTS.items()},
          'benchmark_rows':rows,
          'charge_summary_sha256':sha(charge_summary.read_bytes()) if charge_summary else None,
          'scope':'Common T3 charge factor cancels in response swap ratios; absolute terms use the same stored factor for both models.'}
    output.write_text(json.dumps(info,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--json',type=Path,required=True)
    p.add_argument('--charge-summary',type=Path)
    a=p.parse_args()
    evaluate(a.repo.resolve(),a.inputs,a.json,a.charge_summary)
