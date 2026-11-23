#!/usr/bin/env python3
"""Exact exported-decimal trace diagnostic; no renormalization or solver invocation."""
from __future__ import annotations
from decimal import Decimal as D, localcontext
import argparse
from pathlib import Path
import re
import hashlib
import json

PIN='1395f0c161905341fc89a2978148b1ce7c9ee1dd'
IDS={
 128:('e587c56404db3662229ede647b8bf98e2a1d9874','927ec47c973d9f2196e0f3bdaa5fd8e12d85bba6ec793eceb680a527532eb5c9'),
 129:('924fa20a1d5e9eac2291a85f55ddff4360cd8e25','26ef95e96f79bcd883cc3374e26aac62fdb80be23439e1604a96703468e8dd91'),
 130:('3d674a3c098aefc811df4bba430267dc31f77fed','0b773a800e1af18a4e3a3d2cccce77be888c7ffe7408e4d032fcaf33a0067c6e'),
 131:('972bea0196ed283b6bf198a856035965d9784ab9','2dd36c7aa74a88f8fc979b2423ac2970a605300d0eca973d94cb139bda438bd0'),
 132:('31aec6219cc83931b2572e04fde9865e5ad97495','16c21f2f14b4f066a90179e1a861b465b11003a551938d4a8250003c13b8e740'),
 134:('6e5897f6a54abe144d7ccb6fba0f1a82736c47eb','d2e6b390df90063493abffd9fd63549ee648876e0462ff959535fea13964d7c5'),
 136:('a9dfbc2171f620d6cdd800108abc9ed1b0fdac6a','363c5ba7f547938a8be4def9c05fe135ba05635165d140dbc17f5052c499a2a6')}
TEX='00d6b232355a140d70b64fae5ef01c5a58203866b3137d9bae6a61e8ec21cc38'
VAL={(4,7):'1g7/2',(4,5):'2d5/2',(4,3):'2d3/2',(4,1):'3s1/2',(5,11):'1h11/2'}
PAT=re.compile(r'THE ONE-BODY DENSITY MATRIX FOR 2JO = (\d+) TO =, (\d+)')
Z=54
def hashes(b):
 return hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest(),hashlib.sha256(b).hexdigest()
def quantum(A):return (D(A-2*Z)/2,{129:2,131:4}.get(A,1))
def last_place(x):return D(10)**x.as_tuple().exponent
def parse(raw):
 out={};h=None
 for line in raw.decode('ascii').splitlines():
  if (m:=PAT.fullmatch(line)): h=tuple(map(int,m.groups()));continue
  if h[0]!=0:continue
  dot,i,j,k,l,x=line.split()
  assert dot=='.' and (i,j)==(k,l)
  b=out.setdefault(h[1],{})
  key=(int(i),int(j))
  assert key not in b
  b[key]=D(x)
 assert len(out[0])==16 and set(out[2])==set(VAL) and {0,2}==set(out)
 assert set(out[2])<=set(out[0])
 return out[0],out[2]
def occupation(psi0,psi1,g,t,d):
 c0=(2*D(g)/d).sqrt()
 c1=(6*D(g)/d*t/(t+1)).sqrt()
 S=c0*psi0
 Dn=c1*psi1 if psi1 is not None else D(0)
 u0=c0*last_place(psi0)
 u1=c1*last_place(psi1) if psi1 is not None else D(0)
 return (S-Dn)/2,(S+Dn)/2,u0,u1,c0,c1
def iref(tex,A):
 # These are separate results in sections II.G and III.G of the supplement, not inferred from the density file.
 mark=f'\\subsection*{{$^{{{A}}}$Xe}}'
 bodies=tex.split(mark)
 assert len(bodies)==4
 coeff={k:D(m) for k,m in re.findall(r'F_\{M\}\^\{\(([pn],[pn])\)\}\s*=\s*e\^\{-2\s*y\}\s*\(\s*([.\d]+)',bodies[2].split('\\subsection')[0])}
 pairs=re.findall(r'([^\s&]+)\s*&\s*([^\s&]+)\s*\\\\',bodies[3].split('\\end{tabular}')[0])
 na={lab:[D(p),D(n)] for lab,(p,n) in zip(VAL,pairs[-5:])}
 assert len(coeff)==3 and len(na)==5
 return coeff,na
def one(A, raw, tex):
 t,spin=quantum(A);d=D(spin)*(2*t+1)
 tt=[D(Z),D(A-Z)]
 observed=parse(raw)
 core=[];val=[]
 response,printed=iref(tex,A)
 for k,x in observed[0].items():
  y=observed[1].get(k)
  p,n,u0,u1,w0,w1=occupation(x,y,k[1]+1,t,d)
  record={'orbit':VAL.get(k,'closed core'),'N_major_2j':k,'psi0':x,'psi1':y,'p':p,'n':n,
          'full_ulp_isoscalar':u0,'full_ulp_isovector':u1,'sum_projection':w0,'difference_projection':w1}
  if y is None:
   record['exact_closed_shell_occupation_per_species']=k[1]+1
   core.append(record)
  else:
   record['printed_species_occupations']=printed[k]
   record['raw_minus_printed']=[p-printed[k][0],n-printed[k][1]]
   record['inside_printed_half_last_place']=[abs(a-b)<=last_place(b)/2 for a,b in zip((p,n),printed[k])]
   val.append(record)
 def trace(rows,ex):
  pp,nn,us,uv=(sum(r[k] for r in rows) for k in ('p','n','full_ulp_isoscalar','full_ulp_isovector'))
  out={'observed_p_n':[pp,nn],'expected_p_n':ex,'observed_minus_expected_p_n':[pp-ex[0],nn-ex[1]],
        'observed_minus_expected_p_plus_n':pp+nn-ex[0]-ex[1],'observed_minus_expected_n_minus_p':nn-pp-ex[1]+ex[0],
        'round_to_nearest_species_each_max':(us+uv)/4,'truncate_last_place_species_each_max':(us+uv)/2,
        'round_to_nearest_isoscalar_max':us/2,'round_to_nearest_isovector_max':uv/2}
  return out
 core_t=trace(core,[50,50]);val_t=trace(val,[4,A-104]);all_t=trace(core+val,tt)
 n_p,n_n=all_t['observed_p_n']
 meas={'p,p':n_p*n_p,'p,n':n_p*n_n,'n,n':n_n*n_n}
 ideal={'p,p':tt[0]*tt[0],'p,n':tt[0]*tt[1],'n,n':tt[1]*tt[1]}
 comp={}
 for key,obs in meas.items():
  base=response[key];hd=D(10)**(base.adjusted()-2)/2
  comp[key]={'from_export':obs,'ideal':ideal[key],'difference_fraction':obs/ideal[key]-1,'printed':base,
             'printed_assumed_three_sig_half_step':hd,'export_in_three_sig_rounding_band':abs(obs-base)<=hd,
             'ideal_in_three_sig_rounding_band':abs(ideal[key]-base)<=hd}
 # Fix tau=0 alone leaves the isovector trace unchanged; any fitted multiplicative shift is a diagnostic only.
 f0=(D(A-100)/(sum(val_t['observed_p_n'])))
 f1=D(A-2*Z)/(val_t['observed_p_n'][1]-val_t['observed_p_n'][0])
 return {'twiceJi_plus_one':spin,'T_equals_abs_M_T':t,'all':all_t,'valence':val_t,'core':core_t,
       'inferred_counterfactual_valence_scales_iso0_iso1':[f0,f1],'orbits':core+val,'F_M_y0':comp,
       'valence_only_printed_counts_sum':[sum(x[i] for x in printed.values()) for i in (0,1)],
       'number_sum_violation_is_not_outside_model_space':True}
def run(indir,texpath):
 texraw=texpath.read_bytes();assert hashlib.sha256(texraw).hexdigest()==TEX
 tex=texraw.decode()
 results={}
 with localcontext() as z:
  z.prec=55
  for A,h in IDS.items():
   path=indir/f'{A}Xe_OBDME_SN100PN.dat'
   raw=path.read_bytes();assert hashes(raw)==h
   results[A]={'upstream_Git_blob':h[0],'upstream_SHA256':h[1],**one(A,raw,tex)}
  for A in (129,131):
   r=results[A]
   assert abs(r['all']['observed_minus_expected_p_n'][1])>1000*r['all']['round_to_nearest_species_each_max']
   assert abs(r['valence']['observed_minus_expected_p_plus_n'])>100000*r['valence']['round_to_nearest_isoscalar_max']
   assert abs(r['valence']['observed_minus_expected_n_minus_p'])>100000*r['valence']['round_to_nearest_isovector_max']
  n5=next(x for x in results[129]['orbits'] if x['orbit']=='2d5/2')
  assert abs(n5['n']-6)>10000*(n5['full_ulp_isoscalar']+n5['full_ulp_isovector'])
  assert results[131]['F_M_y0']['n,n']['export_in_three_sig_rounding_band']
  assert not results[131]['F_M_y0']['n,n']['ideal_in_three_sig_rounding_band']
 return {'upstream_author_ref':PIN,'supplement_tex_SHA256':TEX,
         'projection':'Appendix B, Eq. 14, 16-18; exact half-integer Ji from paper and T=abs(M_T).',
         'rounding_bounds':'linear propagated last decimal place (both .5 ULP nearest and 1 ULP truncation); no alteration of density entries',
         'isotopes':results}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('inputs',type=Path);p.add_argument('suptex',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 j=run(a.inputs,a.suptex)
 a.out.write_text(json.dumps(j,indent=2,default=lambda x:str(x) if isinstance(x,D) else None,sort_keys=True)+'\n')
 for A,iso in j['isotopes'].items():
  print(A,'counts',*[f'{n:.11f}' for n in iso['all']['observed_p_n']],
        'full species bound',f"{iso['all']['round_to_nearest_species_each_max']:.5g}",
        'valence species bound',f"{iso['valence']['round_to_nearest_species_each_max']:.5g}")
