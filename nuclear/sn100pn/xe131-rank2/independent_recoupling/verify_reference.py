#!/usr/bin/env python3
"""Replay an independent all-grid AFH-source Xe-131 response comparison."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.request import urlopen
import numpy as np

from afh_recoupling import operators

HASHES = {
    "131Xe_OBDME_SN100PN.dat": "2dd36c7aa74a88f8fc979b2423ac2970a605300d0eca973d94cb139bda438bd0",
    "dmformfactor-V6.m": "dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372",
    "afh_v6_xe131_eight_response_native.csv": "585ab0f483be0d8131566cefc4cb420c6db7e07bebbcd2b887da568cda2f4ffb",
}
RAW = "https://raw.githubusercontent.com/"
URLS = {
    "131Xe_OBDME_SN100PN.dat": RAW + "Raghda-AK/DM-Form-Factor/1395f0c161905341fc89a2978148b1ce7c9ee1dd/131Xe_OBDME_SN100PN.dat",
    "dmformfactor-V6.m": RAW + "lm846802297-bit/LZ-3-body-for-hdm/aac9b33d66f6a7fe1ca5d5d8cc648f5afe52a81b/v6/dmformfactor-V6.m",
    "afh_v6_xe131_eight_response_native.csv": RAW + "lm846802297-bit/radiative-dark-matter-nuclear-currents/b5b053988b1d159283540886f5169e7878eac38f/nuclear/afh/data/native_v6_2026-11-22/native/afh_v6_xe131_eight_response_native.csv",
}

def load(base,name):
    p=base/name
    if p.exists():
        raw=p.read_bytes()
    else:
        with urlopen(URLS[name],timeout=30) as r:
            raw=r.read()
    if hashlib.sha256(raw).hexdigest()!=HASHES[name]:
        raise ValueError(f"source checksum mismatch: {name}")
    if not p.exists():
        p.write_bytes(raw)
    return raw.decode("ascii" if name.endswith(".dat") else "utf-8")


def sn_table(source):
    rows={}
    for line in source.splitlines():
        if line.startswith("THE"):
            two_j,two_t=map(int,re.fullmatch(r"THE ONE-BODY DENSITY MATRIX FOR 2JO = (\d+) TO =, (\d+)", line).groups())
            key=(two_j//2,two_t//2)
            if two_j%2 or two_t not in (0,2) or key in rows:
                raise ValueError("invalid density header")
            rows[key]=[]
        else:
            mark,*arr,v=line.split()
            if mark!=".":
                raise ValueError("invalid source entry")
            rows[key].append((tuple(map(int,arr)),float(v)))
    if {k:len(rows[k]) for k in [(0,0),(0,1),(2,0),(2,1)]} != {(0,0):16,(0,1):5,(2,0):14,(2,1):14}:
        raise ValueError("unexpected density block sizes")
    return rows


def afh_table(source):
    rows={}
    pat=r"sdgXe131\[(\d+)\]=\{[^=]+?=([-+0-9.]+),(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)\};"
    for i,v,j,t,*arr in re.findall(pat,source):
        rows.setdefault((int(j),int(t)),[]).append((tuple(map(int,arr)),float(v)))
    if sum(map(len,rows.values()))!=67:
        raise ValueError("unexpected AFH block size")
    return rows


def amplitude(y,j,t,table):
    terms=[v*np.array(operators(y,*arr,j)) for arr,v in table[j,t]]
    return [math.fsum(c) for c in zip(*terms)]


def form_factors(amp):
    m,p=amp[...,0],amp[...,1]
    return np.stack([m*m,p*p,m*p],axis=-1)


def run(base):
    base.mkdir(exist_ok=True,parents=True)
    sn=sn_table(load(base,"131Xe_OBDME_SN100PN.dat"))
    af=afh_table(load(base,"dmformfactor-V6.m"))
    rows=list(csv.DictReader(load(base,"afh_v6_xe131_eight_response_native.csv").splitlines()))
    q = np.array([float(r["q_GeV"]) for r in rows])
    if len(rows)!=2401 or not np.max(np.abs(q-np.arange(2401)/4000))<2e-15:
        raise ValueError("AFH grid mismatch")
    b=math.sqrt(41.467/(45*131**(-1/3)-25*131**(-2/3)))*5.0677
    y=(q*b/2)**2
    if np.max(np.abs(y-np.array([float(r["y"]) for r in rows])))>1e-13:
        raise ValueError("HO variable mismatch")
    rec = []
    std = math.sqrt(8*math.pi/(4*24))
    iv = math.sqrt(24*math.pi*11.5/(4*24*12.5))
    for rank in (0,2):
        s=std*np.array([amplitude(x,rank,0,sn) for x in y])
        v=iv*np.array([amplitude(x,rank,1,sn) for x in y])
        old=std*np.array([amplitude(x,rank,0,af) for x in y])
        rec.append(((s-v)/2,(s+v)/2,old))
    old_sum=sum(form_factors(r[2]) for r in rec)
    old_native=4*math.pi*np.array([[float(r[key]) for key in ("W_M00_AFH","W_PhiPP00_AFH","W_MPhiPP00_AFH")] for r in rows])
    err=np.abs(old_sum-old_native)
    limit=1e-8 + 1e-9*np.abs(old_native)
    if np.max(err/limit)>1:
        raise ValueError(f"rank-summed B&D AFH native closure failed: max {err.max():.7g}, ratio {(err/limit).max():.7g}")
    (p,n,old2)=rec[1]
    sn0=form_factors(rec[0][0]+rec[0][1])
    sn2=form_factors(p+n)
    old2r=form_factors(old2)
    landmarks={}
    for idx in (0,400,600,627,650,800,1200,2400):
        landmarks[f"{q[idx]:.5f}"] = {
            "SN100PN_rank2_signed_Mp_Mn_Phip_Phin": [p[idx,0],n[idx,0],p[idx,1],n[idx,1]],
            "SN100PN_rank0_equal_pn_response_M_Phi_Int": sn0[idx].tolist(),
            "SN100PN_rank2_equal_pn_response_M_Phi_Int": sn2[idx].tolist(),
            "AFH_BD_rank2_equal_pn_response_M_Phi_Int": old2r[idx].tolist(),
        }
    data = {"inputs_sha256": HASHES,"landmarks":landmarks,"q_nodes":len(rows),
            "status":"AFH-source recoupling verified against immutable native rank-summed B&D tables; independent SN100PN signed rank-two reference",
            "afh_bd_totals_max_native_abs_error_M_Phi_Int":err.max(axis=0).tolist(),
            "afh_bd_totals_max_native_tolerance_utilization":float(np.max(err/limit))}
    return q,np.column_stack([p[:,0],n[:,0],p[:,1],n[:,1]]),data


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--inputs",type=Path,required=True)
    a.add_argument("--json",type=Path,required=True)
    a.add_argument("--signed-csv",type=Path)
    a.add_argument("--candidate-csv",type=Path)
    arg=a.parse_args()
    q,amp,out=run(arg.inputs)
    if arg.signed_csv:
        np.savetxt(arg.signed_csv,np.column_stack((q,amp)),fmt="%.17g",delimiter=",",header="q_GeV,Mp,Mn,Phip,Phin",comments="")
    if arg.candidate_csv:
        rows=np.genfromtxt(arg.candidate_csv,delimiter=",",names=True)
        keys=("J2_amp_Mp","J2_amp_Mn","J2_amp_PhiPPp","J2_amp_PhiPPn")
        a2=np.column_stack([rows[k] for k in keys])
        if not np.array_equal(q,rows["q_GeV"]) or a2.shape != amp.shape or not np.all(np.isfinite(a2)):
            raise ValueError("candidate grid/shape/finiteness failed")
        error=np.abs(a2-amp)
        if np.max(error)>5e-13:
            raise ValueError(f"signed rank-2 independent comparison failed: {error.max()}")
        out["candidate_validation"]={
            "raw_csv_sha256":hashlib.sha256(arg.candidate_csv.read_bytes()).hexdigest(),
            "columns":list(keys),"signed_amplitude_abs_tolerance":5e-13,
            "max_abs_error_each":error.max(axis=0).tolist(),
            "max_abs_error_all":float(error.max()),
        }
    arg.json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    main()
