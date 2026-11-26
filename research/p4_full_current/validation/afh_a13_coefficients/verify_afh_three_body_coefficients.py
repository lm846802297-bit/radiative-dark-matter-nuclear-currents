#!/usr/bin/env python3
from __future__ import annotations
import json, math, platform, sys
import numpy as np

def a13_spatial(T, Q, mN, f1, f2):
    q_afh = np.array([0.0, 0.0, -Q], dtype=float)
    qhat = q_afh / Q
    Tvec = np.asarray(T[1:4], dtype=np.complex128)
    ell_M = -f1 * Tvec
    ell_5 = 1j * (f1 + f2) / (2.0 * mN) * np.cross(q_afh, Tvec)
    q_over_m = q_afh / mN
    r_delta = 0.5 * ((Q * Q / (mN * mN)) * np.vdot(ell_M, ell_M).real - abs(np.dot(q_over_m, ell_M)) ** 2)
    r_sigma = 0.5 * (np.vdot(ell_5, ell_5).real - abs(np.dot(qhat, ell_5)) ** 2)
    r_sigma_delta = float(np.dot(q_over_m, np.real(1j * np.cross(ell_M, np.conjugate(ell_5)))))
    r_sigma_double = abs(np.dot(qhat, ell_5)) ** 2
    ell_E = np.array([0.0, 0.0, Q * (f1 + 2.0 * f2) * T[0] / (4.0 * mN)], dtype=np.complex128)
    r_phi_tilde = 0.5 * ((Q * Q / (mN * mN)) * np.vdot(ell_E, ell_E).real - abs(np.dot(q_over_m, ell_E)) ** 2)
    return r_delta, r_sigma, r_sigma_delta, r_sigma_double, r_phi_tilde

def closed_spatial(T, Q, mN, f1, f2):
    t_perp = abs(T[1]) ** 2 + abs(T[2]) ** 2
    return (f1*f1*Q*Q/(2*mN*mN)*t_perp,
            (f1+f2)**2*Q*Q/(8*mN*mN)*t_perp,
            -f1*(f1+f2)*Q*Q/(2*mN*mN)*t_perp)

def main():
    rng = np.random.default_rng(20261125)
    mN = 0.938918
    n_single = 20000
    max_spatial = {"Delta":0.0,"SigmaPrime":0.0,"SigmaPrimeDelta":0.0}
    max_zero = {"SigmaDoublePrime":0.0,"PhiTildePrime":0.0}
    max_cross_identity = max_charge_reconstruction = max_charge_factor_two = 0.0
    for _ in range(n_single):
        Q=float(rng.uniform(0.02,0.6)); f1=float(rng.uniform(0.1,1.4)); f2=float(rng.uniform(-0.4,0.2))
        T=rng.normal(size=4)+1j*rng.normal(size=4)
        direct=a13_spatial(T,Q,mN,f1,f2); closed=closed_spatial(T,Q,mN,f1,f2)
        for name,a,b in zip(max_spatial,direct[:3],closed):
            rel=abs(a-b)/max(abs(a),abs(b),1.0); max_spatial[name]=max(max_spatial[name],float(rel))
        max_zero["SigmaDoublePrime"]=max(max_zero["SigmaDoublePrime"],float(direct[3]))
        max_zero["PhiTildePrime"]=max(max_zero["PhiTildePrime"],float(direct[4]))
        if f1*(f1+f2)>0:
            identity=-2.0*math.sqrt(max(direct[0]*direct[1],0.0))
            rel=abs(direct[2]-identity)/max(abs(direct[2]),abs(identity),1.0)
            max_cross_identity=max(max_cross_identity,float(rel))
        T0=complex(*rng.normal(size=2)); cM=complex(*rng.normal(size=2)); cPhi=complex(*rng.normal(size=2))
        a=cM*T0; b=cPhi*T0; full=abs(a+b)**2
        rM=abs(a)**2; rPhi=abs(b)**2; rMPhi=2.0*(a*np.conjugate(b)).real
        rel=abs(full-(rM+rPhi+rMPhi))/max(abs(full),abs(rM+rPhi+rMPhi),1.0)
        max_charge_reconstruction=max(max_charge_reconstruction,float(rel))
        explicit_cross=2.0*(cM*np.conjugate(cPhi)).real*abs(T0)**2
        rel2=abs(rMPhi-explicit_cross)/max(abs(rMPhi),abs(explicit_cross),1.0)
        max_charge_factor_two=max(max_charge_factor_two,float(rel2))
    n_gram=2000; max_gram_residual=0.0; max_hermiticity_residual=0.0; min_gram_eigenvalue=math.inf
    for _ in range(n_gram):
        Q=float(rng.uniform(0.02,0.6)); f1=float(rng.uniform(0.1,1.4)); f2=float(rng.uniform(-0.4,0.2))
        X=rng.normal(size=(4,9))+1j*rng.normal(size=(4,9)); L=X@X.conj().T
        max_hermiticity_residual=max(max_hermiticity_residual,float(np.linalg.norm(L-L.conj().T)/np.linalg.norm(L)))
        min_gram_eigenvalue=min(min_gram_eigenvalue,float(np.linalg.eigvalsh(L).min()))
        direct_sum=np.zeros(3)
        for a in range(X.shape[1]): direct_sum += np.asarray(a13_spatial(X[:,a],Q,mN,f1,f2)[:3])
        LT=float((L[1,1]+L[2,2]).real)
        predicted=np.asarray([f1*f1*Q*Q/(2*mN*mN)*LT,(f1+f2)**2*Q*Q/(8*mN*mN)*LT,-f1*(f1+f2)*Q*Q/(2*mN*mN)*LT])
        rel=np.max(np.abs(direct_sum-predicted)/np.maximum(np.maximum(np.abs(direct_sum),np.abs(predicted)),1.0))
        max_gram_residual=max(max_gram_residual,float(rel))
    report={
      "status":"PASS_INDEPENDENT_A13_COEFFICIENT_IDENTITY",
      "scope":"leading-NR one-body Dirac-Pauli vector current; AFH Appendix A Eq. (A13)",
      "seed":20261125,"samples_single_current":n_single,"samples_complex_hermitian_psd_tensors":n_gram,
      "max_relative_residuals":{"spatial_closed_form":max_spatial,"charge_sector_reconstruction":max_charge_reconstruction,"charge_interference_factor_two":max_charge_factor_two,"SigmaPrimeDelta_identity":max_cross_identity,"general_complex_Gram_tensor":max_gram_residual},
      "zero_projection_max_absolute":max_zero,
      "Gram_tensor_checks":{"max_hermiticity_relative_residual":max_hermiticity_residual,"minimum_eigenvalue":min_gram_eigenvalue},
      "python":sys.version,"numpy":np.__version__,"platform":platform.platform()}
    print(json.dumps(report,indent=2))
    checks=[max(max_spatial.values())<2e-14,max(max_zero.values())<2e-14,max_cross_identity<2e-14,max_charge_reconstruction<2e-14,max_charge_factor_two<2e-14,max_gram_residual<2e-14,max_hermiticity_residual<2e-14,min_gram_eigenvalue>-1e-10]
    if not all(checks): raise SystemExit("one or more independent A13 checks failed")

if __name__ == "__main__": main()
