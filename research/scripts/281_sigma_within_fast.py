"""
281_sigma_within_fast.py
========================
Fast 168-case verification of E[sigma_within(v2)]/mean_v2 > E[sigma_within(v0)]/mean_v0.
Same as Script 280 but with adaptive iteration count:
  k<=10: 2000 iterations
  k=11:  1000 iterations
  k=12:   500 iterations
  k=13:   300 iterations
  k=14:   200 iterations
Spectral gap analysis (Script 275) shows |rho2/rho| <= 0.96, so after n iterations
the relative error is <= 0.96^n. At n=200: 0.96^200 ~ 2e-4 (sufficient for 3 sig figs).
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAMBDAS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
KS = list(range(3, 15))
ITERS = {3:2000,4:2000,5:2000,6:2000,7:2000,8:2000,9:2000,10:2000,
         11:1000,12:500,13:300,14:200}

def run_kl(k, lam, n_iter):
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    rho = float(w.max()); t = A/rho; R = (t**2+lam)/(1+t*lam)
    F = (t**4+lam**2)/(1+t**2*lam**2)
    return v, Nl, t, R, F

n_fail = 0; n_total = 0
min_ratio = float('inf'); min_case = None

print(f"{'lam':>5} {'k':>3} {'iters':>5} {'E_s0/mu0':>10} {'E_s2/mu2':>10} {'ratio':>7} {'R':>7} {'F':>6} {'OK':>4}")
sys.stdout.flush()

for lam in LAMBDAS:
    for k in KS:
        n_iter = ITERS[k]
        v, Nl, t, R, F = run_kl(k, lam, n_iter)
        Nl3 = Nl//3
        v0 = v[0::3]; v2 = v[2::3]
        j3 = np.arange(Nl3)
        c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
        c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
        sig_w0 = c0.std(1).mean()
        sig_w2 = c2.std(1).mean()
        mu0 = v0.mean(); mu2 = v2.mean()
        E0 = sig_w0/mu0; E2 = sig_w2/mu2
        ratio = E2/E0
        ok = 'OK' if ratio > 1 else 'FAIL'
        if ratio <= 1: n_fail += 1
        n_total += 1
        if ratio < min_ratio:
            min_ratio = ratio
            min_case = (lam, k, E0, E2, ratio, R, F)
        print(f"{lam:>5.2f} {k:>3} {n_iter:>5}  {E0:>10.7f}  {E2:>10.7f}  {ratio:>7.5f}  {R:>7.5f}  {F:>6.4f}  {ok:>4}")
        sys.stdout.flush()
    print()
    sys.stdout.flush()

print()
print(f"SUMMARY: {n_total} cases, {n_fail} FAIL, {n_total-n_fail} OK")
if min_case:
    lam, k, E0, E2, ratio, R, F = min_case
    print(f"Minimum ratio: {ratio:.6f} at lambda={lam:.2f}, k={k}")
    print(f"  E_sig_w0/mu0={E0:.7f}, E_sig_w2/mu2={E2:.7f}")
    print(f"  R={R:.5f}, F={F:.4f}, ratio>1: {ratio>1}")
print("DONE")
