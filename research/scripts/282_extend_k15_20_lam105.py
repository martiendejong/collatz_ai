"""
282_extend_k15_20_lam105.py
============================
Verify E[sigma_within(v2)]/mean_v2 > E[sigma_within(v0)]/mean_v0 for
k = 15..20, lambda = 1.05 (hardest case, smallest margin in Script 281).

Uses reduced iterations (convergence assured by spectral gap ~ 0.96):
  n_iter=100 gives error ~ 0.96^100 ~ 1.7e-2 (1.7% relative error on eigenvector)
  This is sufficient to determine ratio >> 1.08 (we need ratio > 1, have margin ~ 0.08)

For k=20: N = 3^19 ~ 1.16B -- too large. Cap at k=17 (N ~ 129M).
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 1.05
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)

CASES = [
    (15, 100),
    (16, 80),
    (17, 60),
]

print(f"  k  iters        N   E_s0/mu0   E_s2/mu2  ratio    R      F    sqrt(F)/R")
sys.stdout.flush()

for k, niters in CASES:
    N = 3**(k-1); Nl = N//3
    print(f"  {k}  {niters:>5}  {N:>9}  computing...", flush=True)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N)
    for it in range(niters):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
        if (it+1) % 20 == 0:
            print(f"    iter {it+1}/{niters}...", flush=True)
    rho = float(w.max()); t = A/rho; R = (t**2+lam)/(1+t*lam)
    F = (t**4+lam**2)/(1+t**2*lam**2)
    import math; sqFR = math.sqrt(F)/R
    Nl3 = Nl//3
    v0 = v[0::3]; v2 = v[2::3]
    j3 = np.arange(Nl3)
    c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
    c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
    E0 = c0.std(1).mean() / v0.mean()
    E2 = c2.std(1).mean() / v2.mean()
    ratio = E2/E0
    ok = 'OK' if ratio > 1 else 'FAIL'
    print(f"  {k}  {niters:>5}  {N:>9}  {E0:.7f}  {E2:.7f} {ratio:.5f} {R:.5f} {F:.4f}  {sqFR:.5f}  {ok}")
    sys.stdout.flush()

print("DONE")
