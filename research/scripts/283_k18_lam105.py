"""
283_k18_lam105.py
=================
Quick check: ratio at k=18, lambda=1.05 with 40 iterations.
N = 3^17 = 129M. Memory ~1GB. ~26 seconds/iter at 5 Gops/s.
40 iters -> 0.96^40 ~ 0.20 (20% error on eigenvector -- rough but sufficient to see ratio >> 1).
"""
import sys, numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 1.05
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
k, niters = 18, 40

N = 3**(k-1); Nl = N//3
print(f"k={k}, N={N:,}, niters={niters}", flush=True)
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
    if (it+1) % 10 == 0:
        print(f"  iter {it+1}/{niters}", flush=True)

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
print(f"k={k} iters={niters}: E_s0/mu0={E0:.7f} E_s2/mu2={E2:.7f} ratio={ratio:.5f} R={R:.5f} F={F:.4f} sqrt(F)/R={sqFR:.5f} {'OK' if ratio>1 else 'FAIL'}")
