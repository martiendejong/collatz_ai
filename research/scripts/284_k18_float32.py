"""
284_k18_float32.py
==================
k=18, lambda=1.05 with float32 + int32 to halve memory.
N=129M: float32 needs 493 MiB per array (vs 985 MiB float64).
Total memory ~3.5 GB (feasible on 8GB machine).
float32 gives 7 significant digits -- sufficient for ratio >> 1.
"""
import sys, numpy as np
from math import log2, sqrt

ALPHA = log2(3.0)
lam = 1.05
A = np.float32(lam**-2); B1 = np.float32(lam**(ALPHA-2)); B3 = np.float32(lam**(ALPHA-1))
k, niters = 18, 80

N = 3**(k-1); Nl = N//3; Nl3 = Nl//3
print(f"k={k}, N={N:,}, niters={niters}, using float32/int32", flush=True)

i = np.arange(N, dtype=np.int32)
T4 = ((4*i.astype(np.int64)+2) % N).astype(np.int32)
s_arr = (i // 3).astype(np.int32)
r_arr = (i % 3).astype(np.uint8)
m0 = r_arr == 0; m2 = r_arr == 2
del r_arr
R1 = ((4*s_arr.astype(np.int64)) % Nl).astype(np.int32)
R3 = ((2*s_arr.astype(np.int64)+1) % Nl).astype(np.int32)
del s_arr, i

v = np.ones(N, dtype=np.float32)
for it in range(niters):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A * v[T4]
    w[m2] += B3 * cb[R3[m2]]
    w[m0] += B1 * cb[R1[m0]]
    v = w / w.max()
    if (it+1) % 20 == 0:
        print(f"  iter {it+1}/{niters}", flush=True)

rho = float(w.max()); t = A/rho; R = (t**2+lam)/(1+t*lam); F = (t**4+lam**2)/(1+t**2*lam**2)
j3 = np.arange(Nl3, dtype=np.int32)
v0 = v[0::3]; v2 = v[2::3]
c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1).astype(np.float64)
c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1).astype(np.float64)
E0 = c0.std(1).mean() / float(v0.mean())
E2 = c2.std(1).mean() / float(v2.mean())
ratio = E2/E0
print(f"k={k} niters={niters}: ratio={ratio:.5f} R={R:.5f} F={F:.4f} sqF_R={sqrt(F)/R:.5f} {'OK' if ratio>1 else 'FAIL'}")
print("DONE")
