"""Convergence audit of inc_last at lambda=2: recompute k=14,15,16,17 with
doubled/tripled iterations (fresh, ignore cache) and compare."""
import numpy as np
from math import log2
ALPHA = log2(3.0)
lam = 2.00
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)

def inc_last_fresh(k, iters):
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = ((4*i+2) % N).astype(np.int32)
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = ((4*s_arr) % Nl).astype(np.int32)
    R3 = ((2*s_arr+1) % Nl).astype(np.int32)
    del i, s_arr, r_arr
    v = np.ones(N)
    out = {}
    done = 0
    for target in iters:
        for _ in range(target - done):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
            v = w/w.max()
        done = target
        F = np.log2(v); F -= F.mean()
        prev = None; last = None
        for p in range(k-1):
            M = 3**(p+1)
            cm = F.reshape(N//M, M).mean(axis=0)
            m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
            last = float((m**2).mean()); prev = cm
        out[target] = last
    return out

res = {}
for k, its in [(14, [450, 900]), (15, [300, 700]), (16, [220, 600]), (17, [180, 450])]:
    res[k] = inc_last_fresh(k, its)
    vals = res[k]
    keys = sorted(vals)
    rel = abs(vals[keys[1]]-vals[keys[0]])/vals[keys[1]]
    print(f"k={k}: inc_last @ {keys[0]} it = {vals[keys[0]]:.6e} | @ {keys[1]} it = {vals[keys[1]]:.6e} | rel diff {rel:.2e}", flush=True)

print("\nTR-reeks met geconvergeerde waarden:")
kk = sorted(res)
for a, b in zip(kk, kk[1:]):
    va = res[a][max(res[a])]; vb = res[b][max(res[b])]
    print(f"  TR({a}->{b}) = {vb/va:.5f}")
