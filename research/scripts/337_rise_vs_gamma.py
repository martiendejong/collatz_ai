"""
Hypothesis (frozen form BEFORE measuring): the TR-rise per k-step is proportional
to the gap mass: TR(k->k+1) - TR(k-1->k) ... better: define rise(k) =
TR(k,k+1) - TR(k-1,k). Test: rise(k, lam) = beta * gammabar(lam, k) with
lambda-independent beta. Alternative: beta varies with lam (refuted form).
Measure TR at k=13..17 for lam in {1.30, 1.70, 1.90, 2.00}.
"""
import numpy as np
from math import log2
import os
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

def get_v(lam, k):
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        return np.load(fn)
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = ((4*i+2) % N).astype(np.int32)
    s_arr = i//3; r_arr = (i%3).astype(np.int8)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = ((4*s_arr) % Nl).astype(np.int32)
    R3 = ((2*s_arr+1) % Nl).astype(np.int32)
    del i, s_arr, r_arr
    dt = np.float32 if k >= 16 else np.float64
    v = np.ones(N, dtype=dt)
    its = {13:600, 14:450, 15:320, 16:240, 17:200}[k]
    A_, B1_, B3_ = dt(A), dt(B1), dt(B3)
    for _ in range(its):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A_*v[T4]; w[m2] += B3_*cb[R3[m2]]; w[m0] += B1_*cb[R1[m0]]
        v = w/w.max()
    np.save(fn, v)
    return v

def inc_last(lam, k):
    v = get_v(lam, k)
    N = v.size
    F = np.log2(v.astype(np.float64)); F -= F.mean()
    prev = None; last = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        last = float((m**2).mean()); prev = cm
    return last

def gammabar(lam, k):
    v = get_v(lam, k).astype(np.float64)
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    return float((vbar-cb).mean()/vbar.mean())

for lam in [1.30, 1.70, 1.90, 2.00]:
    incs = {k: inc_last(lam, k) for k in range(13, 18)}
    TR = {k: incs[k+1]/incs[k] for k in range(13, 17)}
    print(f"\nlam={lam}: TR(13-14..16-17) = " + " ".join(f"{TR[k]:.5f}" for k in range(13,17)), flush=True)
    rises = {k: TR[k]-TR[k-1] for k in range(14, 17)}
    gams = {k: gammabar(lam, k) for k in range(14, 17)}
    for k in range(14, 17):
        beta = rises[k]/gams[k] if gams[k] > 0 else float('nan')
        print(f"  k={k}: rise={rises[k]:+.5f}  gammabar={gams[k]:.5f}  beta={beta:+.3f}")
