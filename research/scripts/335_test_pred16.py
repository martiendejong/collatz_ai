"""Test PREDICTIONS #16: TR(2.00, 16->17) = 0.8215 +- 0.0045.
Compute k=17 vector at lambda=2 (N=43M, light), same recipe as Script 333."""
import numpy as np
from math import log2
import os
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)
lam = 2.00
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)

def get_v(k, iters):
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        return np.load(fn)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = ((4*i+2) % N).astype(np.int32)
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = ((4*s_arr) % Nl).astype(np.int32)
    R3 = ((2*s_arr+1) % Nl).astype(np.int32)
    del i, s_arr, r_arr
    v = np.ones(N)
    for _ in range(iters):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    np.save(fn, v)
    return v

def inc_last(k, iters=200):
    v = get_v(k, iters)
    N = v.size
    F = np.log2(v); F -= F.mean()
    prev = None; last = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        last = float((m**2).mean()); prev = cm
    return last

i16 = inc_last(16, 220)
i17 = inc_last(17, 180)
tr = i17/i16
hit = abs(tr - 0.8215) <= 0.0045
print(f"inc_last(16) = {i16:.6e}")
print(f"inc_last(17) = {i17:.6e}")
print(f"TR(16->17) = {tr:.5f}  | voorspeld 0.8215 +- 0.0045 -> {'RAAK' if hit else 'MIS'}")
