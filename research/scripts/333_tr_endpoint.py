"""Decisive creep test in the clean instrument: TR(lam=2, k) for k=9..16.
TR(k) = inc_last(k+1)/inc_last(k). Paper's V-instrument saw +0.003/step creep;
TR at k=12..14 saw none. Extend the series."""
import numpy as np
from math import log2
import os
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)
lam = 2.00
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)

def get_v(k):
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        return np.load(fn)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    it = {9:2500,10:2000,11:1500,12:1000,13:700,14:450,15:300,16:220}[k]
    v = np.ones(N)
    for _ in range(it):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    np.save(fn, v)
    return v

def inc_last(k):
    v = get_v(k)
    N = v.size
    F = np.log2(v); F -= F.mean()
    prev = None; last = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        last = float((m**2).mean()); prev = cm
    return last

vals = {}
for k in range(9, 17):
    vals[k] = inc_last(k)
    print(f"k={k}: inc_last = {vals[k]:.6e}", flush=True)
print("\nTR-reeks (k -> k+1):")
trs = []
for k in range(9, 16):
    tr = vals[k+1]/vals[k]
    trs.append(tr)
    print(f"  TR({k}->{k+1}) = {tr:.5f}")
d = np.diff(trs)
print("drifts:", " ".join(f"{x:+.5f}" for x in d))
print(f"mean drift/step: {d.mean():+.5f}  (paper V-instrument: +0.003)")
