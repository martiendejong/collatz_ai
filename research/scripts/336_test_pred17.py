"""Test prediction #17: TR(2.00, 17->18) in [0.828, 0.838].
k=18 at lambda=2: N = 3^17 = 129M, float64 (~1GB per array), ~350 iterations."""
import numpy as np
from math import log2
import os
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)
lam = 2.00
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)

def inc_last_of(v, k):
    N = v.size
    F = np.log2(v); F -= F.mean()
    prev = None; last = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        m = cm if prev is None else cm - prev[np.arange(M) % (M//3)]
        last = float((m**2).mean()); prev = cm
    return last

k = 18
N = 3**(k-1); Nl = N//3
i = np.arange(N, dtype=np.int64)
T4 = ((4*i+2) % N).astype(np.int32)
s_arr = i // 3
r_arr = (i % 3).astype(np.int8)
m0 = r_arr == 0; m2 = r_arr == 2
R1 = ((4*s_arr) % Nl).astype(np.int32)
R3 = ((2*s_arr+1) % Nl).astype(np.int32)
del i, s_arr, r_arr
v = np.ones(N, dtype=np.float32)
for it in range(300):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = np.float32(A)*v[T4]; w[m2] += np.float32(B3)*cb[R3[m2]]; w[m0] += np.float32(B1)*cb[R1[m0]]
    v = w/w.max()
    if (it+1) % 50 == 0:
        print(f"iter {it+1}/300", flush=True)
i18_a = inc_last_of(v, k)
# convergence check: 100 more iterations
for it in range(100):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = np.float32(A)*v[T4]; w[m2] += np.float32(B3)*cb[R3[m2]]; w[m0] += np.float32(B1)*cb[R1[m0]]
    v = w/w.max()
i18_b = inc_last_of(v, k)
print(f"inc_last(18) @300 = {i18_a:.6e}  @400 = {i18_b:.6e}  rel diff {abs(i18_b-i18_a)/i18_b:.2e}")
i17 = 3.843526e-03
tr = i18_b/i17
hit = 0.828 <= tr <= 0.838
print(f"TR(17->18) = {tr:.5f} | voorspeld [0.828, 0.838] -> {'RAAK' if hit else 'MIS'}")
