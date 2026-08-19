# 397: rho(2, k=19) — endpoint gammabar series one deeper (fork evidence).
# N = 3^18 = 387,420,489. Careful dtypes: float32 vector, int32 indices.
import numpy as np
from math import log2

ALPHA = log2(3.0)
k = 19
N = 3**(k-1)
A = np.float32(0.25); B1 = np.float32(0.75); B3 = np.float32(1.5)
print(f"k={k}: N = {N:,}", flush=True)
i = np.arange(N, dtype=np.int64)
T4 = ((4*i + 2) % N).astype(np.int32)
r8 = (i % 3).astype(np.int8)
s = i // 3
Nl = N//3
R1 = ((4*s) % Nl).astype(np.int32)
R3 = ((2*s + 1) % Nl).astype(np.int32)
del i, s
m0 = r8 == 0; m2 = r8 == 2
del r8
v = np.ones(N, dtype=np.float32)
rho_prev = 0.0
for it in range(3000):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A*v[T4]
    w[m2] += B3*cb[R3[m2]]
    w[m0] += B1*cb[R1[m0]]
    rho = float(w.max())
    v = w/np.float32(rho)
    del w, cb
    if it % 100 == 0:
        print(f"iter {it}: rho = {rho:.8f}", flush=True)
    rho_prev = rho
gb = 4*(1 - rho)/3
print(f"KLAAR k=19: rho(2,19) = {rho:.8f}  gammabar = {gb:.8f}")
print(f"vorige (k=18): gammabar = 0.04634253 -> ratio = {gb/0.04634253:.4f}, c(2)-schatting = {(gb/0.04634253)**2:.4f}")
np.save(r"E:\projects\collatz\research\cache\v_lam2.00_k19.npy", v)
