# 390: background — rho(2, k=18) cold in RAM (float32 + int32 indices):
# extends the gammabar(2,k) series (third c(2) instrument) one depth beyond
# Obs 549 (k<=17). Mass identity: 1 - rho = (3/4) gammabar.
import numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 2.0
k = 18
N = 3**(k-1)
A = np.float32(lam**-2); B1 = np.float32(lam**(ALPHA-2)); B3 = np.float32(lam**(ALPHA-1))
i = np.arange(N, dtype=np.int64)
T4 = ((4*i + 2) % N).astype(np.int32)
r = (i % 3).astype(np.int8)
s = (i // 3).astype(np.int64)
Nl = N//3
R1 = ((4*s) % Nl).astype(np.int32)
R3 = ((2*s + 1) % Nl).astype(np.int32)
del i, s
m0 = r == 0; m2 = r == 2
del r
v = np.ones(N, dtype=np.float32)
rho_prev = 0.0
stable = 0
for it in range(3000):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A*v[T4]
    w[m2] += B3*cb[R3[m2]]
    w[m0] += B1*cb[R1[m0]]
    rho = float(w.max())
    v = w/np.float32(rho)
    if it % 100 == 0:
        print(f"iter {it}: rho = {rho:.8f}", flush=True)
    if abs(rho - rho_prev) < 5e-8:
        stable += 1
        if stable >= 50 and it > 600:
            break
    else:
        stable = 0
    rho_prev = rho
gb = 4*(1 - rho)/3
print(f"KLAAR k=18: rho(2,18) = {rho:.8f}  gammabar = {gb:.8f}")
print(f"vorige (k=17): gammabar = 0.05055271 -> ratio = {gb/0.05055271:.4f}, c(2)-schatting = {(gb/0.05055271)**2:.4f}")
np.save(r"E:\projects\collatz\research\cache\v_lam2.00_k18.npy", v)
