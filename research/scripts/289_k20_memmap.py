"""
289_k20_memmap.py - k=20, lambda=1.05 via disk-backed memmaps.
N = 3^19 = 1,162,261,467. Two float32 memmaps (4.65 GB each) on E:.
RAM use: cb (1.55 GB) + chunk buffers (~1 GB).

Tricks:
- No separate normalization pass: F is 1-homogeneous, so dividing the
  constants A, B1, B3 by the previous iteration's max is equivalent to
  normalizing the vector. Values stay O(1); measured wmax -> rho directly.
- T4 gather via 4 affine segments (j = 4i+2 - s*N is a stride-4 slice).
- cb-additions folded into the same chunk pass (w written exactly once/iter).
- Ping-pong: swap memmap handles instead of copying.

Model prediction (Obs 486): ratio(20) ~ 1.0776 (the trough minimum).
"""
import numpy as np
import os
from math import log2, sqrt

ALPHA = log2(3.0)
lam = 1.05
k = 20
NITERS = 35
CHUNK = 30_000_000

N = 3**(k-1)           # 1,162,261,467
Nl = N // 3
Nl3 = Nl // 3

A  = np.float32(lam**-2)
B1 = np.float32(lam**(ALPHA-2))
B3 = np.float32(lam**(ALPHA-1))

print(f"k={k}, N={N:,}, niters={NITERS}, CHUNK={CHUNK:,}, float32 memmap", flush=True)

d = "E:/temp/collatz_k20"
os.makedirs(d, exist_ok=True)
v = np.memmap(d + "/v.dat", dtype=np.float32, mode="w+", shape=(N,))
w = np.memmap(d + "/w.dat", dtype=np.float32, mode="w+", shape=(N,))
for a in range(0, N, CHUNK):
    v[a:min(a + CHUNK, N)] = 1.0
print("init done", flush=True)

cb = np.empty(Nl, dtype=np.float32)

# segment s holds i with s*N <= 4i+2 < (s+1)*N
seg = [max(0, -(-(s * N - 2) // 4)) for s in range(5)]
seg[4] = N

m_prev = 1.0
rho_est = 0.0
for it in range(NITERS):
    A_eff = np.float32(float(A) / m_prev)
    B1_eff = np.float32(float(B1) / m_prev)
    B3_eff = np.float32(float(B3) / m_prev)
    for a in range(0, Nl, CHUNK):
        b = min(a + CHUNK, Nl)
        cb[a:b] = np.minimum(np.minimum(v[a:b], v[Nl + a:Nl + b]),
                             v[2 * Nl + a:2 * Nl + b])
    wmax = 0.0
    for s in range(4):
        i0, i1 = seg[s], seg[s + 1]
        for a in range(i0, i1, CHUNK):
            b = min(a + CHUNK, i1)
            j0 = 4 * a + 2 - s * N
            buf = A_eff * np.asarray(v[j0: j0 + 4 * (b - a): 4])
            first2 = a + ((2 - a) % 3)
            if first2 < b:
                sv = np.arange((first2 - 2) // 3, (b - 1 - 2) // 3 + 1, dtype=np.int64)
                buf[first2 - a::3] += B3_eff * cb[(2 * sv + 1) % Nl]
                del sv
            first0 = a + ((-a) % 3)
            if first0 < b:
                sv = np.arange(first0 // 3, (b - 1) // 3 + 1, dtype=np.int64)
                buf[first0 - a::3] += B1_eff * cb[(4 * sv) % Nl]
                del sv
            m = float(buf.max())
            if m > wmax:
                wmax = m
            w[a:b] = buf
            del buf
    rho_est = wmax  # with constants scaled by 1/m_prev, wmax -> rho at convergence
    if (it + 1) % 5 == 0 or it == 0:
        print(f"  iter {it+1}/{NITERS} wmax={wmax:.6e} rho_est={rho_est:.6f}", flush=True)
    m_prev = wmax
    v, w = w, v

# analysis: v holds the final vector
sum0 = ssum0 = sum2 = ssum2 = 0.0
CH3 = CHUNK // 3
for a in range(0, Nl3, CH3):
    b = min(a + CH3, Nl3)
    x1 = np.asarray(v[3 * a: 3 * b: 3], dtype=np.float64)
    x2 = np.asarray(v[3 * a + Nl: 3 * b + Nl: 3], dtype=np.float64)
    x3 = np.asarray(v[3 * a + 2 * Nl: 3 * b + 2 * Nl: 3], dtype=np.float64)
    mu = (x1 + x2 + x3) / 3.0
    var = ((x1 - mu) ** 2 + (x2 - mu) ** 2 + (x3 - mu) ** 2) / 3.0
    ssum0 += float(np.sqrt(var).sum())
    sum0 += float(mu.sum())
    y1 = np.asarray(v[3 * a + 2: 3 * b + 2: 3], dtype=np.float64)
    y2 = np.asarray(v[3 * a + 2 + Nl: 3 * b + 2 + Nl: 3], dtype=np.float64)
    y3 = np.asarray(v[3 * a + 2 + 2 * Nl: 3 * b + 2 + 2 * Nl: 3], dtype=np.float64)
    mu2 = (y1 + y2 + y3) / 3.0
    var2 = ((y1 - mu2) ** 2 + (y2 - mu2) ** 2 + (y3 - mu2) ** 2) / 3.0
    ssum2 += float(np.sqrt(var2).sum())
    sum2 += float(mu2.sum())
    del x1, x2, x3, y1, y2, y3, mu, var, mu2, var2

E_s0 = ssum0 / sum0     # = (E[sigma_w]/mu) for r=0
E_s2 = ssum2 / sum2
ratio = E_s2 / E_s0

rho = rho_est
t = float(A) / rho
R_val = (t ** 2 + lam) / (1 + t * lam)
F_val = (t ** 4 + lam ** 2) / (1 + t ** 2 * lam ** 2)
sqFR = sqrt(F_val) / R_val

print(f"\nk={k} niters={NITERS}: E_s0/mu0={E_s0:.7f} E_s2/mu2={E_s2:.7f}")
print(f"  ratio={ratio:.5f} rho={rho:.6f} R={R_val:.5f} F={F_val:.5f} sqrt(F)/R={sqFR:.5f} "
      f"{'OK' if ratio > 1 else 'FAIL!!'}")
print(f"  correction from limit: {sqFR - ratio:.5f}")
print("DONE", flush=True)
