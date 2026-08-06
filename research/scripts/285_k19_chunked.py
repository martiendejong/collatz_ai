"""
285_k19_chunked.py
==================
k=19, lambda=1.05 with chunked T4 computation to avoid the 1.5 GB T4 array.
N = 3^18 = 387M. Memory target: ~5 GB (feasible on 8 GB machine).

Memory layout:
  v, w: N float32 each = 1.548 GB
  cb: N/3 float32 = 0.516 GB
  R1_compact, R3_compact: N/3 int32 each = 0.516 GB each
  T4 computed in CHUNK=10M int32 chunks: 40 MB peak overhead
  Total: ~4.7 GB + OS overhead ~5.2 GB

T4(i) = (4i+2) % N -- computed on-the-fly to save 1.55 GB vs precomputing.

40 iterations -- rough but sufficient (k=17,18 show ratio converges fast in practice).
"""
import numpy as np
from math import log2, sqrt
import sys

ALPHA = log2(3.0)
lam = 1.05
k = 19
niters = 40
CHUNK = 10_000_000  # chunk size for T4 computation

N = 3**(k-1)  # 387,420,489
Nl = N // 3   # 129,140,163
Nl3 = Nl // 3 # 43,046,721

A  = np.float32(lam**-2)
B1 = np.float32(lam**(ALPHA-2))
B3 = np.float32(lam**(ALPHA-1))

print(f"k={k}, N={N:,}, niters={niters}, CHUNK={CHUNK:,}, float32", flush=True)
print(f"  Nl={Nl:,}, Nl3={Nl3:,}", flush=True)
print(f"  Estimated memory: ~5 GB", flush=True)

# Precompute R1 and R3 (compact, for r=0 and r=2 respectively)
# For i with r=0: i=3j, s(i)=j. R1[j] = (4j) % Nl
# For i with r=2: i=3j+2, s(i)=j. R3[j] = (2j+1) % Nl
print("Precomputing R1_compact, R3_compact...", flush=True)
j = np.arange(Nl, dtype=np.int64)
R1_compact = (4 * j % Nl).astype(np.int32)
R3_compact = ((2 * j + 1) % Nl).astype(np.int32)
del j
print("  Done.", flush=True)

# Initialize eigenvector
v = np.ones(N, dtype=np.float32)
w = np.empty(N, dtype=np.float32)

for it in range(niters):
    # Step 1: compute cb = column-wise min
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    # Step 2: compute w = A * v[T4] in chunks (T4 computed on-the-fly)
    for start in range(0, N, CHUNK):
        end = min(start + CHUNK, N)
        i_chunk = np.arange(start, end, dtype=np.int64)
        T4_chunk = ((4 * i_chunk + 2) % N).astype(np.int32)
        w[start:end] = A * v[T4_chunk]
        del i_chunk, T4_chunk

    # Step 3: add B3 * cb[R3] to r=2 positions (i = 2, 5, 8, ...)
    w[2::3] += B3 * cb[R3_compact]
    # Step 4: add B1 * cb[R1] to r=0 positions (i = 0, 3, 6, ...)
    w[0::3] += B1 * cb[R1_compact]

    # Normalize
    wmax = w.max()
    w /= wmax
    v[:] = w

    if (it + 1) % 10 == 0:
        print(f"  iter {it+1}/{niters}, max_before_norm={wmax:.6f}", flush=True)

# Compute Perron eigenvalue and derived quantities
rho = float(wmax)
t = float(A) / rho
R = (t**2 + lam) / (1 + t * lam)
F = (t**4 + lam**2) / (1 + t**2 * lam**2)

# Extract residue-class subvectors
v0 = v[0::3]  # r=0, length Nl
v2 = v[2::3]  # r=2, length Nl

mu0 = float(v0.mean())
mu2 = float(v2.mean())

# Column-triplet analysis
j3 = np.arange(Nl3, dtype=np.int32)
c0 = np.stack([v0[j3], v0[j3 + Nl3], v0[j3 + 2*Nl3]], axis=1).astype(np.float64)
c2 = np.stack([v2[j3], v2[j3 + Nl3], v2[j3 + 2*Nl3]], axis=1).astype(np.float64)

E0 = c0.std(axis=1).mean() / mu0
E2 = c2.std(axis=1).mean() / mu2
ratio = E2 / E0

sqFR = sqrt(F) / R
result = "OK" if ratio > 1 else "FAIL"
print(f"\nk={k} niters={niters}: E_s0/mu0={E0:.7f} E_s2/mu2={E2:.7f}", flush=True)
print(f"  ratio={ratio:.5f} R={R:.5f} F={F:.5f} sqrt(F)/R={sqFR:.5f} {result}", flush=True)
print(f"  correction from limit: {sqFR - ratio:.5f}", flush=True)
print("DONE", flush=True)
