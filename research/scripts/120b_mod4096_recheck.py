"""Re-check mod-4096 spectral gap with 512 samples (vs script 117's 128)."""
import numpy as np, sys
from scipy.sparse.linalg import eigs
import scipy.sparse as sp

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

MOD = 4096; N_SAMP = 512
odd_res = list(range(1, MOD, 2)); N = len(odd_res)
idx = {r: i for i, r in enumerate(odd_res)}
print(f"Building mod-{MOD} chain ({N} states, {N_SAMP} samples)...")
P = np.zeros((N, N), dtype=np.float32)
for i, r in enumerate(odd_res):
    if i % 256 == 0: print(f"  {i}/{N}", end=" ", flush=True)
    K0 = v2(r+1); counts = np.zeros(N, dtype=np.float32); valid = 0
    for k in range(N_SAMP):
        n = r + MOD*k
        if v2(n+1) != K0: continue
        n_out, _, _ = macro_step(n)
        r_out = n_out % MOD
        if r_out % 2 == 1 and r_out in idx:
            counts[idx[r_out]] += 1; valid += 1
    if valid > 0: P[i] = counts/valid
print()
pi = np.ones(N)/N; P64 = P.astype(np.float64)
for _ in range(300):
    pi_new = pi @ P64
    if np.max(np.abs(pi_new-pi)) < 1e-8: break
    pi = pi_new
print(f"Max stationary deviation: {np.max(np.abs(pi-1/N))*N*100:.3f}%")
P_sp = sp.csr_matrix(P64.T)
vals, _ = eigs(P_sp, k=10, which="LM")
vals_sorted = sorted(abs(v) for v in vals)[::-1]
print(f"Top |lambda|: {[f'{v:.6f}' for v in vals_sorted[:10]]}")
print(f"Spectral gap = {1-vals_sorted[1]:.6f}, |lambda2| = {vals_sorted[1]:.6f}")
print(f"Script 117 reported: gap=0.809, |lambda2|=0.191 (N_SAMP=128)")
