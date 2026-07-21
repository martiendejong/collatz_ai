"""Compute mod-2048 spectral gap properly using scipy."""
import numpy as np, sys

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n + 1); m = (n + 1) >> K; x = m * (3**K) - 1; l = v2(x)
    return x >> l, K, l

MOD = 2048
N_SAMP = 256
odd_res = list(range(1, MOD, 2))
N = len(odd_res)
idx = {r: i for i, r in enumerate(odd_res)}

print(f"Building mod-{MOD} chain ({N} states, {N_SAMP} samples)...")
sys.stdout.flush()

P = np.zeros((N, N))
for i, r in enumerate(odd_res):
    if i % 128 == 0: print(f"  {i}/{N}", end=" ", flush=True)
    K0 = v2(r + 1)
    counts = np.zeros(N)
    valid = 0
    for k_iter in range(N_SAMP):
        n = r + MOD * k_iter
        if v2(n + 1) != K0: continue
        n_out, _, _ = macro_step(n)
        r_out = n_out % MOD
        if r_out % 2 == 1 and r_out in idx:
            counts[idx[r_out]] += 1
            valid += 1
    if valid > 0:
        P[i] = counts / valid
print()

pi = np.ones(N) / N
for _ in range(500):
    pi_new = pi @ P
    if np.max(np.abs(pi_new - pi)) < 1e-10: break
    pi = pi_new

uniform = 1.0 / N
max_dev = np.max(np.abs(pi - uniform))
avg_k0 = sum(pi[i] * v2(r + 1) for i, r in enumerate(odd_res))
print(f"Max stationary deviation: {max_dev/uniform*100:.3f}%")
print(f"E[k0] = {avg_k0:.6f}")

try:
    from scipy.sparse.linalg import eigs
    import scipy.sparse as sp
    P_sp = sp.csr_matrix(P.T)
    vals, _ = eigs(P_sp, k=6, which="LM")
    vals_sorted = sorted(vals.real, reverse=True)
    print(f"Top eigenvalues (scipy): {[f'{v:.6f}' for v in vals_sorted[:6]]}")
    print(f"Spectral gap = {1 - vals_sorted[1]:.6f}")
    print(f"lambda2 = {vals_sorted[1]:.6f}")
except ImportError:
    print("scipy not available, using numpy full eigendecomp...")
    vals = np.linalg.eigvals(P)
    vals_sorted = sorted(vals.real, reverse=True)
    print(f"Top 5 (numpy): {[f'{v:.6f}' for v in vals_sorted[:5]]}")
    print(f"Spectral gap = {1 - vals_sorted[1]:.6f}")
    print(f"lambda2 = {vals_sorted[1]:.6f}")

# Also collect j-distribution at mod-2048
from collections import Counter
j_dist = Counter()
for r in odd_res:
    K = v2(r + 1)
    m_red = (r + 1) >> K
    prod = m_red * (3**K)
    l0 = v2(prod - 1)
    j = 11 - K - l0
    j_dist[min(j, 12)] += 1

print("\nj-distribution at mod-2048:")
for j in sorted(j_dist):
    print(f"  j={j}: {j_dist[j]} elements")
