"""
117_mod4096_gap.py
===================
Fifth data point for the Expander Conjecture: mod-4096 spectral gap.

Existing series (lambda2 = second eigenvalue):
  mod-256  (128 states):  gap=0.938, lambda2=0.062
  mod-512  (256 states):  gap=0.913, lambda2=0.087
  mod-1024 (512 states):  gap=0.886, lambda2=0.114
  mod-2048 (1024 states): gap=0.840, lambda2=0.160

If lambda2 ~ N^0.44: at N=2048 => lambda2 ~ 0.160 * 2^0.44 = 0.218, gap=0.782
If exponent decreasing toward 0: lambda2 approaching constant, gap approaching constant > 0.
"""
import numpy as np
import sys

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n + 1); m = (n + 1) >> K; x = m * (3**K) - 1; l = v2(x)
    return x >> l, K, l

MOD = 4096
N_SAMP = 128
odd_res = list(range(1, MOD, 2))
N = len(odd_res)
idx = {r: i for i, r in enumerate(odd_res)}

print(f"Building mod-{MOD} chain ({N} states, {N_SAMP} samples/state)...")
sys.stdout.flush()

import numpy as np
P = np.zeros((N, N), dtype=np.float32)
for i, r in enumerate(odd_res):
    if i % 256 == 0: print(f"  {i}/{N}", end=" ", flush=True)
    K0 = v2(r + 1)
    counts = np.zeros(N, dtype=np.float32)
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

# Stationary distribution
pi = np.ones(N, dtype=np.float64) / N
P64 = P.astype(np.float64)
for _ in range(300):
    pi_new = pi @ P64
    if np.max(np.abs(pi_new - pi)) < 1e-9: break
    pi = pi_new

uniform = 1.0 / N
max_dev = float(np.max(np.abs(pi - uniform)))
avg_k0 = sum(pi[i] * v2(r + 1) for i, r in enumerate(odd_res))
print(f"Max stationary deviation: {max_dev/uniform*100:.3f}%")
print(f"E[k0] = {avg_k0:.6f}")

# Spectral gap via scipy sparse eigensolver
try:
    from scipy.sparse.linalg import eigs
    import scipy.sparse as sp
    P_sp = sp.csr_matrix(P64.T)
    vals, _ = eigs(P_sp, k=8, which="LM")
    vals_sorted = sorted(abs(v) for v in vals)[::-1]
    print(f"Top eigenvalues (scipy): {[f'{v:.6f}' for v in vals_sorted[:8]]}")
    gap = 1 - vals_sorted[1]
    lambda2 = vals_sorted[1]
    print(f"Spectral gap = {gap:.6f}")
    print(f"lambda2 = {lambda2:.6f}")
except ImportError:
    print("scipy not available, using numpy full eigendecomp (slow)...")
    vals = np.linalg.eigvals(P64)
    vals_abs = sorted(abs(v) for v in vals)[::-1]
    print(f"Top 8: {[f'{v:.6f}' for v in vals_abs[:8]]}")
    print(f"Spectral gap = {1 - vals_abs[1]:.6f}")
    print(f"lambda2 = {vals_abs[1]:.6f}")

# Summary table with prediction
print()
print("=" * 60)
print("EXPANDER CONJECTURE STATUS")
print("=" * 60)
data = [
    (128,  0.938189, 0.061811),
    (256,  0.912523, 0.087477),
    (512,  0.885971, 0.114029),
    (1024, 0.839642, 0.160358),
]
data.append((N, gap, lambda2))
import math
print(f"{'States':>8} {'Gap':>10} {'Lambda2':>10} {'alpha':>8}")
print("-" * 40)
prev_l = None
for i, (n, g, l) in enumerate(data):
    if prev_l is not None and prev_l > 0:
        alpha = math.log(l / prev_l) / math.log(n / data[i-1][0])
    else:
        alpha = float('nan')
    print(f"{n:>8} {g:>10.6f} {l:>10.6f} {alpha:>8.4f}")
    prev_l = l

# Fit power law to last 4 points (excluding first since it might be outlier)
if len(data) >= 4:
    ns = np.log([d[0] for d in data[1:]])
    ls = np.log([d[2] for d in data[1:]])
    slope, intercept = np.polyfit(ns, ls, 1)
    print(f"\nFitted lambda2 ~ N^{slope:.4f} (from last {len(data)-1} points)")
    print(f"=> {'gap -> 0' if slope > 0 else 'gap -> constant'}")
