"""
121_gap_survey.py
==================
Standardized spectral gap survey: fixed total computation budget ~524k transitions.
N_SAMP = 524288 // N_states for each modulus.
Uses |lambda| (spectral radius) for gap computation.
"""
import numpy as np, math
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

BUDGET = 524288
results = []

for MOD in [256, 512, 1024, 2048, 4096, 8192, 16384]:
    odd_res = list(range(1, MOD, 2)); N = len(odd_res)
    NS = max(64, BUDGET // N)
    idx = {r: i for i, r in enumerate(odd_res)}
    P = np.zeros((N, N), dtype=np.float32)
    for i, r in enumerate(odd_res):
        K0 = v2(r+1); counts = np.zeros(N, dtype=np.float32); valid = 0
        for k in range(NS):
            n = r + MOD*k
            if v2(n+1) != K0: continue
            n_out, _, _ = macro_step(n)
            r_out = n_out % MOD
            if r_out % 2 == 1 and r_out in idx:
                counts[idx[r_out]] += 1; valid += 1
        if valid > 0: P[i] = counts / valid
    P64 = P.astype(np.float64)
    pi = np.ones(N)/N
    for _ in range(300):
        pi_new = pi @ P64
        if np.max(np.abs(pi_new - pi)) < 1e-9: break
        pi = pi_new
    dev = np.max(np.abs(pi - 1/N)) * N * 100
    P_sp = sp.csr_matrix(P64.T)
    k_eig = min(10, N-1)
    vals, _ = eigs(P_sp, k=k_eig, which="LM")
    vals_abs = sorted(abs(v) for v in vals)[::-1]
    l2 = vals_abs[1]
    gap = 1 - l2
    results.append((N, MOD, NS, gap, l2, dev))
    print(f"mod-{MOD:6d} ({N:5d} states, NS={NS:5d}): gap={gap:.6f} |l2|={l2:.6f} dev={dev:.2f}%")

print()
print("=" * 70)
print("EXPANDER CONJECTURE SURVEY (fixed-budget, |lambda| spectral radius)")
print("=" * 70)
print(f"{'States':>8} {'Mod':>8} {'NS':>6} {'Gap':>10} {'|l2|':>8} {'alpha':>8}")
print("-" * 55)
for i, (N, MOD, NS, g, l, dev) in enumerate(results):
    if i == 0:
        print(f"{N:>8} {MOD:>8} {NS:>6} {g:>10.6f} {l:>8.6f} {'—':>8}")
    else:
        alpha = math.log(l/results[i-1][4]) / math.log(N/results[i-1][0])
        print(f"{N:>8} {MOD:>8} {NS:>6} {g:>10.6f} {l:>8.6f} {alpha:>8.4f}")

# Power law fit (all points)
ns = np.log([r[0] for r in results])
ls = np.log([r[4] for r in results])
slope, intercept = np.polyfit(ns, ls, 1)
print(f"\nFit |l2| ~ N^alpha (all {len(results)} points): alpha = {slope:.4f}")
print(f"Extrapolated: gap -> 0 at N ~= {int(np.exp(-intercept/slope)):,} states")
if slope < 0.5:
    print("Note: alpha < 0.5 suggests lambda2 might saturate (gap might stay > 0)")
