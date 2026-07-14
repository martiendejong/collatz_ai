"""E32: the EXACT Krasikov-Lagarias 2003 system L_k^NT(lambda) (arXiv math/0205002).
Classes [3^k] = {m mod 3^k, m = 2 mod 3}. alpha = log2(3). Constraints:
 (L1) m=2 mod 9:  c^m <= c^{4m} L^-2 + cbar^{(4m-2)/3} L^(a-2)
 (L2) m=5 mod 9:  c^m <= c^{4m} L^-2
 (L3) m=8 mod 9:  c^m <= c^{4m} L^-2 + cbar^{(2m-1)/3} L^(a-1)
 (L4) cbar^t = min over the three refinements of t in [3^k]
Feasible with c>0 <=> rho(F_L) >= 1; gamma = log2(lambda). K-L: k=11 -> 0.84.
"""
import sys, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
out = {}
ALPHA = math.log2(3)

def solve_k(k, iters=400):
    N = 3 ** (k - 1)                     # size of [3^k]
    M = 3 ** k
    Nc = 3 ** (k - 2)                    # size of [3^(k-1)]
    Mc = 3 ** (k - 1)
    m = 3 * np.arange(N, dtype=np.int64) + 2          # class reps
    i4 = (((4 * m) % M) - 2) // 3                     # index of 4m
    mod9 = m % 9
    is2 = mod9 == 2
    is8 = mod9 == 8
    # branch targets in [3^(k-1)]
    t2 = np.zeros(N, dtype=np.int64); t8 = np.zeros(N, dtype=np.int64)
    t2[is2] = (((4 * m[is2] - 2) // 3) % Mc)
    t8[is8] = (((2 * m[is8] - 1) // 3) % Mc)
    # refinement indices in [3^k] for a class t in [3^(k-1)]: t, t+3^(k-1), t+2*3^(k-1)
    def ref_idx(t):
        return ((t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3)
    r2 = ref_idx(t2); r8 = ref_idx(t8)

    def min_ratio(lam):
        w0 = lam ** -2.0
        w2 = lam ** (ALPHA - 2.0)
        w8 = lam ** (ALPHA - 1.0)
        c = np.ones(N)
        for it in range(iters):
            cb2 = np.minimum(np.minimum(c[r2[0]], c[r2[1]]), c[r2[2]])
            cb8 = np.minimum(np.minimum(c[r8[0]], c[r8[1]]), c[r8[2]])
            f = w0 * c[i4]
            f = f + np.where(is2, w2 * cb2, 0.0) + np.where(is8, w8 * cb8, 0.0)
            r = f.max()
            c = f / r
        # feasibility indicator: growth factor of the iteration
        cb2 = np.minimum(np.minimum(c[r2[0]], c[r2[1]]), c[r2[2]])
        cb8 = np.minimum(np.minimum(c[r8[0]], c[r8[1]]), c[r8[2]])
        f = w0 * c[i4] + np.where(is2, w2 * cb2, 0.0) + np.where(is8, w8 * cb8, 0.0)
        return (f / c).min()

    lo, hi = 1.5, 2.0
    for _ in range(30):
        mid = (lo + hi) / 2
        if min_ratio(mid) >= 1.0: lo = mid
        else: hi = mid
    return lo, math.log2(lo)

for k in (5, 7, 9, 11):
    lam, gamma = solve_k(k)
    out[k] = dict(lam=round(lam, 5), gamma=round(gamma, 5))
    print(f"k={k}: lambda = {lam:.5f}, gamma = log2(lambda) = {gamma:.5f}", flush=True)
print(json.dumps(out, indent=1))
