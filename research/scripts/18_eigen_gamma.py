"""E28/G1b: class-coupled eigenproblem for the density exponent gamma(K).
Backward-tree node-counting N_s(y) ~ u_s y^gamma over classes s = (parity b, t mod 3^K):
  (A_g u)_s = 2^(-g) * u_(0, 2t)                       [doubling child, always]
            + [b=0 and t=1 mod 3] * 3^g * min_j u_(1, (t-1)/3 + j*3^(K-1))   [third child]
gamma(K) = the g with Perron radius rho(A_g) = 1 (A_g monotone homogeneous; power iteration).
"""
import sys, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
out = {}

def gamma_of_K(K, iters=3000, tol=1e-10):
    M = 3 ** K
    # precompute transitions
    dbl = (2 * np.arange(M)) % M               # child1 ternary
    branch_mask = (np.arange(M) % 3 == 1)
    base = (np.arange(M) - 1) // 3             # (t-1)/3 for t = 1 mod 3 (valid where mask)
    offs = 3 ** (K - 1)

    def rho(g):
        w1 = 2.0 ** (-g)
        w2 = 3.0 ** g
        u0 = np.ones(M); u1 = np.ones(M)       # u0 = even classes, u1 = odd classes
        r_prev = 1.0
        for it in range(iters):
            # child1 of (b,t) is (0, 2t) for both parities
            c1 = u0[dbl]
            # child2 only from even, t=1 mod 3 -> (1, base + j*offs), min over j
            m3 = np.minimum(np.minimum(u1[(base) % M], u1[(base + offs) % M]),
                            u1[(base + 2 * offs) % M])
            n0 = w1 * c1 + np.where(branch_mask, w2 * m3, 0.0)
            n1 = w1 * c1
            r = max(n0.max(), n1.max())
            u0, u1 = n0 / r, n1 / r
            if it > 50 and abs(r - r_prev) < tol:
                return r
            r_prev = r
        return r_prev

    lo, hi = 0.3, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if rho(mid) > 1: hi = mid
        else: lo = mid
    return (lo + hi) / 2

for K in range(3, 10):
    g = gamma_of_K(K)
    out[K] = round(g, 5)
    print(f"K={K} (mod 3^{K}, {2*3**K} classes): gamma = {g:.5f}", flush=True)
print(json.dumps(out, indent=1))
