"""E30: K-L amortized class-coupled system.
u_c = sum_{leaves l of block tree from c} theta_l^gamma * sum_{c' in refine(l)} u_{c'}
gamma* = gamma with Perron rho = 1. Leaves harvested at info-exhaustion (o=K, theta<1)
or horizon (d=0 remaining, theta<1). Refinement sum v(partial) = sum of u over trit extensions.
One backward DP computes B_gamma u for all roots at once.
"""
import sys, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
out = {}
L23 = math.log2(3)

def solve(K, D, power_iters=60):
    # state packing: (parity p, o, r mod 3^(K-o))
    offs = []; tot = 0
    sizes = []
    for o in range(K + 1):
        offs.append(tot); sizes.append(3 ** (K - o)); tot += 2 * 3 ** (K - o)
    def idx(p, o, r): return offs[o] + p * sizes[o] + r
    dbl_to = np.empty(tot, dtype=np.int64)
    br_to = np.full(tot, -1, dtype=np.int64)
    o_of = np.empty(tot, dtype=np.int64)
    for o in range(K + 1):
        M = sizes[o]
        for p in (0, 1):
            for r in range(M):
                i = idx(p, o, r)
                o_of[i] = o
                dbl_to[i] = idx(0, o, (2 * r) % M)
                if p == 0 and o < K and r % 3 == 1:
                    br_to[i] = idx(1, o + 1, (r - 1) // 3)

    full = 2 * 3 ** K          # full-info classes = states with o=0
    def refinement_v(u):
        """v(state) = MIN of u over full-class refinements (leaf = ONE integer, class unknown
        among refinements -> rigorous lower bound takes worst case ONCE per block)."""
        v = np.empty(tot)
        v[offs[0]:offs[0] + full] = u
        for o in range(1, K + 1):
            Mo = sizes[o]; Mp = sizes[o - 1]
            for p in (0, 1):
                src = v[offs[o - 1] + p * Mp: offs[o - 1] + (p + 1) * Mp].reshape(3, Mo)
                v[offs[o] + p * Mo: offs[o] + (p + 1) * Mo] = src.min(axis=0)
        return v

    def apply_B(gamma, u):
        w1 = 2.0 ** (-gamma); w2 = 3.0 ** gamma
        v = refinement_v(u)
        f = np.zeros(tot)
        # horizon harvest (d=0): states with theta<1: D_taken = D, e = D - o
        hz = (D - o_of) > o_of * L23
        f[hz] = v[hz]
        for d in range(1, D + 1):
            nf = w1 * f[dbl_to]
            has_br = br_to >= 0
            nf[has_br] += w2 * f[br_to[has_br]]
            # exhaustion harvest: arriving at o=K with d steps remaining, steps_taken = D-d
            # e = (D-d) - K; harvest iff e > K*L23
            if (D - d) - K > K * L23:
                ex = o_of == K
                nf[ex] = v[ex]          # harvest overrides (leaf: stop propagating)
            f = nf
        return f[offs[0]:offs[0] + full]   # value at full-info roots

    def rho(gamma):
        u = np.ones(full)
        r = 1.0
        for it in range(power_iters):
            nu = apply_B(gamma, u)
            r = nu.max()
            if r <= 0: return 0.0
            u = nu / r
        return r

    lo, hi = 0.5, 1.0
    for _ in range(22):
        mid = (lo + hi) / 2
        if rho(mid) >= 1.0: lo = mid
        else: hi = mid
    return lo

for K, D in [(4, 40), (5, 50), (6, 60), (7, 70), (8, 80)]:
    g = solve(K, D)
    out[f"K={K},D={D}"] = round(g, 4)
    print(f"K={K} D={D}: amortized gamma = {g:.4f}", flush=True)
print(json.dumps(out, indent=1))
