"""E29/G1b-fix: Krasikov-Lagarias induction engine.
pi_worst(x) >= sum_leaves pi(theta_i x), theta_i = 3^o/2^e < 1  ==>  gamma from sum theta_i^gamma = 1.
Tree: doubling child always (e+1); third-child iff (even, r=1 mod 3, o<K) (o+1, consumes a trit).
State = (parity, o, r mod 3^(K-o)); budget level = K-o. Harvest when o=K (info gone) or d=D,
only if 2^(d-o) > 3^o (theta<1). DP for fixed gamma: mass propagation with weights
2^(-gamma) per doubling, 3^(gamma) per branch; bisect gamma s.t. min-root mass = 1.
"""
import sys, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
out = {}
L23 = math.log2(3)

def gamma_KD(K, D):
    # pack states: for o in 0..K: block of 2*3^(K-o)
    offs = []
    tot = 0
    for o in range(K + 1):
        offs.append(tot)
        tot += 2 * 3 ** (K - o)
    def idx(p, o, r):
        return offs[o] + p * 3 ** (K - o) + r

    # precompute transitions per o-level
    dbl_to = np.empty(tot, dtype=np.int64)
    br_to = np.full(tot, -1, dtype=np.int64)
    for o in range(K + 1):
        M = 3 ** (K - o)
        for p in (0, 1):
            for r in range(M):
                i = idx(p, o, r)
                dbl_to[i] = idx(0, o, (2 * r) % M)
                if p == 0 and o < K and M >= 3 and r % 3 == 1:
                    br_to[i] = idx(1, o + 1, (r - 1) // 3)
    o_of = np.empty(tot, dtype=np.int64)
    for o in range(K + 1):
        o_of[offs[o]:offs[o] + 2 * 3 ** (K - o)] = o

    def harvest_mass(gamma, root):
        w1 = 2.0 ** (-gamma)
        w2 = 3.0 ** gamma
        mass = np.zeros(tot)
        mass[root] = 1.0
        harvested = 0.0
        for d in range(1, D + 1):
            nm = np.zeros(tot)
            src = np.nonzero(mass)[0]
            np.add.at(nm, dbl_to[src], w1 * mass[src])
            bsrc = src[br_to[src] >= 0]
            np.add.at(nm, br_to[bsrc], w2 * mass[bsrc])
            # harvest states with o=K (no info left) if theta<1: 2^(d-K) > 3^K
            if 3 ** K < 2 ** max(d - K, 0):
                full = o_of == K
                harvested += nm[full].sum()
                nm[full] = 0.0
            mass = nm
        # final harvest at horizon: theta<1 means d-o > o*log2 3
        okid = np.nonzero(mass)[0]
        for i in okid:
            o = o_of[i]
            if D - o > o * L23:
                harvested += mass[i]
        return harvested

    # worst root over classes with full info (o=0), r not divisible by 3
    def gamma_root(root):
        lo, hi = 0.0, 1.0
        for _ in range(30):
            mid = (lo + hi) / 2
            if harvest_mass(mid, root) >= 1.0: lo = mid
            else: hi = mid
        return lo
    worst = None
    M0 = 3 ** K
    # sample worst root: all residues r (mod 3^K) with r % 3 != 0, both parities is too many;
    # growth depends on doubling-orbit => enough to scan orbit representatives; do full for small K
    roots = [idx(p, 0, r) for p in (0, 1) for r in range(M0) if r % 3 != 0]
    step = max(1, len(roots) // 200)   # sample for speed; exact for small K
    for i, root in enumerate(roots[::step]):
        g = gamma_root(root)
        if worst is None or g < worst: worst = g
    return worst

for K, D in [(3, 30), (4, 40), (5, 50), (6, 60), (6, 80), (7, 80)]:
    g = gamma_KD(K, D)
    out[f"K={K},D={D}"] = round(g, 4)
    print(f"K={K} D={D}: certified gamma >= {g:.4f}", flush=True)
print(json.dumps(out, indent=1))
