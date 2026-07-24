"""
188_feed_chains.py
==================
THE OPEN LEMMA, REDUCED: flow decay of feed-dominated chains.

Structure of the reduction (research log Obs 358):
  - Intra-triple variation propagates DOWN the feed edges with weight
    phi(m) = feed-share of class m's Perron equation; backbone edges mix
    along the single N-cycle m -> 4m with geometric weights (lam^-2/rho)^j.
  - Weak contraction therefore requires CHAINS of consecutive
    feed-dominated classes: phi > 1-eps at m, at its feed target, etc.
  - Feed domination phi(m) ~ 1 requires the backbone side v(4m) to be
    desert-suppressed. Desert suppression is exponential in desert depth
    and desert density decays 3-adically -- so chain flow should decay
    GEOMETRICALLY in chain length g. If it does (uniformly in k), the
    L2/flow form of the Open Lemma follows by summation.

Measurement: at the feasibility edge (k = 12..15), for each class the
maximal g such that phi > 1-eps holds along g consecutive feed steps
starting there; report the FLOW carried by {g >= j} for j = 1..8,
for eps in {0.2, 0.1, 0.05}. Prediction (density model / Open Lemma):
geometric decay in j with ratio bounded away from 1, stable in k.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

def make_maps(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    return N, Nl, T4, (r == 0), R1, (r == 2), R3

def edge_vector(k):
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    lo_l, hi_l = 1.5, 1.999
    v = np.ones(N, dtype=np.float64)
    for _ in range(36):
        lam = 0.5 * (lo_l + hi_l)
        A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
        w = v.copy()
        for _ in range(60):
            cb = np.minimum(np.minimum(w[:Nl], w[Nl:2 * Nl]), w[2 * Nl:])
            w2 = A * w[T4]
            w2[m1] += B1 * cb[R1[m1]]
            w2[m3] += B3 * cb[R3[m3]]
            g = w2.max()
            w = w2 / g
        if g >= 1.0:
            lo_l, v = lam, w
        else:
            hi_l = lam
    lam = lo_l
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    for _ in range(300):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v = w2 / w2.max()
    return lam, v

print("feed-chain flow decay at the feasibility edge")
print("k   eps    flow(g>=1) g>=2     g>=3     g>=4     g>=5     g>=6     ratio")
for k in (12, 13, 14, 15):
    lam, v = edge_vector(k)
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    rho_v = A * v[T4]
    feed = np.zeros(N)
    feed[m1] = B1 * cb[R1[m1]]
    feed[m3] = B3 * cb[R3[m3]]
    tot = rho_v + feed
    phi = feed / tot
    # feed-target class (level k-1) lifted to level-k argmin lift for chain
    # continuation: follow the min lift (the one the min actually selects)
    tgt = np.full(N, -1, dtype=np.int64)
    stack = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    amin = stack.argmin(axis=0)              # which lift the min picks
    r1full = np.zeros(N, dtype=np.int64)
    r1full[m1] = R1[m1]
    r1full[m3] = R3[m3]
    has = m1 | m3
    tgt[has] = r1full[has] + amin[r1full[has]] * Nl
    for eps in (0.2, 0.1, 0.05):
        dom = (phi > 1 - eps)
        # chain length from m: g = number of consecutive feed-dominated
        # classes along the selected-min feed path starting at m
        glen = dom.astype(np.int64).copy()
        alive = dom & (tgt >= 0)
        pos = np.where(alive, tgt, 0)
        for _ in range(10):
            alive = alive & dom[pos]
            if not alive.any():
                break
            glen[alive] += 1
            nxt_ok = alive & (tgt[pos] >= 0)
            pos = np.where(nxt_ok, tgt[pos], pos)
            alive = nxt_ok
        W = v / v.sum()
        flows = [float(W[glen >= j].sum()) for j in range(1, 7)]
        rat = (flows[4] / flows[2]) ** 0.5 if flows[2] > 0 else float("nan")
        print(f"{k:<3} {eps:<5}  " + "  ".join(f"{f:.5f}" for f in flows)
              + f"   {rat:.3f}")
