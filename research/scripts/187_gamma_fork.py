"""
187_gamma_fork.py
=================
DECIDING STATISTICS FOR THE GAMMA FORK (PREDICTIONS.md, R2250):
  CEILING model: gamma_inf = H(1/log2 3) = 0.9507 -- the LP saturates at the
                 multifractal worst local exponent (= dim D, Obs 336!).
  DENSITY model: gamma_inf = 1.

gamma(k) itself cannot discriminate before k ~ 27. But the models differ NOW
in WHERE the triple spread (the min-loss 1-q) lives:
  - CEILING needs a POSITIVE-FLOW fraction of triples to keep O(1) spread
    (the min couples to a fat set of thin backward directions).
  - DENSITY predicts the spread flees to flow-starved desert classes:
    flow-weighted spread -> 0 geometrically while worst-case spread may
    persist on a measure->0 set (cf. Obs 327: sup fails, flow works).

Measurements at the feasibility edge, k = 10..15:
 M1: q(k) = 3*sum(cbar)/sum(c) and flow-weighted triple loss
     Lflow(k) = sum_t w_t * log2(mean_t/min_t) / sum_t w_t, w_t = triple flow.
 M2: spread-mass curve: fraction of FLOW in triples with loss > eps,
     eps in {0.01, 0.05, 0.2} -- ceiling: stabilizes; density: -> 0.
 M3: desert localization: are the worst-loss triples the 3-adic deserts?
     Correlate triple loss with desert depth d3(m) = v3(m+1) of the triple's
     base class (desert theorem, Obs 319-320).
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

def edge_vector(k, iters=400):
    """Perron vector at (near) the feasibility edge lambda*(k)."""
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    lo_l, hi_l = 1.5, 1.999
    v = np.ones(N, dtype=np.float64)
    for _ in range(40):                      # bisection on rho(lambda) = 1
        lam = 0.5 * (lo_l + hi_l)
        A = lam ** -2.0
        B1 = lam ** (ALPHA - 2.0)
        B3 = lam ** (ALPHA - 1.0)
        w = v.copy()
        g = 1.0
        for _ in range(60):
            cb = np.minimum(np.minimum(w[:Nl], w[Nl:2 * Nl]), w[2 * Nl:])
            w2 = A * w[T4]
            w2[m1] += B1 * cb[R1[m1]]
            w2[m3] += B3 * cb[R3[m3]]
            g = w2.max()
            w = w2 / g
        if g >= 1.0:
            lo_l = lam
            v = w                            # keep best feasible-side vector
        else:
            hi_l = lam
    # polish at lo_l
    lam = lo_l
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    for _ in range(iters):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v = w2 / w2.max()
    return lam, v, Nl

print("k   gamma*   q(k)      Lflow      P(flow: loss>.01) >.05    >.2     "
      "corr(loss,desert)")
for k in range(10, 16):
    lam, v, Nl = edge_vector(k)
    t = v.reshape(3, Nl)                     # triples: lifts r, r+Nl, r+2Nl
    tmin = t.min(axis=0)
    tmean = t.mean(axis=0)
    tsum = t.sum(axis=0)                     # triple flow weight
    loss = np.log2(tmean / tmin)
    W = tsum / tsum.sum()
    q = 3.0 * tmin.sum() / v.sum()
    Lflow = float((W * loss).sum())
    p1 = float(W[loss > 0.01].sum())
    p5 = float(W[loss > 0.05].sum())
    p20 = float(W[loss > 0.2].sum())
    # desert depth of the triple base class r -> m = 3r+2 mod 3^(k-1); depth
    # d3 = v3(m+1) (desert theorem invariant)
    r = np.arange(Nl, dtype=np.int64)
    m = 3 * r + 2
    x = m + 1
    d3 = np.zeros(Nl, dtype=np.int64)
    xx = x.copy()
    for _ in range(k):
        div = (xx % 3 == 0) & (xx > 0)
        if not div.any():
            break
        d3[div] += 1
        xx = np.where(div, xx // 3, xx)
    corr = float(np.corrcoef(loss, d3)[0, 1])
    print(f"{k:<3} {log2(lam):.5f}  {q:.5f}  {Lflow:.6f}   {p1:.4f}"
          f"          {p5:.4f}  {p20:.5f}   {corr:+.3f}")
print()
print("CEILING predicts the flow-fractions stabilize; DENSITY predicts")
print("geometric decay (spread flees to flow-starved deserts, corr > 0).")
