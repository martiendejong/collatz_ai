"""
197_shell_tilt.py
=================
The last missing constant of the corrected Task-4 route (Obs 391):
the REALIZED TILT per domination shell.

Corrected Lemma-D route:  W{chain >= g} <= C'_tilt * env^g.
This script measures, for each chain depth g:
    count_g = counting mass of {chain >= g}      (uniform measure)
    W_g     = flow mass of {chain >= g}          (Perron weights)
    tilt_g  = W_g / count_g                      (realized tilt factor)
    env-normalized shell: W_g / env^g            (must stay bounded)
Prediction (Prop tilt / density-beats-tilt): tilt_g bounded, flat-ish in
g and k -- the flow measure does not concentrate on domination chains
faster than the counting envelope decays. Any growth of W_g/env^g in g
would refute the corrected route.
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


for k in (11, 13, 15):
    lam, v = edge_vector(k)
    N, Nl, T4, m1, R1, m3, R3 = make_maps(k)
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    env = (B1 + B3) / 3.0
    W = v / v.sum()
    F = np.log2(v)
    G = F[T4] - F

    stack = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    amin = stack.argmin(axis=0)
    r1full = np.zeros(N, dtype=np.int64)
    r1full[m1] = R1[m1]
    r1full[m3] = R3[m3]
    has = m1 | m3
    tgt = np.full(N, -1, dtype=np.int64)
    tgt[has] = r1full[has] + amin[r1full[has]] * Nl

    print(f"\nk={k}  lam={lam:.6f}  env={env:.4f}")
    for eps in (0.05, 0.10):
        t0 = -log2(eps * lam ** 2)
        dom = G <= -t0
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
        print(f"    eps={eps:.2f} (t0={t0:.2f})")
        print("      g   count_g     W_g       tilt=W/cnt  W_g/env^g")
        for g in range(1, 7):
            sel = glen >= g
            cnt = float(sel.sum()) / N
            wg = float(W[sel].sum())
            if cnt == 0:
                break
            print(f"      {g}  {cnt:9.6f}  {wg:9.6f}   {wg/cnt:8.3f}   "
                  f"{wg/env**g:8.4f}")
