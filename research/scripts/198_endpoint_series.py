"""
198_endpoint_series.py
======================
THE ENDPOINT ROUTE (Obs 395): extend the endpoint-variance series to
k = 16, 17.

Motivation: the tower endpoint variance Var_count(Xt_{k-1}) measured at
k = 11..15 (script 194) decays geometrically IN k at rate
0.829/0.833/0.839/0.836 ~ 0.835 -- numerically identical to the
attenuation constant kappa_deep = 0.839 +- 0.002 from the independent
linearized instrument. Since (Samuelson + edge-rate, both proved)
    1 - gamma_k <= 4.917 * CV_w,top(k),
and CV_w,top^2 tracks the endpoint variance up to bounded factors,
GEOMETRIC ENDPOINT DECAY ALONE implies gamma -> 1 -- no Summation
Theorem, no chain recursion, no multi-scale uniformity.

This script extends the series two depths:
  k=16: short warm-started bisection for lam* in [1.846, 1.862].
  k=17: fixed lam = 1.86168 (the exact-certified value, cert_k17), power
        iteration only.
Checkpoints every 50 iterations (the k=20 lesson: never run long
computations without checkpoints). Output flushed line by line.

Reported per k: lam, Var_count(Xt_end), per-triple linear CV stats,
flow-weighted CV_w,top, q_k, 1-q vs sqrt2*CV_w margin, and the ratio to
the previous depth.
"""
import os
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
OUT = os.path.join(os.path.dirname(__file__), "..", "certificates")


def make_maps(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    return N, Nl, T4, (r == 0), R1, (r == 2), R3


def power_iter(v, n_iter, lam, maps, ckpt=None):
    N, Nl, T4, m1, R1, m3, R3 = maps
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    g = 1.0
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        g = w2.max()
        v = w2 / g
        if ckpt and (it + 1) % 50 == 0:
            np.save(ckpt, v)
            print(f"      ckpt iter {it+1}, growth {g:.8f}", flush=True)
    return v, g


def endpoint_stats(k, lam, v):
    N = 3 ** (k - 1)
    Nl = N // 3
    T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    mean = T.mean(axis=0)
    mn = T.min(axis=0)
    sig = T.std(axis=0)
    # endpoint tower increment: Xt_end(i) = log2 v(i) - log2 tripmean
    X = np.log2(T) - np.log2(mean)[None, :]
    var_end = float(np.var(X))
    cv = sig / mean
    q_k = float(3.0 * mn.sum() / v.sum())
    cvw = float(np.sqrt((mean * cv ** 2).sum() / mean.sum()))
    lin = float(np.sqrt(2.0) * sig.sum() / mean.sum())
    print(f"  k={k}  lam={lam:.6f}", flush=True)
    print(f"    Var_count(Xt_end) = {var_end:.6f}", flush=True)
    print(f"    1-q = {1-q_k:.6f}   sqrt2*E_w[cv] = {lin:.6f}   "
          f"sqrt2*CV_w = {np.sqrt(2)*cvw:.6f}", flush=True)
    print(f"    Samuelson margin = {lin/(1-q_k):.3f}   "
          f"CV_w,top = {cvw:.6f}", flush=True)
    return var_end


print("endpoint series extension (k=16, 17)", flush=True)
prev_var = 0.002696          # k=15 value from script 194 run (Var at p=14)

# ---- k = 16: warm-started bisection ---------------------------------------
k = 16
maps = make_maps(k)
lo_l, hi_l = 1.846, 1.862
v = np.ones(3 ** (k - 1), dtype=np.float64)
for step in range(14):
    lam = 0.5 * (lo_l + hi_l)
    w, g = power_iter(v.copy(), 50, lam, maps)
    if g >= 1.0:
        lo_l, v = lam, w
    else:
        hi_l = lam
    print(f"    bisect step {step}: lam={lam:.6f} growth={g:.8f}", flush=True)
lam16 = lo_l
v, g = power_iter(v, 250, lam16, maps,
                  ckpt=os.path.join(OUT, "k16_eig_198.npy"))
np.save(os.path.join(OUT, "k16_eig_198.npy"), v)
var16 = endpoint_stats(k, lam16, v)
print(f"    ratio vs k=15: {var16/prev_var:.4f}", flush=True)

# ---- k = 17: certified lambda, power iteration only -----------------------
k = 17
maps = make_maps(k)
lam17 = 1.86168              # exact-certified edge (cert_k17, R30-47)
v = np.ones(3 ** (k - 1), dtype=np.float64)
v, g = power_iter(v, 350, lam17, maps,
                  ckpt=os.path.join(OUT, "k17_eig_198.npy"))
np.save(os.path.join(OUT, "k17_eig_198.npy"), v)
print(f"    final growth {g:.8f} (should be ~1: certified edge)", flush=True)
var17 = endpoint_stats(k, lam17, v)
print(f"    ratio vs k=16: {var17/var16:.4f}", flush=True)
print("done", flush=True)
