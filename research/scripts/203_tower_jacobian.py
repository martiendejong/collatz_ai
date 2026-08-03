"""
203_tower_jacobian.py
=====================
The cross-depth mechanism, factored exactly (sub-question B, final
form). Push the top-triple projection P_W through the depth-(k+1)
fixed-point equation: backbone maps top triples to top triples and the
feed target's top digit tracks the source's top digit, so EXACTLY

    x = (1/g) [ A * x o T4  +  b * y o R ],

with x = P_W v^(k+1) (top deviation of the profile) and y = P'_W cb
(top deviation of the min field, one level down). Hence the cross-depth
ratio factors as

    d_lin = ( f1 * f2 * f3 )^2,
    f1 = |x_{k+1}| / |y|      resolvent gain (within-system, operator)
    f2 = |y| / |u|            min/mean top-deviation ratio (q-machinery)
    f3 = |u| / |x_k|          tower correspondence (cross-system)

where u = P'_W (mean of lifts) and x_k = P_W v^(k). All profiles are
mean-normalized so relative deviations are comparable across systems.
THE question: which factor carries the creep. A flat f1 and f2 with
creeping f3 would locate the entire fork in the block-equation
correspondence (measured corr 0.989-0.995 at k=13 -- its slow approach
to 1 would BE the creep, and its limit existence would BE Conjecture G).

Measured at lam = 1.70, pairs (k, k+1) for k = 13..17.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)


def perron(k, n_iter=300):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1m = ((4 * s) % Nl)[r == 0]
    R3m = ((2 * s + 1) % Nl)[r == 2]
    m1, m3 = (r == 0), (r == 2)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1m]
        w2[m3] += B3 * cb[R3m]
        v = w2 / w2.max()
    return v / v.mean()          # mean-normalized profile


def top_dev(field):
    """top-triple deviation of a field (triples spaced len/3)."""
    n = field.size // 3
    m = (field[:n] + field[n:2 * n] + field[2 * n:]) / 3.0
    out = np.empty_like(field)
    out[:n] = field[:n] - m
    out[n:2 * n] = field[n:2 * n] - m
    out[2 * n:] = field[2 * n:] - m
    return out, m


profiles = {k: perron(k) for k in range(13, 19)}
print(f"tower-jacobian factorisation at lam = {LAM}", flush=True)
print("  k->k+1   f1(resolvent)  f2(min/mean)  f3(tower)   "
      "d_lin     d_log(ref)", flush=True)
DLOG = {13: 0.7590, 14: 0.7662, 15: 0.7690, 16: 0.7719, 17: 0.7753}
for k in range(13, 18):
    vk1 = profiles[k + 1]
    vk = profiles[k]
    Nl = vk1.size // 3
    x1, _ = top_dev(vk1)
    lifts = np.stack([vk1[:Nl], vk1[Nl:2 * Nl], vk1[2 * Nl:]])
    cb = lifts.min(axis=0)
    mn = lifts.mean(axis=0)
    y, _ = top_dev(cb)
    u, _ = top_dev(mn)
    xk, _ = top_dev(vk)
    n1 = float(np.linalg.norm(x1)) / np.sqrt(x1.size)
    ny = float(np.linalg.norm(y)) / np.sqrt(y.size)
    nu = float(np.linalg.norm(u)) / np.sqrt(u.size)
    nk = float(np.linalg.norm(xk)) / np.sqrt(xk.size)
    f1, f2, f3 = n1 / ny, ny / nu, nu / nk
    dlin = (n1 / nk) ** 2
    print(f"  {k}->{k+1}   {f1:.5f}       {f2:.5f}      {f3:.5f}    "
          f"{dlin:.5f}   {DLOG[k]:.4f}", flush=True)
print("done", flush=True)
