"""
204_smoothing_mechanism.py
==========================
WHY does the min-smoothing weaken (Obs 407: f2 = 0.861 -> 0.882,
+0.005/depth)? The smoothing lives where the ARGMIN SWITCHES within a
top-triple (the min clips tops); on constant-argmin triples the min
follows a single lift and smooths nothing. Two candidate sources for
the creep, both measured here per depth (systems k+1 = 14..18,
lam = 1.70):

  (a) the switching rate falls with depth, or
  (b) the three lifts' top-deviations correlate increasingly
      (homogenization: the min co-moves with everything).

Also measured: conditional f2 on same-argmin vs switching triples, the
per-lift deviation norms (sqrt3-type factor vs the mean field), and the
share-weighted reconstruction of f2^2 (sanity: parts must reassemble
the whole).
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
    return v / v.mean()


def top_dev(field):
    n = field.size // 3
    m = (field[:n] + field[n:2 * n] + field[2 * n:]) / 3.0
    out = np.empty_like(field)
    out[:n] = field[:n] - m
    out[n:2 * n] = field[n:2 * n] - m
    out[2 * n:] = field[2 * n:] - m
    return out


print(f"smoothing mechanism at lam = {LAM}", flush=True)
print("  sys   switch%  corr(lifts)  |dev_lift|/|u|  f2_same  f2_switch"
      "   f2_all", flush=True)
for kp1 in (14, 15, 16, 17, 18):
    v = perron(kp1)
    N = v.size
    Nl = N // 3
    L0, L1, L2 = v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]
    cb = np.minimum(np.minimum(L0, L1), L2)
    mn = (L0 + L1 + L2) / 3.0
    amin = np.stack([L0, L1, L2]).argmin(axis=0)

    y = top_dev(cb)
    u = top_dev(mn)
    d0, d1, d2 = top_dev(L0.copy()), top_dev(L1.copy()), top_dev(L2.copy())

    # pairwise correlation of lift top-deviations
    def corr(a, b):
        return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    cmean = (corr(d0, d1) + corr(d0, d2) + corr(d1, d2)) / 3.0
    liftnorm = (np.linalg.norm(d0) + np.linalg.norm(d1)
                + np.linalg.norm(d2)) / 3.0 / np.linalg.norm(u)

    # switching triples in the base space
    Nl3 = Nl // 3
    a0, a1, a2 = amin[:Nl3], amin[Nl3:2 * Nl3], amin[2 * Nl3:]
    same = (a0 == a1) & (a1 == a2)
    switch_rate = 1.0 - float(same.mean())
    mask_same = np.concatenate([same, same, same])
    ns2 = float(np.sum(y[mask_same] ** 2)), float(np.sum(u[mask_same] ** 2))
    nw2 = (float(np.sum(y[~mask_same] ** 2)),
           float(np.sum(u[~mask_same] ** 2)))
    f2_same = (ns2[0] / ns2[1]) ** 0.5 if ns2[1] > 0 else float("nan")
    f2_switch = (nw2[0] / nw2[1]) ** 0.5 if nw2[1] > 0 else float("nan")
    f2_all = float(np.linalg.norm(y) / np.linalg.norm(u))
    print(f"  {kp1:2d}   {switch_rate*100:6.3f}  {cmean:9.5f}   "
          f"{liftnorm:10.5f}    {f2_same:.5f}  {f2_switch:.5f}   "
          f"{f2_all:.5f}", flush=True)
print("done", flush=True)
