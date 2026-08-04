"""
208_cv_convergence.py
=====================
Obs 413 -- direct CV convergence at frozen lam.

The floor->f2 translation (Obs 412) rests on CV(v^(k)) staying bounded
away from zero as k -> inf.  Two routes:

  (A) operator level: sigma_W/rho = 0.755 flat (Obs 405) -- measured
  (B) Perron-vector level: CV_k(lam) itself -- this script

For a primitive non-uniform operator (varying row weights b(i): B1 vs B3)
the Perron eigenvector is non-constant by Perron-Frobenius, so CV > 0 for
every finite k.  The question is whether CV_k -> CV_inf > 0 as k -> inf.

We measure CV_k = std(v^(k)) / mean(v^(k)) at lam = 1.70 for k = 12..18,
and track whether it converges to a positive limit.  The SHAPE factor of
f2 is essentially SHAPE ~ 1 - c * CV, so CV_inf > 0 <=> SHAPE_inf < 1 <=>
f2_inf < SCALE_inf < 1 (Conjecture G via the floor-to-f2 chain).

Second output: the fraction of eigenvector mass in the top quartile vs
the bottom quartile -- the "richness concentration" that drives the
rich-rough coupling.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)


def perron(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl   = N // 3
    R1   = (4 * s) % Nl
    R3   = (2 * s + 1) % Nl
    m1, m3 = (r == 0), (r == 2)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb   = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w2   = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v    = w2 / w2.max()
    v = v / v.mean()
    return v


print(f"CV convergence at fixed lam = {LAM}", flush=True)
print("  k   N          CV(v)      top-Q/mean  bot-Q/mean  top/bot", flush=True)

for k in range(12, 19):
    v = perron(k)
    cv = float(np.std(v) / np.mean(v))
    q25 = float(np.quantile(v, 0.25))
    q75 = float(np.quantile(v, 0.75))
    mean = float(np.mean(v))
    top_ratio = q75 / mean
    bot_ratio = q25 / mean
    print(f"  {k:2d}  {v.size:10d}  {cv:.6f}   {top_ratio:.4f}      {bot_ratio:.4f}      "
          f"{top_ratio/bot_ratio:.4f}", flush=True)

print("done", flush=True)
