"""
206_prefix_association.py
=========================
The shared-prefix association test (the FKG/EPW proof candidate for the
rich-rough coupling, Obs 409).

In a multiplicative cascade every top-triple member factors as
V = P * W (shared prefix x member factor), so
   log(level) = log P + log mean(W),
   log(spread) = log P + log spread(W),
and the shared log P term gives Cov(log L, log S) >= Var(log P) > 0 --
positive association WITH an explicit lower bound (EPW 1967 genre), and
persistence of the coupling = prefix variance bounded below (measured:
Var(F) ~ 1.4 flat; provable-adjacent via C_inj > 0).

The test of how multiplicative our tower really is: regress
log(spread_u) on log(level) per top-triple. Pure prefix-multiplicativity
predicts slope beta = 1; R^2 measures the prefix share; the trend of
(beta, corr, R^2) across k tells whether the coupling's strength is
stable (pro-G through provable structure) or eroding.

Systems 14..18, lam = 1.70. Spread taken of the mean field over the
triple (the object entering f2's denominator); levels are triple means.
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


print(f"shared-prefix association test at lam = {LAM}", flush=True)
print("  sys   corr(logL,logS)  slope    R2      Var(logL)  Var(logS)",
      flush=True)
for kp1 in (14, 15, 16, 17, 18):
    v = perron(kp1)
    Nl = v.size // 3
    mn = (v[:Nl] + v[Nl:2 * Nl] + v[2 * Nl:]) / 3.0
    Nl3 = Nl // 3
    M = np.stack([mn[:Nl3], mn[Nl3:2 * Nl3], mn[2 * Nl3:]])
    level = M.mean(axis=0)
    spread = M.std(axis=0)
    ok = spread > 0
    x = np.log(level[ok])
    y = np.log(spread[ok])
    x -= x.mean()
    y -= y.mean()
    vx = float(np.mean(x * x))
    vy = float(np.mean(y * y))
    cxy = float(np.mean(x * y))
    corr = cxy / (vx * vy) ** 0.5
    slope = cxy / vx
    r2 = corr ** 2
    print(f"  {kp1:2d}   {corr:12.5f}   {slope:.5f}  {r2:.5f}  "
          f"{vx:.5f}    {vy:.5f}", flush=True)
print("done", flush=True)
