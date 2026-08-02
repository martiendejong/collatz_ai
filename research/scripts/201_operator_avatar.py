"""
201_operator_avatar.py
======================
The operator avatar of Conjecture G: measure the L2 spectral radius of
the compressed linearized difference operator across depth, at frozen
lam = 1.70.

Objects (per depth k, frozen lam):
  v*      Perron profile (300 iterations), growth rho(k).
  sel(r)  = argmin_t v*(r + t*Nl)  (frozen argmins).
  L       delta -> A*delta(T4 i) + b(i)*delta(R(i) + sel(R(i))*Nl)
          (the linearization of the K-L operator at v*).
  P_W     projection onto W = {zero mean on every top-digit triple}.
  M       = P_W o L o P_W;  sigma_W(k) = |dominant eigenvalue| of M,
          estimated by the late-iteration norm growth of the power
          method (robust to complex/negative dominant pairs).

Payload: sigma_W(k)/rho(k) across k = 12..17 -- the candidate operator
form of the endpoint contraction d. REGISTERED (loose): if the operator
route is right, sigma_W/rho is flat in k and < 1; its value identifies
which measured constant it embodies (d ~ 0.77, kappa_deep ~ 0.84, or
the within-tower ratio ~ 0.70); a k-creep in sigma_W/rho would mirror
the g-column and locate the creep at operator level.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)


def build(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    Rfull = np.where(r == 0, (4 * s) % Nl, (2 * s + 1) % Nl)
    bfull = np.where(r == 0, B1, np.where(r == 2, B3, 0.0))
    has = r != 1
    return N, Nl, T4, Rfull, bfull, has


def perron(k, n_iter=300):
    N, Nl, T4, Rfull, bfull, has = build(k)
    v = np.ones(N, dtype=np.float64)
    g = 1.0
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[has] += bfull[has] * cb[Rfull[has]]
        g = w2.max()
        v = w2 / g
    # growth measured as Rayleigh-type ratio at the fixed point
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    w2 = A * v[T4]
    w2[has] += bfull[has] * cb[Rfull[has]]
    rho = float(w2.sum() / v.sum())
    stack = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    sel = stack.argmin(axis=0).astype(np.int64)
    return v, rho, sel


def sigma_W(k, sel, n_iter=500, seed=1):
    N, Nl, T4, Rfull, bfull, has = build(k)
    tgt = Rfull + sel[Rfull] * Nl          # selected-lift feed index
    rng = np.random.default_rng(seed)
    d = rng.standard_normal(N)

    def PW(x):
        m = (x[:Nl] + x[Nl:2 * Nl] + x[2 * Nl:]) / 3.0
        y = x.copy()
        y[:Nl] -= m
        y[Nl:2 * Nl] -= m
        y[2 * Nl:] -= m
        return y

    d = PW(d)
    d /= np.linalg.norm(d)
    rates = []
    for it in range(n_iter):
        y = A * d[T4]
        y[has] += bfull[has] * d[tgt[has]]
        y = PW(y)
        nrm = np.linalg.norm(y)
        rates.append(nrm)
        d = y / nrm
    # late-window geometric mean of per-step norm growth
    tail = np.array(rates[-100:])
    return float(np.exp(np.mean(np.log(tail))))


print(f"operator avatar at frozen lam = {LAM}", flush=True)
print("  k    rho       sigma_W    sigma_W/rho", flush=True)
for k in (12, 13, 14, 15, 16, 17):
    v, rho, sel = perron(k)
    sw = sigma_W(k, sel)
    print(f"  {k:2d}  {rho:.6f}  {sw:.6f}   {sw/rho:.5f}", flush=True)
print("done", flush=True)
