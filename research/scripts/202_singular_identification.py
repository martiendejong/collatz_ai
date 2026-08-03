"""
202_singular_identification.py
==============================
Sub-question (B) of Obs 405: WHICH operator quantity bounds the realized
endpoint ratio d_k?

Candidate: d is a one-shot cross-depth ratio, so for a NON-NORMAL M the
right bound is the operator norm s_max(M) = largest singular value, not
the spectral radius sigma_W (measured flat 0.755*rho). Three numbers per
depth, frozen lam = 1.70:

  sigma_W/rho   spectral radius of M = P_W L P_W (known, script 201)
  s_max/rho     largest singular value of M (power method on M^T M;
                adjointness machine-checked per depth)
  r_real/rho    one-step Rayleigh growth of M on the REALIZED difference
                field x = P_W v (the actual Perron profile's top-scale
                deviation -- note P_W v IS the linear endpoint field)

Readings: if s_max/rho >= d_k everywhere and is flat/convergent, the
identification inequality "d <= s_max/rho" has the right direction and
shape for a proof; where r_real sits between sigma and s_max locates
the realization inside the operator's non-normal range.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)
D_OWN = {12: None, 13: 0.7535, 14: 0.7590, 15: 0.7662, 16: 0.7690,
         17: 0.7719}          # d_{k-1}? (profile column, for reference)


def build(k):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    inv4 = pow(4, -1, N)
    T4inv = ((i - 2) * inv4) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    Rfull = np.where(r == 0, (4 * s) % Nl, (2 * s + 1) % Nl)
    bfull = np.where(r == 0, B1, np.where(r == 2, B3, 0.0))
    return N, Nl, T4, T4inv, Rfull, bfull, (r != 1)


def perron(k, n_iter=300):
    N, Nl, T4, T4inv, Rfull, bfull, has = build(k)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[has] += bfull[has] * cb[Rfull[has]]
        v = w2 / w2.max()
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    w2 = A * v[T4]
    w2[has] += bfull[has] * cb[Rfull[has]]
    rho = float(w2.sum() / v.sum())
    stack = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    sel = stack.argmin(axis=0).astype(np.int64)
    return v, rho, sel


def make_ops(k, sel):
    N, Nl, T4, T4inv, Rfull, bfull, has = build(k)
    tgt = Rfull + sel[Rfull] * Nl

    def PW(x):
        m = (x[:Nl] + x[Nl:2 * Nl] + x[2 * Nl:]) / 3.0
        y = x.copy()
        y[:Nl] -= m
        y[Nl:2 * Nl] -= m
        y[2 * Nl:] -= m
        return y

    def L(x):
        y = A * x[T4]
        y[has] += bfull[has] * x[tgt[has]]
        return y

    bw = np.zeros(N)
    bw[has] = bfull[has]

    def LT(y):
        z = A * y[T4inv]
        z += np.bincount(tgt[has], weights=bw[has] * y[has], minlength=N)
        return z

    def M(x):
        return PW(L(PW(x)))

    def MT(y):
        return PW(LT(PW(y)))

    return PW, L, LT, M, MT, N


print(f"singular identification at frozen lam = {LAM}", flush=True)
print("  k    rho      adjoint_err  sigma_ref  s_max/rho  r_real/rho",
      flush=True)
for k in (12, 13, 14, 15, 16, 17):
    v, rho, sel = perron(k)
    PW, L, LT, M, MT, N = make_ops(k, sel)
    rng = np.random.default_rng(7)
    x0, y0 = rng.standard_normal(N), rng.standard_normal(N)
    aerr = abs(float(np.dot(L(x0), y0) - np.dot(x0, LT(y0)))) / (
        np.linalg.norm(x0) * np.linalg.norm(y0))
    # largest singular value: power method on M^T M
    z = PW(rng.standard_normal(N))
    z /= np.linalg.norm(z)
    smax2 = 0.0
    for _ in range(60):
        w = MT(M(z))
        smax2 = float(np.linalg.norm(w))
        z = w / smax2
    smax = smax2 ** 0.5
    # realized one-step Rayleigh growth on x = P_W v
    x = PW(v)
    r_real = float(np.linalg.norm(M(x)) / np.linalg.norm(x))
    print(f"  {k:2d}  {rho:.5f}  {aerr:.2e}   0.755      "
          f"{smax/rho:.5f}    {r_real/rho:.5f}", flush=True)
print("done", flush=True)
