"""
200_frozen_lambda_curve.py
==========================
THE FROZEN-LAMBDA CURVE: d_k(lambda) on a grid INCLUDING lambda = 2.

The Open Lemma in its cleanest form (Obs 399/400): the frozen-lambda
limit contraction stays below 1 uniformly as lambda -> 2. The own-edge
series could never reach lambda = 2 (lam_k -> 2 entangled with k); the
frozen grid CAN: at fixed lambda -- even 2.0 exactly -- the depth-k
system is merely subcritical (rho < 1) and its Perron profile is
well-defined.

Measured here: Var_end(k, lam) for k = 13..17 and
lam in {1.70, 1.80, 1.85, 1.90, 1.95, 2.00}; then
d_k(lam) = Var_end(k+1, lam)/Var_end(k, lam) for the four pairs.

Readings:
  - d_k(lam) flat in k at fixed lam  => the frozen-lambda limit exists
    and is measured directly; the lambda = 2 column is the endpoint
    contraction itself.
  - d_k(2.0) flat below 1  => direct evidence for the Open Lemma at the
    actual endpoint (gamma -> 1 reading).
  - d_k(2.0) drifting up with k  => the ceiling mechanism, located.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
GRID = (1.70, 1.80, 1.85, 1.90, 1.95, 2.00)
KS = (13, 14, 15, 16, 17)


def perron(k, lam, n_iter=300):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1m = ((4 * s) % Nl)[r == 0]
    R3m = ((2 * s + 1) % Nl)[r == 2]
    m1, m3 = (r == 0), (r == 2)
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1m]
        w2[m3] += B3 * cb[R3m]
        v = w2 / w2.max()
    return v


def var_end(k, v):
    Nl = 3 ** (k - 2)
    T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    X = np.log2(T) - np.log2(T.mean(axis=0))[None, :]
    return float(np.var(X))


print("frozen-lambda curve: Var_end(k, lam) and d_k(lam)", flush=True)
for lam in GRID:
    vs = {}
    for k in KS:
        vs[k] = var_end(k, perron(k, lam))
    ds = [vs[KS[j + 1]] / vs[KS[j]] for j in range(len(KS) - 1)]
    print(f"  lam={lam:.2f}  Var: " +
          " ".join(f"{vs[k]:.6f}" for k in KS), flush=True)
    print(f"           d_13/d_14/d_15/d_16: " +
          " ".join(f"{d:.4f}" for d in ds), flush=True)
print("done", flush=True)
