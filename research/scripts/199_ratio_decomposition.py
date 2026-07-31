"""
199_ratio_decomposition.py
==========================
Decompose the endpoint-ratio creep (Obs 396/397) into a pure-depth
factor and a lambda-ladder factor.

Var_end(k) is a function of BOTH the depth k and the edge parameter
lam*(k). Exactly:
    r_k = Var(k+1, lam_{k+1}) / Var(k, lam_k)
        = [Var(k+1, lam_k) / Var(k, lam_k)]     (depth factor d_k)
        x [Var(k+1, lam_{k+1}) / Var(k+1, lam_k)]  (lambda factor l_k)
The cross terms Var(k+1, lam_k) are computed here (subcritical Perron
vectors: rho < 1, profile well-defined). If d_k is flat below 1 and the
creep lives in l_k, the fork reduces to the lam*-ladder geometry --
bounded (lam <= 2) and analytically charted territory (S1 monotonicity,
Thm 19 edge rates). If d_k itself creeps, the depth mechanism is the
open content.

Own-edge values (194/198b/198c): Var(13)=0.003850, Var(14)=0.003230,
Var(15)=0.002696, Var(16)=0.002274, Var(17)=0.001932 at
lam = 1.818824 / 1.830772 / 1.841968 / 1.852192 / 1.861680.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

OWN = {13: (1.818824, 0.003850), 14: (1.830772, 0.003230),
       15: (1.841968, 0.002696), 16: (1.852192, 0.002274),
       17: (1.861680, 0.001932)}


def perron(k, lam, n_iter=250):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    m1, m3 = (r == 0), (r == 2)
    R1m, R3m = R1[m1], R3[m3]
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


print("ratio decomposition r_k = depth d_k x lambda l_k", flush=True)
print("  k   Var(k+1,lam_k)  d_k     l_k     r_k(own)", flush=True)
for k in (13, 14, 15, 16):
    lam_k, var_k = OWN[k]
    _, var_k1 = OWN[k + 1]
    v = perron(k + 1, lam_k)
    cross = var_end(k + 1, v)
    d = cross / var_k
    l = var_k1 / cross
    print(f"  {k}  {cross:.6f}       {d:.4f}  {l:.4f}  {var_k1/var_k:.4f}",
          flush=True)
print("done", flush=True)
