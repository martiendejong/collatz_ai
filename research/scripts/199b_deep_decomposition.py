"""
199b_deep_decomposition.py
==========================
Extend the ratio decomposition (script 199) to the deep end: d_17 and
d_18 -- the two deepest pure-depth contraction factors reachable on
this machine.

Cross terms computed sequentially (RAM):
  Var(18, lam_17 = 1.861680)   -- 129M classes, ~5 GB
  Var(19, lam_18 = 1.870749)   -- 387M classes, chunked, ~9 GB
Then, with own-edge values Var(17) = 0.001932, Var(18) = 0.001651,
Var(19) = 0.001406:
  d_17 = Var(18, lam_17)/Var(17)   l_17 = Var(18, lam_18)/Var(18, lam_17)
  d_18 = Var(19, lam_18)/Var(18)   l_18 = Var(19, lam_19)/Var(19, lam_18)
Deciding question (rem:decomp of density_one.tex): does d_k keep
creeping (0.788/0.787/0.797/0.805 at k=13..16) or saturate below 1?
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)


def var_end_of(v, Nl):
    T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    X = np.log2(T) - np.log2(T.mean(axis=0))[None, :]
    return float(np.var(X))


def perron_inmem(k, lam, n_iter=250):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = ((4 * i + 2) % N).astype(np.int32)
    s = (i // 3).astype(np.int64)
    r = (i % 3).astype(np.int8)
    del i
    Nl = N // 3
    m1, m3 = (r == 0), (r == 2)
    del r
    R1m = ((4 * s[m1]) % Nl).astype(np.int32)
    R3m = ((2 * s[m3] + 1) % Nl).astype(np.int32)
    del s
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    v = np.ones(N, dtype=np.float64)
    g = 1.0
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1m]
        w2[m3] += B3 * cb[R3m]
        g = w2.max()
        v = w2 / g
        if (it + 1) % 100 == 0:
            print(f"    k={k} iter {it+1} growth {g:.8f}", flush=True)
    return v, g, Nl


def perron_chunked(k, lam, n_iter=250, chunk=3 ** 16):
    N = 3 ** (k - 1)
    Nl = N // 3
    i64 = np.arange(N, dtype=np.int64)
    T4 = ((4 * i64 + 2) % N).astype(np.int32)
    del i64
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    v = np.ones(N, dtype=np.float64)
    w2 = np.empty(N, dtype=np.float64)
    g = 1.0
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        for a in range(0, N, chunk):
            b = a + chunk
            idx = np.arange(a, b, dtype=np.int64)
            w2[a:b] = A * v[T4[a:b]]
            r = (idx % 3).astype(np.int8)
            s = idx // 3
            m1 = r == 0
            m3 = r == 2
            w2[a:b][m1] += B1 * cb[(4 * s[m1]) % Nl]
            w2[a:b][m3] += B3 * cb[(2 * s[m3] + 1) % Nl]
            del idx, r, s, m1, m3
        g = w2.max()
        np.divide(w2, g, out=v)
        if (it + 1) % 100 == 0:
            print(f"    k={k} iter {it+1} growth {g:.8f}", flush=True)
    return v, g, Nl


VAR17, VAR18, VAR19 = 0.001932, 0.001651, 0.001406

print("deep decomposition: d_17, d_18", flush=True)
v, g, Nl = perron_inmem(18, 1.861680)
cross18 = var_end_of(v, Nl)
del v
d17 = cross18 / VAR17
l17 = VAR18 / cross18
print(f"  Var(18, lam_17) = {cross18:.6f}  growth {g:.8f}", flush=True)
print(f"  d_17 = {d17:.4f}   l_17 = {l17:.4f}   r_17(own) = "
      f"{VAR18/VAR17:.4f}", flush=True)

v, g, Nl = perron_chunked(19, 1.870749)
cross19 = var_end_of(v, Nl)
del v
d18 = cross19 / VAR18
l18 = VAR19 / cross19
print(f"  Var(19, lam_18) = {cross19:.6f}  growth {g:.8f}", flush=True)
print(f"  d_18 = {d18:.4f}   l_18 = {l18:.4f}   r_18(own) = "
      f"{VAR19/VAR18:.4f}", flush=True)
print(f"  d-series: 0.7878/0.7873/0.7967/0.8049/{d17:.4f}/{d18:.4f}",
      flush=True)
print("done", flush=True)
