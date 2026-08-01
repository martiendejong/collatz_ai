"""
200b_g_sequence.py
==================
The lambda-free sequence g(k): extend the frozen-lambda column at
lam = 1.70 to k = 18, 19.

Script 200 found the d-matrix approximately additively separable,
d_k(lam) ~ f(lam) + g(k): lambda-span per k is ~constant (0.064) and
k-drift per lambda is ~constant (0.010) -- so the Open Lemma's hard
"uniform as lam -> 2" clause dissolves empirically, and the open core
reduces to the convergence of ONE lambda-independent sequence g(k),
measurable at ANY lambda. We measure it at lam = 1.70 (deep
subcritical: cheapest, cleanest spectral gap).

Known column (script 200): Var_end(k, 1.70) =
0.001976/0.001494/0.001126/0.000855/0.000655 at k = 13..17,
d = 0.7560/0.7535/0.7590/0.7662.

REGISTERED PREDICTION (separability + own-edge g-increments):
d_17(1.70) = 0.771 +- 0.004, d_18(1.70) = 0.775 +- 0.006, increments
shrinking if g saturates. Persistent non-shrinking increments would
keep the fork open in its final lambda-free form.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)
VAR17 = 0.000655
CHUNK = 3 ** 16


def var_end_of(v, Nl):
    T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
    X = np.log2(T) - np.log2(T.mean(axis=0))[None, :]
    return float(np.var(X))


def perron_inmem(k, n_iter=300):
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
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1m]
        w2[m3] += B3 * cb[R3m]
        v = w2 / w2.max()
    return v, Nl


def perron_chunked(k, n_iter=300):
    N = 3 ** (k - 1)
    Nl = N // 3
    i64 = np.arange(N, dtype=np.int64)
    T4 = ((4 * i64 + 2) % N).astype(np.int32)
    del i64
    v = np.ones(N, dtype=np.float64)
    w2 = np.empty(N, dtype=np.float64)
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        for a in range(0, N, CHUNK):
            b = a + CHUNK
            idx = np.arange(a, b, dtype=np.int64)
            w2[a:b] = A * v[T4[a:b]]
            r = (idx % 3).astype(np.int8)
            s = idx // 3
            m1 = r == 0
            m3 = r == 2
            w2[a:b][m1] += B1 * cb[(4 * s[m1]) % Nl]
            w2[a:b][m3] += B3 * cb[(2 * s[m3] + 1) % Nl]
            del idx, r, s, m1, m3
        np.divide(w2, w2.max(), out=v)
        if (it + 1) % 100 == 0:
            print(f"    k={k} iter {it+1}", flush=True)
    return v, Nl


print(f"g-sequence at frozen lam = {LAM}", flush=True)
v, Nl = perron_inmem(18)
var18 = var_end_of(v, Nl)
del v
d17 = var18 / VAR17
print(f"  Var_end(18) = {var18:.6f}   d_17(1.70) = {d17:.4f}", flush=True)

v, Nl = perron_chunked(19)
var19 = var_end_of(v, Nl)
del v
d18 = var19 / var18
print(f"  Var_end(19) = {var19:.6f}   d_18(1.70) = {d18:.4f}", flush=True)
print(f"  column d(1.70): 0.7560/0.7535/0.7590/0.7662/{d17:.4f}/{d18:.4f}",
      flush=True)
print("done", flush=True)
