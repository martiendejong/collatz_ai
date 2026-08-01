"""
200c_g_k20.py
=============
Seventh point of the lambda-free g-column: Var_end(20, lam = 1.70)
=> d_19(1.70). Registered expectation: increments have plateaued at
~+0.003 (Obs 402); a shrinking increment (< +0.002) supports g
saturating at ~0.777 (=> d_inf(2.0) ~ 0.841 = kappa_deep); a persistent
+0.003 keeps the slow-drift fork alive.

float32 chunked (validated pipeline: reproduces recorded 1-q to 4e-5),
no stored feed maps, ~11 GB peak, checkpoints every 25 iterations,
200 iterations (deep subcritical: fast convergence).
Known: Var_end(19, 1.70) = 0.000389 (script 200b).
"""
import os
import numpy as np
from math import log2

ALPHA = log2(3.0)
OUT = os.path.join(os.path.dirname(__file__), "..", "certificates")
LAM = 1.70
A = np.float32(LAM ** -2.0)
B1 = np.float32(LAM ** (ALPHA - 2.0))
B3 = np.float32(LAM ** (ALPHA - 1.0))
k = 20
N = 3 ** (k - 1)
Nl = N // 3
CHUNK = 3 ** 16

v = np.ones(N, dtype=np.float32)
w2 = np.empty(N, dtype=np.float32)
for it in range(200):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    for a in range(0, N, CHUNK):
        b = a + CHUNK
        idx = np.arange(a, b, dtype=np.int64)
        t4 = (4 * idx + 2) % N
        w2[a:b] = A * v[t4]
        r = (idx % 3).astype(np.int8)
        s = idx // 3
        m1 = r == 0
        m3 = r == 2
        w2[a:b][m1] += B1 * cb[(4 * s[m1]) % Nl]
        w2[a:b][m3] += B3 * cb[(2 * s[m3] + 1) % Nl]
        del idx, t4, r, s, m1, m3
    g = float(w2.max())
    np.divide(w2, np.float32(g), out=v)
    if (it + 1) % 25 == 0:
        np.save(os.path.join(OUT, "k20_lam170_200c.npy"), v)
        print(f"  ckpt iter {it+1}, growth {g:.8f}", flush=True)

sum_x = sum_x2 = 0.0
for a in range(0, Nl, CHUNK):
    b = min(a + CHUNK, Nl)
    v1 = np.asarray(v[a:b], dtype=np.float64)
    v2 = np.asarray(v[Nl + a:Nl + b], dtype=np.float64)
    v3 = np.asarray(v[2 * Nl + a:2 * Nl + b], dtype=np.float64)
    lm = np.log2((v1 + v2 + v3) / 3.0)
    for vt in (v1, v2, v3):
        x = np.log2(vt) - lm
        sum_x += float(x.sum())
        sum_x2 += float((x * x).sum())
n = 3 * Nl
var20 = sum_x2 / n - (sum_x / n) ** 2
d19 = var20 / 0.000389
print(f"Var_end(20, 1.70) = {var20:.6f}", flush=True)
print(f"d_19(1.70) = {d19:.4f}", flush=True)
print(f"column: 0.7560/0.7535/0.7590/0.7662/0.7690/0.7719/{d19:.4f}",
      flush=True)
print("done", flush=True)
