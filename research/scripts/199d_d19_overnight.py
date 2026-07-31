"""
199d_d19_overnight.py
=====================
The single most informative remaining number: d_19 = the pure-depth
contraction at the deepest reachable pair.

    d_19 = Var(20, lam_19) / Var(19, lam_19),
with Var(19, lam_19) = 0.001406 (own edge, script 198d) and the cross
term Var(20, lam_19 = 1.878186) computed here: one k=20 Perron run
(1.16e9 classes) at the frozen k=19 edge.

PREDICTION (Obs 399, falsifiable): the d-increments decay geometrically
(0.0094/0.0082/0.0056/0.0041, ratio ~0.74), extrapolating to
d_19 = 0.818 +- 0.004 and d_inf ~ 0.826 < 1 (the saturation reading of
the fork). d_19 > ~0.825 would weaken the saturation fit; d_19 in band
confirms it at the deepest measurable point.

Engineering: float32 (validated: the saved k20 float32 vector reproduces
the recorded q to 4e-5), no stored index maps (chunkwise arithmetic),
~11 GB peak, checkpoints every 25 iterations, 175 iterations (the k=20
convergence lesson: 60 suffice at 50; margin taken).
"""
import os
import numpy as np
from math import log2

ALPHA = log2(3.0)
OUT = os.path.join(os.path.dirname(__file__), "..", "certificates")

k = 20
N = 3 ** (k - 1)
Nl = N // 3
LAM = 1.878186          # frozen k=19 edge
A = np.float32(LAM ** -2.0)
B1 = np.float32(LAM ** (ALPHA - 2.0))
B3 = np.float32(LAM ** (ALPHA - 1.0))
CHUNK = 3 ** 16

v = np.ones(N, dtype=np.float32)
w2 = np.empty(N, dtype=np.float32)
g = 1.0
for it in range(175):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    for a in range(0, N, CHUNK):
        b = a + CHUNK
        idx = np.arange(a, b, dtype=np.int64)
        t4 = ((4 * idx + 2) % N)
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
        np.save(os.path.join(OUT, "k20_lam19_199d.npy"), v)
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
var_cross = sum_x2 / n - (sum_x / n) ** 2
d19 = var_cross / 0.001406
l19 = 0.001205 / var_cross
print(f"Var(20, lam_19) = {var_cross:.6f}  growth {g:.8f}", flush=True)
print(f"d_19 = {d19:.4f}   l_19 = {l19:.4f}   r_19(own) = "
      f"{0.001205/0.001406:.4f}", flush=True)
print(f"d-series: 0.7878/0.7873/0.7967/0.8049/0.8105/0.8146/{d19:.4f}",
      flush=True)
print("done", flush=True)
