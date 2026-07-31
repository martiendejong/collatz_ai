"""
198d_k19.py
===========
Endpoint series, depth k = 19 (387M classes) -- the eighth ratio of the
fork sequence (Obs 397).

Runs at the published edge lam* = 2^0.90934 = 1.87823 (the certified
k=19 record's edge gamma; R30-47) -- no bisection, 250 power iterations,
checkpoints every 50 (k=20 lesson). Memory-lean layout for a 31GB/14GB-
free machine: stored arrays are v, w2, T4 (int32), cb only (~8.8 GB);
branch masks and feed maps are recomputed chunkwise from index
arithmetic. Final growth is reported -- a deviation much beyond ~1e-4
would flag the edge value itself.
"""
import os
import numpy as np
from math import log2

ALPHA = log2(3.0)
OUT = os.path.join(os.path.dirname(__file__), "..", "certificates")

k = 19
N = 3 ** (k - 1)
Nl = N // 3
LAM = 2.0 ** 0.90934
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)
CHUNK = 3 ** 16          # 43M, divides N

i64 = np.arange(N, dtype=np.int64)
T4 = ((4 * i64 + 2) % N).astype(np.int32)
del i64

v = np.ones(N, dtype=np.float64)
w2 = np.empty(N, dtype=np.float64)
g = 1.0
for it in range(250):
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
    g = w2.max()
    np.divide(w2, g, out=v)
    if (it + 1) % 50 == 0:
        np.save(os.path.join(OUT, "k19_eig_198d.npy"), v)
        print(f"  ckpt iter {it+1}, growth {g:.8f}", flush=True)

T = np.stack([v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]])
mean = T.mean(axis=0)
mn = T.min(axis=0)
sig = T.std(axis=0)
X = np.log2(T) - np.log2(mean)[None, :]
var_end = float(np.var(X))
cv = sig / mean
q_k = float(3.0 * mn.sum() / v.sum())
cvw = float(np.sqrt((mean * cv ** 2).sum() / mean.sum()))
lin = float(np.sqrt(2.0) * sig.sum() / mean.sum())
print(f"k=19  lam={LAM:.6f}  growth={g:.8f}  gamma={log2(LAM):.5f}",
      flush=True)
print(f"  Var_count(Xt_end) = {var_end:.6f}", flush=True)
print(f"  1-q = {1-q_k:.6f}  sqrt2*E_w[cv] = {lin:.6f}  "
      f"sqrt2*CV_w = {np.sqrt(2)*cvw:.6f}  margin = {lin/(1-q_k):.3f}",
      flush=True)
print(f"  ratio vs k=18 (0.001651): {var_end/0.001651:.4f}", flush=True)
print("done", flush=True)
