"""
198c_k18.py
===========
The deciding measurement for the endpoint-rate creep question (Obs 396):
extend the endpoint series to k = 18 (129M classes, in-memory).

The series Var_count(Xt_end): ratios 0.829/0.833/0.839/0.836/0.843/0.850
at k = 11..17 -- upward creep +0.005/depth. Saturation below 1 vs slow
drift toward 1 is undecidable on six points; k=18 adds the seventh ratio.
If the next ratio drops back toward ~0.84, saturation; if it continues
to ~0.855+, drift stays on the table and Endpoint Decay is in danger.

Method: cold-start bisection (the 198-run-1 lesson), 14 steps x 60 iters
in [1.863, 1.877] (lam* interpolation ~1.8695), then 200 polish
iterations with checkpoints every 50 (the k=20 lesson). int32 index
maps to halve memory (~5 GB peak).
"""
import os
import numpy as np
from math import log2

ALPHA = log2(3.0)
OUT = os.path.join(os.path.dirname(__file__), "..", "certificates")

k = 18
N = 3 ** (k - 1)
i = np.arange(N, dtype=np.int64)
T4 = ((4 * i + 2) % N).astype(np.int32)
s = (i // 3).astype(np.int32)
r = (i % 3).astype(np.int8)
del i
Nl = N // 3
R1 = ((4 * s.astype(np.int64)) % Nl).astype(np.int32)
R3 = ((2 * s.astype(np.int64) + 1) % Nl).astype(np.int32)
del s
m1, m3 = (r == 0), (r == 2)
del r
R1m, R3m = R1[m1], R3[m3]
del R1, R3

def sweep(v, lam, n_iter, ckpt=None):
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    g = 1.0
    for it in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1m]
        w2[m3] += B3 * cb[R3m]
        g = w2.max()
        v = w2 / g
        if ckpt and (it + 1) % 50 == 0:
            np.save(ckpt, v)
            print(f"    ckpt iter {it+1}, growth {g:.8f}", flush=True)
    return v, g

lo_l, hi_l = 1.863, 1.877
v = np.ones(N, dtype=np.float64)
for step in range(14):
    lam = 0.5 * (lo_l + hi_l)
    w, g = sweep(np.ones(N, dtype=np.float64), lam, 60)
    print(f"  step {step}: lam={lam:.6f} growth={g:.8f}", flush=True)
    if g >= 1.0:
        lo_l, v = lam, w
    else:
        hi_l = lam

lam18 = lo_l
v, g = sweep(v, lam18, 200, ckpt=os.path.join(OUT, "k18_eig_198c.npy"))
np.save(os.path.join(OUT, "k18_eig_198c.npy"), v)

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
print(f"k=18  lam={lam18:.6f}  growth={g:.8f}  gamma={log2(lam18):.5f}",
      flush=True)
print(f"  Var_count(Xt_end) = {var_end:.6f}", flush=True)
print(f"  1-q = {1-q_k:.6f}  sqrt2*E_w[cv] = {lin:.6f}  "
      f"sqrt2*CV_w = {np.sqrt(2)*cvw:.6f}  margin = {lin/(1-q_k):.3f}",
      flush=True)
print(f"  ratio vs k=17 (0.001932): {var_end/0.001932:.4f}", flush=True)
print("done", flush=True)
