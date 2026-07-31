"""
198b_k16_redo.py
================
Redo of the k=16 endpoint point of script 198. The warm-started
bisection there carried the eigenvector across lambda trials with only
50 iterations per trial; the growth estimate went stale and the bracket
turned the wrong way at step 6, landing at lam = 1.8465 with converged
growth 1.0018 (i.e. visibly below the true edge; interpolation of the
lam* series puts lam*_16 near 1.852). Lesson logged: bisection trials
must COLD-START, exactly as the original edge_vector() does.

This script: cold-start bisection, 18 steps x 60 iters in the narrowed
bracket [1.848, 1.856], then 250 polish iterations, then the endpoint
statistics.
"""
import os
import numpy as np
from math import log2

ALPHA = log2(3.0)
OUT = os.path.join(os.path.dirname(__file__), "..", "certificates")

k = 16
N = 3 ** (k - 1)
i = np.arange(N, dtype=np.int64)
T4 = (4 * i + 2) % N
s, r = np.divmod(i, 3)
Nl = N // 3
R1 = (4 * s) % Nl
R3 = (2 * s + 1) % Nl
m1, m3 = (r == 0), (r == 2)

lo_l, hi_l = 1.848, 1.856
v = np.ones(N, dtype=np.float64)
for step in range(18):
    lam = 0.5 * (lo_l + hi_l)
    A, B1, B3 = lam ** -2.0, lam ** (ALPHA - 2.0), lam ** (ALPHA - 1.0)
    w = np.ones(N, dtype=np.float64)          # COLD start per trial
    for _ in range(60):
        cb = np.minimum(np.minimum(w[:Nl], w[Nl:2 * Nl]), w[2 * Nl:])
        w2 = A * w[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        g = w2.max()
        w = w2 / g
    print(f"  step {step}: lam={lam:.6f} growth={g:.8f}", flush=True)
    if g >= 1.0:
        lo_l, v = lam, w
    else:
        hi_l = lam

lam16 = lo_l
A, B1, B3 = lam16 ** -2.0, lam16 ** (ALPHA - 2.0), lam16 ** (ALPHA - 1.0)
for it in range(250):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
    w2 = A * v[T4]
    w2[m1] += B1 * cb[R1[m1]]
    w2[m3] += B3 * cb[R3[m3]]
    g = w2.max()
    v = w2 / g
    if (it + 1) % 50 == 0:
        np.save(os.path.join(OUT, "k16_eig_198b.npy"), v)
        print(f"    polish iter {it+1}, growth {g:.8f}", flush=True)

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
print(f"k=16  lam={lam16:.6f}  growth={g:.8f}", flush=True)
print(f"  Var_count(Xt_end) = {var_end:.6f}", flush=True)
print(f"  1-q = {1-q_k:.6f}  sqrt2*E_w[cv] = {lin:.6f}  "
      f"sqrt2*CV_w = {np.sqrt(2)*cvw:.6f}  margin = {lin/(1-q_k):.3f}",
      flush=True)
print(f"  ratio vs k=15 (0.002696): {var_end/0.002696:.4f}", flush=True)
print(f"  ratio k=17 (0.001932) vs this: {0.001932/var_end:.4f}", flush=True)
print("done", flush=True)
