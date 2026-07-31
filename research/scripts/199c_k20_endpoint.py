"""
199c_k20_endpoint.py
====================
The ninth point of the endpoint series -- k = 20 -- computed from the
SAVED polished eigenvector of the certified record run (15-jul,
certificates/k20_polished.npy, lam = 1.885, float32, 1.16e9 classes).
No new eigenvector computation needed.

This directly tests the Obs 399 prediction from the d/l decomposition:
r_k should DECLINE from 0.852 toward ~0.83 (saturation reading). A
value near 0.85 or above would push back toward the drift reading.
"""
import os
import numpy as np
from math import log2

PATH = os.path.join(os.path.dirname(__file__), "..", "certificates",
                    "k20_polished.npy")
LAM = 1.885
CHUNK = 3 ** 15          # ~14M

v = np.load(PATH, mmap_mode="r")
N = v.shape[0]
Nl = N // 3
print(f"k=20 loaded: N={N}  dtype={v.dtype}  lam={LAM}", flush=True)

sum_x = 0.0
sum_x2 = 0.0
sum_min = 0.0
sum_all = 0.0
sum_mean = 0.0
sum_sig = 0.0
sum_mcv2 = 0.0
for a in range(0, Nl, CHUNK):
    b = min(a + CHUNK, Nl)
    v1 = np.asarray(v[a:b], dtype=np.float64)
    v2 = np.asarray(v[Nl + a:Nl + b], dtype=np.float64)
    v3 = np.asarray(v[2 * Nl + a:2 * Nl + b], dtype=np.float64)
    mean = (v1 + v2 + v3) / 3.0
    mn = np.minimum(np.minimum(v1, v2), v3)
    lm = np.log2(mean)
    for vt in (v1, v2, v3):
        x = np.log2(vt) - lm
        sum_x += float(x.sum())
        sum_x2 += float((x * x).sum())
    sig2 = ((v1 - mean) ** 2 + (v2 - mean) ** 2 + (v3 - mean) ** 2) / 3.0
    sig = np.sqrt(sig2)
    sum_min += float(mn.sum())
    sum_all += float(v1.sum() + v2.sum() + v3.sum())
    sum_mean += float(mean.sum())
    sum_sig += float(sig.sum())
    sum_mcv2 += float((sig2 / mean).sum())

n = 3 * Nl
var_end = sum_x2 / n - (sum_x / n) ** 2
q_k = 3.0 * sum_min / sum_all
cvw = (sum_mcv2 / sum_mean) ** 0.5
lin = 2.0 ** 0.5 * sum_sig / sum_mean
print(f"  Var_count(Xt_end) = {var_end:.6f}", flush=True)
print(f"  1-q = {1-q_k:.6f}  sqrt2*E_w[cv] = {lin:.6f}  "
      f"sqrt2*CV_w = {2.0**0.5*cvw:.6f}  margin = {lin/(1-q_k):.3f}",
      flush=True)
print(f"  ratio vs k=19 (0.001406): {var_end/0.001406:.4f}", flush=True)
print("done", flush=True)
