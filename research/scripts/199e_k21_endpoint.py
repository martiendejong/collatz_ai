"""
199e_k21_endpoint.py
====================
The tenth point of the endpoint series -- k = 21 -- computed from the
SAVED exact-integer certificate of the record run (25-jul,
research/k21/cert_k21.npy, int64, 3.49e9 classes, the object behind
pi(x) >= x^0.9184).

Values are ~1e7 so integer rounding is a ~1e-7 relative perturbation --
negligible for profile statistics. Caveat recorded: this is the
certified feasible solution, not the fully converged Perron vector; at
k=20 the analogous difference (polished vs certificate) was below the
reported precision.

Outputs: Var_count(Xt_end)(21), 1-q(21) (arbitrates the min-loss curve
value 0.02647 at lam = 1.89015), CV_w, Samuelson margin, and
r_20 = Var(21)/Var(20) -- the tenth ratio of the fork series.
"""
import os
import numpy as np

PATH = os.path.join(os.path.dirname(__file__), "..", "k21",
                    "cert_k21.npy")
CHUNK = 3 ** 15

v = np.lib.format.open_memmap(PATH, mode="r")
N = v.shape[0]
Nl = N // 3
print(f"k=21 certificate loaded: N={N}  dtype={v.dtype}", flush=True)

sum_x = sum_x2 = 0.0
sum_min = sum_all = 0.0
sum_mean = sum_sig = sum_mcv2 = 0.0
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
    sum_min += float(mn.sum())
    sum_all += float(v1.sum() + v2.sum() + v3.sum())
    sum_mean += float(mean.sum())
    sum_sig += float(np.sqrt(sig2).sum())
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
print(f"  min-loss curve at lam=1.89015 predicts 1-q = 0.02647",
      flush=True)
print(f"  r_20 = Var(21)/Var(20) = {var_end/0.001205:.4f}", flush=True)
print("done", flush=True)
