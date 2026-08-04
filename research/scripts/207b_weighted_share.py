"""
207b_weighted_share.py
======================
Obs 412 -- refinement of the common-cause test (Obs 411).

In 207_common_cause.py the proxy for the shared cause was the raw D3-count
along the type walk.  Under the multiplicative K-L model the Perron
eigenvector value at a base triple is the product of the branch weights
encountered along the ancestral type walk.  The TRUE shared cause in log
space is therefore:

    Z_log = sum_{g=1}^{G} log(b_{step_g})   where b in {B1, B3}

The raw D3-count is a coarse 0/1 binarization of this (log(B3) vs log(B1)).
The log-weight version uses the actual magnitudes and hence should absorb more
of the common-cause structure.

Prediction: partial_corr(logL, logS | Z_log) << 0.66 (the 207 residual);
ideally closer to zero, confirming the multiplicative common-cause skeleton.

Setup identical to 207 except the walk accumulates log-weight instead of
a D3 count.  Systems 15..17, lam=1.70, G=8.
"""
import numpy as np
from math import log2, log

ALPHA = log2(3.0)
LAM   = 1.70
A     = LAM ** -2.0
B1    = LAM ** (ALPHA - 2.0)
B3    = LAM ** (ALPHA - 1.0)
LOG_B1 = log(B1)   # negative  (~-0.220)
LOG_B3 = log(B3)   # positive  (~+0.310)
G = 8


def perron_and_logweight(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl   = N // 3
    R1   = (4 * s) % Nl
    R3   = (2 * s + 1) % Nl
    m1, m3 = (r == 0), (r == 2)

    # Perron eigenvector
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb   = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w2   = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v    = w2 / w2.max()
    v = v / v.mean()

    # Ancestral log-weight along the type walk, G steps.
    # Each step descends one level (Freshness Lemma: modulus /= 3).
    pos      = np.arange(Nl, dtype=np.int64)
    logw     = np.zeros(Nl, dtype=np.float64)
    alive    = np.ones(Nl, dtype=bool)
    mod      = Nl
    for _ in range(G):
        t    = pos % 3
        d1   = alive & (t == 0)
        d3   = alive & (t == 2)
        logw[d1] += LOG_B1
        logw[d3] += LOG_B3
        alive = alive & (t != 1)          # D2 absorbs
        sp   = pos // 3
        mod //= 3
        nxt  = np.where(t == 0, (4 * sp) % mod, (2 * sp + 1) % mod)
        pos  = np.where(alive, nxt, pos)

    return v, logw


print(f"weighted-share common-cause test  lam={LAM}  G={G}", flush=True)
print("  sys  corr(L,S)  partial|logW   partial_raw   "
      "corr(L,logW)  corr(S,logW)", flush=True)

PARTIAL_RAW = {15: 0.65995, 16: 0.65851, 17: 0.65776}   # from 207

for kp1 in (15, 16, 17):
    v, logw = perron_and_logweight(kp1)
    Nl  = v.size // 3
    mn  = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:]) / 3.0
    Nl3 = Nl // 3
    M   = np.stack([mn[:Nl3], mn[Nl3:2*Nl3], mn[2*Nl3:]])
    level  = M.mean(axis=0)
    spread = M.std(axis=0)
    W      = np.stack([logw[:Nl3], logw[Nl3:2*Nl3], logw[2*Nl3:]]).mean(axis=0)

    ok  = spread > 0
    x   = np.log(level[ok]);  x -= x.mean()
    y   = np.log(spread[ok]); y -= y.mean()
    z   = W[ok].astype(np.float64); z -= z.mean()

    def c(a, b):
        return float(np.mean(a*b) / np.sqrt(np.mean(a*a)*np.mean(b*b)))

    cxy, cxz, cyz = c(x, y), c(x, z), c(y, z)
    partial = (cxy - cxz*cyz) / np.sqrt((1 - cxz**2)*(1 - cyz**2))

    print(f"  {kp1:2d}  {cxy:.5f}    {partial:.5f}        "
          f"{PARTIAL_RAW[kp1]:.5f}        {cxz:.5f}       {cyz:.5f}",
          flush=True)

print("done", flush=True)
