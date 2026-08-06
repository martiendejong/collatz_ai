"""
292_slack_k15_17.py
===================
Slack ratio s2/s0 (Obs 490) at lambda=1.05 for k=15,16,17 (float64, RAM).
Trend matters: (3b) tail needs lim sup s2/s0 < lambda = 1.05.
Through k=14: s2/s0 rises 0.55 -> 0.74. Where does it head?
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 1.05
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
ITERS = {15: 150, 16: 100, 17: 80}

for k in [15, 16, 17]:
    N = 3**(k-1); Nl = N//3; Nl3 = Nl//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr == 0; m2 = r_arr == 2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    del i, s_arr, r_arr

    v = np.ones(N)
    rho = 1.0
    for _ in range(ITERS[k]):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max())
        w /= rho
        v = w
    del T4, R1, R3, m0, m2

    t = A/rho
    R = (t**2+lam)/(1+t*lam)
    W = lam*(1-t**3)/(1+t*lam)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    del v

    m = np.arange(Nl3)
    e = np.array([0, 1, 2])
    s_idx = m[:, None] + e[None, :]*Nl3

    def cols(vec, idx):
        return np.stack([vec[idx[:, 0]], vec[idx[:, 1]], vec[idx[:, 2]]], axis=1)

    def gap(c):
        return c.mean(axis=1) - c.min(axis=1)

    tgt0 = (4*s_idx) % Nl
    S0 = (t*gap(cols(v2, tgt0)) + t*lam**ALPHA*gap(cols(cb, tgt0))
          - gap(cols(v0, s_idx)))
    tgt2p = (4*s_idx+3) % Nl
    tgt2c = (2*s_idx+1) % Nl
    S2 = (t*gap(cols(v1, tgt2p)) + t*lam**(ALPHA+1)*gap(cols(cb, tgt2c))
          - gap(cols(v2, s_idx)))
    s0 = float(S0.mean()); s2 = float(S2.mean())
    print(f"k={k}: s2/s0={s2/s0:.5f} margin={lam-s2/s0:+.5f} rho={rho:.6f}", flush=True)
    del cb, v0, v1, v2, S0, S2, s_idx, tgt0, tgt2p, tgt2c
print("DONE")
