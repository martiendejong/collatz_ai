"""
295_lambda_edge.py
==================
The lambda edge: the density-one chain uses lambda in (1, 2], but the verification
grid starts at 1.05. Check lambda = 1.01, 1.02, 1.03 (k up to 12):
  - exact criterion c2/c0 < R and margin behavior as lambda -> 1+
  - slack ratio s2/s0 vs lambda (required < lambda, so the cushion shrinks to 1!)
  - identity mu2/mu0 = R sanity
At lambda=1 the operator may become degenerate/symmetric; the margin may vanish.
This probes whether a second unbounded direction (lambda -> 1+) hides in the chain.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
ITERS = {6: 3000, 8: 2500, 10: 1500, 12: 800}

for lam in [1.01, 1.02, 1.03, 1.05]:
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    for k in [6, 8, 10, 12]:
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
            w_ = A*v[T4]
            w_[m2] += B3*cb[R3[m2]]
            w_[m0] += B1*cb[R1[m0]]
            rho = float(w_.max()); w_ /= rho; v = w_
        del T4, R1, R3, m0, m2
        t = A/rho
        R = (t**2+lam)/(1+t*lam)
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
        mu0 = v0.mean(); mu2 = v2.mean()
        c0 = cb[0::3].mean(); c2 = cb[2::3].mean()
        g0 = mu0-c0; g2 = mu2-c2

        m = np.arange(Nl3); e = np.array([0,1,2])
        s_idx = m[:,None] + e[None,:]*Nl3
        def cols(vec, idx): return np.stack([vec[idx[:,0]], vec[idx[:,1]], vec[idx[:,2]]], axis=1)
        def gap(c): return c.mean(axis=1) - c.min(axis=1)
        w_cb = t*lam**ALPHA
        tgt0 = (4*s_idx) % Nl
        S0 = t*gap(cols(v2, tgt0)) + w_cb*gap(cols(cb, tgt0)) - gap(cols(v0, s_idx))
        S2 = (t*gap(cols(v1, (4*s_idx+3) % Nl)) + lam*w_cb*gap(cols(cb, (2*s_idx+1) % Nl))
              - gap(cols(v2, s_idx)))
        s0 = float(S0.mean()); s2 = float(S2.mean())

        print(f"lam={lam:.2f} k={k:2d}: R={R:.6f} c2/c0={c2/c0:.6f} margin={R-c2/c0:+.3e} "
              f"g2/(R*g0)={g2/(R*g0):.5f} s2/s0={s2/s0:.5f} slack_margin={lam-s2/s0:+.5f} "
              f"id={mu2/mu0-R:+.1e}", flush=True)
print("DONE")
