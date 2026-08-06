"""
291_slack_all_cases.py
======================
Measure the slack ratio s2/s0 (Obs 490) for all 144 cases (k=3..14 x 12 lambdas).
Step (3b) is EXACTLY equivalent to lam*s0 > s2, i.e. s2/s0 < lam.

Also verifies the exact identity (1-t^3)(g2 - R*g0) = W*(s0 - s2/lam) per case.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAMBDAS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
KS = list(range(3, 15))
ITERS = {3:2000, 4:2000, 5:2000, 6:2000, 7:2000, 8:1500, 9:1500, 10:1000,
         11:1000, 12:500, 13:300, 14:200}

results = []
for lam in LAMBDAS:
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    for k in KS:
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

        t = A/rho
        R = (t**2+lam)/(1+t*lam)
        W = lam*(1-t**3)/(1+t*lam)
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

        m = np.arange(Nl3)
        e = np.array([0, 1, 2])
        s_idx = m[:, None] + e[None, :]*Nl3

        def cols(vec, idx):
            return np.stack([vec[idx[:, 0]], vec[idx[:, 1]], vec[idx[:, 2]]], axis=1)

        def gap(c):
            return c.mean(axis=1) - c.min(axis=1)

        tgt0 = (4*s_idx) % Nl
        col0 = cols(v0, s_idx)
        pass0 = cols(v2, tgt0)
        cb0 = cols(cb, tgt0)
        S0 = t*gap(pass0) + t*lam**ALPHA*gap(cb0) - gap(col0)

        tgt2p = (4*s_idx+3) % Nl
        tgt2c = (2*s_idx+1) % Nl
        col2 = cols(v2, s_idx)
        pass2 = cols(v1, tgt2p)
        cb2 = cols(cb, tgt2c)
        S2 = t*gap(pass2) + t*lam**(ALPHA+1)*gap(cb2) - gap(col2)

        s0 = float(S0.mean()); s2 = float(S2.mean())
        g0 = float(gap(col0).mean()); g2 = float(gap(col2).mean())
        lhs = (1-t**3)*(g2 - R*g0)
        rhs = W*(s0 - s2/lam)
        idratio = lhs/rhs if rhs != 0 else float("nan")
        ok = "OK" if lam*s0 > s2 else "FAIL!!"
        results.append((lam, k, s2/s0, lam - s2/s0, idratio))
        print(f"lam={lam:.2f} k={k:2d}: s2/s0={s2/s0:.5f} margin(lam-s2/s0)={lam-s2/s0:+.5f} "
              f"id={idratio:.8f} {ok}", flush=True)
        del v, w, cb, col0, col2, pass0, pass2, cb0, cb2, S0, S2, T4, R1, R3, m0, m2

print("\n=== SUMMARY ===")
worst = max(results, key=lambda r: r[2])
print(f"max s2/s0 = {worst[2]:.5f} at lam={worst[0]:.2f} k={worst[1]}")
print(f"min margin (lam - s2/s0) = {min(r[3] for r in results):+.5f}")
bad_id = [r for r in results if abs(r[4]-1) > 1e-6]
print(f"identity violations (>1e-6): {len(bad_id)}")
print("DONE")
