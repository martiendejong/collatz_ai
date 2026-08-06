"""
288_large_k_all_lam.py
======================
Verify ratio > 1 for k=15..18 at lambda=1.10, 1.20, 1.30, 1.40.
(lambda=1.05 already verified in Scripts 282/284.)
(For lambda>=1.50, R > 1 so the condition c2/c0 < R is trivially satisfied.)

Uses float32 throughout for speed. 40 iterations each.
Output: ratio = E[sw2]/mu2 / (E[sw0]/mu0), min should be > 1.
"""
import numpy as np
from math import log2, sqrt

ALPHA = log2(3.0)
LAMBDAS = [1.10, 1.20, 1.30, 1.40]
KS = [15, 16, 17, 18]
NITERS = 40

results = []
for k in KS:
    N = 3**(k-1); Nl = N//3; Nl3 = Nl//3
    print(f"\n--- k={k}, N={N:,} ---", flush=True)

    for lam in LAMBDAS:
        A  = np.float32(lam**-2)
        B1 = np.float32(lam**(ALPHA-2))
        B3 = np.float32(lam**(ALPHA-1))

        # Precompute index arrays
        i = np.arange(N, dtype=np.int64)
        T4 = ((4*i+2) % N).astype(np.int32)
        s_arr, r_arr = np.divmod(i, 3)
        m0 = (r_arr == 0)
        m2 = (r_arr == 2)
        R1 = ((4*s_arr) % Nl).astype(np.int32)
        R3 = ((2*s_arr+1) % Nl).astype(np.int32)
        del i, s_arr, r_arr

        v = np.ones(N, dtype=np.float32)
        rho = np.float32(1.0)
        for it in range(NITERS):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A * v[T4]
            w[m2] += B3 * cb[R3[m2]]
            w[m0] += B1 * cb[R1[m0]]
            wmax = float(w.max())
            rho = np.float32(wmax)
            w /= np.float32(wmax)
            v = w
            if (it+1) % 10 == 0:
                print(f"  lam={lam:.2f} k={k} iter {it+1}/{NITERS} wmax={wmax:.6f}", flush=True)

        v0 = v[0::3]; v2 = v[2::3]
        mu0 = float(v0.mean()); mu2 = float(v2.mean())
        j3 = np.arange(Nl3, dtype=np.int64)
        c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1).astype(np.float64)
        c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1).astype(np.float64)

        sw0 = c0.std(axis=1); sw2 = c2.std(axis=1)
        E_s0 = sw0.mean() / mu0
        E_s2 = sw2.mean() / mu2
        ratio = E_s2 / E_s0

        t = float(A) / float(rho)
        R_val = (t**2 + lam) / (1 + t*lam)
        sqFR = sqrt((t**4 + lam**2)/(1 + t**2*lam**2)) / R_val

        results.append((lam, k, ratio, sqFR, E_s0, E_s2))
        print(f"  lam={lam:.2f} k={k}: ratio={ratio:.6f} sqFR={sqFR:.5f} E_s0={E_s0:.7f} E_s2={E_s2:.7f}", flush=True)

        del v, w, c0, c2, sw0, sw2, j3, T4, R1, R3, m0, m2

print("\n=== SUMMARY ===")
for lam, k, ratio, sqFR, E_s0, E_s2 in results:
    status = "OK" if ratio > 1.0 else "FAIL!!"
    print(f"lam={lam:.2f} k={k:2d}: ratio={ratio:.6f} sqFR={sqFR:.5f} {status}")
print(f"min ratio = {min(r[2] for r in results):.6f}")
print("DONE")
