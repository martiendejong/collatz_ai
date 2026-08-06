"""
286_corrLtoL_bound.py
=====================
Compute corrLtoL = (E[sig_w2]/mu2 / sqrt(E[sig_w2^2]/mu2^2)) /
                   (E[sig_w0]/mu0 / sqrt(E[sig_w0^2]/mu0^2))
for all 144 cases (k=3..14, 12 lambda values).

Also compute bound B = R / sqrt(F) (the threshold below which ratio < 1).
Show corrLtoL > B for all cases (i.e., the L1/L2 correction does not undermine ratio > 1).

Uses the pre-computed eigenvectors from Script 281 (re-runs the iteration).
Fast: 144 cases, reuses code from Script 281.
"""
import numpy as np
from math import log2, sqrt
import sys

ALPHA = log2(3.0)
LAMBDAS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
KS = list(range(3, 15))
ITERS = {3:2000,4:2000,5:2000,6:2000,7:2000,8:1500,9:1500,10:1000,11:1000,12:500,13:300,14:200}

results = []
for lam in LAMBDAS:
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    for k in KS:
        niters = ITERS[k]
        N = 3**(k-1); Nl = N//3; Nl3 = Nl//3
        i = np.arange(N, dtype=np.int64)
        T4 = (4*i+2) % N
        s_arr, r_arr = np.divmod(i, 3)
        m0 = r_arr == 0; m2 = r_arr == 2
        R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
        del i, s_arr, r_arr

        v = np.ones(N)
        for _ in range(niters):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A*v[T4]
            w[m2] += B3*cb[R3[m2]]
            w[m0] += B1*cb[R1[m0]]
            w /= w.max()
            v = w

        rho = float(w.max()); t = A/rho
        R = (t**2+lam)/(1+t*lam)
        F = (t**4+lam**2)/(1+t**2*lam**2)
        sqFR = sqrt(F)/R
        B_bound = R/sqrt(F)  # = 1/sqFR, threshold for corrLtoL

        v0 = v[0::3]; v2 = v[2::3]
        mu0 = v0.mean(); mu2 = v2.mean()
        j3 = np.arange(Nl3, dtype=np.int64)
        c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
        c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

        # L1 components (mean of std)
        sw0 = c0.std(axis=1); sw2 = c2.std(axis=1)
        E_L1_0 = sw0.mean() / mu0; E_L1_2 = sw2.mean() / mu2
        ratio_L1 = E_L1_2 / E_L1_0

        # L2 components (sqrt of mean-of-sq)
        E_L2_0 = sqrt(np.mean(sw0**2)) / mu0
        E_L2_2 = sqrt(np.mean(sw2**2)) / mu2
        ratio_L2 = E_L2_2 / E_L2_0  # should be close to sqFR

        # corrLtoL factor
        f0 = E_L1_0 / E_L2_0  # L1/L2 for r=0
        f2 = E_L1_2 / E_L2_2  # L1/L2 for r=2
        corrLtoL = f2 / f0

        margin = corrLtoL - B_bound
        status = "OK" if corrLtoL > B_bound else "FAIL"
        results.append((lam, k, ratio_L1, sqFR, corrLtoL, B_bound, margin, status))
        print(f"lam={lam:.2f} k={k:2d}: ratio={ratio_L1:.5f} sqFR={sqFR:.5f} "
              f"corrLtoL={corrLtoL:.5f} B={B_bound:.5f} margin={margin:.5f} {status}", flush=True)
        del v, w, c0, c2, sw0, sw2, j3, T4, R1, R3, m0, m2

fails = [r for r in results if r[7]=="FAIL"]
print(f"\n=== SUMMARY ===")
print(f"Total: {len(results)} cases, {len(fails)} FAIL")
print(f"min corrLtoL = {min(r[4] for r in results):.5f}")
print(f"min B_bound  = {min(r[5] for r in results):.5f}")
print(f"min margin   = {min(r[6] for r in results):.5f}")
print(f"min ratio_L1 = {min(r[2] for r in results):.5f}")
print("DONE")
