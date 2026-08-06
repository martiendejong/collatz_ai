"""
287_corrLtoL_fixed.py
=====================
FIXED version of Script 286: correctly tracks rho (Perron eigenvalue).
Bug in 286: rho = float(w.max()) was evaluated AFTER w /= w.max(), giving rho=1.0 always.
Fix: save wmax = w.max() before normalization, set rho = wmax.

Computes for all 144 cases (k=3..14, lambda=1.05..2.00):
  ratio_L1  = E[sw2]/mu2 / (E[sw0]/mu0)  -- same as Script 281/286, should be > 1
  ratio_L2  = sqrt(E[sw2^2])/mu2 / sqrt(E[sw0^2]/mu0)  -- approx sqrt(F)/R
  corrLtoL  = ratio_L1 / ratio_L2  -- should be > 1/ratio_L2
  sqFR      = sqrt(F)/R with CORRECT t = A/rho
  B_bound   = 1/sqFR  -- CORRECT threshold for corrLtoL
  margin    = corrLtoL - B_bound (positive means OK)
"""
import numpy as np
from math import log2, sqrt

ALPHA = log2(3.0)
LAMBDAS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
KS = list(range(3, 15))
ITERS = {3:2000, 4:2000, 5:2000, 6:2000, 7:2000, 8:1500, 9:1500, 10:1000,
         11:1000, 12:500, 13:300, 14:200}

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
        rho = 1.0
        for _ in range(niters):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A * v[T4]
            w[m2] += B3 * cb[R3[m2]]
            w[m0] += B1 * cb[R1[m0]]
            wmax = float(w.max())
            rho = wmax          # FIX: capture before normalization
            w /= wmax
            v = w

        # Correct t and derived quantities
        t = A / rho
        R_val = (t**2 + lam) / (1 + t * lam)
        F_val = (t**4 + lam**2) / (1 + t**2 * lam**2)
        sqFR = sqrt(F_val) / R_val
        B_bound = 1.0 / sqFR   # = R/sqrt(F), threshold for corrLtoL

        v0 = v[0::3]; v2 = v[2::3]
        mu0 = v0.mean(); mu2 = v2.mean()
        j3 = np.arange(Nl3, dtype=np.int64)
        c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
        c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)

        sw0 = c0.std(axis=1); sw2 = c2.std(axis=1)
        E_L1_0 = sw0.mean() / mu0
        E_L1_2 = sw2.mean() / mu2
        ratio_L1 = E_L1_2 / E_L1_0

        E_L2_0 = sqrt(float(np.mean(sw0**2))) / mu0
        E_L2_2 = sqrt(float(np.mean(sw2**2))) / mu2
        ratio_L2 = E_L2_2 / E_L2_0

        f0 = E_L1_0 / E_L2_0   # L1/L2 ratio for r=0
        f2 = E_L1_2 / E_L2_2   # L1/L2 ratio for r=2
        corrLtoL = f2 / f0

        margin = corrLtoL - B_bound
        status = "OK" if corrLtoL > B_bound else "FAIL"
        results.append((lam, k, ratio_L1, sqFR, corrLtoL, B_bound, margin, status, rho, t))
        print(f"lam={lam:.2f} k={k:2d}: ratio={ratio_L1:.5f} sqFR={sqFR:.5f} "
              f"corrLtoL={corrLtoL:.5f} B={B_bound:.5f} margin={margin:+.6f} {status} "
              f"rho={rho:.5f} t={t:.5f}", flush=True)
        del v, w, c0, c2, sw0, sw2, j3, T4, R1, R3, m0, m2

fails = [r for r in results if r[7] == "FAIL"]
print(f"\n=== SUMMARY ===")
print(f"Total: {len(results)} cases, {len(fails)} FAIL")
print(f"min ratio_L1 = {min(r[2] for r in results):.5f}")
print(f"min sqFR     = {min(r[3] for r in results):.5f}")
print(f"min corrLtoL = {min(r[4] for r in results):.5f}")
print(f"min B_bound  = {min(r[5] for r in results):.5f}")
print(f"min margin   = {min(r[6] for r in results):+.6f}")
if fails:
    print("FAILS (corrLtoL < B_bound):")
    for r in fails:
        print(f"  lam={r[0]:.2f} k={r[1]:2d}: corrLtoL={r[4]:.5f} B={r[5]:.5f} margin={r[6]:+.6f}")
else:
    print("All 144 cases: corrLtoL > B_bound (PASS)")
print("DONE")
