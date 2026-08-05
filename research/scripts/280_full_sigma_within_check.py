"""
280_full_sigma_within_check.py
==============================
Full 168-case verification of the KEY CONDITION for step (3b):
  E[sigma_within(v2)]/mean_v2 > E[sigma_within(v0)]/mean_v0

For all k in {3..14} and lambda in {1.05, 1.10, ..., 2.00}.
168 = 12 lambda values x 14 k values.

Also computes:
- E[sigma^2_within(v2)]/E[sigma^2_within(v0)]: L2 ratio (should be > F > R^2)
- F = (t^4+lam^2)/(1+t^2*lam^2): Obs 471 formula
- Minimum ratio across all cases.
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)

LAMBDAS = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 1.95, 2.00]
KS = list(range(3, 15))

def run_kl(k, lam, n_iter=2000):
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    rho = float(w.max()); t = A/rho; R = (t**2+lam)/(1+t*lam)
    F = (t**4+lam**2)/(1+t**2*lam**2)
    return v, Nl, t, R, F

results = []
n_fail = 0
n_total = 0
min_ratio = float('inf')
min_case = None

print(f"{'lam':>5} {'k':>3} {'E_s0/mu0':>10} {'E_s2/mu2':>10} {'ratio':>7} {'R':>7} {'F':>7} {'L2ratio':>8} {'OK':>4}")
sys.stdout.flush()

for lam in LAMBDAS:
    for k in KS:
        v, Nl, t, R, F = run_kl(k, lam)
        Nl3 = Nl//3
        v0 = v[0::3]; v2 = v[2::3]
        j3 = np.arange(Nl3)
        c0 = np.stack([v0[j3], v0[j3+Nl3], v0[j3+2*Nl3]], axis=1)
        c2 = np.stack([v2[j3], v2[j3+Nl3], v2[j3+2*Nl3]], axis=1)
        sig_w0 = c0.std(1).mean()
        sig_w2 = c2.std(1).mean()
        sig2_w0 = c0.var(1).mean()
        sig2_w2 = c2.var(1).mean()
        mu0 = v0.mean(); mu2 = v2.mean()
        E0 = sig_w0/mu0; E2 = sig_w2/mu2
        ratio = E2/E0
        L2ratio = sig2_w2/sig2_w0 if sig2_w0 > 0 else 999.0
        ok = 'OK' if ratio > 1 else 'FAIL'
        if ratio <= 1:
            n_fail += 1
        n_total += 1
        if ratio < min_ratio:
            min_ratio = ratio
            min_case = (lam, k, E0, E2, ratio, R, F, L2ratio)
        print(f"{lam:>5.2f} {k:>3}  {E0:>10.7f}  {E2:>10.7f}  {ratio:>7.5f}  {R:>7.5f}  {F:>7.4f}  {L2ratio:>8.4f}  {ok:>4}")
        sys.stdout.flush()
    print()

print()
print(f"SUMMARY: {n_total} cases, {n_fail} FAIL, {n_total-n_fail} OK")
if min_case:
    lam, k, E0, E2, ratio, R, F, L2ratio = min_case
    print(f"Minimum ratio: {ratio:.6f} at lambda={lam:.2f}, k={k}")
    print(f"  E_sig_w0/mu0={E0:.7f}, E_sig_w2/mu2={E2:.7f}")
    print(f"  R={R:.5f}, F={F:.4f}, L2ratio={L2ratio:.4f}")
    print(f"  ratio > 1: {ratio > 1}")
print()
print("CONCLUSION: E[sigma_within_v2]/mean_v2 > E[sigma_within_v0]/mean_v0")
print("  for ALL tested cases? -> FAIL count:", n_fail)
print("done")
