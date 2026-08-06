"""
307_gamma_extrapolation.py
==========================
Fast gamma(k) estimation via geometric extrapolation of rho_k(lambda)
(minutes instead of days; uses the Obs 502b discovery that rho converges
geometrically in k with stable rate).

Method:
  1. Compute rho_k(lambda) for k = 9..16 on a lambda-grid near the edge.
  2. Per lambda: fit rho_k = rho_inf - C*r^k on the deep k's; extrapolate to target k.
  3. Solve rho_k(lambda*) = 1 -> gamma(k) = log2(lambda*).
  4. VALIDATE against certified records gamma(17)=0.8953, gamma(19)=0.9069,
     gamma(20)=0.9146, gamma(21)=0.9184 (certified = lower bounds with margin
     ~2e-4 in min-ratio, so true lambda* slightly above certified lambda).
  5. Output gamma(22..28) and compare with the (2/3)^(1/6)-law and CEILING fork.

NOTE: this yields ESTIMATES (prediction tests), not certificates. Records still
need full exact-integer verification.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
KS = list(range(9, 17))
ITERS = {9: 3000, 10: 2500, 11: 2000, 12: 1500, 13: 1000, 14: 700, 15: 450, 16: 300}
LAMS = [1.85, 1.87, 1.89, 1.91, 1.93]

def rho_k(lam, k, niters):
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
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
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max())
        w /= rho
        v = w
    # CW gap for convergence check
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A*v[T4]
    w[m2] += B3*cb[R3[m2]]
    w[m0] += B1*cb[R1[m0]]
    r = w/v
    return float(r.max()), float(r.max()/r.min()-1)

print("computing rho_k(lambda) grid...")
RHO = {}
for lam in LAMS:
    for k in KS:
        rho, gap = rho_k(lam, k, ITERS[k])
        RHO[(lam, k)] = rho
        print(f"  lam={lam:.2f} k={k}: rho={rho:.9f} cwgap={gap:.1e}", flush=True)

def extrapolate(lam, ktarget):
    """fit rho_k = rho_inf - C*r^k on k=12..16, extrapolate"""
    kk = np.array([12, 13, 14, 15, 16], dtype=float)
    rr = np.array([RHO[(lam, int(x))] for x in kk])
    # solve 3-param geometric via ratios of successive differences
    d = np.diff(rr)
    r_est = np.mean(d[1:]/d[:-1])
    # rho_inf from last point + geometric tail
    C_last = d[-1]/(r_est**kk[-2]*(r_est-1)) if r_est != 1 else 0
    rho_inf = rr[-1] + d[-1]*r_est/(1-r_est)
    rho_t = rho_inf - (rho_inf - rr[-1])*r_est**(ktarget-kk[-1])
    return rho_t, r_est, rho_inf

def gamma_est(ktarget):
    vals = []
    for lam in LAMS:
        rho_t, r_est, _ = extrapolate(lam, ktarget)
        vals.append((lam, rho_t))
    # solve rho_t(lambda) = 1 by linear interpolation on the grid
    ls = np.array([x[0] for x in vals]); rs = np.array([x[1] for x in vals])
    # rho decreases in lambda beyond edge? find bracketing pair
    for j in range(len(ls)-1):
        if (rs[j]-1)*(rs[j+1]-1) <= 0:
            lam_star = ls[j] + (1-rs[j])*(ls[j+1]-ls[j])/(rs[j+1]-rs[j])
            return log2(lam_star), lam_star
    # extrapolate linearly if not bracketed
    sl = (rs[-1]-rs[0])/(ls[-1]-ls[0])
    lam_star = ls[0] + (1-rs[0])/sl
    return log2(lam_star), lam_star

print("\n=== VALIDATION against certified records (certified are lower bounds) ===")
cert = {17: 0.8953, 19: 0.9069, 20: 0.9146, 21: 0.9184}
for kt in [17, 19, 20, 21]:
    g, ls = gamma_est(kt)
    print(f"k={kt}: extrapolated gamma={g:.4f} lambda*={ls:.4f} | certified {cert[kt]} (diff {g-cert[kt]:+.4f})")

print("\n=== PREDICTIONS k=22..28 ===")
print("law (2/3)^(1/6): gamma_law(k) = 1 - 0.0816*(2/3)^((k-21)/6)")
for kt in range(22, 29):
    g, ls = gamma_est(kt)
    glaw = 1 - 0.0816*(2/3)**((kt-21)/6)
    print(f"k={kt}: extrapolated gamma={g:.4f} | law={glaw:.4f} (diff {g-glaw:+.4f})")
print("\nDONE")
