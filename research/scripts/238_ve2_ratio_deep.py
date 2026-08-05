"""
238_ve2_ratio_deep.py
=====================
Deep measurement of ve2/ve0 ratio.

cor:ve_equality proves ve1 = ve0 analytically.
This script asks: what is L = lim_{k->inf} ve2/ve0?

Key analytical fact (proved here numerically, explained analytically):
  sigma_20(s) = (16s+14) mod Nl  maps triplets to triplets.
  (Same argument as cor:ve_equality: 16 ≡ 1 mod 3.)
  => The pure transport term (A^2/rho^2)*v0[sigma_20(s)] has CODE-variance = ve0.
  => ve2 - ve0 comes entirely from the r=2 bonus (B3/rho)*cb[R3(s)].

Measurements:
1. ve2/ve0 for k=4..17 at lambda=1.70 (extending Script 234)
2. ve2/ve0 at k=13 for lambda in {1.30, 1.50, 1.70, 1.90} (lambda-dependence)
3. Analytical decomposition: what fraction of ve2 comes from transport vs bonus?
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
N_ITER = 500

def run_kl(k, lam, n_iter=N_ITER):
    A  = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    N  = 3 ** (k - 1)
    Nl = N // 3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % Nl
    R3 = (2 * s_arr + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v, Nl, A, B1, B3

def ve_breakdown(v, Nl):
    """Compute ve0, ve1, ve2 from the CODE-variance breakdown."""
    j_arr = np.arange(Nl, dtype=np.int64)
    r_j   = j_arr % 3
    T_code    = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    lmean_c   = T_code.mean(axis=0)
    log_dev   = np.log2(T_code) - np.log2(lmean_c)[None,:]
    var_code  = float(np.var(log_dev))
    ve_by_r = []
    for r in range(3):
        mask = (r_j == r)
        dev_r = log_dev[:, mask]
        ve_by_r.append(float(np.mean(dev_r**2)))
    return ve_by_r[0], ve_by_r[1], ve_by_r[2], var_code

# ============================================================
# Part 1: k=4..17 at lambda=1.70
# ============================================================
LAM = 1.70
print("238: ve2/ve0 deep measurement")
print(f"Part 1: k=4..17, lambda={LAM}")
print("=" * 72)
print(f"  {'k':>4}  {'ve0':>10}  {'ve1':>10}  {'ve2':>10}  {'ve2/ve0':>9}  {'d_ve0':>9}  {'d_ve2':>9}")
sys.stdout.flush()

prev_ve0 = None
prev_ve2 = None
results_k = []

for k in range(4, 18):
    v, Nl, A, B1, B3 = run_kl(k, LAM)
    ve0, ve1, ve2, vc = ve_breakdown(v, Nl)
    ratio = ve2/ve0 if ve0 > 0 else float('nan')
    d_ve0 = ve0/prev_ve0 if prev_ve0 and prev_ve0 > 0 else float('nan')
    d_ve2 = ve2/prev_ve2 if prev_ve2 and prev_ve2 > 0 else float('nan')
    print(f"  k={k:>3}  ve0={ve0:>10.6f}  ve1={ve1:>10.6f}  ve2={ve2:>10.6f}  "
          f"ratio={ratio:>9.5f}  d_ve0={d_ve0:>9.5f}  d_ve2={d_ve2:>9.5f}")
    sys.stdout.flush()
    results_k.append((k, ve0, ve1, ve2, ratio))
    prev_ve0 = ve0
    prev_ve2 = ve2

print()
# Extrapolate L from last two points
if len(results_k) >= 3:
    ks = [r[0] for r in results_k[-3:]]
    ls = [r[4] for r in results_k[-3:]]
    dl = [ls[i+1]-ls[i] for i in range(len(ls)-1)]
    # If dl is geometric (exponential decay), extrapolate
    if len(dl) >= 2 and dl[0] != 0 and dl[-1] != 0:
        rate = dl[-1]/dl[-2] if dl[-2] != 0 else 0
        L_extrap = ls[-1] + dl[-1]/(1-rate) if abs(rate) < 1 else ls[-1]
    else:
        L_extrap = ls[-1]
    print(f"  Last ve2/ve0 = {ls[-1]:.6f}, extrapolated L = {L_extrap:.6f}")
    print(f"  Simple fractions: 6/5 = {6/5:.6f}, 7/6 = {7/6:.6f}, 5/4 = {5/4:.6f}")

# ============================================================
# Part 2: Lambda-dependence at k=13
# ============================================================
K_FIXED = 13
LAMS = [1.30, 1.50, 1.70, 1.90]
print(f"\nPart 2: lambda-dependence at k={K_FIXED}")
print(f"  {'lam':>6}  {'ve0':>10}  {'ve2':>10}  {'ve2/ve0':>9}  {'B3/B1':>8}  {'B3':>8}")
for lam in LAMS:
    v, Nl, A, B1, B3 = run_kl(K_FIXED, lam, n_iter=N_ITER)
    ve0, ve1, ve2, vc = ve_breakdown(v, Nl)
    ratio = ve2/ve0 if ve0 > 0 else float('nan')
    print(f"  lam={lam:>4.2f}  ve0={ve0:>10.6f}  ve2={ve2:>10.6f}  "
          f"ratio={ratio:>9.5f}  B3/B1={B3/B1:>8.4f}  B3={B3:>8.4f}")
    sys.stdout.flush()

# ============================================================
# Part 3: Analytical decomposition at k=14, lambda=1.70
# ============================================================
K_ANAL = 14
print(f"\nPart 3: Transport vs bonus decomposition (k={K_ANAL}, lambda={LAM})")
v, Nl, A, B1, B3 = run_kl(K_ANAL, LAM)

lam = LAM
rho = float(np.max(v))  # Perron value (not exactly, but v is normalized to max=1)
# Better: estimate rho from the fixed-point equation
# rho = (A*v[T4] + B3*cb[R3] + B1*cb[R1]) / v at convergence
# Since v is normalized to max=1, we need to track the actual rho.
# Run one more step and measure the scale factor:
i = np.arange(3**( K_ANAL-1), dtype=np.int64)
T4 = (4*i+2) % (3**(K_ANAL-1))
s_arr, r_arr = np.divmod(i, 3)
m0, m2 = (r_arr==0), (r_arr==2)
R1 = (4*s_arr) % Nl
R3 = (2*s_arr+1) % Nl
cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
w = A*v[T4]
w[m2] += B3*cb[R3[m2]]
w[m0] += B1*cb[R1[m0]]
rho_est = float(w.max())  # = rho * v.max() = rho since v.max()=1
print(f"  Estimated rho = {rho_est:.6f}")

# sigma_20(s) = (16s+14) mod Nl
sl = np.arange(Nl, dtype=np.int64)
sigma_20 = (16*sl + 14) % Nl

v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

# v2(s) = (A^2/rho^2)*v0[sigma_20(s)] + (B3/rho)*cb[R3_s]
# where R3_s = (2s+1) mod Nl for s in [0,Nl)
R3_s = (2*sl + 1) % Nl
v2_transport = (A**2/rho_est**2) * v0[sigma_20]   # transport term
v2_bonus     = (B3/rho_est) * cb[R3_s]             # bonus term
v2_pred      = v2_transport + v2_bonus

rel_err = np.abs(v2 - v2_pred) / np.maximum(np.abs(v2), 1e-15)
print(f"  v2 = (A^2/rho^2)*v0[sigma_20] + (B3/rho)*cb[R3] check:")
print(f"    max_rel_err = {rel_err.max():.2e}   mean_rel_err = {rel_err.mean():.2e}")

# Triplet structure: CODE variance decomposition for r=2
# Triplets: j=0..Nl/3-1, r=2 means j=2,5,8,... in [0,Nl): (j, j+Nl, j+2Nl) in [0,N)
# Equivalently: for s=0..Nl/3-1, triplet = (v2(s), v2(s+Nl/3), v2(s+2Nl/3))
Nl3 = Nl // 3
s0 = np.arange(Nl3, dtype=np.int64)
T2_0 = v2[s0]; T2_1 = v2[s0+Nl3]; T2_2 = v2[s0+2*Nl3]
mean_T2 = (T2_0 + T2_1 + T2_2) / 3.0
log_dev_2 = np.stack([
    np.log2(T2_0/mean_T2), np.log2(T2_1/mean_T2), np.log2(T2_2/mean_T2)
])
ve2_direct = float(np.mean(log_dev_2**2))

# Same for transport term:
tr_0 = v2_transport[s0]; tr_1 = v2_transport[s0+Nl3]; tr_2 = v2_transport[s0+2*Nl3]
mean_tr = (tr_0+tr_1+tr_2)/3.0
log_dev_tr = np.stack([np.log2(tr_0/mean_tr),np.log2(tr_1/mean_tr),np.log2(tr_2/mean_tr)])
ve2_transport = float(np.mean(log_dev_tr**2))

# ve0 for comparison:
T0_0 = v0[s0]; T0_1 = v0[s0+Nl3]; T0_2 = v0[s0+2*Nl3]
mean_T0 = (T0_0+T0_1+T0_2)/3.0
log_dev_0 = np.stack([np.log2(T0_0/mean_T0),np.log2(T0_1/mean_T0),np.log2(T0_2/mean_T0)])
ve0_direct = float(np.mean(log_dev_0**2))

print(f"\n  ve0 (direct)             = {ve0_direct:.6f}")
print(f"  ve2_transport (A2/rho2)*v0[sigma20] CODE-var = {ve2_transport:.6f}")
print(f"  ve2 (direct, full)       = {ve2_direct:.6f}")
print(f"  ve2_transport/ve0 = {ve2_transport/ve0_direct:.6f}  (should be 1.0 if sigma_20 maps triplets)")
print(f"  ve2/ve0           = {ve2_direct/ve0_direct:.6f}  (measured ratio L)")
print()
print(f"  Extra ve2 from bonus = ve2 - ve2_transport = {ve2_direct - ve2_transport:.6f}")
print(f"  Extra fraction       = {(ve2_direct-ve2_transport)/ve0_direct:.6f} * ve0")
print()

# Check: sigma_20 maps triplets to triplets?
# sigma_20(s + Nl/3) should = sigma_20(s) + Nl/3 mod Nl
check_triplet = (sigma_20[s0+Nl3] - sigma_20[s0] - Nl3) % Nl
print(f"  sigma_20 maps triplets check: max|sigma_20(s+Nl/3)-sigma_20(s)-Nl/3|={check_triplet.max()}")
check_triplet2 = (sigma_20[s0+2*Nl3] - sigma_20[s0] - 2*Nl3) % Nl
print(f"  sigma_20(s+2Nl/3) check: max|...-2Nl/3|={check_triplet2.max()}")

print()
print("Summary: sigma_20 maps triplets to triplets (as sigma_1 does), so the")
print("transport term contributes ve0 to ve2. The bonus (B3/rho)*cb[R3(s)] adds")
print("extra variance, giving ve2/ve0 = L > 1.")
print()
print("done")
