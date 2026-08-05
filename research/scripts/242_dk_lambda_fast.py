"""
242_dk_lambda_fast.py
=====================
Fast lambda-scan for d_k = ve0(k+1)/ve0(k) using k=8..12 only.
Complements Script 241 (which uses k=12..15, slower).

Key question: is d_k < 1 for ALL lambda in relevant range?
Script 241 (deep, slow) confirmed lam=1.30: d~0.57, lam=1.40: d~0.63.
This script gives all lambda values quickly at shallower depth.

Also: does d_k depend on k at these shallow depths? Shows the trend.
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
    return v, Nl

def ve0(v, Nl):
    """CODE-variance of r=0 component via block-triplet structure."""
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    v0 = v[0::3]
    a0 = v0[s0]; a1 = v0[s0+Nl3]; a2 = v0[s0+2*Nl3]
    mean_a = (a0+a1+a2)/3.0
    log_dev = np.stack([np.log2(a0/mean_a), np.log2(a1/mean_a), np.log2(a2/mean_a)])
    return float(np.mean(log_dev**2))

def Vk_full(v, Nl):
    """V_k = variance of log2(v / top-triple-mean) over all components."""
    Nl3 = Nl // 3
    T = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    lm = T.mean(axis=0)
    log_dev = np.log2(T) - np.log2(lm)[None,:]
    return float(np.var(log_dev))

LAMS = [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]
K_RANGE = [8, 9, 10, 11, 12]  # max N = 3^12 = 531k, very fast

print("242: Fast lambda-scan of d_k = ve0(k+1)/ve0(k)")
print(f"Lambda scan: {LAMS}")
print(f"k-range: {K_RANGE}")
print("="*90)

# Header
hdr = f"{'lam':>6}  " + "  ".join(f"{'ve0('+str(k)+')':>10}" for k in K_RANGE)
hdr += "  " + "  ".join(f"{'d_k='+str(k):>8}" for k in K_RANGE[:-1])
print(hdr)
sys.stdout.flush()

all_results = {}
for lam in LAMS:
    ves = {}
    for k in K_RANGE:
        v, Nl = run_kl(k, lam)
        ves[k] = ve0(v, Nl)

    dk_vals = [ves[k+1]/ves[k] for k in K_RANGE[:-1]]
    all_results[lam] = (ves, dk_vals)

    ve_str = "  ".join(f"{ves[k]:>10.6f}" for k in K_RANGE)
    dk_str = "  ".join(f"{d:>8.5f}" for d in dk_vals)
    print(f"lam={lam:.2f}  {ve_str}  {dk_str}")
    sys.stdout.flush()

print()
print("Summary: d_k at k=11 (ve0(12)/ve0(11)) as function of lambda:")
for lam, (ves, dk_vals) in all_results.items():
    d11 = ves[12]/ves[11]
    print(f"  lam={lam:.2f}: d_11={d11:.5f}  A={lam**-2:.4f}  B1={lam**(ALPHA-2):.4f}  B3/B1=lam={lam:.4f}")
sys.stdout.flush()

print()
print("d_k < 1 for all lambda?", all(
    dk < 1.0
    for lam, (ves, dk_vals) in all_results.items()
    for dk in dk_vals
))

print()
print("V_k (full CODE-variance) at k=11, 12 for each lambda:")
print(f"{'lam':>6}  {'Vk(11)':>10}  {'Vk(12)':>10}  {'dVk':>8}")
for lam in LAMS:
    v11, Nl11 = run_kl(11, lam)
    v12, Nl12 = run_kl(12, lam)
    Vk11 = Vk_full(v11, Nl11)
    Vk12 = Vk_full(v12, Nl12)
    print(f"lam={lam:.2f}  {Vk11:>10.6f}  {Vk12:>10.6f}  {Vk12/Vk11:>8.5f}")
sys.stdout.flush()

print()
print("done")
