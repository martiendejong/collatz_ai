"""
241_dk_lambda_scan.py
=====================
Lambda-dependence of d_k = ve0(k+1)/ve0(k) (CODE-variance decay rate).

Key question: is d_k < 1 for ALL lambda in the relevant range?
If so, V_k -> 0 for all lambda, not just lambda=1.70.

Script 238 Part 2 measured ve0 at fixed k=13 for lambda in {1.30,1.50,1.70,1.90},
but not d_k. This script measures d_k for a fine grid of lambda values.

Also measure: the LIMIT d_inf ~ d_k for large k (extrapolated from k=12,13,14).

Expected range for lambda: the Perron eigenvalue of the 3x+1 map is
  lambda = (3/2)^{1/alpha} where alpha = log2(3) = 1.585
Actually from the K-L analysis, the natural lambda is around 1.5-2.0.
The convergence of d_k to d_inf(lambda) is the key function to understand.
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
    v0 = v[0::3]  # interleaved r=0
    a0 = v0[s0]; a1 = v0[s0+Nl3]; a2 = v0[s0+2*Nl3]
    mean_a = (a0+a1+a2)/3.0
    log_dev = np.stack([np.log2(a0/mean_a), np.log2(a1/mean_a), np.log2(a2/mean_a)])
    return float(np.mean(log_dev**2))

def Vk(v, Nl):
    """Full CODE-variance V_k = (ve0 + ve1 + ve2) / 3."""
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    T = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    lm = T.mean(axis=0)
    log_dev = np.log2(T) - np.log2(lm)[None,:]
    return float(np.var(log_dev))

# Lambda scan
LAMS = [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]
K_RANGE = [12, 13, 14]  # Measure d_k = ve0(k+1)/ve0(k) at k=12,13

print("241: Lambda-dependence of d_k = ve0(k+1)/ve0(k)")
print(f"Lambda scan: {LAMS}")
print(f"k-range: {K_RANGE} + one extra for ratio")
print("="*80)
print(f"{'lam':>6}  " + "  ".join(f"{'ve0('+str(k)+')':>12}" for k in K_RANGE+[K_RANGE[-1]+1])
      + "  " + "  ".join(f"{'d_k='+str(k):>10}" for k in K_RANGE))
sys.stdout.flush()

all_results = {}
for lam in LAMS:
    ves = {}
    for k in K_RANGE + [K_RANGE[-1]+1]:
        if k not in ves:
            v, Nl = run_kl(k, lam)
            ves[k] = ve0(v, Nl)

    dk_vals = [ves[k+1]/ves[k] for k in K_RANGE]
    all_results[lam] = (ves, dk_vals)

    ve_str = "  ".join(f"{ves[k]:>12.6f}" for k in K_RANGE+[K_RANGE[-1]+1])
    dk_str = "  ".join(f"{d:>10.5f}" for d in dk_vals)
    print(f"lam={lam:.2f}  {ve_str}  {dk_str}")
    sys.stdout.flush()

print()
print("Summary: d_k at k=13 as function of lambda (= ve0(14)/ve0(13)):")
for lam, (ves, dk_vals) in all_results.items():
    print(f"  lam={lam:.2f}: d_13={dk_vals[-1]:.5f}  B3/B1={lam:.4f}  A={lam**-2:.4f}  B1={lam**(log2(3)-2):.4f}")
    sys.stdout.flush()

print()
# Also measure V_k (full CODE-variance) for comparison
print("V_k = (ve0+ve1+ve2)/3 at k=13, k=14 for each lambda:")
print(f"{'lam':>6}  {'Vk(13)':>10}  {'Vk(14)':>10}  {'dVk=Vk14/Vk13':>15}")
for lam in LAMS:
    v13, Nl13 = run_kl(13, lam)
    v14, Nl14 = run_kl(14, lam)
    Vk13 = Vk(v13, Nl13)
    Vk14 = Vk(v14, Nl14)
    print(f"lam={lam:.2f}  Vk13={Vk13:.6f}  Vk14={Vk14:.6f}  dVk={Vk14/Vk13:.5f}")
    sys.stdout.flush()

print()
print("done")
