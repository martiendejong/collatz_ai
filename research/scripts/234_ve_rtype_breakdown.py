"""
234_ve_rtype_breakdown.py
=========================
Verifies cor:ve_equality and measures ve2_CODE / ve0_CODE ratio.

var_end_CODE = (1/3)(ve0_CODE + ve1_CODE + ve2_CODE).
Cor cor:ve_equality proved: ve1_CODE = ve0_CODE analytically.
This script checks that numerically and measures ve2_CODE / ve0_CODE.

For CODE variance: triplet at j is (v[j], v[j+Nl], v[j+2Nl]).
All three have r = j%3. Break by r: j in [0,Nl) with r=j%3.
  ve0_CODE = mean squared deviation for j%3==0 (Nl/3 indices, 3 deviations each)
  ve1_CODE = mean squared deviation for j%3==1
  ve2_CODE = mean squared deviation for j%3==2
Check: (ve0+ve1+ve2)/3 == var_end_CODE.
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)
N_ITER = 500

print("234: per-r-type CODE variance breakdown (lambda=1.70)")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}  N_iter={N_ITER}")
print("=" * 70)
sys.stdout.flush()

for k in range(4, 14):
    N  = 3 ** (k - 1)
    Nl = N // 3

    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % Nl
    R3 = (2 * s_arr + 1) % Nl

    v = np.ones(N, dtype=np.float64)
    for it in range(N_ITER):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()

    # CODE variance breakdown
    # j in [0, Nl): r = j%3, triplet = (v[j], v[j+Nl], v[j+2*Nl])
    j_arr = np.arange(Nl, dtype=np.int64)
    r_j   = j_arr % 3

    T_code    = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])  # shape 3 x Nl
    lmean_c   = T_code.mean(axis=0)                         # shape Nl
    log_dev   = np.log2(T_code) - np.log2(lmean_c)[None,:] # shape 3 x Nl

    # Total (check)
    var_end_code = float(np.var(log_dev))

    # Per r-type: average over the Nl/3 j's with r=j%3=r, and all 3 rows of T_code
    ve_by_r = []
    for r in range(3):
        mask = (r_j == r)   # Nl/3 True entries
        dev_r = log_dev[:, mask]   # shape 3 x (Nl/3)
        ve_by_r.append(float(np.mean(dev_r**2)))

    ve0, ve1, ve2 = ve_by_r
    check = (ve0 + ve1 + ve2) / 3.0
    ratio_10 = ve1 / ve0 if ve0 > 0 else float('nan')
    ratio_20 = ve2 / ve0 if ve0 > 0 else float('nan')

    print(f"k={k:2d}  ve0={ve0:.6f}  ve1={ve1:.6f}  ve2={ve2:.6f}  "
          f"ve1/ve0={ratio_10:.6f}  ve2/ve0={ratio_20:.6f}  "
          f"check={(ve0+ve1+ve2)/3:.6f}  var_code={var_end_code:.6f}")
    sys.stdout.flush()

print()
print("Expected: ve1/ve0 = 1.000000 (proved, cor:ve_equality)")
print("Measured: ve2/ve0 = ?")
print("done")
