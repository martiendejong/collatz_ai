"""
233_small_k_varend.py
=====================
var_end en d_k voor KLEINE k (k=3..12) bij lambda=1.70.
Doel: volledig beeld van d_k van k=3 tot k=12 (daarna k=13..19 bekend).
Geeft ook: de sigma1-autocorrelatie rho1 voor kleine k.
Vergelijkt var_end_CODE (Scripts 229/231) en var_end_BLOCK (type-scheiding).
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)
N_ITER = 500  # meer iteraties voor kleine k (sneller)

print("233: var_end en d_k voor k=3..12, lambda=1.70")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}")
print("=" * 70)
sys.stdout.flush()

varend_prev = None
K_RANGE = range(3, 13)

for k in K_RANGE:
    N  = 3 ** (k - 1)
    Nl = N // 3
    if Nl == 0: Nl = 1

    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % max(Nl, 1)
    R3 = (2 * s_arr + 1) % max(Nl, 1)

    v = np.ones(N, dtype=np.float64)
    for it in range(N_ITER):
        if Nl >= 1:
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:]) if N >= 3 else np.ones(1)
        else:
            cb = np.ones(1)
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()

    # var_end_CODE (Scripts 229/231 quantity)
    if Nl >= 1:
        T_code = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
        lmean_code = T_code.mean(axis=0)
        var_end_code = float(np.var(np.log2(T_code) - np.log2(lmean_code)[None,:]))
    else:
        var_end_code = 0.0

    # var_end_BLOCK (within-triplet type variation)
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    m_s = (v0 + v1 + v2) / 3.0
    X0b = np.log2(v0 / m_s); X1b = np.log2(v1 / m_s); X2b = np.log2(v2 / m_s)
    var_end_block = float(np.var(np.stack([X0b, X1b, X2b])))

    # sigma1 autocorrelation
    if Nl >= 1:
        sl = np.arange(Nl, dtype=np.int64)
        sigma1 = (4 * sl + 2) % Nl
        f0 = np.log2(v0)
        f0_shifted = f0[sigma1]
        f0_c  = f0 - f0.mean()
        f0s_c = f0_shifted - f0_shifted.mean()
        var_f0 = float(np.var(f0))
        rho1 = float(np.mean(f0_c * f0s_c)) / var_f0 if var_f0 > 1e-15 else 0.0
    else:
        rho1 = 0.0

    d_k = var_end_code / varend_prev if varend_prev and varend_prev > 0 else float('nan')

    print(f"k={k:2d}  N={N:8d}  ve_CODE={var_end_code:.6f}  d_k-1={d_k:.6f}  rho1={rho1:.6f}  ve_BLOCK={var_end_block:.6f}")
    sys.stdout.flush()
    varend_prev = var_end_code

print()
print("Bekende waarden k=13..19 (Scripts 229-231):")
known = [(13,0.001976,0.7560), (14,0.001494,0.7535), (15,0.001126,0.7590),
         (16,0.000855,0.7662), (17,0.000655,0.7692), (18,0.000504,0.7723), (19,0.000389,0.7753)]
for kk, ve, dk in known:
    print(f"  k={kk}: ve_CODE={ve:.6f}  d_k={dk:.4f}")
print()
print("done")
