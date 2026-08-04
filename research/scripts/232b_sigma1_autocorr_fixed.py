"""
232b_sigma1_autocorr_fixed.py
==============================
Bug fix for Script 232: v[:Nl] in interleaved ordering contains all three
r-types (v[j] has r=j%3, s=j//3). Correct extraction:
  v0 = v[0::3]  (r=0 nodes indexed by s=0..Nl-1)
  v1 = v[1::3]  (r=1 nodes indexed by s=0..Nl-1)
  v2 = v[2::3]  (r=2 nodes indexed by s=0..Nl-1)

Then v1[s] = (A/rho)*v0[sigma1(s)] exactly at the fixed point.
So rho1 = Cov[f0(s), f0(sigma1(s))] / Var[f0] = Cov[f0(s), f1(s)] / Var[f0].

Also computes two versions of var_end:
  var_end_BLOCK: Var over (r,s) of [log v_r(s) - log mean_r(s)]
                 where mean_r(s) = (v0[s]+v1[s]+v2[s])/3
                 (within-triplet TYPE variation at same s-position)
  var_end_CODE:  Var over j of [log v[j] - log mean(v[j],v[j+Nl],v[j+2Nl])]
                 (what Scripts 229/231 compute -- within-r-type SPATIAL variation)

This clarifies which var_end is the 'correct' one for the manuscript.
"""
import sys
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)
N_ITER = 300
RHO_TRUE = 1.055823   # known from Script 231

print("232b: sigma1-autocorrelatie gecorrigeerd (lambda=1.70)")
print(f"     A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}  N_iter={N_ITER}")
print("=" * 70)
sys.stdout.flush()

for k in [13, 14, 15, 16]:
    N  = 3 ** (k - 1)
    Nl = N // 3
    print(f"\nk={k}: N={N:,}  Nl={Nl:,}")
    sys.stdout.flush()

    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % Nl
    R3 = (2 * s_arr + 1) % Nl

    sl = np.arange(Nl, dtype=np.int64)
    sigma1 = (4 * sl + 2) % Nl   # sigma1 on s-coordinates

    v = np.ones(N, dtype=np.float64)
    for it in range(N_ITER):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()

    # CORRECT r-type extraction (interleaved ordering)
    v0 = v[0::3]   # r=0 nodes: v[0], v[3], v[6], ... length Nl
    v1 = v[1::3]   # r=1 nodes: v[1], v[4], v[7], ...
    v2 = v[2::3]   # r=2 nodes: v[2], v[5], v[8], ...

    # Verify r=1 equation: v1[s] = (A/rho)*v0[sigma1[s]]
    v1_pred = (A / RHO_TRUE) * v0[sigma1]
    rel_err = np.abs(v1 - v1_pred) / np.maximum(v1, 1e-15)
    print(f"  v1 = (A/rho)*v0[sigma1] check: max_rel_err = {rel_err.max():.2e}  mean_rel_err = {rel_err.mean():.2e}")

    f0 = np.log2(v0)
    f1 = np.log2(v1)
    f2 = np.log2(v2)

    # sigma1 autocorrelation (correct)
    f0_shifted = f0[sigma1]   # f0 at sigma1(s), should equal f1 + const
    f0_c   = f0 - f0.mean()
    f0s_c  = f0_shifted - f0_shifted.mean()
    f1_c   = f1 - f1.mean()

    var_f0 = float(np.var(f0))
    rho1   = float(np.mean(f0_c * f0s_c)) / var_f0
    rho_cross = float(np.mean(f0_c * f1_c)) / var_f0   # should == rho1

    mixing_gain = float(np.var(f0_shifted - f0)) / var_f0  # = 2*(1-rho1)

    print(f"  rho1       = Cov[f0(s), f0(sigma1(s))] / Var[f0] = {rho1:.6f}")
    print(f"  rho_cross  = Cov[f0(s), f1(s)] / Var[f0]         = {rho_cross:.6f}  (should == rho1)")
    print(f"  mixing_gain = Var[f0(sigma1)-f0] / Var[f0]        = {mixing_gain:.6f}  (should == 2*(1-rho1)={2*(1-rho1):.6f})")

    # VAR_END: two versions
    # 1. BLOCK version (within-triplet TYPE variation at same s)
    m_s_block = (v0 + v1 + v2) / 3.0
    X0b = np.log2(v0 / m_s_block)
    X1b = np.log2(v1 / m_s_block)
    X2b = np.log2(v2 / m_s_block)
    var_end_block = float(np.var(np.stack([X0b, X1b, X2b])))
    ve0b, ve1b, ve2b = float(np.var(X0b)), float(np.var(X1b)), float(np.var(X2b))

    # 2. CODE version (what Scripts 229/231 compute)
    T_code = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    lmean_code = T_code.mean(axis=0)
    var_end_code = float(np.var(np.log2(T_code) - np.log2(lmean_code)[None,:]))

    print(f"  var_end_BLOCK = {var_end_block:.8f}  (within-triplet TYPE variation)")
    print(f"    ve0b={ve0b:.6f}  ve1b={ve1b:.6f}  ve2b={ve2b:.6f}")
    print(f"  var_end_CODE  = {var_end_code:.8f}  (Scripts 229/231 quantity)")
    print(f"  ratio BLOCK/CODE = {var_end_block/var_end_code:.4f}")

    sys.stdout.flush()

print()
print("done")
