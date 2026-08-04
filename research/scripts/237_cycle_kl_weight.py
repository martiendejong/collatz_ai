"""
237_cycle_kl_weight.py
======================
SP3A: Perron weight v^(k)(c) at candidate cycle residue classes.
SP3B: sigma1 anti-correlation (rho1 < 0) as quantitative constraint on cycle structures.

A Collatz cycle with k_cyc odd steps has elements that, reduced mod 3^{K-1},
land at specific K-L positions. The K-L Perron weight at those positions
measures how "compatible" the cycle is with the global density distribution.

SP3A: For small K (K=4..10), compute v^(K) and evaluate it at:
  - The trivial cycle position (n=1 -> K-L index 1, which has r=1, s=0)
  - Random high/low weight positions (for comparison)
  - Candidate cycle residue classes from SP1A

SP3B: sigma1 autocorrelation rho1 < 0 as cycle constraint.
  A cycle element at K-L position s with r=1 satisfies:
    v(r=1, s) = (A/rho) * v(r=0, sigma1(s))
  If s is a cycle candidate (high weight), then sigma1(s) must also be high weight.
  But rho1 < 0 means sigma1 ANTI-correlates: high-weight s -> lower-weight sigma1(s).
  Quantify: what fraction of positions s can have BOTH v1(s) and v0(sigma1(s)) above median?
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

print("237: K-L Perron weight at cycle candidates (SP3A + SP3B)")
print(f"     lambda={LAM}  A={A:.6f}  B1={B1:.6f}  B3={B3:.6f}")
print("=" * 72)
sys.stdout.flush()

# ============================================================
# SP3A: Perron weight at specific positions
# ============================================================
print("\n--- SP3A: K-L Perron weight at cycle candidate positions ---")
print()

# Trivial cycle: n=1.
# In K-L indexing at depth K: the position of 1.
# 1 mod 3 = 1 (r=1), 1 mod 9 = 1 (r=1, s=0), 1 mod 27 = 1, etc.
# K-L index i = 3*s + r = 3*0 + 1 = 1.
# So the trivial cycle (n=1) corresponds to K-L index i=1 for all K.

TRIVIAL_IDX = 1  # K-L index for n=1 (r=1, s=0)

for K in range(4, 12):
    N  = 3 ** (K - 1)
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

    # Weight at trivial cycle position (i=1, r=1, s=0)
    w_trivial = v[TRIVIAL_IDX]
    w_max = v.max()  # always 1 (normalized)
    w_mean = v.mean()
    w_median = float(np.median(v))
    w_p99 = float(np.percentile(v, 99))

    # Percentile rank of the trivial position
    pct_rank = float(np.mean(v < w_trivial)) * 100.0

    # r-type of trivial position
    r_trivial = TRIVIAL_IDX % 3  # =1
    s_trivial = TRIVIAL_IDX // 3  # =0

    # sigma1 image of s=0: sigma1(0) = (4*0+2) % Nl = 2 % Nl
    sigma1_s0 = (4 * 0 + 2) % Nl
    w_sigma1_s0 = v[3 * sigma1_s0 + 0]  # r=0 node at sigma1(s=0)

    print(f"K={K:2d}  N={N:8d}  v[1]={w_trivial:.6f}  mean={w_mean:.6f}  median={w_median:.6f}  "
          f"p99={w_p99:.6f}  rank={pct_rank:.1f}%  "
          f"v0[sigma1(0)]={w_sigma1_s0:.6f}")
    sys.stdout.flush()

print()
print("Interpretation:")
print("  v[1] = K-L Perron weight at position 1 (trivial cycle element).")
print("  'rank' = percentile rank of v[1] among all N positions.")
print("  If v[1] is low-percentile, the trivial cycle sits in a 'sparse' region of the Perron vector.")
print()

# ============================================================
# SP3B: sigma1 anti-correlation constraint
# ============================================================
print("--- SP3B: sigma1 anti-correlation as cycle constraint ---")
print()
print("For a cycle element with r=1 at position s:")
print("  v(r=1, s) = (A/rho) * v(r=0, sigma1(s))   [K-L r=1 equation]")
print("A cycle needs high v(r=1, s) (to be a dominant node in the density).")
print("This forces v(r=0, sigma1(s)) to be high as well.")
print("But rho1 = Corr[log v0(s), log v0(sigma1(s))] < 0 (anti-correlation).")
print("=> High v0(s) tends to imply LOWER v0(sigma1(s)).")
print("=> A cycle element at high-weight s is surrounded by LOWER-weight sigma1(s).")
print()
print("Correct joint fraction: fraction of s where BOTH v0(s) > median(v0) AND v0(sigma1(s)) > median(v0).")
print("Note: v1(s) = (A/rho)*v0(sigma1(s)), so asking about v1(s) is equivalent to asking about")
print("v0(sigma1(s)). The meaningful anti-correlation is in the (s, sigma1(s)) pair of v0 values.")
print("With rho1 < 0: high v0(s) -> lower v0(sigma1(s)), so joint-high fraction < 0.25.")
print()
print("Relevance: a Collatz cycle element at r=0 position s generates r=1 element at sigma1^{-1}(s).")
print("For the cycle to persist, the r=0 backbone must maintain high weight at both s AND sigma1(s).")
print("The anti-correlation of sigma1 suppresses this joint requirement.")
print()

print(f"{'K':>4}  {'rho1':>8}  {'f(s)&f(s1)>med':>16}  {'expected_ind':>13}  {'frac/0.25':>10}")

for K in range(4, 14):
    N  = 3 ** (K - 1)
    Nl = N // 3

    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr == 0), (r_arr == 2)
    R1 = (4 * s_arr) % Nl
    R3 = (2 * s_arr + 1) % Nl

    sl = np.arange(Nl, dtype=np.int64)
    sigma1 = (4 * sl + 2) % Nl

    v = np.ones(N, dtype=np.float64)
    for it in range(N_ITER):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()

    # r-type extraction (interleaved)
    v0 = v[0::3]  # r=0 nodes, indexed by s
    v1 = v[1::3]  # r=1 nodes

    # sigma1 autocorrelation (on r=0 backbone)
    f0 = np.log2(v0)
    f0_shifted = f0[sigma1]
    f0_c  = f0 - f0.mean()
    f0s_c = f0_shifted - f0_shifted.mean()
    var_f0 = float(np.var(f0))
    rho1 = float(np.mean(f0_c * f0s_c)) / var_f0 if var_f0 > 1e-15 else 0.0

    # CORRECT joint fraction: both v0(s) and v0(sigma1(s)) are above median of v0
    med_v0 = float(np.median(v0))
    high_v0_s      = v0        > med_v0   # s where v0(s) is above median
    high_v0_sig1s  = v0[sigma1] > med_v0  # s where v0(sigma1(s)) is above median

    frac_joint = float(np.mean(high_v0_s & high_v0_sig1s))
    expected_ind = 0.25  # if rho1=0
    ratio = frac_joint / expected_ind  # < 1 means anti-correlation reduces joint prob

    # Also: r=1 weight percentile rank vs r=0
    pct_v1_vs_all = float(np.mean(v1.mean() > v0))  # is mean(v1) above median(v0)?

    print(f"K={K:>3}  rho1={rho1:>8.4f}  frac_joint={frac_joint:>14.4f}  "
          f"expected_ind={expected_ind:>13.4f}  ratio={ratio:>10.4f}")
    sys.stdout.flush()

print()
print("Interpretation:")
print("  ratio < 1: sigma1 anti-correlation reduces joint-high probability.")
print("  A cycle element needs v0 high at BOTH s AND sigma1(s).")
print("  With rho1 < 0, this is less likely than independence predicts.")
print()
print("  For a k-element cycle traversing k different s-positions, and assuming")
print("  consecutive cycle positions are related by sigma1 (worst case for r=1 elements),")
print("  the probability all k positions are jointly high is:")
print("  P(all k high) ~ frac_joint^k  [geometric suppression].")
print()

import math

# Use K=10 for illustration
K_demo = 10
N  = 3 ** (K_demo - 1)
Nl = N // 3
i  = np.arange(N, dtype=np.int64)
T4 = (4 * i + 2) % N
s_arr, r_arr = np.divmod(i, 3)
m0, m2 = (r_arr == 0), (r_arr == 2)
R1 = (4 * s_arr) % Nl
R3 = (2 * s_arr + 1) % Nl
sl = np.arange(Nl, dtype=np.int64)
sigma1_K10 = (4 * sl + 2) % Nl
v = np.ones(N, dtype=np.float64)
for it in range(N_ITER):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w  = A * v[T4]
    w[m2] += B3 * cb[R3[m2]]
    w[m0] += B1 * cb[R1[m0]]
    v = w / w.max()
v0 = v[0::3]
med_v0 = float(np.median(v0))
high_v0_s = v0 > med_v0
high_v0_sig1s = v0[sigma1_K10] > med_v0
frac_joint_K10 = float(np.mean(high_v0_s & high_v0_sig1s))

f0 = np.log2(v0)
f0s = f0[sigma1_K10]
f0_c = f0 - f0.mean(); f0s_c = f0s - f0s.mean()
rho1_K10 = float(np.mean(f0_c * f0s_c)) / float(np.var(f0))

print(f"At K={K_demo}: frac_joint = {frac_joint_K10:.4f}, rho1 = {rho1_K10:.4f}.")
print(f"  For k=35000 (minimum cycle length, Simons-de Weger):")
print(f"  P(all jointly high) ~ {frac_joint_K10:.4f}^35000")
if frac_joint_K10 > 0:
    log10_p = 35000 * math.log10(frac_joint_K10)
    print(f"  = 10^{log10_p:.0f}  [astronomically small if < 0.25]")
print()
print()
print("Note: This is not a rigorous proof (cycle elements are not independent).")
print("But it quantifies the K-L tension: cycles require a highly improbable configuration")
print("of the Perron eigenvector under the anti-correlation structure of sigma1.")
print()
print("done")
