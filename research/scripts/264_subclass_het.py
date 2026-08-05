"""
264_subclass_het.py
===================
STRUCTURAL PROOF PATH for c2/c0 <= R = mean_v2/mean_v0.

KEY INSIGHT (from K-L forcing analysis):
  v2[s] cb-forcing for s ≡ 0 mod 3: uses c1 = (A/rho)*c0   (SMALL)
  v2[s] cb-forcing for s ≡ 1 mod 3: uses c0                  (MEDIUM)
  v2[s] cb-forcing for s ≡ 2 mod 3: uses c2                  (LARGE, self-referential)

  v0[s] cb-forcing for s ≡ 0 mod 3: uses c0  (its own column-min)
  v0[s] cb-forcing for s ≡ 1 mod 3: uses c2  (cross-class)
  v0[s] cb-forcing for s ≡ 2 mod 3: uses c2  (cross-class)

v2 has MORE between-sub-class forcing variation (c1 << c0 << c2 range),
while v0 has LESS variation (c0, c2 range only). This creates larger
relative variance in v2 => lower min-to-mean ratio for v2 => c2/c0 < R.

VERIFICATION:
  1. Check cb-classes feeding into each v0/v2 sub-class
  2. Compute sub-class means and column-minimum means
  3. Verify the between-sub-class variance is larger for v2
  4. Verify c2/mean_v2 < c0/mean_v0
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 700 + 100*max(0, k-8)
    A  = lam**-2.0; B1 = lam**(ALPHA-2.0); B3 = lam**(ALPHA-1.0)
    N  = 3**(k-1); Nl = N//3
    i  = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0_mask, m2_mask = (r_arr==0), (r_arr==2)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A*v[T4]
        w[m2_mask] += B3*cb[R3[m2_mask]]
        w[m0_mask] += B1*cb[R1[m0_mask]]
        v = w/w.max()
    return v, Nl, A, B1, B3

def analyze(k, lam):
    v, Nl, A, B1, B3 = run_kl(k, lam)
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    j  = np.arange(Nl, dtype=np.int64)
    s  = np.arange(Nl, dtype=np.int64)

    # rho
    v1 = v[1::3]; v0 = v[0::3]
    sigma1 = (4*s+2) % Nl
    rho = A / float(np.mean(v1 / v0[sigma1]))

    # v0/v2 sub-class means
    v0_arr = v[0::3]  # r=0 nodes
    v2_arr = v[2::3]  # r=2 nodes

    # R3 mapping: for v2 nodes, R3(s) = (2s+1) mod Nl
    # Class of R3(s) = (2s+1)%3: s=0->1, s=1->0, s=2->2
    R3_s = (2*s+1) % Nl
    R3_class = R3_s % 3  # class of cb used by v2[s]

    # R1 mapping: for v0 nodes, R1(s) = (4s)%Nl
    # Class of R1(s) = s%3 (since 4s%3 = s%3)
    R1_s = (4*s) % Nl
    R1_class = R1_s % 3  # class of cb used by v0[s]

    # Verify cb classes
    # For v2: should be R3_class = (2s+1)%3
    # For v0: should be R1_class = s%3

    # Sub-class means of v0 and v2
    sc0_v0 = float(np.mean(v0_arr[s%3==0]))
    sc1_v0 = float(np.mean(v0_arr[s%3==1]))
    sc2_v0 = float(np.mean(v0_arr[s%3==2]))

    sc0_v2 = float(np.mean(v2_arr[s%3==0]))
    sc1_v2 = float(np.mean(v2_arr[s%3==1]))
    sc2_v2 = float(np.mean(v2_arr[s%3==2]))

    mean_v0 = float(np.mean(v0_arr))
    mean_v2 = float(np.mean(v2_arr))

    # Column minimums for v0 and v2 type
    # cb_v0_col[j] = min(v0[j], v0[j+Nl//3], v0[j+2*Nl//3]) for j in [0, Nl//3)
    Nl3 = Nl//3
    j3 = np.arange(Nl3, dtype=np.int64)
    cb_v0_col = np.minimum(np.minimum(v0_arr[j3], v0_arr[j3+Nl3]), v0_arr[j3+2*Nl3])
    cb_v2_col = np.minimum(np.minimum(v2_arr[j3], v2_arr[j3+Nl3]), v2_arr[j3+2*Nl3])

    c0 = float(np.mean(cb_v0_col))  # = mean of v0 column-minimums
    c2 = float(np.mean(cb_v2_col))  # = mean of v2 column-minimums
    R = mean_v2/mean_v0

    # cb classes for v2 per sub-class
    # For s≡0: R3(s)≡1 mod3 (uses c1 = (A/rho)*c0)
    # For s≡1: R3(s)≡0 mod3 (uses c0)
    # For s≡2: R3(s)≡2 mod3 (uses c2, self-referential)
    cb_arr = np.array([float(np.mean(cb[j%3==0])),   # c0
                       float(np.mean(cb[j%3==1])),   # c1
                       float(np.mean(cb[j%3==2]))])  # c2

    # Forcing cb class per v2 sub-class
    cb_for_v2_sc0 = cb_arr[1]  # class 1 = c1
    cb_for_v2_sc1 = cb_arr[0]  # class 0 = c0
    cb_for_v2_sc2 = cb_arr[2]  # class 2 = c2

    # Forcing cb class per v0 sub-class (R1(s) class = s%3)
    cb_for_v0_sc0 = cb_arr[0]  # class 0 = c0
    cb_for_v0_sc1 = cb_arr[1]  # class 1 = c1
    cb_for_v0_sc2 = cb_arr[2]  # class 2 = c2

    # Between-sub-class CV (coefficient of variation) of cb forcing
    v2_forcing = np.array([cb_for_v2_sc0, cb_for_v2_sc1, cb_for_v2_sc2])
    v0_forcing = np.array([cb_for_v0_sc0, cb_for_v0_sc1, cb_for_v0_sc2])
    cv_v2 = float(np.std(v2_forcing)/np.mean(v2_forcing))
    cv_v0 = float(np.std(v0_forcing)/np.mean(v0_forcing))

    # min-to-mean ratio for v0 and v2
    m2m_v0 = c0/mean_v0
    m2m_v2 = c2/mean_v2

    return {
        'c0': c0, 'c2': c2,
        'mean_v0': mean_v0, 'mean_v2': mean_v2,
        'R': R, 'c2_c0': c2/c0,
        'rho': rho, 'A': A,
        'm2m_v0': m2m_v0, 'm2m_v2': m2m_v2,
        'cv_v2': cv_v2, 'cv_v0': cv_v0,
        'sc0_v2': sc0_v2, 'sc1_v2': sc1_v2, 'sc2_v2': sc2_v2,
        'sc0_v0': sc0_v0, 'sc1_v0': sc1_v0, 'sc2_v0': sc2_v0,
        'cb_v2': v2_forcing, 'cb_v0': v0_forcing,
    }

print("264: Sub-class heterogeneity -- structural proof of c2/c0 <= R")
print("="*75)
print()

# Detailed analysis at lam=1.70, k=10
lam, k = 1.70, 10
d = analyze(k, lam)
print(f"Detailed: lam={lam}, k={k}")
print(f"  c2/c0 = {d['c2_c0']:.6f}  vs  R = mean_v2/mean_v0 = {d['R']:.6f}")
print(f"  c2/c0 <= R: {d['c2_c0'] <= d['R']}")
print()
print(f"  min-to-mean ratio: m2m_v0 = c0/mean_v0 = {d['m2m_v0']:.6f}")
print(f"  min-to-mean ratio: m2m_v2 = c2/mean_v2 = {d['m2m_v2']:.6f}")
print(f"  m2m_v2 < m2m_v0: {d['m2m_v2'] < d['m2m_v0']}  (KEY: v2 more compressed)")
print()
print(f"  cb forcing per v2 sub-class: sc0={d['cb_v2'][0]:.5f}(c1), sc1={d['cb_v2'][1]:.5f}(c0), sc2={d['cb_v2'][2]:.5f}(c2)")
print(f"  cb forcing per v0 sub-class: sc0={d['cb_v0'][0]:.5f}(c0), sc1={d['cb_v0'][1]:.5f}(c1), sc2={d['cb_v0'][2]:.5f}(c2)")
print(f"  CV of cb forcing for v2: {d['cv_v2']:.5f}")
print(f"  CV of cb forcing for v0: {d['cv_v0']:.5f}")
print(f"  CV(v2) > CV(v0): {d['cv_v2'] > d['cv_v0']}  (KEY: v2 more heterogeneous forcing)")
print()
print(f"  v2 sub-class means: sc0={d['sc0_v2']:.5f}, sc1={d['sc1_v2']:.5f}, sc2={d['sc2_v2']:.5f}")
print(f"  v0 sub-class means: sc0={d['sc0_v0']:.5f}, sc1={d['sc1_v0']:.5f}, sc2={d['sc2_v0']:.5f}")
print(f"  CV of v2 sub-class means: {np.std([d['sc0_v2'],d['sc1_v2'],d['sc2_v2']])/np.mean([d['sc0_v2'],d['sc1_v2'],d['sc2_v2']]):.5f}")
print(f"  CV of v0 sub-class means: {np.std([d['sc0_v0'],d['sc1_v0'],d['sc2_v0']])/np.mean([d['sc0_v0'],d['sc1_v0'],d['sc2_v0']]):.5f}")
print()

# Scan
print(f"\nScan (k=8): m2m_v2 < m2m_v0 and CV(v2)>CV(v0) for all lambda?")
print(f"{'lam':>6}  {'m2m_v0':>8}  {'m2m_v2':>8}  {'m2m_v2<v0':>10}  {'CV_v2':>8}  {'CV_v0':>8}  {'CV(v2)>v0':>10}  {'c2/c0<=R':>10}")
k = 8
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    d = analyze(k, lam)
    ok_m2m = d['m2m_v2'] < d['m2m_v0']
    ok_cv  = d['cv_v2'] > d['cv_v0']
    ok_r   = d['c2_c0'] <= d['R']
    print(f"lam={lam:.2f}  {d['m2m_v0']:>8.5f}  {d['m2m_v2']:>8.5f}  {str(ok_m2m):>10}  {d['cv_v2']:>8.5f}  {d['cv_v0']:>8.5f}  {str(ok_cv):>10}  {str(ok_r):>10}")
    sys.stdout.flush()

print(f"""
STRUCTURAL ARGUMENT FOR c2/c0 <= R:

For v2 sub-classes, the cb forcing is:
  sub-class 0 (s≡0): R3(s)%3 = 1 => uses c1 = (A/rho)*c0   [SMALLEST]
  sub-class 1 (s≡1): R3(s)%3 = 0 => uses c0               [MEDIUM]
  sub-class 2 (s≡2): R3(s)%3 = 2 => uses c2               [LARGEST]
  Range: c1 to c2. Ratio c2/c1 = c2/c0 * rho/A >> 1 (large range!)

For v0 sub-classes, the cb forcing is:
  sub-class 0 (s≡0): R1(s)%3 = 0 => uses c0               [MEDIUM]
  sub-class 1 (s≡1): R1(s)%3 = 1 => uses c1               [SMALLEST]
  sub-class 2 (s≡2): R1(s)%3 = 2 => uses c2               [LARGEST]
  Same range: c1 to c2. Same pattern!

Wait -- same pattern for v0 and v2?? Then the between-sub-class forcing
is the same and does NOT explain the gap c2/mean_v2 < c0/mean_v0.

=> The difference must come from the T4 TERM (A or A^2/rho), not the cb term.
=> For v2: T4 term is (A^2/rho)*v0 (smaller coefficient).
=> For v0: T4 term is A*v2 (larger, since v2 > v0 and A > A^2/rho).

CORRECTED ARGUMENT: The T4 term for v2 is SMALLER relative to cb than for v0.
This means v2's forcing is MORE dominated by the cb term (more uniform across s),
while v0 has stronger T4 influence (more v2-heterogeneity injected).

Hmm, need to re-examine. See numerical CV values above.
""")

print("done")
