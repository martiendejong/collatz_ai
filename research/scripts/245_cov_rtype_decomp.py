"""
245_cov_rtype_decomp.py
=======================
Decompose Cov(u_v2_sig0, u_cb_sig0) by the r-type of sigma0(j) = 4j mod Nl.

Analytical mechanism (from analysis):
  sigma0(j) = 4j%Nl has r-type (4j)%3 = j%3 (since 4 equiv 1 mod 3).
  For j equiv 2 mod 3 (r-type=2):
    v2_at_sigma0[j] = v2[4j%Nl] is ONE OF THE THREE values in the min defining cb[4j%Nl].
    => When v2[4j%Nl] is large (not the min), cb excludes it => cb smaller.
    => Structural anti-correlation: large u_v2 => small u_cb at same triplet position.

  For j equiv 0 or 1 mod 3 (r-type=0 or 1):
    cb[4j%Nl] uses v0 or v1, NOT v2 => no direct coupling => Cov could be positive or neutral.

  PREDICTION: Cov(j equiv 2) < 0 (strong negative)
               Cov(j equiv 0) and Cov(j equiv 1) positive or neutral
               Overall Cov < 0 because j equiv 2 contribution dominates.

This script verifies this decomposition at k=12, lambda=1.70 and also shows lambda dependence.
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
    return v, Nl, A, B1

def code_var_nl(arr, Nl):
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    a0 = arr[s0]; a1 = arr[s0+Nl3]; a2 = arr[s0+2*Nl3]
    mean_a = (a0 + a1 + a2) / 3.0
    log_dev = np.stack([np.log2(a0/mean_a), np.log2(a1/mean_a), np.log2(a2/mean_a)])
    return float(np.mean(log_dev**2)), log_dev  # log_dev shape: (3, Nl3)

def decomp_cov(v, Nl, A, B1):
    """Decompose Cov by r-type of sigma0(j) = 4j%Nl."""
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    sigma0 = (4 * np.arange(Nl, dtype=np.int64)) % Nl

    v2_interleaved = v[2::3]
    v2_at_sigma0 = v2_interleaved[sigma0]
    cb_at_sigma0 = cb[sigma0]

    # CODE-variance and log-devs
    ve0_direct, ld_v0 = code_var_nl(v[0::3], Nl)
    ve2, ld_v2_sig0 = code_var_nl(v2_at_sigma0, Nl)
    ve_cb, ld_cb_sig0 = code_var_nl(cb_at_sigma0, Nl)
    f = A * v2_at_sigma0 + B1 * cb_at_sigma0
    ve_f, _ = code_var_nl(f, Nl)

    # Total Cov
    cov_total = float(np.mean(ld_v2_sig0 * ld_cb_sig0))

    # Decompose by r-type of sigma0(j)
    # sigma0(j) = 4j%Nl has r-type (4j)%3 = j%3
    # The log_dev matrix has shape (3, Nl3) where axis 0 indexes triplet position (0,1,2)
    # and axis 1 indexes the triplet group j=0..Nl3-1.
    # For the r-type of sigma0(j): j in [0, Nl3) has j%3 in {0,1,2}

    Nl3 = Nl // 3
    j_idx = np.arange(Nl3, dtype=np.int64)  # j ranges over [0, Nl3)
    r_of_j = j_idx % 3  # r-type of sigma0(j) = j%3

    # ld_v2_sig0[:, j] for j in {j: j%3 == r} for each r
    cov_by_r = {}
    frac_by_r = {}
    for r in range(3):
        mask = (r_of_j == r)  # shape (Nl3,), True for triplet groups with r-type r
        frac = float(np.sum(mask)) / Nl3  # fraction of triplet groups with this r
        if np.sum(mask) == 0:
            cov_by_r[r] = 0.0
            frac_by_r[r] = 0.0
            continue
        # ld matrices for this r subset: shape (3, count_r)
        ld_v2_r = ld_v2_sig0[:, mask]
        ld_cb_r = ld_cb_sig0[:, mask]
        cov_r = float(np.mean(ld_v2_r * ld_cb_r))
        cov_by_r[r] = cov_r
        frac_by_r[r] = frac

    # Verify: cov_total should equal sum of fraction-weighted cov_r
    cov_reconstructed = sum(frac_by_r[r] * cov_by_r[r] for r in range(3))

    return {
        've0': ve0_direct, 've2': ve2, 've_cb': ve_cb, 've_f': ve_f,
        'cov_total': cov_total, 'cov_by_r': cov_by_r, 'frac_by_r': frac_by_r,
        'cov_reconstructed': cov_reconstructed
    }

print("245: Cov decomposition by r-type of sigma0(j)")
print("Analytical prediction: Cov(j equiv 2) < 0, Cov(j equiv 0,1) >= 0")
print("="*80)

# Main analysis at k=12, lambda=1.70
K = 12
for lam in [1.70, 1.30, 2.00]:
    v, Nl, A, B1 = run_kl(K, lam)
    res = decomp_cov(v, Nl, A, B1)

    print(f"\nlam={lam:.2f}, k={K}:")
    print(f"  ve0={res['ve0']:.6f}  ve2={res['ve2']:.6f}  ve_cb={res['ve_cb']:.6f}")
    print(f"  L=ve2/ve0={res['ve2']/res['ve0']:.4f}  ve_f=ve0? {abs(res['ve_f']-res['ve0'])/res['ve0']:.1e}")
    print(f"  Cov TOTAL = {res['cov_total']:+.8f}")
    print(f"  Cov by r-type:")
    for r in range(3):
        print(f"    r={r}: frac={res['frac_by_r'][r]:.4f}  Cov={res['cov_by_r'][r]:+.8f}  "
              f"contribution={res['frac_by_r'][r]*res['cov_by_r'][r]:+.8f}")
    print(f"  Reconstructed = {res['cov_reconstructed']:+.8f}  "
          f"(should match Cov TOTAL, err={abs(res['cov_reconstructed']-res['cov_total']):.1e})")
    print(f"  PREDICTION CHECK: Cov(r=2) < 0? {res['cov_by_r'][2] < 0}")
    sys.stdout.flush()

print()
print("="*80)
print("Lambda scan at k=12: Cov by r-type")
print(f"{'lam':>6}  {'cov_r0':>12}  {'cov_r1':>12}  {'cov_r2':>12}  {'cov_tot':>12}")
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    v, Nl, A, B1 = run_kl(K, lam)
    res = decomp_cov(v, Nl, A, B1)
    print(f"lam={lam:.2f}  {res['cov_by_r'][0]:>12.8f}  {res['cov_by_r'][1]:>12.8f}  "
          f"{res['cov_by_r'][2]:>12.8f}  {res['cov_total']:>12.8f}")
    sys.stdout.flush()

print()
print("done")
