"""
244_cov_block_lambda_scan.py
============================
Corrected lambda-scan of Cov(u_v2, u_cb) using the BLOCK minimum cb.

BUG IN SCRIPT 243: used cross-type cb = min(v0, v1, v2) at same s-position.
The K-L formula uses the BLOCK cb:
  cb[j] = min(v[j], v[j+Nl], v[j+2Nl])  for j in [0, Nl)
= CODE-triplet minimum within each r-type.

Script 240 used the correct block cb and found Cov(u_v2, u_cb) < 0 at lambda=1.70.
This script verifies Cov < 0 for ALL lambda in [1.30, 2.00].

Also: verify the ve0=CODE-var(f) formula holds for all lambda.
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
        # BLOCK minimum: cb[j] = min(v[j], v[j+Nl], v[j+2Nl]) for j in [0, Nl)
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v, Nl, A, B1

def code_var_nl(arr, Nl):
    """CODE-variance of Nl-dim array via block-triplet structure."""
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    a0 = arr[s0]; a1 = arr[s0+Nl3]; a2 = arr[s0+2*Nl3]
    mean_a = (a0 + a1 + a2) / 3.0
    log_dev = np.stack([np.log2(a0/mean_a), np.log2(a1/mean_a), np.log2(a2/mean_a)])
    return float(np.mean(log_dev**2)), log_dev

LAMS = [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]
K = 12  # depth (N = 3^11 = 177k, fast)

print("244: Lambda-scan of Cov(u_v2_sig0, u_cb_sig0) using BLOCK cb")
print(f"Depth k={K}, lambda scan: {LAMS}")
print("Block cb: cb[j] = min(v[j], v[j+Nl], v[j+2Nl])")
print("="*95)
print(f"{'lam':>6}  {'ve0':>10}  {'ve2':>10}  {'ve_cb':>10}  {'cov_term':>12}  {'cov<0?':>6}  {'L':>8}  {'ve_f=ve0?':>10}")
sys.stdout.flush()

results = {}
for lam in LAMS:
    v, Nl, A, B1 = run_kl(K, lam)

    # Block cb
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    # sigma0 = R1 = 4s mod Nl
    sigma0 = (4 * np.arange(Nl, dtype=np.int64)) % Nl

    # Interleaved r=0 and r=2 components
    v0_interleaved = v[0::3]  # v0[s] = v[3s]
    v2_interleaved = v[2::3]  # v2[s] = v[3s+2]

    # Pullback through sigma0
    v2_at_sigma0 = v2_interleaved[sigma0]
    cb_at_sigma0 = cb[sigma0]

    # Blend f = A*v2 + B1*cb (evaluated at sigma0)
    f = A * v2_at_sigma0 + B1 * cb_at_sigma0

    # CODE-variances
    ve0_direct, ld_v0 = code_var_nl(v0_interleaved, Nl)
    ve2, ld_v2_sig0 = code_var_nl(v2_at_sigma0, Nl)  # = CODE-var(v2) since sigma0 maps triplets
    ve_cb, ld_cb_sig0 = code_var_nl(cb_at_sigma0, Nl)  # = CODE-var(block_cb)
    ve_f, ld_f = code_var_nl(f, Nl)

    # Covariance: E[ld_v2 * ld_cb] (means ~= 0 for log-devs)
    cov_term = float(np.mean(ld_v2_sig0 * ld_cb_sig0))

    # Weights
    mean_v2 = float(np.mean(v2_at_sigma0))
    mean_cb = float(np.mean(cb_at_sigma0))
    w2 = A * mean_v2 / (A * mean_v2 + B1 * mean_cb)

    L = ve2 / ve0_direct
    ve_f_check = abs(ve_f - ve0_direct) / max(ve0_direct, 1e-20)

    results[lam] = dict(ve0=ve0_direct, ve2=ve2, ve_cb=ve_cb, cov=cov_term, w2=w2, L=L)

    print(f"lam={lam:.2f}  {ve0_direct:>10.6f}  {ve2:>10.6f}  {ve_cb:>10.6f}  "
          f"{cov_term:>12.8f}  {'YES' if cov_term<0 else 'NO':>6}  {L:>8.4f}  "
          f"{'OK' if ve_f_check < 1e-10 else 'FAIL':>10}")
    sys.stdout.flush()

print()
print(f"Cov(u_v2, u_cb) < 0 for ALL lambda: {all(results[lam]['cov'] < 0 for lam in LAMS)}")
print(f"L = ve2/ve0 > 1 for ALL lambda: {all(results[lam]['L'] > 1 for lam in LAMS)}")
print(f"ve_f = ve0 for ALL lambda: (checking above)")

print()
print("Summary table: cov and L vs lambda:")
for lam in LAMS:
    r = results[lam]
    print(f"  lam={lam:.2f}: cov={r['cov']:+.8f}  L={r['L']:.5f}  w2={r['w2']:.4f}  "
          f"ve_cb/ve2={r['ve_cb']/r['ve2']:.4f}")

# Depth scan at lambda=1.70 and lambda=2.00
print()
print("Depth scan: cov vs k at lam=1.70 and lam=2.00")
print(f"{'k':>4}  {'cov(1.70)':>14}  {'cov(2.00)':>14}  {'L(1.70)':>8}  {'L(2.00)':>8}")
for k in [5, 6, 7, 8, 9, 10, 11, 12]:
    row = {}
    for lam in [1.70, 2.00]:
        v, Nl, A, B1 = run_kl(k, lam)
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        sigma0 = (4 * np.arange(Nl, dtype=np.int64)) % Nl
        v2_interleaved = v[2::3]
        v2_at_sigma0 = v2_interleaved[sigma0]
        cb_at_sigma0 = cb[sigma0]
        ve0_d, ld_v0 = code_var_nl(v[0::3], Nl)
        ve2_d, ld_v2_sig0 = code_var_nl(v2_at_sigma0, Nl)
        ve_cb_d, ld_cb_sig0 = code_var_nl(cb_at_sigma0, Nl)
        cov_term = float(np.mean(ld_v2_sig0 * ld_cb_sig0))
        row[lam] = (cov_term, ve2_d / ve0_d)
    print(f"k={k:>2}  {row[1.70][0]:>14.8f}  {row[2.00][0]:>14.8f}  {row[1.70][1]:>8.4f}  {row[2.00][1]:>8.4f}")
    sys.stdout.flush()

print()
print("done")
