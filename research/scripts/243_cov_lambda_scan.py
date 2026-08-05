"""
243_cov_lambda_scan.py
======================
Lambda-scan of the anti-correlation Cov(u2, u_cb) and related quantities.

Background: Script 240 showed (at lambda=1.70):
  ve0 = CODE-var(A*v2 + B1*cb) EXACTLY (sigma0 = R1)
  Cov(u2, u_cb) < 0 at ALL depths k=4..15
  This anti-correlation makes ve0 < ve2 (L > 1) and implies the blend COMPRESSES.

Key question: does Cov(u2, u_cb) < 0 hold for ALL lambda?
If yes, the L > 1 mechanism (and hence the compression to ve0 < ve2) is universal.

Also measure: weight w2 = A*mean_v2 / (A*mean_v2 + B1*mean_cb) vs lambda.
And: does the "near-analytical L>1 argument" (Script 240 rem:ve0_blend) hold universally?
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
    # Compute rho = max of updated field before normalization
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    ww = A * v[T4]
    ww[m2] += B3 * cb[R3[m2]]
    ww[m0] += B1 * cb[R1[m0]]
    rho = ww.max()
    return v, Nl, rho, A, B1, B3

def code_cov(v, Nl):
    """CODE-variance of ve0, ve2, ve_cb and Cov(u2, u_cb), u2, u_cb."""
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]
    cb = np.minimum(np.minimum(v0, v1), v2)

    # CODE-variances via block-triplet structure
    def triplet_stats(x):
        a0 = x[s0]; a1 = x[s0+Nl3]; a2 = x[s0+2*Nl3]
        mean_a = (a0+a1+a2)/3.0
        u0 = np.log2(a0/mean_a)
        u1 = np.log2(a1/mean_a)
        u2 = np.log2(a2/mean_a)
        log_dev = np.stack([u0, u1, u2])
        ve = float(np.mean(log_dev**2))
        return ve, np.concatenate([u0, u1, u2])

    ve0, u0_all = triplet_stats(v0)
    ve2, u2_all = triplet_stats(v2)
    ve_cb, ucb_all = triplet_stats(cb)

    # Cov(u2, u_cb) over all (s, r) pairs
    cov_u2_ucb = float(np.mean(u2_all * ucb_all))
    mean_u2 = float(np.mean(u2_all))
    mean_ucb = float(np.mean(ucb_all))
    cov_u2_ucb -= mean_u2 * mean_ucb  # true covariance

    # Weights
    mean_v2 = float(np.mean(v2))
    mean_cb = float(np.mean(cb))

    return ve0, ve2, ve_cb, cov_u2_ucb, mean_v2, mean_cb

LAMS = [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]
K = 12  # depth (N = 3^11 = 177k, fast)

print("243: Lambda-scan of anti-correlation Cov(u2, u_cb)")
print(f"Depth k={K}, lambda scan: {LAMS}")
print(f"N = 3^{K-1} = {3**(K-1)}")
print("="*100)
print(f"{'lam':>6}  {'ve0':>10}  {'ve2':>10}  {'ve_cb':>10}  {'Cov(u2,ucb)':>14}  {'w2':>6}  {'L=ve2/ve0':>10}  {'cov<0?':>6}")
sys.stdout.flush()

results = {}
for lam in LAMS:
    v, Nl, rho, A, B1, B3 = run_kl(K, lam)
    ve0, ve2, ve_cb, cov, mean_v2, mean_cb = code_cov(v, Nl)

    w2 = A * mean_v2 / (A * mean_v2 + B1 * mean_cb)
    L = ve2 / ve0

    results[lam] = dict(ve0=ve0, ve2=ve2, ve_cb=ve_cb, cov=cov, w2=w2, L=L)
    print(f"lam={lam:.2f}  {ve0:>10.6f}  {ve2:>10.6f}  {ve_cb:>10.6f}  {cov:>14.6f}  {w2:>6.4f}  {L:>10.5f}  {'YES' if cov<0 else 'NO':>6}")
    sys.stdout.flush()

print()
print("Cov(u2, u_cb) < 0 universally:", all(results[lam]['cov'] < 0 for lam in LAMS))
print("L = ve2/ve0 > 1 universally:", all(results[lam]['L'] > 1 for lam in LAMS))

print()
print("Near-analytical L>1 check: is ve_cb/ve2 < (1+w2)/w_cb?")
print(f"{'lam':>6}  {'ve_cb/ve2':>10}  {'(1+w2)/wcb':>12}  {'satisfied?':>10}")
for lam in LAMS:
    r = results[lam]
    ratio = r['ve_cb'] / r['ve2']
    w2 = r['w2']
    threshold = (1 + w2) / (1 - w2)
    print(f"lam={lam:.2f}  {ratio:>10.4f}  {threshold:>12.4f}  {'YES' if ratio < threshold else 'NO':>10}")
sys.stdout.flush()

# Also do depth scan at fixed lambda=1.70 and lambda=2.00 to see cov trend
print()
print("Depth scan at lam=1.70 and lam=2.00 (Cov vs k):")
print(f"{'k':>4}  {'Cov(1.70)':>12}  {'Cov(2.00)':>12}  {'L(1.70)':>8}  {'L(2.00)':>8}")
for k in [6, 8, 10, 12]:
    row = {}
    for lam in [1.70, 2.00]:
        v, Nl, rho, A, B1, B3 = run_kl(k, lam)
        ve0, ve2, ve_cb, cov, mean_v2, mean_cb = code_cov(v, Nl)
        row[lam] = (cov, ve2/ve0)
    print(f"k={k:>2}  {row[1.70][0]:>12.6f}  {row[2.00][0]:>12.6f}  {row[1.70][1]:>8.4f}  {row[2.00][1]:>8.4f}")
    sys.stdout.flush()

print()
print("done")
