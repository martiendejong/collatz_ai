"""
258_logspace_anticorr.py
========================
Check whether the LOG-SPACE anti-correlation Corr(log v2[2s+1], log v2[s])
is stable with k (unlike the LINEAR-space Corr which weakens).

FROM Script 256: Corr(v2[2s+1], v2[s]) -> 0 as k increases (linear space weakens).
FROM Script 257b: d_k -> 0.756 < 1 (stays bounded below 1).

RESOLUTION: d_k depends on LOG-SPACE correlations, not linear-space.
If Corr(log v2[2s+1], log v2[s]) is STABLE with k, that explains why d_k is constant.

Also: directly compare Corr(log v2[s], log cb[(2s+1)]) with d_k.
The CODE-variance formula:
  d_k = ve0(k+1) / ve0(k)
where ve0 = within-CODE-triplet variance of log(v).

The Cov(log v2[s], log cb[(2s+1)]) enters the K-L operator formula for d_k.
"""
import numpy as np
from math import log2, log
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = 600 + 100 * max(0, k - 8)
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

def compute_ve0(v, Nl):
    """Correct CODE-variance = column-triplet variance of log(v)."""
    lv = np.log(v)
    lv_col = np.column_stack([lv[:Nl], lv[Nl:2*Nl], lv[2*Nl:]])
    return float(np.mean(np.var(lv_col, axis=1, ddof=0)))

def analyze_logspace(k, lam):
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    s = np.arange(Nl, dtype=np.int64)

    # Log-space values
    lv2 = np.log(v2)
    lcb = np.log(cb)

    # Corr(log v2[s], log v2[2s+1]) -- log-space version of doubling anti-corr
    max_s = min((Nl - 9) // 9, Nl // 3)
    s0 = np.arange(max(max_s, 2), dtype=np.int64)

    c_lin = float(np.corrcoef(v2[s0], v2[(2*s0+1)%Nl])[0,1])
    c_log = float(np.corrcoef(lv2[s0], lv2[(2*s0+1)%Nl])[0,1])

    # Corr(log v2[s], log cb[(2s+1)]) -- the key covariance for CODE-variance
    R3_s0 = (2*s0 + 1) % Nl
    c_v2_cb_log = float(np.corrcoef(lv2[s0], lcb[R3_s0])[0,1])

    # Cov(log v2[s], log cb[(2s+1)])
    cov_v2_cb_log = float(np.cov(lv2[s0], lcb[R3_s0])[0,1])

    # LOG-SPACE class means by s mod 3
    lv2_classes = [float(np.mean(lv2[s%3==r])) for r in range(3)]
    lv2_std = [float(np.std(lv2[s%3==r])) for r in range(3)]
    # Log-space class means (a0_log etc.)
    la0, la1, la2 = lv2_classes

    # Between-class contribution in log-space
    mean_la = np.mean(lv2_classes)
    var_btwn_la = np.mean(np.array(lv2_classes)**2) - mean_la**2
    cov_btwn_la = (la0*la1 + la1*la0 + la2*la2)/3 - mean_la**2
    c_btwn_log = cov_btwn_la / var_btwn_la if var_btwn_la > 0 else 0.0

    # CODE-variance
    ve0 = compute_ve0(v, Nl)

    return c_lin, c_log, c_v2_cb_log, cov_v2_cb_log, lv2_classes, lv2_std, c_btwn_log, ve0

print("258: Log-space anti-correlation stability with k")
print("="*70)

print(f"\nDepth scan lam=1.70:")
print(f"{'k':>4}  {'Corr_lin':>9}  {'Corr_log':>9}  {'Corr_log_cb':>12}  "
      f"{'logCls0':>8}  {'logCls1':>8}  {'logCls2':>8}  {'btwn_log':>9}  {'ve0':>9}")
prev_ve0 = None
lam = 1.70
for k in range(6, 15):
    v, Nl = run_kl(k, lam)
    c_lin, c_log, c_lvc, cov_lvc, la, la_std, c_btwn, ve0 = analyze_logspace(k, lam)
    print(f"k={k:>2}  {c_lin:>9.4f}  {c_log:>9.4f}  {c_lvc:>12.4f}  "
          f"{la[0]:>8.3f}  {la[1]:>8.3f}  {la[2]:>8.3f}  {c_btwn:>9.4f}  {ve0:>9.6f}")
    sys.stdout.flush()
    prev_ve0 = ve0

# LOG-SPACE ANALYSIS at k=8 in detail
print(f"\n\nDetailed log-space analysis k=8, lam=1.70:")
k, lam = 8, 1.70
v, Nl = run_kl(k, lam)
v2 = v[2::3]
cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
s = np.arange(Nl, dtype=np.int64)
lv2 = np.log(v2)
lcb = np.log(cb)

# Log-space class means and CVs
print(f"  Log-space class means (lv2 by s mod 3):")
for r in range(3):
    mask = (s % 3 == r)
    m = float(np.mean(lv2[mask]))
    sd = float(np.std(lv2[mask]))
    print(f"    log v2[s=={r} mod 3]: mean={m:.3f}  std={sd:.4f}  CV_log={sd/abs(m):.4f}")

# Linear-space class means for comparison
print(f"  Linear-space class means (v2 by s mod 3):")
for r in range(3):
    mask = (s % 3 == r)
    m = float(np.mean(v2[mask]))
    sd = float(np.std(v2[mask]))
    print(f"    v2[s=={r} mod 3]: mean={m:.4f}  std={sd:.5f}  CV={sd/m:.4f}")

# Log-space class mean ratio
la = [float(np.mean(lv2[s%3==r])) for r in range(3)]
print(f"\n  Log-space class mean differences: la0-la1 = {la[0]-la[1]:.4f}  (< 0 means la0 < la1?)")
print(f"  la0={la[0]:.3f}, la1={la[1]:.3f}, la2={la[2]:.3f}")
print(f"  exp(la0)/exp(la1) = {np.exp(la[0]-la[1]):.4f} = a0/a1 in linear space")

# Log-space doubling anti-corr stability
max_s = (Nl-9)//9
s0 = np.arange(max_s, dtype=np.int64)
print(f"\n  Log-space doubling anti-corr:")
c_log = float(np.corrcoef(lv2[s0], lv2[(2*s0+1)%Nl])[0,1])
c_lin = float(np.corrcoef(v2[s0], v2[(2*s0+1)%Nl])[0,1])
print(f"    Corr(lv2[2s+1], lv2[s]) = {c_log:+.4f}  (log-space)")
print(f"    Corr(v2[2s+1], v2[s])   = {c_lin:+.4f}  (linear-space)")

# Log-space class means doubling map prediction
mean_la = np.mean(la)
var_btwn = np.mean(np.array(la)**2) - mean_la**2
cov_btwn = (la[0]*la[1]+la[1]*la[0]+la[2]*la[2])/3 - mean_la**2
c_btwn = cov_btwn / var_btwn
print(f"\n  Log-space between-class Corr = {c_btwn:+.4f}")
print(f"  Fraction of log-space Corr = {c_btwn/c_log:.4f}")

print("\n\nLambda scan k=12: log-space vs linear-space anti-corr")
print(f"{'lam':>6}  {'Corr_lin':>9}  {'Corr_log':>9}  {'btwn_log':>9}  "
      f"{'d_12':>7}  {'la0-la1':>8}")
k = 12
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v, Nl = run_kl(k, lam)
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)
    lv2 = np.log(v2)
    max_s = (Nl-9)//9
    s0 = np.arange(max_s, dtype=np.int64)
    c_lin = float(np.corrcoef(v2[s0], v2[(2*s0+1)%Nl])[0,1])
    c_log = float(np.corrcoef(lv2[s0], lv2[(2*s0+1)%Nl])[0,1])
    la = [float(np.mean(lv2[s%3==r])) for r in range(3)]
    mean_la = np.mean(la); var_la = np.mean(np.array(la)**2)-mean_la**2
    cov_la = (la[0]*la[1]+la[1]*la[0]+la[2]*la[2])/3 - mean_la**2
    c_btwn = cov_la/var_la if var_la > 0 else 0.0
    ve0 = compute_ve0(v, Nl)
    v13, Nl13 = run_kl(k+1, lam)
    ve0_13 = compute_ve0(v13, Nl13)
    dk = ve0_13 / ve0
    print(f"lam={lam:.2f}  {c_lin:>9.4f}  {c_log:>9.4f}  {c_btwn:>9.4f}  "
          f"{dk:>7.4f}  {la[0]-la[1]:>8.4f}")
    sys.stdout.flush()

print("\ndone")
