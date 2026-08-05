"""
257b_dk_correct.py
==================
Corrected d_k computation.

BUG IN 257: used lv.reshape(Nl, 3) which groups {v[3s], v[3s+1], v[3s+2]}
(fixed s, varying r-type). This is the WRONG triplet.

CORRECT: CODE-triplet for j in [0, Nl) is {v[j], v[j+Nl], v[j+2Nl]}
(fixed r-type position, three s-positions in the same CODE-block).
This is the COLUMN-TRIPLET: same j-offset, varying CODE-block.

The cb block-minimum is ALSO the column-triplet:
  cb[j] = min(v[j], v[j+Nl], v[j+2Nl])

So CODE-variance = within-column-triplet variance of log(v).

RESULT EXPECTED: d_k = ve0(k+1)/ve0(k) < 1 (Conjecture G).
"""
import numpy as np
from math import log2
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

def compute_ve0_correct(v, Nl):
    """CODE-variance = within-column-triplet variance of log(v)."""
    lv = np.log(v)
    # Column triplet j: {lv[j], lv[j+Nl], lv[j+2Nl]} for j in [0, Nl)
    lv_col = np.column_stack([lv[:Nl], lv[Nl:2*Nl], lv[2*Nl:]])  # shape (Nl, 3)
    trip_var = np.var(lv_col, axis=1, ddof=0)
    return float(np.mean(trip_var))

def compute_ve0_wrong(v, Nl):
    """WRONG: within-row-triplet (fixed s, varying r)."""
    lv = np.log(v)
    lv_row = lv.reshape(Nl, 3)  # groups {v[3s], v[3s+1], v[3s+2]}
    trip_var = np.var(lv_row, axis=1, ddof=0)
    return float(np.mean(trip_var))

print("257b: Corrected d_k = ve0(k+1)/ve0(k) with correct CODE-variance")
print("="*70)

# === COMPARE CORRECT vs WRONG ===
print(f"\nComparison at lam=1.70:")
print(f"{'k':>4}  {'ve0_correct':>12}  {'ve0_wrong':>12}  {'d_k_correct':>12}  {'d_k_wrong':>10}")
prev_ve0_c = None
prev_ve0_w = None
lam = 1.70
for k in range(5, 15):
    v, Nl = run_kl(k, lam)
    ve0c = compute_ve0_correct(v, Nl)
    ve0w = compute_ve0_wrong(v, Nl)
    dk_c = ve0c / prev_ve0_c if prev_ve0_c is not None else float('nan')
    dk_w = ve0w / prev_ve0_w if prev_ve0_w is not None else float('nan')
    print(f"k={k:>2}  {ve0c:>12.6f}  {ve0w:>12.6f}  "
          f"{dk_c:>12.6f}  {dk_w:>10.6f}")
    prev_ve0_c = ve0c
    prev_ve0_w = ve0w
    sys.stdout.flush()

# === LAMBDA SCAN with correct d_k ===
print(f"\n\nLambda scan k=12->13 (correct formula):")
print(f"{'lam':>6}  {'ve0(12)':>10}  {'ve0(13)':>10}  {'d_12':>9}")
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    v12, Nl12 = run_kl(12, lam)
    v13, Nl13 = run_kl(13, lam)
    ve0_12 = compute_ve0_correct(v12, Nl12)
    ve0_13 = compute_ve0_correct(v13, Nl13)
    dk = ve0_13 / ve0_12
    print(f"lam={lam:.2f}  {ve0_12:>10.6f}  {ve0_13:>10.6f}  {dk:>9.6f}")
    sys.stdout.flush()

print("\ndone")
