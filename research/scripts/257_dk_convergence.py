"""
257_dk_convergence.py
=====================
Directly measure d_k = ve0(k+1)/ve0(k) for k=5..14.
The CENTRAL QUESTION: does d_k converge to d_inf < 1 (strong Conj G),
or does d_k -> 1 (weak form: d_k < 1 for all finite k but margin vanishes)?

d_k is the CODE-variance ratio of the K-L Perron eigenvector:
  ve0(k) = CODE-variance at depth k = within-triplet variance of log(v^(k))

From earlier measurements (Scripts 201, 240, etc.):
  d_k ~= 0.740-0.742 for k=6..11 at lam=1.70.
  Was this EXACTLY constant or drifting?

If d_k -> d_inf < 1: strong Conjecture G, and the CODE-variance decays geometrically.
If d_k -> 1: weak form, the margin 1-d_k vanishes.

METHOD: Run K-L iteration to high convergence at depth k and k+1, measure ve0 both.
ve0(k) = Var(log v_k | within triplet) = E[(log v_k[i] - log v_k[i*])^2] for triplet pairs.

The CODE-variance is the WITHIN-s-CODE-TRIPLET variance of log(v^(k)).
A s-CODE-TRIPLET at depth k is {(s,0), (s,1), (s,2)} for s in [0, N_l).
ve0(k) = (1/N_l) sum_s Var(log v^(k)[3s+r] for r=0,1,2)

This is computed from the K-L eigenvector at depth k.
"""
import numpy as np
from math import log2, log
import sys

ALPHA = log2(3.0)

def run_kl(k, lam, n_iter=None):
    """Run K-L iteration, return converged v and Nl."""
    if n_iter is None:
        n_iter = 500 + 100 * max(0, k - 8)
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
    """Compute CODE-variance = within-s-triplet variance of log(v)."""
    # s-CODE-triplet at index i = (s, r): {(s,0), (s,1), (s,2)} for each s
    # In v: v[3s+0], v[3s+1], v[3s+2] for s = 0..Nl-1
    lv = np.log(v)
    # Reshape to (Nl, 3) — each row is one s-triplet
    lv_3 = lv.reshape(Nl, 3)
    # Within-triplet variance for each s: Var({lv[3s+r] : r=0,1,2})
    trip_var = np.var(lv_3, axis=1, ddof=0)
    # ve0 = mean over all s-triplets
    ve0 = float(np.mean(trip_var))
    return ve0

print("257: d_k = ve0(k+1)/ve0(k) convergence test")
print("="*70)

# === MAIN DEPTH SCAN ===
print(f"\nDepth scan: d_k at lam=1.70 and lam=2.00")
print(f"{'k':>4}  {'ve0_k (1.70)':>13}  {'d_k (1.70)':>11}  {'ve0_k (2.00)':>13}  {'d_k (2.00)':>11}")

lams = [1.70, 2.00]
prev_ve0 = {lam: None for lam in lams}

for k in range(5, 15):
    ve0s = {}
    for lam in lams:
        v, Nl = run_kl(k, lam)
        ve0s[lam] = compute_ve0(v, Nl)
        sys.stdout.flush()

    row = f"k={k:>2}  "
    for lam in lams:
        row += f"{ve0s[lam]:>13.6f}  "
        if prev_ve0[lam] is not None:
            dk = ve0s[lam] / prev_ve0[lam]
            row += f"{dk:>11.6f}  "
        else:
            row += f"{'---':>11}  "

    print(row)
    for lam in lams:
        prev_ve0[lam] = ve0s[lam]
    sys.stdout.flush()

# === LAMBDA SCAN AT k=12 AND k=13 ===
print(f"\n\nLambda scan: d_12 = ve0(12+1)/ve0(12)")
print(f"{'lam':>6}  {'ve0(12)':>13}  {'ve0(13)':>13}  {'d_12':>9}")
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    v12, Nl12 = run_kl(12, lam)
    v13, Nl13 = run_kl(13, lam)
    ve0_12 = compute_ve0(v12, Nl12)
    ve0_13 = compute_ve0(v13, Nl13)
    dk = ve0_13 / ve0_12
    print(f"lam={lam:.2f}  {ve0_12:>13.6f}  {ve0_13:>13.6f}  {dk:>9.6f}")
    sys.stdout.flush()

# === CONVERGENCE ANALYSIS ===
print(f"\n\nConvergence analysis: d_k as k -> inf (lam=1.70)")
print("Computing d_k for k=6..14 with 2*convergence iterations:")
prev_ve0 = None
lam = 1.70
for k in range(6, 15):
    v, Nl = run_kl(k, lam, n_iter=1200)
    ve0 = compute_ve0(v, Nl)
    if prev_ve0 is not None:
        dk = ve0 / prev_ve0
        print(f"k={k:>2}: ve0={ve0:.6f}  d_{k-1}={dk:.6f}")
    else:
        print(f"k={k:>2}: ve0={ve0:.6f}  (reference)")
    prev_ve0 = ve0
    sys.stdout.flush()

print("\ndone")
