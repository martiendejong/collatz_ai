"""
256_depth_class_means.py
========================
Depth scan: does the class-mean asymmetry a0/a1 shrink or grow with k?
Does Corr(v2[2s+1], v2[s]) vanish as k -> inf, or stays bounded away from 0?

From Script 252 depth scan at lam=1.70:
  Corr(v2[3s+2], v2[s]) weakens: -0.879 (k=5) to -0.162 (k=11).
  But that is Corr at one SPECIFIC scale (3s+2 vs s), not the OVERALL Corr.

The overall Corr(v2[2s+1], v2[s]) over all s might stabilize.

KEY QUESTION: does the ANTI-CORRELATION converge to a nonzero negative limit as k -> inf?

Implications:
  - If Corr -> 0: the mechanism weakens and may not guarantee d_k < 1 for large k.
  - If Corr -> c < 0: robust mechanism, Conjecture G holds for all k.
  - If Corr -> 0 like 1/k or 1/sqrt(k): need more work.

Tests:
1. Depth scan k=5..14 for Corr(v2[2s+1], v2[s]) at lam=1.70
2. Depth scan for a0/a1 and CV0/CV1 (class-mean and CV ratios)
3. Check whether the BETWEEN-CLASS contribution decays faster than the total
4. Lambda scan to see if there's a lambda where Corr converges to 0
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)
N_ITER_BASE = 600

def run_kl(k, lam, n_iter=None):
    if n_iter is None:
        n_iter = N_ITER_BASE + 100 * max(0, k - 8)
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

def compute_stats(v, Nl):
    v2 = v[2::3]
    s = np.arange(Nl, dtype=np.int64)

    # Class means and CVs
    a = np.array([float(np.mean(v2[s%3==r])) for r in range(3)])
    std = np.array([float(np.std(v2[s%3==r])) for r in range(3)])
    cv = std / a

    # Between-class corr (analytic)
    mean_a = a.mean()
    var_btwn = np.mean(a**2) - mean_a**2
    cov_btwn = (a[0]*a[1] + a[1]*a[0] + a[2]*a[2])/3 - mean_a**2
    corr_btwn = cov_btwn / var_btwn if var_btwn > 0 else 0.0

    # Overall doubling anti-corr (using small subset to avoid boundary effects)
    max_s = min((Nl - 9) // 9, Nl // 3)
    s0 = np.arange(max(max_s, 1), dtype=np.int64)
    corr_actual = float(np.corrcoef(v2[s0], v2[(2*s0+1)%Nl])[0,1])

    # Within class 2->2 corr
    s2 = s0[s0 % 3 == 2]
    if len(s2) > 2:
        corr_22 = float(np.corrcoef(v2[s2], v2[(2*s2+1)%Nl])[0,1])
    else:
        corr_22 = float('nan')

    # Corr(v2[3s+2], v2[s]) (for comparison with Script 252 depth scan)
    s_p = np.arange(min(Nl//9, max_s), dtype=np.int64)
    if len(s_p) > 2:
        corr_01 = float(np.corrcoef(v2[3*s_p+2], v2[s_p])[0,1])
    else:
        corr_01 = float('nan')

    return a, cv, corr_actual, corr_btwn, corr_22, corr_01

print("256: Depth scan -- class-mean asymmetry and anti-correlation vs k")
print("="*70)

# === DEPTH SCAN lam=1.70 ===
print(f"\nDepth scan lam=1.70:")
print(f"{'k':>4}  {'Nl':>6}  {'a0':>7}  {'a1':>7}  {'a0/a1':>7}  "
      f"{'CV0':>6}  {'CV1':>6}  {'corr_act':>9}  {'corr_btwn':>10}  {'corr_01':>8}  {'corr_22':>8}")
lam = 1.70
for k in range(5, 15):
    v, Nl = run_kl(k, lam)
    a, cv, c_act, c_btwn, c22, c01 = compute_stats(v, Nl)
    print(f"k={k:>2}  {Nl:>6}  {a[0]:>7.4f}  {a[1]:>7.4f}  {a[0]/a[1]:>7.4f}  "
          f"{cv[0]:>6.4f}  {cv[1]:>6.4f}  {c_act:>9.4f}  {c_btwn:>10.4f}  "
          f"{c01:>8.4f}  {c22:>8.4f}")
    sys.stdout.flush()

# === LAMBDA SCAN k=12 ===
print(f"\n\nLambda scan k=12:")
print(f"{'lam':>6}  {'a0':>7}  {'a1':>7}  {'a0/a1':>7}  "
      f"{'CV0':>6}  {'CV1':>6}  {'corr_act':>9}  {'corr_btwn':>10}")
k = 12
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    v, Nl = run_kl(k, lam)
    a, cv, c_act, c_btwn, c22, c01 = compute_stats(v, Nl)
    print(f"lam={lam:.2f}  {a[0]:>7.4f}  {a[1]:>7.4f}  {a[0]/a[1]:>7.4f}  "
          f"{cv[0]:>6.4f}  {cv[1]:>6.4f}  {c_act:>9.4f}  {c_btwn:>10.4f}")
    sys.stdout.flush()

# === CONVERGENCE RATE OF Corr vs k ===
print(f"\n\nConvergence of Corr(v2[2s+1], v2[s]) with k (lam=1.70):")
print(f"{'k':>4}  {'corr':>9}  {'ratio':>8}")
prev = None
lam = 1.70
for k in range(6, 15):
    v, Nl = run_kl(k, lam)
    a, cv, c_act, c_btwn, c22, c01 = compute_stats(v, Nl)
    ratio = c_act / prev if prev is not None else float('nan')
    print(f"k={k:>2}  {c_act:>9.4f}  {ratio:>8.4f}")
    prev = c_act
    sys.stdout.flush()

print("\ndone")
