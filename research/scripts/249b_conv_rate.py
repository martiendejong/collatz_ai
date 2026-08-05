"""
249b_conv_rate.py
=================
Measure the actual convergence rate of power iteration toward the Perron eigenvector,
by tracking ||v_n - v_inf|| / ||v_{n-1} - v_inf||.

This gives the true |lambda_2/lambda_1| without needing the left eigenvector.

Compare to sqrt(d_k) where d_k = ve0(k+1)/ve0(k).
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)

def run_kl_with_history(k, lam, n_warmup=800, n_measure=200):
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

    def apply_F(v):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        return w / w.max()

    # Warmup to get v_inf
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_warmup):
        v = apply_F(v)
    v_inf = v.copy()

    # Measure convergence: perturb and track
    # Actually: run fresh from ones and track convergence to v_inf
    v2 = np.ones(N, dtype=np.float64)
    # Additional warmup (v2 and v_inf both converge, but v2 from scratch)
    # Use a perturbed start to see convergence
    rng = np.random.default_rng(0)
    v2 = v_inf + 0.01 * rng.standard_normal(N)
    v2 = v2 / v2.max()

    ratios = []
    prev_diff = None
    for _ in range(n_measure):
        v2_new = apply_F(v2)
        diff = float(np.linalg.norm(v2_new - v_inf))
        if prev_diff is not None and prev_diff > 1e-15 and diff > 1e-15:
            ratios.append(diff / prev_diff)
        prev_diff = diff
        v2 = v2_new

    return v_inf, Nl, ratios

def code_var_ve0(v, Nl):
    v0 = v[0::3]
    Nl3 = Nl // 3
    a0 = v0[:Nl3]; a1 = v0[Nl3:2*Nl3]; a2 = v0[2*Nl3:]
    mean_a = (a0 + a1 + a2) / 3.0
    ld = np.log2(np.stack([a0, a1, a2]) / mean_a)
    return float(np.mean(ld**2))

print("249b: Convergence rate of power iteration = |lambda_2/lambda_1|")
print("Compare sqrt(d_k) from CODE-variance to conv_rate")
print("="*65)

# Depth scan at lam=1.70
print(f"\nDepth scan, lam=1.70:")
print(f"{'k':>4}  {'conv_rate':>12}  {'d_k':>10}  {'sqrt(d_k)':>10}  {'match?':>8}")
LAM = 1.70
prev_ve0 = None
for k in range(6, 12):
    v_inf, Nl, ratios = run_kl_with_history(k, LAM, n_warmup=600, n_measure=150)
    if len(ratios) > 20:
        conv_rate = float(np.median(np.array(ratios[10:])))
    else:
        conv_rate = float('nan')

    ve0 = code_var_ve0(v_inf, Nl)
    d_k = ve0 / prev_ve0 if prev_ve0 is not None else float('nan')
    prev_ve0 = ve0

    print(f"k={k:>2}  {conv_rate:>12.6f}  {d_k:>10.6f}  {d_k**0.5:>10.6f}  "
          f"{'OK' if abs(conv_rate - d_k**0.5) < 0.05 else 'DIFF':>8}")
    sys.stdout.flush()

# Lambda scan at k=8
print(f"\nLambda scan, k=8:")
print(f"{'lam':>6}  {'conv_rate':>12}  {'d_k(8->9)':>12}  {'sqrt(d_k)':>10}  {'ratio':>8}")
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    v8, Nl8, ratios8 = run_kl_with_history(8, lam, n_warmup=600, n_measure=150)
    v9, Nl9, _ = run_kl_with_history(9, lam, n_warmup=600, n_measure=50)
    ve0_8 = code_var_ve0(v8, Nl8)
    ve0_9 = code_var_ve0(v9, Nl9)
    d_k = ve0_9 / ve0_8
    conv_rate = float(np.median(np.array(ratios8[10:]))) if len(ratios8) > 10 else float('nan')
    ratio = conv_rate / d_k**0.5 if d_k > 0 else float('nan')
    print(f"lam={lam:.2f}  {conv_rate:>12.6f}  {d_k:>12.6f}  {d_k**0.5:>10.6f}  {ratio:>8.4f}")
    sys.stdout.flush()

print()
print("done")
