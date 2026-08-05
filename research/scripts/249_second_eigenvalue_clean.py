"""
249_second_eigenvalue_clean.py
==============================
Clean computation of the second eigenvalue of the linearized K-L operator
at the Perron eigenvector, via deflated power iteration.

The linearized (Frechet derivative) K-L map at v* is:
  L[delta_v](i) = A * delta_v[T4(i)] + B * delta_v[argmin_i]
where argmin_i is the position achieving the min in cb[s(i)].

The argmin-frozen linearization: argmin is fixed at the argmin of v*.

Method:
  1. Converge to Perron eigenvector v* (power iteration)
  2. Deflate: project out the Perron component
  3. Apply linearized map repeatedly, track convergence ratio
  => |lambda_2 / lambda_1| = second eigenvalue ratio

Also: directly measure CODE-variance decay d_k = ve0(k+1)/ve0(k)
and compare d_k vs (lambda_2/lambda_1)^2.

Connection: CODE-variance is a quadratic functional of deviations from v*,
so if linear deviations decay at rate r = |lambda_2/lambda_1|,
then CODE-variance decays at rate r^2 = (lambda_2/lambda_1)^2.
Thus d_k -> (second eigenvalue ratio)^2 as k -> inf.
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)
N_ITER_CONVERGE = 800

def run_kl(k, lam, n_iter=N_ITER_CONVERGE):
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
        w = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    return v, Nl, T4, m0, m2, R1, R3, A, B1, B3

def apply_linearized(delta_v, v_star, Nl, T4, m0, m2, R1, R3, A, B1, B3):
    """Apply the argmin-frozen linearized K-L map to delta_v at v_star."""
    # Argmin of block min: for each s in [0,Nl), argmin is the index j in {s, s+Nl, s+2Nl}
    # where v_star[j] is minimized.
    v_blocks = np.stack([v_star[:Nl], v_star[Nl:2*Nl], v_star[2*Nl:]], axis=1)  # (Nl, 3)
    argmin_blocks = np.argmin(v_blocks, axis=1)  # (Nl,) values in {0,1,2}
    argmin_idx = np.arange(Nl, dtype=np.int64) + argmin_blocks * Nl  # (Nl,) indices into v

    # Linearized cb: delta_cb[s] = delta_v[argmin_idx[s]]
    delta_cb = delta_v[argmin_idx]  # (Nl,)

    # Apply linearized map
    delta_w = A * delta_v[T4]
    delta_w[m2] += B3 * delta_cb[R3[m2]]
    delta_w[m0] += B1 * delta_cb[R1[m0]]
    return delta_w

def deflate(delta_v, v_star):
    """Project out Perron component: delta_v -= <delta_v, v_star> / <v_star, v_star> * v_star."""
    coeff = float(np.dot(delta_v, v_star)) / float(np.dot(v_star, v_star))
    return delta_v - coeff * v_star


def code_var_ve0(v, Nl):
    """CODE-variance ve0 of the r=0 component."""
    v0 = v[0::3]  # length Nl
    Nl3 = Nl // 3
    a0 = v0[:Nl3]; a1 = v0[Nl3:2*Nl3]; a2 = v0[2*Nl3:]
    mean_a = (a0 + a1 + a2) / 3.0
    ld = np.log2(np.stack([a0, a1, a2]) / mean_a)
    return float(np.mean(ld**2))

print("249: Second eigenvalue of linearized K-L operator via deflated power iteration")
print("Also: compare d_k = ve0(k+1)/ve0(k) vs (lambda2/lambda1)^2")
print("="*75)

LAM = 1.70

for k in [6, 7, 8, 9]:
    print(f"\nk={k}, lam={LAM}:")
    sys.stdout.flush()

    v_star, Nl, T4, m0, m2, R1, R3, A, B1, B3 = run_kl(k, LAM)

    # Second eigenvalue via deflated power iteration
    N = len(v_star)
    rng = np.random.default_rng(42)
    delta = rng.standard_normal(N)
    delta = deflate(delta, v_star)
    delta /= float(np.linalg.norm(delta))

    ratios = []
    for it in range(300):
        delta_new = apply_linearized(delta, v_star, Nl, T4, m0, m2, R1, R3, A, B1, B3)
        delta_new = deflate(delta_new, v_star)
        norm_new = float(np.linalg.norm(delta_new))
        if norm_new < 1e-15:
            break
        ratios.append(norm_new)
        delta = delta_new / norm_new

    # Convergence ratio: ratios[n+1]/ratios[n] -> |lambda_2|/lambda_1 * (norm factors)
    # Actually, ratios[n] = ||L^n delta|| / ||L^{n-1} delta|| -> |lambda_2| / lambda_1
    if len(ratios) > 50:
        r_arr = np.array(ratios[50:])  # skip transient
        lam2_ratio = float(np.median(r_arr))
        print(f"  |lambda_2/lambda_1| = {lam2_ratio:.6f}  (from {len(r_arr)} ratios)")
        print(f"  Predicted d_k = (lambda_2/lambda_1)^2 = {lam2_ratio**2:.6f}")

    # Measure actual d_k = ve0(k+1)/ve0(k)
    ve0_k = code_var_ve0(v_star, Nl)
    v_kp1, Nl_kp1, *_ = run_kl(k+1, LAM)
    ve0_kp1 = code_var_ve0(v_kp1, Nl_kp1)
    d_k = ve0_kp1 / ve0_k
    print(f"  Actual d_k = ve0(k+1)/ve0(k) = {d_k:.6f}")
    print(f"  Ratio d_k / (lambda2)^2 = {d_k / lam2_ratio**2:.4f}")
    sys.stdout.flush()

print()
print("Lambda scan at k=8:")
print(f"{'lam':>6}  {'|lam2|':>10}  {'lam2^2':>10}  {'d_k':>10}  {'ratio':>8}")
for lam in [1.30, 1.50, 1.70, 1.90, 2.00]:
    v_star, Nl, T4, m0, m2, R1, R3, A, B1, B3 = run_kl(8, lam)
    N = len(v_star)
    rng = np.random.default_rng(42)
    delta = rng.standard_normal(N)
    delta = deflate(delta, v_star)
    delta /= float(np.linalg.norm(delta))
    ratios = []
    for it in range(200):
        delta_new = apply_linearized(delta, v_star, Nl, T4, m0, m2, R1, R3, A, B1, B3)
        delta_new = deflate(delta_new, v_star)
        norm_new = float(np.linalg.norm(delta_new))
        if norm_new < 1e-15: break
        ratios.append(norm_new)
        delta = delta_new / norm_new
    lam2_ratio = float(np.median(np.array(ratios[50:]))) if len(ratios) > 50 else float('nan')
    ve0_k = code_var_ve0(v_star, Nl)
    v_kp1, Nl_kp1, *_ = run_kl(9, lam)
    ve0_kp1 = code_var_ve0(v_kp1, Nl_kp1)
    d_k = ve0_kp1 / ve0_k
    print(f"lam={lam:.2f}  {lam2_ratio:>10.6f}  {lam2_ratio**2:>10.6f}  {d_k:>10.6f}  {d_k/lam2_ratio**2:>8.4f}")
    sys.stdout.flush()

print()
print("done")
