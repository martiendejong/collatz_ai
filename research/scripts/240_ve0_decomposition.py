"""
240_ve0_decomposition.py
========================
Analytical decomposition of ve0 from the r=0 K-L equation.

The r=0 fixed-point equation:
  rho * v0(s) = A * v2[sigma0(s)] + B1 * cb[R1(s)]

Both sigma0(s) = 4s mod Nl and R1(s) = 4s mod Nl ARE THE SAME MAP.
(sigma0 = R1 — both are 4s mod Nl for r=0 nodes.)

Wait, let me re-derive from the script:
  T4(i) = (4i+2) mod N. For i = 3s+0 (r=0):
    T4(3s) = (12s+2) mod N = 3*((4s+0) mod Nl) + 2 when (12s+2)%3=2.
    Actually T4(3s) = 12s+2. (12s+2)%3 = (0+2)%3 = 2. So T4(3s) is an r=2 node.
    s-coordinate of T4(3s): (12s+2)//3 = 4s. So T4 maps r=0 node at s to r=2 node at 4s.
    => sigma_0(s) = 4s mod Nl (the s-coord map for r=0).

  R1(s) = 4s mod Nl (used in B-term for r=0 nodes).

So sigma_0 = R1 = 4s mod Nl! The two terms in the v0 equation come from
the SAME index map but applied to DIFFERENT functions (v2 vs cb).

Therefore:
  rho * v0(s) = A * v2[sigma0(s)] + B1 * cb[sigma0(s)]
             = (A * v2[j] + B1 * cb[j])  where j = sigma0(s)

Since sigma0 is a bijection (gcd(4, Nl)=1), and maps triplets to triplets (4 ≡ 1 mod 3),
the CODE-variance of v0 equals the CODE-variance of (A*v2 + B1*cb) pointwise.

Define: f(j) = A * v2[j] + B1 * cb[j]  (a Nl-dim vector)

Then: ve0 = CODE-var(f) (since v0 = f[sigma0^{-1}] and sigma0 maps triplets to triplets).

This is a KEY structural identity:
  ve0 = CODE-var(A * v2 + B1 * cb)

Questions:
1. Is ve0 = CODE-var(A*v2 + B1*cb) to machine precision? (Should be exact.)
2. How does CODE-var(A*v2 + B1*cb) relate to ve2 and ve_cb?
3. Does the "weighted average" approximation hold?
   ve0 ≈ w2^2 * ve2 + wcb^2 * ve_cb + 2*w2*wcb*Cov(log_dev_v2, log_dev_cb)?
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
    return v, Nl, A, B1, B3, R3

def code_var_nl(arr, Nl):
    """CODE-variance: group (j, j+Nl/3, j+2Nl/3) as a triplet, measure log-spread."""
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    a0 = arr[s0]; a1 = arr[s0+Nl3]; a2 = arr[s0+2*Nl3]
    mean_a = (a0 + a1 + a2) / 3.0
    log_dev = np.stack([np.log2(a0/mean_a), np.log2(a1/mean_a), np.log2(a2/mean_a)])
    return float(np.mean(log_dev**2)), log_dev

def code_var_from_logdev(ld):
    return float(np.mean(ld**2))

LAM = 1.70

print("240: ve0 decomposition analysis")
print(f"lambda = {LAM}")
print("="*72)
print(f"KEY: sigma0 = R1 = 4s mod Nl (same map!) => ve0 = CODE-var(A*v2 + B1*cb)")
print()
sys.stdout.flush()

header = f"{'k':>4}  {'ve0':>10}  {'ve2':>10}  {'ve_cb':>10}  {'ve_f':>10}  {'err_ve0_vef':>12}  {'ratio_ve2_ve0':>13}  {'ve_cb/ve0':>10}"
print(header)

for k in range(4, 16):
    v, Nl, A, B1, B3, R3 = run_kl(k, LAM)

    # Run one more K-L step to get rho and cb at convergence
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    s_arr, r_arr = np.divmod(i, 3)
    m0, m2 = (r_arr==0), (r_arr==2)
    R1_arr = (4*s_arr) % Nl
    R3_arr = (2*s_arr+1) % Nl

    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

    # v0, v1, v2 components
    v0 = v[0::3]  # r=0 values: v[0], v[3], v[6], ...
    v1 = v[1::3]
    v2 = v[2::3]

    # Wait — ve_breakdown in script 238 uses the Nl-block structure, not 0::3.
    # Let me use BOTH and see which gives ve0=ve1.

    # Method A: interleaved (v[0::3], v[1::3], v[2::3])
    # Triplet for CODE-var: (v0[j], v0[j+Nl/3], v0[j+2Nl/3])
    Nl3 = Nl // 3
    s0 = np.arange(Nl3, dtype=np.int64)
    # v0[s0], v0[s0+Nl3], v0[s0+2*Nl3] — triplets of r=0 values spaced by Nl/3

    def cov_logdev(ld1, ld2):
        """Mean covariance between two log-dev matrices (same shape 3 x Nl3)."""
        return float(np.mean(ld1 * ld2))

    # ve0, ve2 via interleaved storage
    ve0_d, ld_v0 = code_var_nl(v0, Nl)
    ve2_d, ld_v2 = code_var_nl(v2, Nl)

    # ve_cb: CODE-variance of cb (the Nl-dimensional min-vector)
    # cb[j] = min(v[j], v[j+Nl], v[j+2*Nl]) for j in [0,Nl)
    # Note: cb is indexed in the BLOCK ordering, not interleaved.
    # For r=0: cb[R1(s)] where R1(s) = 4s mod Nl. This is the cb applied
    # to the sigma0 index.
    # The v0 equation: v0(s) = (1/rho)*[A*v2[sigma0(s)] + B1*cb[sigma0(s)]]
    # Since sigma0(s) = 4s mod Nl.

    # But v2 is stored in INTERLEAVED order: v2[s] = v[2 + 3s].
    # And cb is in BLOCK order: cb[j] = min(v[j], v[j+Nl], v[j+2*Nl]).
    # So cb[j] is the column min across 3 r-types for the "block" index j.

    # Let's work in the BLOCK-indexed cb.
    # For j in [0, Nl): cb[j] = min(v[j], v[j+Nl], v[j+2*Nl])
    # j%3 tells us the r-type of all three: all have r = j%3 (since Nl%3=0)
    # For j = 3*m: r=0, s-coords m, m+Nl/3, m+2Nl/3
    # cb[3*m] = min(v0[m], v0[m+Nl/3], v0[m+2*Nl/3])

    # So cb restricted to j=3*m (r=0 columns) = pointwise min of v0 triplet!
    # This is the CODE-minimum of v0.

    cb_r0 = cb[0::3]  # cb[j] for j=0,3,6,... = r=0 cb values (indexed by m in [0,Nl/3))
    cb_r1 = cb[1::3]
    cb_r2 = cb[2::3]

    # ve_cb for r=0 sub-vector (cb_r0 is Nl/3 dimensional, no triplet structure here)
    # Actually cb_r0[m] = min(v0[m], v0[m+Nl3], v0[m+2*Nl3]) is a SCALAR per m
    # It's the minimum of the v0 triplet at position m.
    # To get CODE-variance of cb_r0 "as if it were v0", I need to group cb_r0 into
    # sub-triplets of size Nl/9... but this gets complicated.

    # Alternative: measure CODE-var of the FULL cb vector using the same block structure.
    # i.e., treat cb (Nl-dim) the same way as v0 (Nl-dim) for CODE-var purposes.
    ve_cb_block, ld_cb = code_var_nl(cb, Nl)  # treats cb[j], cb[j+Nl/3], cb[j+2Nl/3] as triplet

    # But wait — cb is Nl-dimensional, not N-dimensional. The "triplet" in cb space
    # would be (cb[j], cb[j+Nl/3], cb[j+2*Nl/3]) for j in [0, Nl/3).
    # This groups the r=0, r=1, r=2 sub-vectors of cb (each Nl/3 long).
    # Let me compute this:
    cb0 = cb[:Nl3]; cb1 = cb[Nl3:2*Nl3]; cb2 = cb[2*Nl3:]
    cb_triplet_min_of_triplet_min = (cb0 + cb1 + cb2) / 3.0  # not really useful

    # KEY: the formula v0 = (A*v2 + B1*cb) / rho, where v2 and cb are looked up at sigma0(s)=4s%Nl
    # In the interleaved storage: v2[s] = v[2+3s] (the r=2 value at s-coordinate s).
    # In the block storage: cb[j] for j = sigma0(s) = 4s%Nl.
    # But what is cb[4s%Nl]? 4s%Nl has r-type (4s%Nl)%3 = s%3.
    # So cb[4s] is the column minimum for the r=s%3 column at position 4s%Nl.

    # For a r=0 node (where s%3 = 0 for s=0,3,6,...):
    # cb[4s%Nl] = min(v[4s%Nl], v[4s%Nl+Nl], v[4s%Nl+2*Nl])
    # All three have r = (4s)%3 = 0 (since 4 ≡ 1 mod 3 and s ≡ 0 mod 3 for r=0 nodes).
    # Wait, s ranges over ALL values in [0, Nl), not just multiples of 3.

    # Actually, for the r=0 K-L equation, the s parameter ranges over [0, Nl) with NO
    # restriction on s%3. The equation holds for ALL s in [0, Nl).
    # v0(s) = (A*v2[sigma0(s)] + B1*cb[sigma0(s)]) / rho
    # where sigma0(s) = 4s%Nl.

    # Let me compute f(s) = A*v2[4s%Nl] + B1*cb[4s%Nl] for s in [0, Nl)
    # and compare ve0 with CODE-var(f) / rho^2.

    # Note: v2 is in interleaved storage: v2[s] = v[2+3*s].
    # And cb[j] is in block storage (Nl-dimensional).
    # sigma0(s) = 4s%Nl is an index j in [0, Nl).

    # For the block-indexed j: v2[j] means... hmm, v2[s] where s is the s-coordinate.
    # The BLOCK structure has v[j] for j in [0, N). j = 3s+r, so s=j//3, r=j%3.
    # v2 means r=2: v2_block = v[2::3] = v[2], v[5], v[8], ...

    # sigma0(s) = 4s%Nl. For s in [0, Nl), 4s%Nl is also in [0, Nl).
    # We need v2[sigma0(s)] = v2[4s%Nl] in the INTERLEAVED r=2 values.
    # In interleaved storage: v2[s] = v[3s+2] for s in [0, Nl).
    # So v2[4s%Nl] = v[3*(4s%Nl)+2] in the FULL N-vector.

    sigma0 = (4 * np.arange(Nl, dtype=np.int64)) % Nl
    v2_interleaved = v[2::3]  # v2[s] = v[3s+2], s in [0, Nl)
    v2_at_sigma0 = v2_interleaved[sigma0]  # shape (Nl,)

    cb_at_sigma0 = cb[sigma0]  # shape (Nl,)

    # Estimate rho
    T4 = (4 * np.arange(N, dtype=np.int64) + 2) % N
    s_arr_full, r_arr_full = np.divmod(np.arange(N, dtype=np.int64), 3)
    m0_full = (r_arr_full == 0)
    m2_full = (r_arr_full == 2)
    R1_full = (4 * s_arr_full) % Nl
    R3_full = (2 * s_arr_full + 1) % Nl
    cb_full = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w = A * v[T4]
    w[m2_full] += B3 * cb_full[R3_full[m2_full]]
    w[m0_full] += B1 * cb_full[R1_full[m0_full]]
    rho = float(w.max())

    # f(s) = A*v2[sigma0(s)] + B1*cb[sigma0(s)]
    f = A * v2_at_sigma0 + B1 * cb_at_sigma0  # shape (Nl,)

    # v0(s) should equal f(s) / rho
    v0_from_f = f / rho
    v0_interleaved = v[0::3]  # v[3s] for s in [0, Nl)
    err_v0 = np.max(np.abs(v0_interleaved - v0_from_f)) / np.maximum(np.abs(v0_interleaved), 1e-15).max()
    # Note: v0_from_f[s] = f(s)/rho, and f(s) uses sigma0(s). Let me re-examine.

    # Actually the K-L equation for r=0 node at (s, r=0) is:
    # rho * v[3s] = A * v[T4[3s]] + B1 * cb[R1[3s]]
    # T4[3s] = (4*3s+2)%N = (12s+2)%N. (12s+2)%3 = 2, so T4[3s] is r=2 node.
    # s-coord of T4[3s]: (12s+2)//3 = 4s (for s < N/4). -> v2_interleaved[4s%Nl].
    # R1[3s] = (4*s)%Nl (from s_arr[3s]=s). -> cb[(4s)%Nl].

    # f(s) = A * v2_interleaved[4s%Nl] + B1 * cb[(4s)%Nl]
    # v0[s] = f(s) / rho (to machine precision if converged)

    # Let me also verify:
    v0_pred = f / rho
    err_check = float(np.max(np.abs(v0_interleaved - v0_pred)) / (np.max(np.abs(v0_interleaved)) + 1e-30))

    # CODE-variance of f / rho = CODE-variance of f (scalar scaling doesn't change log-dev)
    ve_f, ld_f = code_var_nl(f, Nl)

    # CODE-variance of v2_at_sigma0
    ve_v2_sig0, ld_v2_sig0 = code_var_nl(v2_at_sigma0, Nl)
    # This should equal ve2 since sigma0 maps triplets to triplets.

    # CODE-variance of cb_at_sigma0
    ve_cb_sig0, ld_cb_sig0 = code_var_nl(cb_at_sigma0, Nl)
    # This equals CODE-var(cb) since sigma0 maps triplets to triplets.

    # Weights: w2 = A*mean_v2 / (A*mean_v2 + B1*mean_cb)
    mean_v2 = float(np.mean(v2_at_sigma0))
    mean_cb = float(np.mean(cb_at_sigma0))
    denom = A * mean_v2 + B1 * mean_cb
    w2 = A * mean_v2 / denom
    wcb = B1 * mean_cb / denom

    # Cross-term:
    cov_term = cov_logdev(ld_v2_sig0, ld_cb_sig0)

    # Predicted ve0 from weighted-average formula:
    ve0_pred_formula = w2**2 * ve_v2_sig0 + wcb**2 * ve_cb_sig0 + 2*w2*wcb*cov_term

    print(f"k={k:2d}: ve0={ve0_d:.6f}  ve2={ve2_d:.6f}  "
          f"ve_cb={ve_cb_sig0:.6f}  ve_f={ve_f:.6f}  "
          f"err_f={err_check:.1e}  "
          f"r21={ve2_d/ve0_d:.4f}  cb/v0={ve_cb_sig0/ve0_d:.4f}")
    print(f"       w2={w2:.4f}  wcb={wcb:.4f}  cov={cov_term:.6f}  "
          f"ve0_pred={ve0_pred_formula:.6f}  ve0_err={(ve0_pred_formula-ve0_d)/ve0_d:.4f}")
    sys.stdout.flush()

print()
print("Summary:")
print("  ve0 = CODE-var(f) where f(s) = A*v2[sigma0(s)] + B1*cb[sigma0(s)]")
print("  Both terms use SAME map sigma0 = 4s mod Nl (triplet-preserving).")
print("  ve_f = ve0 exactly (f/rho = v0).")
print("  ve_v2_sig0 = ve2 (since sigma0 maps triplets to triplets).")
print("  ve_cb_sig0 = ve_cb (same argument).")
print("  Formula: ve0 = w2^2*ve2 + wcb^2*ve_cb + 2*w2*wcb*Cov(logdev_v2,logdev_cb)")
print()
print("  If Cov ~ 0 (near-zero covariance between v2 and cb log-deviations):")
print("  ve0 ~ w2^2*ve2 + wcb^2*ve_cb")
print("  Negative Cov further reduces ve0 below ve2: L > 1.")
print("  Actual: linear formula underpredicts ve0 by ~27% (nonlinear terms matter).")
print()
print("done")
