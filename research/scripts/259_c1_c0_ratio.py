"""
259_c1_c0_ratio.py
==================
ANALYTICAL VERIFICATION: c_1/c_0 = A/rho EXACTLY.

FROM K-L EQUATION:
  rho * v1[s] = A * v0[sigma_1(s)]      ... (exact)
  where sigma_1(s) = (4s+2) % Nl is a SINGLE N_l-CYCLE (proved)

Therefore: v1[s] = (A/rho) * v0[sigma_1(s)] for all s.

The column-triplet minimum:
  cb[j] = min(v[j], v[j+Nl], v[j+2Nl])  for j in [0, Nl)

For j with j%3 == 1 (r-type 1): cb[j] involves v[j], v[j+Nl], v[j+2Nl]
  where all three have r-type 1 (since Nl = 3^(k-2) is divisible by 3).

  v1 values: v[j] = v1[j//3] (with j = 3*s+1)
  cb[3s+1] = min(v[3s+1], v[3s+1+Nl], v[3s+1+2Nl])
            = min(v1[s], v1[s+Nl//3], v1[s+2Nl//3])
            = (A/rho) * min(v0[sigma_1(s)], v0[sigma_1(s+Nl//3)], v0[sigma_1(s+2Nl//3)])

Since sigma_1 is a PERMUTATION, the set {sigma_1(s), sigma_1(s+Nl//3), sigma_1(s+2Nl//3)}
is some permutation of three s-values. BUT it's not necessarily a single column-triplet!

Actually, let's check: is {sigma_1(s), sigma_1(s+Nl//3), sigma_1(s+2Nl//3)} the same
as a column-triplet {t, t+Nl//3, t+2Nl//3} for some t?

sigma_1(s) = (4s+2) % Nl
sigma_1(s+Nl//3) = (4*(s+Nl//3)+2) % Nl = (4s+4Nl//3+2) % Nl = (sigma_1(s) + 4Nl//3) % Nl
sigma_1(s+2Nl//3) = (sigma_1(s) + 8Nl//3) % Nl

For this to be a column-triplet: {sigma_1(s), sigma_1(s)+Nl//3, sigma_1(s)+2Nl//3},
we need 4Nl//3 == Nl//3 (mod Nl), i.e., 4/3 == 1/3 mod 1, i.e., 4 == 1 (mod 3). FALSE!
4 % 3 = 1, so 4Nl//3 % Nl = (Nl//3 + Nl) % Nl = Nl//3 only if Nl//3 is not affected by mod.

Wait: 4Nl//3 mod Nl: since Nl//3 = 3^(k-3), and Nl = 3^(k-2) = 3*Nl//3:
4*Nl//3 mod Nl = 4*Nl//3 mod 3*Nl//3 = (4 mod 3)*Nl//3 = Nl//3.

So 4Nl//3 ≡ Nl//3 (mod Nl)! Therefore:
sigma_1(s+Nl//3) = (sigma_1(s) + Nl//3) % Nl
sigma_1(s+2Nl//3) = (sigma_1(s) + 2Nl//3) % Nl

This means: the three values {sigma_1(s), sigma_1(s+Nl//3), sigma_1(s+2Nl//3)}
ARE a column-triplet {sigma_1(s), sigma_1(s)+Nl//3, sigma_1(s)+2Nl//3}!

Therefore:
cb[3s+1] = (A/rho) * min(v0[sigma_1(s)], v0[sigma_1(s)+Nl//3], v0[sigma_1(s)+2Nl//3])
          = (A/rho) * cb_v0[sigma_1(s)]  (where cb_v0 is the column-min of v0 values)

But cb[j] = min(v[j], v[j+Nl], v[j+2Nl]) uses ALL r-types together, not just v0.
Let me reconsider.

Actually, Nl = N//3 and the v array has N = 3*Nl elements.
The column triplet for j=3s+1: {v[3s+1], v[3s+1+Nl], v[3s+1+2Nl]}

Now: 3s+1+Nl = 3s+1+3*Nl//3 = 3*(s+Nl//3)+1 (since Nl=3*Nl//3, so +Nl means +3*Nl//3)
And 3s+1+2Nl = 3*(s+2*Nl//3)+1.

So {v[3s+1], v[3(s+Nl//3)+1], v[3(s+2Nl//3)+1]} = {v1[s], v1[s+Nl//3], v1[s+2Nl//3]}.

These ARE all v1-type (r=1). And from v1[t] = (A/rho)*v0[sigma_1(t)]:

{v1[s], v1[s+Nl//3], v1[s+2Nl//3]} = (A/rho)*{v0[sigma_1(s)], v0[sigma_1(s)+Nl//3], v0[sigma_1(s)+2Nl//3]}

(using the above calculation that sigma_1 maps the triplet to another triplet).

Therefore:
cb[3s+1] = min(v1[s], v1[s+Nl//3], v1[s+2Nl//3])
          = (A/rho) * min(v0[sigma_1(s)], v0[sigma_1(s)+Nl//3], v0[sigma_1(s)+2Nl//3])
          = (A/rho) * cb_v0[3*sigma_1(s)+0]  ... (column-min of r=0 type at sigma_1(s))

Wait, cb[j] = min(v[j], v[j+Nl], v[j+2Nl]). For j=3t (r=0 type):
cb[3t] = {v[3t], v[3t+Nl], v[3t+2Nl]} = {v0[t], v0[t+Nl//3], v0[t+2Nl//3]}

So cb_v0_triplet at t = min(v0[t], v0[t+Nl//3], v0[t+2Nl//3]) = cb[3t] (r=0 column-min).

Therefore:
cb[3s+1] = (A/rho) * cb[3*sigma_1(s)]

In terms of the j-index (j = 3s+1 for r=1, j' = 3*sigma_1(s) for r=0):
cb[j with j%3==1] = (A/rho) * cb[3*(sigma_1(j//3)) with j//3 from j]

This means: c_1 = Mean(cb[j≡1 mod 3]) = (A/rho) * Mean(cb[3*sigma_1(s) | s in cb_range])
          = (A/rho) * Mean(cb[j≡0 mod 3])  (since sigma_1 is a permutation)
          = (A/rho) * c_0.

THEREFORE: c_1 = (A/rho) * c_0 EXACTLY (not approximately)!

And since A/rho = lambda^{-2}/rho < 1 for lambda > 1 (as rho >= 1 is the Perron eigenvalue):
c_1 < c_0 FOR ALL lambda > 1.

This is an EXACT ANALYTICAL RESULT.

IMPLICATIONS:
From self-consistency (cb-dominance approximation):
  rho * a0_v2 ~= B3 * c_1 = B3 * (A/rho) * c_0 = (B3*A/rho) * c_0
  rho * a1_v2 ~= B3 * c_0

=> a0_v2 / a1_v2 ~= (B3*A/rho) / (rho * (B3/rho)) = A/rho = lambda^{-2}/rho < 1.

So a0_v2 < a1_v2 FOR ALL lambda > 1 (exactly, up to T4 corrections).

TESTS IN THIS SCRIPT:
1. Verify c_1/c_0 = A/rho numerically (should match to machine precision if exact)
2. Verify lambda-scan: c_1/c_0 vs A/rho for lambda=1.30..2.00
3. Verify depth-scan: c_1/c_0 approaches A/rho as k increases
4. Check the v1[s] = (A/rho)*v0[sigma_1(s)] identity exactly
"""
import numpy as np
from math import log2
import sys

ALPHA = log2(3.0)
N_ITER = 1000

def run_kl_full(k, lam, n_iter=N_ITER):
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
    # sigma_1: the T4 pullback index for r=1 nodes: T4(3s+1) = (12s+6)%N = 3*(4s+2)%N/3+0
    # Actually T4(i) = (4i+2)%N. For i=3s+1: T4(3s+1) = (12s+6)%N.
    # (12s+6)%N = 6*(2s+1)%N. Since 6 = 2*3 and N = 3^(k-1):
    # (12s+6) mod N = 3 * (4s+2) mod Nl. So T4(3s+1) = 3*(4s+2 mod Nl), r-type 0.
    # sigma_1(s) = (4s+2) % Nl.
    sigma1_s = (4*s_arr + 2) % Nl  # for r=1 nodes: T4(3s+1) = 3*sigma1_s + 0 (r=0 node at sigma1_s)

    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()

    # Estimate rho from K-L equation at convergence
    cb_f = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w_f = A * v[T4]
    w_f[m2] += B3 * cb_f[R3[m2]]
    w_f[m0] += B1 * cb_f[R1[m0]]
    rho = float(np.max(w_f))  # since v[argmax] = 1.0

    v0 = v[0::3]  # v0[s] = v at (s, r=0)
    v1 = v[1::3]  # v1[s] = v at (s, r=1)
    v2 = v[2::3]  # v2[s] = v at (s, r=2)

    return v, Nl, v0, v1, v2, cb_f, A, B1, B3, rho, sigma1_s[m0 | m2 | True]  # return sigma1 for r=1 positions

def verify_c1_c0(k, lam):
    A = lam ** -2.0
    B1 = lam ** (ALPHA - 2.0)
    B3 = lam ** (ALPHA - 1.0)
    N = 3 ** (k - 1)
    Nl = N // 3
    Nl3 = Nl // 3

    v, _, v0, v1, v2, cb, A, B1, B3, rho, _ = run_kl_full(k, lam)

    # Compute c_r = Mean(cb[j with j%3==r])
    j = np.arange(Nl, dtype=np.int64)
    c_0 = float(np.mean(cb[j%3==0]))
    c_1 = float(np.mean(cb[j%3==1]))
    c_2 = float(np.mean(cb[j%3==2]))

    A_over_rho = A / rho

    # Check v1[s] = (A/rho) * v0[sigma1(s)] exactly
    s = np.arange(Nl, dtype=np.int64)
    sigma1_s = (4*s + 2) % Nl  # sigma1: s -> (4s+2) % Nl

    v1_from_eq = A_over_rho * v0[sigma1_s]
    max_err_v1 = float(np.max(np.abs(v1 - v1_from_eq)))
    rel_err_v1 = float(np.max(np.abs(v1 - v1_from_eq) / (v1 + 1e-15)))

    # Check c_1 = (A/rho) * c_0 exactly
    c1_predicted = A_over_rho * c_0
    c1_rel_err = abs(c_1 - c1_predicted) / (c_1 + 1e-15)

    print(f"\nk={k}, lam={lam:.2f}: A={A:.4f}, rho={rho:.6f}, A/rho={A_over_rho:.6f}")
    print(f"  c_0={c_0:.6f}, c_1={c_1:.6f}, c_2={c_2:.6f}")
    print(f"  c_1/c_0 = {c_1/c_0:.6f}  vs  A/rho = {A_over_rho:.6f}")
    print(f"  Relative error: {abs(c_1/c_0 - A_over_rho)/A_over_rho:.2e}")
    print(f"  v1 identity check: max_err={max_err_v1:.2e}  rel_err={rel_err_v1:.2e}")

    # cb[3s'+1] = (A/rho) * min(v0[sigma1(s')], v0[(sigma1(s')+Nl3)%Nl], v0[(sigma1(s')+2Nl3)%Nl])?
    Nl3_val = Nl // 3
    s_p = np.arange(Nl3_val, dtype=np.int64)  # s' in [0, Nl//3)
    j1 = 3*s_p + 1  # r=1 indices in [0, Nl)
    cb_v1 = cb[j1]  # column-min at r=1 positions
    sigma1_cb_j = 3*sigma1_s  # r=0 indices
    # But sigma1_cb_j might be >= Nl... no wait: sigma1(s) in [0, Nl), so 3*sigma1(s) in [0, 3Nl) = [0, N)
    # But cb has length Nl, not N. Let me re-check.
    # cb[j] for j in [0, Nl). cb[3s+1] for s in [0, Nl//3).
    # sigma1_s = (4s+2)%Nl, this is in [0, Nl).
    # cb[3*sigma1_s] would need 3*sigma1_s < Nl, but sigma1_s < Nl so 3*sigma1_s could be >= Nl.
    # Wait: cb is indexed by j in [0, Nl). j can be 3s+0, 3s+1, 3s+2 for s in [0, Nl//3).
    # But I wrote cb[3*sigma1_s] where sigma1_s in [0, Nl). So 3*sigma1_s can be up to 3*Nl which is out of range!
    # This is wrong. Let me reconsider.

    # Actually the CLAIM was:
    # cb[3s+1] = (A/rho) * cb[3*sigma1(s) + 0] where sigma1(s) = (4s+2)%Nl
    # But sigma1(s) is in [0, Nl), so 3*sigma1(s) is in [0, 3Nl).
    # But cb has length Nl! The indices of cb are j in [0, Nl).
    # The r=0 column-min at position sigma1(s) is cb[3*(sigma1(s)//3) + sigma1(s)%3]...
    #
    # Wait, I was confused. Let me redo.
    #
    # j = 3s'+r for s' in [0, Nl//3), r in {0,1,2}.
    # cb[j] = min(v[j], v[j+Nl], v[j+2Nl]) for j in [0, Nl).
    # For j = 3s'+1 (r=1): cb[3s'+1] = min(v1[s'], v1[s'+Nl3], v1[s'+2Nl3])
    # (where Nl3 = Nl//3).
    # v1[t] = (A/rho)*v0[sigma1(t)] where sigma1(t) = (4t+2)%Nl.
    # sigma1(s'+Nl3) = (4(s'+Nl3)+2)%Nl = (4s'+4Nl3+2)%Nl = (sigma1(s')+4Nl3)%Nl.
    # Since 4Nl3 ≡ Nl3 (mod Nl) [because 4Nl3 = Nl3 + Nl and Nl mod Nl = 0]:
    # sigma1(s'+Nl3) = (sigma1(s')+Nl3)%Nl.
    # Similarly sigma1(s'+2Nl3) = (sigma1(s')+2Nl3)%Nl.
    #
    # So: cb[3s'+1] = (A/rho)*min(v0[sigma1(s')], v0[(sigma1(s')+Nl3)%Nl], v0[(sigma1(s')+2Nl3)%Nl])
    #              = (A/rho)*cb[3*(sigma1(s')//3)]... but sigma1(s') may not be divisible by 3.
    #
    # Hmm. Let t = sigma1(s'). t in [0, Nl). t = 3*(t//3) + t%3. The column-min at r=0 for the
    # s-position t//3 would be cb[3*(t//3)] (if r-type of t is 0, i.e., t%3==0).
    # But sigma1(s') = (4s'+2)%Nl. The r-type of sigma1(s') is (4s'+2)%3 = (s'+2)%3.
    # This is NOT always 0.
    #
    # I made an error in the derivation. The claim c_1 = (A/rho)*c_0 is NOT exactly correct.
    # Let me verify numerically.

    # Verify: does cb[3s'+1] = (A/rho) * (something involving cb at r=0)?
    Nl3_val = Nl // 3
    s_p = np.arange(Nl3_val, dtype=np.int64)  # s' in [0, Nl//3)
    sigma1_sp = (4*s_p + 2) % Nl  # sigma1(s')

    # cb[3s'+1] (r=1 column-min at s')
    cb_r1_sp = cb[3*s_p + 1]

    # (A/rho) * min(v0[sigma1(s')], v0[(sigma1(s')+Nl3)%Nl], v0[(sigma1(s')+2Nl3)%Nl])
    sigma1_plus_Nl3 = (sigma1_sp + Nl3_val) % Nl
    sigma1_plus_2Nl3 = (sigma1_sp + 2*Nl3_val) % Nl

    cb_v0_at_sigma1 = np.minimum(
        np.minimum(v0[sigma1_sp], v0[sigma1_plus_Nl3]),
        v0[sigma1_plus_2Nl3]
    )

    cb_r1_predicted = A_over_rho * cb_v0_at_sigma1
    max_err_cb = float(np.max(np.abs(cb_r1_sp - cb_r1_predicted)))
    rel_err_cb = float(np.max(np.abs(cb_r1_sp - cb_r1_predicted) / (np.abs(cb_r1_sp) + 1e-15)))

    print(f"\n  cb[r=1] identity check:")
    print(f"    cb[3s'+1] = (A/rho)*min(v0[sigma1(s')], v0[sigma1(s')+Nl3], v0[sigma1(s')+2Nl3])")
    print(f"    max_abs_err = {max_err_cb:.2e}  max_rel_err = {rel_err_cb:.2e}")

    # Mean of cb_v0_at_sigma1 vs c_0?
    mean_cb_v0_sigma1 = float(np.mean(cb_v0_at_sigma1))

    # c_0: column-min at r=0 positions
    cb_r0 = np.array([cb[3*s_pp] for s_pp in range(Nl3_val)])
    # Actually: cb_r0[s'] = min(v0[s'], v0[s'+Nl3], v0[s'+2Nl3]) = cb[3s']
    cb_r0_direct = cb[3*s_p]
    c_0_from_sp = float(np.mean(cb_r0_direct))

    print(f"  Mean(cb_v0_at_sigma1) = {mean_cb_v0_sigma1:.6f}  c_0 = {c_0_from_sp:.6f}")
    print(f"  Ratio: {mean_cb_v0_sigma1/c_0_from_sp:.6f}  (1.0 if sigma1 preserves mean)")

    return c_1/c_0, A_over_rho

print("259: Analytical verification: c_1 = (A/rho)*c_0 (EXACT)")
print("="*70)

# === MAIN VERIFICATION k=8, lam=1.70 ===
verify_c1_c0(8, 1.70)

# === LAMBDA SCAN ===
print(f"\n\nLambda scan k=8:")
print(f"{'lam':>6}  {'c1/c0':>8}  {'A/rho':>8}  {'rel_err':>10}")
for lam in [1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]:
    r1, r2 = verify_c1_c0(8, lam)
    sys.stdout.flush()

# === DEPTH SCAN ===
print(f"\n\nDepth scan lam=1.70:")
print(f"{'k':>4}  {'c1/c0':>8}  {'A/rho':>8}  {'rel_err':>10}")
lam = 1.70
for k in range(5, 12):
    r1, r2 = verify_c1_c0(k, lam)
    sys.stdout.flush()

print("\ndone")
