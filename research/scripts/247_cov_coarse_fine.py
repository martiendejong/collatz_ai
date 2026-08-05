"""
247_cov_coarse_fine.py
======================
Investigate the coarse-to-fine structure of the Cov(ld_v2_sig0, ld_cb_sig0) < 0.

For j equiv 2 mod 3, the CODE-triplet of v2_at_sigma0 at group g uses:
  v2_at_sigma0 values: {v2[j*], v2[j*+Nl/3], v2[j*+2Nl/3]}  where j* = sigma0(g_base)
  cb_at_sigma0 values: {cb[j*], cb[j*+Nl/3], cb[j*+2Nl/3]}

For j* equiv 2 mod 3, cb[j*] = min(v2[j*//3], v2[j*//3+Nl/3], v2[j*//3+2Nl/3]).
These are the COARSE CODE-triplet of v2 at base j*//3.

Meanwhile, v2_at_sigma0 slot 0 = v2[j*] = v2_interleaved[j*].
And j* = 3*(j*//3) + 2 = 3q+2 where q = j*//3.

Question: Is v2_interleaved[3q+2] one of the three values in cb[3q+2]=min(v2[q],v2[q+Nl/3],v2[q+2Nl/3])?

NO - as shown analytically, v2_interleaved[3q+2] is a DIFFERENT position.
v2_interleaved[3q+2] = v[3*(3q+2)+2] = v[9q+8] (fine grid)
v2[q] = v[3*3q+2] = v[9q+2] (different position)

But there IS a FINER-LEVEL CODE-triplet relationship:
  The Nl-dim v2_interleaved array has CODE-triplets at spacing Nl//3.
  The COARSE triplet base q uses positions {q, q+Nl//3, q+2Nl//3}.
  Position 3q+2 is in a DIFFERENT CODE-triplet at base (3q+2)%(Nl//3).

This script checks: what is the actual relationship between:
  1. v2[3q+2] and cb[3q+2] = min(v2[q], v2[q+Nl/3], v2[q+2Nl/3])
  2. The coarse CODE-triplet structure

Key test: at what level does the anti-correlation arise?
"""
import numpy as np
from math import log2
import sys

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
    return v, Nl

K = 8
LAM = 1.70

v, Nl = run_kl(K, LAM)
v2 = v[2::3]  # shape (Nl,): v2[s] = v[3s+2]

# Block cb: cb[j] = min(v[j], v[j+Nl], v[j+2Nl]) for j in [0, Nl)
cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])

# For j equiv 2 mod 3: cb[j] = min(v2[j//3], v2[j//3+Nl/3], v2[j//3+2Nl/3])
# Verify this:
Nl3 = Nl // 3
j2_base = np.arange(Nl3, dtype=np.int64) * 3 + 2  # j values equiv 2 mod 3: 2,5,8,...
q = j2_base // 3  # = 0,1,...,Nl3-1
cb_formula = np.minimum(np.minimum(v2[q], v2[q+Nl3]), v2[q+2*Nl3])
print(f"k={K}, lam={LAM}: cb formula check (j equiv 2):")
print(f"  max|cb[3q+2] - min(v2[q],v2[q+Nl3],v2[q+2Nl3])| = {np.max(np.abs(cb[j2_base] - cb_formula)):.2e}")
print()

# Now look at the coarse-to-fine relationship:
# v2 CODE-triplet at base q (in Nl-dim space with Nl3 spacing):
#   {v2[q], v2[q+Nl3], v2[q+2Nl3]} = the three values in cb[3q+2]
# v2_interleaved[3q+2] = v2[3q+2] -- a FINE-level position
# Is 3q+2 related to {q, q+Nl3, q+2Nl3}?
# 3q+2 vs q: difference = 2q+2 (large for large q)
# 3q+2 vs q+Nl3: difference = 2q+2-Nl3
# 3q+2 vs q+2Nl3: difference = 2q+2-2Nl3

print("Coarse-to-fine relationship: v2[3q+2] vs COARSE triplet {v2[q], v2[q+Nl3], v2[q+2Nl3]}")
print("These are DIFFERENT positions in v2_interleaved.")
print()

# Key question: does v2[3q+2] predict the RELATIVE position within the coarse triplet?
# Compute: for each q, compare:
#   rank of v2[3q+2] among the fine-level CODE-triplet it belongs to
#   rank of v2[q] in its coarse triplet {v2[q],v2[q+Nl3],v2[q+2Nl3]}
# The fine CODE-triplet of v2 at position 3q+2: belongs to group (3q+2)%(Nl3)
# This varies by q, so let's just look at the scatter

# sigma0 map and CODE-triplets
sigma0 = (4 * np.arange(Nl, dtype=np.int64)) % Nl
v2_sig0 = v2[sigma0]  # v2 pulled back through sigma0

# CODE-triplet of v2_sig0 (groups of Nl3 elements):
v2_sig0_t0 = v2_sig0[:Nl3]
v2_sig0_t1 = v2_sig0[Nl3:2*Nl3]
v2_sig0_t2 = v2_sig0[2*Nl3:]
mean_v2_sig0 = (v2_sig0_t0 + v2_sig0_t1 + v2_sig0_t2) / 3.0
ld_v2_sig0_0 = np.log2(v2_sig0_t0 / mean_v2_sig0)
ld_v2_sig0_1 = np.log2(v2_sig0_t1 / mean_v2_sig0)
ld_v2_sig0_2 = np.log2(v2_sig0_t2 / mean_v2_sig0)

# CODE-triplet of cb_sig0
cb_sig0 = cb[sigma0]
cb_sig0_t0 = cb_sig0[:Nl3]
cb_sig0_t1 = cb_sig0[Nl3:2*Nl3]
cb_sig0_t2 = cb_sig0[2*Nl3:]
mean_cb_sig0 = (cb_sig0_t0 + cb_sig0_t1 + cb_sig0_t2) / 3.0
ld_cb_sig0_0 = np.log2(cb_sig0_t0 / mean_cb_sig0)
ld_cb_sig0_1 = np.log2(cb_sig0_t1 / mean_cb_sig0)
ld_cb_sig0_2 = np.log2(cb_sig0_t2 / mean_cb_sig0)

# Total Cov
cov_total = (np.mean(ld_v2_sig0_0*ld_cb_sig0_0) + np.mean(ld_v2_sig0_1*ld_cb_sig0_1) +
             np.mean(ld_v2_sig0_2*ld_cb_sig0_2)) / 3.0
print(f"Total Cov(ld_v2_sig0, ld_cb_sig0) = {cov_total:.6e}")

# Now decompose by r-type of sigma0 (= group type of j):
# j equiv r mod 3 → group has r-type r
# Slot 0 of CODE-triplet: j in [0, Nl3), corresponds to groups g=0..Nl3-1
# Slot 0 group with r-type 2: groups g where g%3==2, i.e., g ∈ {2,5,8,...}
g_all = np.arange(Nl3, dtype=np.int64)
for r in range(3):
    mask = (g_all % 3 == r)
    # Slot 0 Cov for this r-type group
    cov_slot0_r = np.mean(ld_v2_sig0_0[mask] * ld_cb_sig0_0[mask])
    cov_slot1_r = np.mean(ld_v2_sig0_1[mask] * ld_cb_sig0_1[mask])
    cov_slot2_r = np.mean(ld_v2_sig0_2[mask] * ld_cb_sig0_2[mask])
    cov_r = (cov_slot0_r + cov_slot1_r + cov_slot2_r) / 3.0 * (np.sum(mask)/Nl3)
    print(f"r={r}: slot0={cov_slot0_r:+.4e}  slot1={cov_slot1_r:+.4e}  slot2={cov_slot2_r:+.4e}  "
          f"total_contrib={cov_r:+.4e}")

print()

# For r=2 groups (g%3==2): COARSE triplet analysis
# Group g with g%3==2: base = g*3+2 in the j-space... no
# j_base = g (j runs from 0 to Nl3-1 = Nl//3-1), sigma0(j) = 4j%Nl
# For g%3==2 (r=2 groups): the three slots are j = g, g+Nl3, g+2Nl3

# Let's look at what cb[sigma0(g)] relates to:
# sigma0(g) for g equiv 2 mod 3: 4g%Nl equiv 2 mod 3
# cb[4g%Nl] = min(v2[(4g%Nl)//3], v2[(4g%Nl)//3+Nl3], v2[(4g%Nl)//3+2Nl3])
# v2_sig0[g] = v2[4g%Nl]
# Is v2[4g%Nl] one of the three values in cb[4g%Nl]?
# v2[4g%Nl] = v2_interleaved[4g%Nl] (at index 4g%Nl in Nl-dim space)
# cb[4g%Nl] uses indices {(4g%Nl)//3, (4g%Nl)//3+Nl3, (4g%Nl)//3+2Nl3}
# (4g%Nl)//3 ≠ 4g%Nl in general!

mask_r2 = (g_all % 3 == 2)
g_r2 = g_all[mask_r2]
j_star = (4 * g_r2) % Nl  # sigma0(g) for g equiv 2
q_coarse = j_star // 3    # coarse base
# Check: is j_star one of {q, q+Nl3, q+2Nl3}?
in_triplet = (j_star == q_coarse) | (j_star == q_coarse + Nl3) | (j_star == q_coarse + 2*Nl3)
frac_in = np.sum(in_triplet) / len(g_r2)
print(f"r=2 groups: frac where v2[j*] IN coarse triplet of cb[j*]: {frac_in:.4f}  "
      f"({np.sum(in_triplet)} / {len(g_r2)})")
print(f"(expected ~1/Nl3 = 1/{Nl3} = {1/Nl3:.4f} for large k)")

# So the 'same values' claim only holds for ~1/Nl3 fraction
# What about a LARGER structural relationship?
# v2[j*] where j* = 3q+2 (q = j*//3):
# The FINE-level position 3q+2 in v2_interleaved
# The COARSE triplet positions: {q, q+Nl3, q+2Nl3}
# These are at VERY DIFFERENT locations in v2_interleaved space

# Cov breakdown: slot 0 vs slot 1 vs slot 2 for r=2 groups
print()
print(f"r=2 group Cov breakdown by slot:")
cov_0 = np.mean(ld_v2_sig0_0[mask_r2] * ld_cb_sig0_0[mask_r2])
cov_1 = np.mean(ld_v2_sig0_1[mask_r2] * ld_cb_sig0_1[mask_r2])
cov_2 = np.mean(ld_v2_sig0_2[mask_r2] * ld_cb_sig0_2[mask_r2])
print(f"  slot 0: Cov = {cov_0:+.6e}")
print(f"  slot 1: Cov = {cov_1:+.6e}")
print(f"  slot 2: Cov = {cov_2:+.6e}")

# If the 'same values' argument were right, slot 0 should be most negative
# (since slot 0 = base position g, sigma0(g) = 4g%Nl, and cb[4g%Nl] is
#  the min of a triplet that MIGHT include v2[4g%Nl] for some g)

print()
print("done")
