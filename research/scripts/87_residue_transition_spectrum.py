"""
87_residue_transition_spectrum.py
===================================
The fundamental 128x128 transition matrix of the Collatz macro-step
on odd residues mod 256.

T[a][b] = P(output ≡ b mod 256 | input ≡ a mod 256)

This is EXACT (one full period: 256 odd m-values per source residue a).

Key questions:
1. What is the spectrum of T? Spectral gap determines mixing speed.
2. What is the stationary distribution of T? (Should be related to the
   density of elements at each odd residue in Collatz orbits)
3. What is the restriction of T to BSet? (15x15 sub-matrix of the full
   booster transition, regardless of hop length h -- this is the effective
   Markov chain on boosters)
4. Why is P(h=1) ≈ 15/128 and not exactly uniform?
   => Understand which residues send outputs preferentially to BSet vs non-BSet.
5. What does the "BSet-to-BSet via non-BSet" transition look like?
   => 2-step transfer operator restricted to BSet.

Also: compute the expected k0 of "random BSet hits" for comparison with
the h=1 and h=2 observed destination distributions.
"""
import sys
import numpy as np
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
BList = sorted(BSet)
k0_of = {r: v2(r+1) for r in BList}

# All 128 odd residues mod 256
ODD_RES = sorted([r for r in range(256) if r % 2 == 1])  # 1,3,5,...,255
res_idx = {r: i for i, r in enumerate(ODD_RES)}
N_RES = len(ODD_RES)  # 128

print(f"Computing 128x128 transition matrix on odd residues mod 256...")
print(f"Using 256 odd m-values per source (exact one-period computation).\n")

# For each odd residue a, we assume n ≡ a mod 256 (with n odd).
# The macro-step: k = v2(n+1), output = (3^k * m - 1) / 2^l
# But k = v2(n+1) depends on n, not just n mod 256!
# n+1 is even (since n is odd). n mod 256 = a determines only the low 8 bits of n+1.
# v2(n+1) = v2(a+1) if the carry from the low 8 bits doesn't propagate further.
# But n+1 can have higher powers of 2 depending on higher bits of n.
#
# IMPORTANT: T[a][b] is the transition probability AVERAGED over all n ≡ a mod 256.
# Since v2(n+1) depends on higher bits, we must average over those higher bits.
#
# For a given odd residue a:
# - n+1 ≡ a+1 mod 256
# - If a+1 ≡ 0 mod 2^j but ≢ 0 mod 2^{j+1}: k=j with probability 1 (if j < 8)
#   UNLESS a+1 ≡ 0 mod 256 (i.e., a=255), in which case k >= 8 with P(k=j) = 1/2^{j-8} for j>=8.
#
# Key insight: k = v2(n+1) depends on the FULL n, not just n mod 256.
# For a ≡ r mod 2^{k0+1} with exactly v2(a+1) = k0 (and k0 < 8):
#   k = k0 with probability 1 (since a+1 has exactly k0 trailing zeros mod 256,
#   and higher bits of n are "random" odd multiples of 2^k0)
#   ACTUALLY: n+1 = 2^k0 * m where m = (n+1)/2^k0. m might be even! If m is even,
#   then v2(n+1) > k0. So k = k0 requires m to be ODD. P(m odd) depends on higher bits.
#
# Correct approach: for n ≡ a mod 256, n+1 ≡ a+1 mod 256.
# Write a+1 = 2^k0 * u where u is odd and k0 = v2(a+1) (with k0 ≤ 7, since a is odd so a+1 is even).
# Then n+1 = 2^k0 * (u + t * 256/2^k0 * 2) for some integer t.
# m = (n+1)/2^k0 = u + t * (256/2^k0) * 2 ... hmm this gets complicated.
#
# PRACTICAL APPROACH: For each source residue a, use 256 odd m-values to compute
# the output distribution. Here m = (n+1)/2^{v2(n+1)}, but we treat m as ranging
# over all odd values uniformly (the "high-n limit" assumption).
#
# This is exactly what scripts 84-86 did for BOOSTER sources. Now we do it for ALL 128 sources.
#
# For source residue a with k0 = v2(a+1):
# The macro-step output = (3^k0 * m - 1) / 2^{v2(3^k0 * m - 1)} for m odd.
# Using 256 odd m values (one period for output mod 256).

M_PERIOD = 256  # 256 odd m values

# Build T matrix
T = np.zeros((N_RES, N_RES))

for a in ODD_RES:
    i = res_idx[a]
    k0 = v2(a + 1)  # the k value for n ≡ a mod 256
    pow3k = 3**k0

    output_counts = Counter()
    for m_idx in range(M_PERIOD):
        m = 2 * m_idx + 1  # m = 1, 3, ..., 511
        val = pow3k * m - 1
        l = v2(val)
        out = val >> l
        r_out = out % 256
        output_counts[r_out] += 1

    # Fill row i of T
    total = M_PERIOD
    for b, cnt in output_counts.items():
        if b in res_idx:
            j = res_idx[b]
            T[i][j] = cnt / total
        # (if b not in ODD_RES: this can't happen since out is always odd)

print(f"T matrix computed: {N_RES}x{N_RES}")
print(f"Row sums (should be ≈1.0): min={T.sum(axis=1).min():.6f}, max={T.sum(axis=1).max():.6f}")
print()

# =====================================================================
# SPECTRAL ANALYSIS
# =====================================================================
print("=== SPECTRAL ANALYSIS OF T ===\n")

eigenvalues = np.linalg.eigvals(T)
eigenvalues_real = np.sort(np.abs(eigenvalues))[::-1]  # Sort by magnitude

print(f"Top 10 eigenvalues by magnitude:")
for i, ev in enumerate(eigenvalues_real[:10]):
    print(f"  λ_{i+1} = {ev:.8f}")

print(f"\nSpectral gap: λ_1 - λ_2 = {eigenvalues_real[0] - eigenvalues_real[1]:.8f}")
print(f"Second eigenvalue λ_2 = {eigenvalues_real[1]:.8f}")
print(f"Mixing time estimate: τ ≈ 1/gap = {1/(eigenvalues_real[0]-eigenvalues_real[1]):.2f} steps")

# =====================================================================
# STATIONARY DISTRIBUTION
# =====================================================================
print("\n=== STATIONARY DISTRIBUTION ===\n")

# Power iteration
pi = np.ones(N_RES) / N_RES
for _ in range(100000):
    pi_new = pi @ T
    pi_new /= pi_new.sum()
    if np.max(np.abs(pi_new - pi)) < 1e-12:
        pi = pi_new
        break
    pi = pi_new

print("Top 20 stationary distribution weights:")
top_idx = np.argsort(pi)[::-1][:20]
for rank, idx in enumerate(top_idx):
    r = ODD_RES[idx]
    is_bset = " [BSet]" if r in BSet else ""
    k0 = v2(r+1)
    print(f"  rank {rank+1}: r={r:3d} (k0={k0}) π={pi[idx]:.6f}{is_bset}")

# BSet total stationary weight
bset_pi = sum(pi[res_idx[r]] for r in BSet)
print(f"\nTotal stationary weight on BSet: {bset_pi:.6f} ({100*bset_pi:.3f}%)")
print(f"Expected if uniform: {15/128:.6f} ({100*15/128:.3f}%)")

# Avg k0 under stationary distribution (for BSet elements only)
avg_k0_stat = sum(pi[res_idx[r]] * k0_of[r] for r in BSet) / bset_pi
print(f"Avg k0 of BSet elements under stationary: {avg_k0_stat:.4f}")
print(f"Avg k0 of BSet elements uniform: {sum(k0_of[r] for r in BSet)/15:.4f}")

# =====================================================================
# THE "RANDOM BSet HIT" DISTRIBUTION
# =====================================================================
print("\n=== EXPECTED k0 FOR RANDOM BSet HITS ===\n")
print("For random odd n, P(v2(n+1)=k0) = 1/2^k0 for k0=1,2,3,...")
print("Expected k0 of a random BSet hit, weighted by P(k0):\n")

# Weight of each k0 type in BSet: n_elements(k0) * P(v2(r+1)=k0)
k0_groups = {}
for r in BList:
    k0 = k0_of[r]
    if k0 not in k0_groups:
        k0_groups[k0] = 0
    k0_groups[k0] += 1

total_weight = sum(cnt / 2**k0 for k0, cnt in k0_groups.items())
print(f"{'k0':>5} {'n_elements':>12} {'P(k0)=1/2^k0':>15} {'weight':>10} {'fraction':>10}")
for k0 in sorted(k0_groups.keys()):
    cnt = k0_groups[k0]
    pk0 = 1 / 2**k0
    w = cnt * pk0
    frac = w / total_weight
    print(f"  k0={k0}: {cnt:12d}  {pk0:.6f}       {w:.6f}   {frac:.4f}")

expected_k0_random = sum(k0 * cnt / 2**k0 for k0, cnt in k0_groups.items()) / total_weight
print(f"\nExpected k0 for random BSet hit: {expected_k0_random:.4f}")
print(f"Actual h=1 destination avg k0: ~4.0-4.3 (from script 84/86)")
print(f"Actual h=2 destination avg k0: ~2.4-3.2 (from script 85/86)")
print(f"Actual h=3 destination avg k0: ~1.9-3.7 (from script 86)")
print(f"True random limit: {expected_k0_random:.4f}")

# =====================================================================
# T RESTRICTED TO BSet: BSet -> BSet via 1-step
# =====================================================================
print("\n=== T RESTRICTED TO BSet (15x15 h=1 TRANSITION MATRIX) ===\n")

T_bset = np.zeros((15, 15))
bset_idx = {r: i for i, r in enumerate(BList)}

for i, r_src in enumerate(BList):
    row_src = res_idx[r_src]
    total_h1 = sum(T[row_src][res_idx[r_dst]] for r_dst in BList)
    for j, r_dst in enumerate(BList):
        if total_h1 > 0:
            T_bset[i][j] = T[row_src][res_idx[r_dst]] / total_h1

print(f"Note: T_bset[i][j] = fraction of h=1 transitions from r_i to r_j")
print(f"(normalized so rows sum to 1)\n")

# Stationary dist of T_bset
pi_bset = np.ones(15) / 15
for _ in range(100000):
    pi_new = pi_bset @ T_bset
    pi_new /= pi_new.sum()
    if np.max(np.abs(pi_new - pi_bset)) < 1e-12:
        pi_bset = pi_new
        break
    pi_bset = pi_new

print("Stationary distribution of BSet h=1 chain (per booster):")
for i, r in enumerate(BList):
    k0 = k0_of[r]
    print(f"  r={r:3d} (k0={k0}): π={pi_bset[i]:.6f} ({100*pi_bset[i]:.3f}%)")

avg_k0_bset_stat = sum(pi_bset[i] * k0_of[BList[i]] for i in range(15))
print(f"\nAvg k0 under BSet stationary dist: {avg_k0_bset_stat:.4f}")

# =====================================================================
# 2-STEP TRANSFER: BSet --(2 steps)--> BSet
# T^2 restricted to BSet gives P(h=1 at step 1 AND h=1 at step 2 | start at r)
# vs T * T restricted to BSet gives P(output in BSet after 2 steps)
# =====================================================================
print("\n=== P(h<=2) FROM SPECTRAL PERSPECTIVE ===\n")

# The probability of hitting BSet in exactly 2 steps starting from BSet element r:
# = sum over non-BSet r' of P(r -> r' in 1 step) * P(r' -> BSet in 1 step)
# This uses only the transition matrix T.

print("P(h=2) from transition matrix T (vs exact computation):")
for r_src in BList:
    i = res_idx[r_src]
    k0 = k0_of[r_src]
    # P(h=2) = sum over non-BSet states r' of T[r_src][r'] * P(r' in BSet in 1 step)
    p_h2 = 0.0
    for r_mid in ODD_RES:
        if r_mid in BSet:
            continue
        j_mid = res_idx[r_mid]
        # P(r_mid -> BSet in 1 step) = sum over BSet of T[r_mid][r_dst]
        p_bset_from_mid = sum(T[j_mid][res_idx[r_dst]] for r_dst in BSet)
        p_h2 += T[i][j_mid] * p_bset_from_mid
    print(f"  r={r_src:3d} k0={k0}: P(h=2) from T = {100*p_h2:.4f}%  "
          f"(vs exact script 85: varies per booster)")

print("\n[Done]")
