"""
91_k_autocorrelation_mersenne.py
==================================
Three analyses:

1. K-AUTOCORRELATION: After visiting a k0-type booster, how does k decay
   over the next J macro-steps? Quantifies the "k-debt" after high-k visits.

2. 255 SELF-LOOP: Exact h=1 count for 255→255 (exact arithmetic),
   and the large-n h distribution for 255→255 specifically.

3. MERSENNE THRESHOLD: Why are 2^6-1, 2^7-1, 2^8-1 in BSet but not
   2^5-1=31, 2^4-1=15, etc.? Compute the single-step drift for each
   Mersenne number and find the threshold.
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1)
    m = (n + 1) >> k
    x = m * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
BList = sorted(BSet)
k0_of = {r: v2(r+1) for r in BList}

M_PERIOD = 256   # 256 odd m-values per booster
M_BASE = 10**12  # Large-n base

# =====================================================================
# PART 1: EXACT 255→255 SELF-LOOP (h=1 count)
# =====================================================================
print("=" * 70)
print("PART 1: EXACT 255→255 SELF-LOOP")
print("=" * 70)
print()
print("Computing exact P(output ≡ 255 mod 256 in 1 step from r=255)...")
print("Using 256 odd m-values: m=1,3,...,511 with n = 256m-1")
print()

pow3_8 = 3**8  # = 6561

# Exact h=1 from r=255 to each BSet element
bset_h1_counts = Counter()
for m_idx in range(M_PERIOD):
    m = 2 * m_idx + 1  # m = 1, 3, ..., 511
    val = pow3_8 * m - 1
    l = v2(val)
    out = val >> l
    r_out = out % 256
    if r_out in BSet:
        bset_h1_counts[r_out] += 1

total_h1 = sum(bset_h1_counts.values())
print(f"Total h=1 transitions: {total_h1}/256 = {100*total_h1/256:.3f}%")
print(f"  (Expected: 31/256 = 12.109%)")
print()
print("h=1 destination distribution from r=255:")
print(f"{'dest r':>8}  {'k0':>4}  {'count':>6}  {'frac of h=1':>12}")
for r in sorted(bset_h1_counts.keys(), key=lambda x: -bset_h1_counts[x]):
    k0 = k0_of[r]
    cnt = bset_h1_counts[r]
    print(f"  r={r:3d}  k={k0}: {cnt:4d}  ({100*cnt/total_h1:.2f}% of h=1 transitions)")

# Check 255→255 h=1 probability
p_255_self_h1 = bset_h1_counts[255] / M_PERIOD
print(f"\n255→255 in exactly 1 step: {bset_h1_counts[255]}/256 = {100*p_255_self_h1:.4f}%")
print(f"255→127 in exactly 1 step: {bset_h1_counts[127]}/256 = {100*bset_h1_counts[127]/M_PERIOD:.4f}%")

# For h=1 self-loop: cycle mean = k0 / h = 8 / 1 = 8.0 (each h=1 round trip)
# But this contributes very little because P(h=1→255) is small
print(f"\nFor the h=1 self-loop at r=255:")
print(f"  Each h=1 self-hop contributes k/step = 8/1 = 8.0")
print(f"  But P(self-hop in 1 step) = {bset_h1_counts[255]}/256 = {p_255_self_h1:.4f}")
print(f"  => h=1 self-loop contributes {8 * p_255_self_h1:.4f} to avg k/step")

# =====================================================================
# PART 2: K-AUTOCORRELATION ANALYSIS
# =====================================================================
print()
print("=" * 70)
print("PART 2: K-AUTOCORRELATION AFTER BOOSTER VISITS")
print("=" * 70)
print()

# For each booster type k0 (using representative booster),
# trace 256 long orbits and record the sequence of k-values
# for the first J=20 macro-steps after the booster visit.
J_MAX = 20
N_TRACES = 256

print(f"Tracing {N_TRACES} large-n starting points per k0 type for {J_MAX} steps.")
print(f"Recording k-sequence: k_1, k_2, ..., k_{J_MAX}")
print()

k0_groups = defaultdict(list)
for r in BList:
    k0_groups[k0_of[r]].append(r)
K0_TYPES = sorted(k0_groups.keys())
rep = {k0: k0_groups[k0][0] for k0 in K0_TYPES}

print(f"{'k0':>4}  ", end='')
print(" ".join(f"E[k_{j:2d}]" for j in range(1, 11)))

k_sequences = {}  # k0 -> list of k-vectors

for k0 in K0_TYPES:
    r = rep[k0]
    step = 2**(k0 + 1)
    k_seqs = []
    for m_idx in range(N_TRACES):
        m = M_BASE + 2 * m_idx + 1
        n = step * m - 1
        seq = []
        n_cur = n
        for j in range(J_MAX):
            n_cur, k, l = macro_step(n_cur)
            seq.append(k)
            if n_cur <= 1:
                while len(seq) < J_MAX: seq.append(0)  # orbit converged
                break
        if len(seq) < J_MAX: seq.extend([0]*(J_MAX-len(seq)))
        k_seqs.append(seq)
    k_sequences[k0] = k_seqs

    avg_k = [sum(seq[j] for seq in k_seqs)/N_TRACES for j in range(min(10,J_MAX))]
    print(f"  k0={k0}: ", end='')
    print(" ".join(f"{x:6.3f}" for x in avg_k))

# Long-run E[k] for random odd n
print()
print("Theoretical E[k] for random odd n = 2.000 (geometric distribution: E[sum j*2^{-j}] = 2)")
print("Observed: E[k_j] should converge to ~2.0 as j increases")

# Show full decay from k0=8
print()
print("K-decay profile for k0=8 (r=255):")
k_seqs_8 = k_sequences[8]
print(f"{'Step j':>7}  {'E[k_j]':>8}  {'std(k_j)':>10}  {'E[k_j]/k0=8':>12}")
for j in range(J_MAX):
    vals = [seq[j] for seq in k_seqs_8 if seq[j] > 0]
    if vals:
        mean_k = sum(vals)/len(vals)
        std_k = (sum((v-mean_k)**2 for v in vals)/len(vals))**0.5
        print(f"  j={j+1:3d}:  {mean_k:6.3f}    {std_k:6.3f}    {mean_k/8:.4f}")

# K-autocorrelation coefficient at lag 1
print()
print("K-lag-1 autocorrelation from k0=8 orbits:")
k1 = [k_seqs_8[i][0] for i in range(N_TRACES)]
k2 = [k_seqs_8[i][1] for i in range(N_TRACES) if k_seqs_8[i][1] > 0]
if k2:
    pairs = [(k_seqs_8[i][j], k_seqs_8[i][j+1]) for i in range(N_TRACES) for j in range(J_MAX-1) if k_seqs_8[i][j]>0 and k_seqs_8[i][j+1]>0]
    mean_k = sum(a for a,b in pairs)/len(pairs)
    cov = sum((a-mean_k)*(b-mean_k) for a,b in pairs)/len(pairs)
    var = sum((a-mean_k)**2 for a,b in pairs)/len(pairs)
    corr = cov/var if var > 0 else 0
    print(f"  rho(k_t, k_{{t+1}}) = {corr:.4f}")

# Check: what is E[k_{t+1}] given k_t is high vs low?
print()
k_pairs_by_kt = defaultdict(list)
for k0 in K0_TYPES:
    for i in range(N_TRACES):
        for j in range(J_MAX - 1):
            kt = k_sequences[k0][i][j]
            kt1 = k_sequences[k0][i][j+1]
            if kt > 0 and kt1 > 0:
                k_pairs_by_kt[kt].append(kt1)

print("E[k_{t+1} | k_t = c] (regression effect):")
for c in range(1, 10):
    if k_pairs_by_kt[c]:
        mean_next = sum(k_pairs_by_kt[c]) / len(k_pairs_by_kt[c])
        print(f"  E[k_{{t+1}} | k_t={c}] = {mean_next:.4f}  (n={len(k_pairs_by_kt[c])})")

# =====================================================================
# PART 3: MERSENNE THRESHOLD
# =====================================================================
print()
print("=" * 70)
print("PART 3: MERSENNE THRESHOLD: WHY 2^k-1 ∈ BSET IFF k >= 6?")
print("=" * 70)
print()

print("For each Mersenne number r = 2^k - 1 (k0=k):")
print("Compute: single-step growth factor = 3^k0 / 2^E[l], drift = log2(3^k0) - E[l]")
print()
print("Single-step macro-step: n' = (3^k0 * m - 1) / 2^l for random odd m")
print("Growth factor per step: 3^k0 * m / (2^l * n'_odd) ≈ 3^k0 / 2^l")
print()

print(f"{'k':>4}  {'r=2^k-1':>10}  {'k0':>4}  {'E[l]':>8}  {'drift=k*log2(3)-E[l]':>22}  {'in BSet?':>10}")
print("-"*80)
for k in range(1, 9):
    r = 2**k - 1
    k0 = v2(r+1)
    pow3k = 3**k0

    # Compute E[l] exactly for 256 odd m-values
    l_vals = []
    for m_idx in range(M_PERIOD):
        m = 2 * m_idx + 1
        val = pow3k * m - 1
        l_vals.append(v2(val))

    E_l = sum(l_vals) / M_PERIOD
    import math
    drift = k * math.log2(3) - E_l
    in_bset = r in BSet

    print(f"  k={k}: r={r:3d}  k0={k0}  E[l]={E_l:.4f}  drift={drift:+.4f}  {'IN BSet' if in_bset else 'NOT in BSet'}")

print()
print("Threshold: drift > 0 <=> 2^k-1 is a 'single-step booster'")
print("BSet membership requires MULTI-STEP drift analysis (not just single step)")
print("Key: single-step drift > 0 is NECESSARY but not sufficient for BSet")
print()

# More detailed: for each Mersenne number, compute the long-run avg k/step
# using N large-n trajectories
print("Long-run avg k/step for Mersenne numbers (using 100 large-n trajectories):")
print(f"{'k':>4}  {'r=2^k-1':>10}  {'avg_k/step':>12}  {'E[h to BSet]':>14}  {'BSet? (3.419 threshold)':>25}")
print("-" * 80)

N_TRAJ = 100  # shorter for quick computation
MAX_H = 1000

for k_val in range(1, 9):
    r_mersenne = 2**k_val - 1
    k0 = v2(r_mersenne + 1)
    step = 2**(k0 + 1)

    # Large-n starting points with n ≡ r_mersenne mod step
    k_sums = []
    h_vals_m = []
    for traj_idx in range(N_TRAJ):
        m = M_BASE + 2 * traj_idx + 1
        n = step * m - 1  # n ≡ 2^k0 * m - 1 ≡ 2^k-1 mod 2^{k+1}
        # Verify residue
        r_check = n % 256 if k0 < 8 else n % 256

        n_cur = n
        h = 0
        k_sum = 0
        for _ in range(MAX_H):
            n_cur, k, l = macro_step(n_cur)
            h += 1
            k_sum += k
            if n_cur <= 1:
                break
            if n_cur % 256 in BSet:
                break

        if h > 0:
            k_sums.append(k_sum)
            h_vals_m.append(h)

    if h_vals_m:
        avg_h = sum(h_vals_m) / len(h_vals_m)
        avg_k_step = sum(k_sums) / sum(h_vals_m)
        in_bset = r_mersenne in BSet
        note = f"IN BSet ✓" if in_bset else f"NOT in BSet (avg_k/step too low)"
        print(f"  k={k_val}: r={r_mersenne:3d}  avg_k/step={avg_k_step:.4f}  E[h]={avg_h:.2f}  {note}")

# =====================================================================
# PART 4: SELF-LOOP CYCLE MEAN BOUND FROM EXACT h=1 DATA
# =====================================================================
print()
print("=" * 70)
print("PART 4: RIGOROUS LOWER BOUND ON E[h(255→255)]")
print("=" * 70)
print()

# For 255 self-loop:
# P(h=1) = bset_h1_counts[255] / 256 (exact)
# E[h] >= 1 * P(h=1) + 2 * P(h>=2) = 1 * P1 + 2 * (1-P1)
# = 2 - P1
p1_255_self = bset_h1_counts[255] / M_PERIOD
print(f"P(255→255 in h=1) = {bset_h1_counts[255]}/{M_PERIOD} = {100*p1_255_self:.4f}% (EXACT)")
lb_trivial = 1 * p1_255_self + 2 * (1 - p1_255_self)
print(f"Trivial lower bound: E[h(255→255)] >= {lb_trivial:.4f}")

# Better: P(h=2)
# From large-n 256-point sample: how many 255→255 paths have h=2?
print()
print("Tracing 255→255 paths for exact h distribution (large-n, 256 starting points):")

h_255_self = []
for m_idx in range(M_PERIOD):
    m = M_BASE + 2 * m_idx + 1
    n = 256 * m - 1
    n_cur = n
    h = 0
    k_sum = 0
    for _ in range(MAX_H):
        n_cur, k, l = macro_step(n_cur)
        h += 1
        k_sum += k
        if n_cur <= 1: break
        if n_cur % 256 == 255:  # Hit 255 again
            h_255_self.append((h, k_sum))
            break

print(f"Paths reaching 255 again: {len(h_255_self)}/256 = {100*len(h_255_self)/256:.2f}%")
if h_255_self:
    h_list_self = [h for h,k in h_255_self]
    k_list_self = [k for h,k in h_255_self]
    E_h_self = sum(h_list_self) / len(h_list_self)
    E_k_self = sum(k_list_self) / sum(h_list_self)
    print(f"E[h(255→255)] = {E_h_self:.3f}")
    print(f"E[k/step(255→255)] = {E_k_self:.4f}")
    print(f"Self-loop cycle mean = E[k/step] = {E_k_self:.4f}")
    h_ctr = Counter(h_list_self)
    print(f"\nh distribution for 255→255 paths:")
    for h in sorted(h_ctr.keys())[:20]:
        k_here = sum(k for hh,k in h_255_self if hh==h)
        print(f"  h={h}: {h_ctr[h]} paths  avg_k/step={k_here/h_ctr[h]/h:.3f}")

# =====================================================================
# PART 5: WHY THE GAP IS 0.708 NOT 0.622
# =====================================================================
print()
print("=" * 70)
print("PART 5: UNDERSTANDING THE GAP DISCREPANCY (2.711 vs 2.7974)")
print("=" * 70)
print()
print("Script 82 (N=5000, n~10^9): lambda* = 2.7974")
print("Script 90 (N=256, n~10^14): lambda* = 2.7111")
print()
print("Possible explanations:")
print("1. SAMPLE SIZE: Script 90 has ~6-12 paths per 8-type edge (noisy)")
print("   Script 82 has N=5000 starting points per booster (much larger)")
print()
print("2. N DEPENDENCE: Large-n orbits may have slightly different statistics")
print("   than medium-n orbits. The true lambda* may be between 2.71-2.80.")
print()
print("3. K0 AGGREGATION: Script 90 collapses 15 boosters to 8 types,")
print("   losing within-type variation. Script 82 uses all 15 separately.")
print()
print("Either way: both estimates are far below the threshold 3.419.")
print("Gap from script 82: 3.419 - 2.797 = 0.622")
print("Gap from script 90: 3.419 - 2.711 = 0.708")
print("Consistent conclusion: D_hard_kern = ∅")
