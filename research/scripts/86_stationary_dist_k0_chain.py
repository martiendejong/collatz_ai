"""
86_stationary_dist_k0_chain.py
================================
Two analyses:

1. Exact P(h=3) and k_avg_dest(h=3): verify the h-dependent k-destination drift
   conjecture (k_avg_dest decreases as h increases).

2. Stationary distribution of the 8-state k0-type Markov chain:
   - h=1 transition matrix P[k0_src][k0_dst] (from script 84 exact data)
   - Stationary distribution π (power iteration)
   - Avg k under stationary distribution
   - Compare to empirical booster visit frequencies from script 82
"""
import sys
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
M_PERIOD = 256  # 256 odd m values

# k0 value for each booster
k0_of = {r: v2(r+1) for r in BList}
# Inverse: k0 -> list of boosters with that k0
k0_groups = defaultdict(list)
for r in BList:
    k0_groups[k0_of[r]].append(r)
K0_TYPES = sorted(k0_groups.keys())  # [1,2,3,4,5,6,7,8]

print("=== K0 GROUPS ===")
for k0 in K0_TYPES:
    rs = k0_groups[k0]
    print(f"  k0={k0}: {rs}")

# =====================================================================
# PART 1: Exact P(h=1), P(h=2), P(h=3) and destination k-averages
# Using ONE representative booster per k0 type (they are identical)
# =====================================================================
print("\n=== EXACT P(h=1), P(h=2), P(h=3) AND k_avg_dest BY k0 TYPE ===\n")

rep = {k0: k0_groups[k0][0] for k0 in K0_TYPES}  # representative booster per k0

type_results = {}

for k0 in K0_TYPES:
    r = rep[k0]
    pow3k = 3**k0

    h1_dest = Counter()
    h2_dest = Counter()
    h3_dest = Counter()
    h_gt3 = 0

    for m_idx in range(M_PERIOD):
        m = 2 * m_idx + 1

        # Step 1 from booster r
        val1 = pow3k * m - 1
        l1 = v2(val1)
        out1 = val1 >> l1
        r1 = out1 % 256

        if out1 <= 1:
            continue  # converged

        if r1 in BSet:
            h1_dest[r1] += 1
            continue

        # Step 2 from out1
        k2 = v2(out1 + 1)
        m2 = (out1 + 1) >> k2
        val2 = m2 * (3**k2) - 1
        l2 = v2(val2)
        out2 = val2 >> l2
        r2 = out2 % 256

        if out2 <= 1:
            continue

        if r2 in BSet:
            h2_dest[r2] += 1
            continue

        # Step 3 from out2
        k3 = v2(out2 + 1)
        m3 = (out2 + 1) >> k3
        val3 = m3 * (3**k3) - 1
        l3 = v2(val3)
        out3 = val3 >> l3
        r3 = out3 % 256

        if out3 <= 1:
            continue

        if r3 in BSet:
            h3_dest[r3] += 1
        else:
            h_gt3 += 1

    ph1 = sum(h1_dest.values()) / M_PERIOD
    ph2 = sum(h2_dest.values()) / M_PERIOD
    ph3 = sum(h3_dest.values()) / M_PERIOD
    ph_gt3 = h_gt3 / M_PERIOD

    ka1 = sum(k0_of[r2]*c for r2,c in h1_dest.items()) / max(1, sum(h1_dest.values()))
    ka2 = sum(k0_of[r2]*c for r2,c in h2_dest.items()) / max(1, sum(h2_dest.values()))
    ka3 = sum(k0_of[r2]*c for r2,c in h3_dest.items()) / max(1, sum(h3_dest.values()))

    type_results[k0] = {
        'h1': h1_dest, 'h2': h2_dest, 'h3': h3_dest,
        'ph1': ph1, 'ph2': ph2, 'ph3': ph3, 'ph_gt3': ph_gt3,
        'ka1': ka1, 'ka2': ka2, 'ka3': ka3
    }

    print(f"k0={k0} (r={r}):")
    print(f"  P(h=1)={100*ph1:.3f}%  k_avg_dest(h=1)={ka1:.3f}")
    print(f"  P(h=2)={100*ph2:.3f}%  k_avg_dest(h=2)={ka2:.3f}")
    print(f"  P(h=3)={100*ph3:.3f}%  k_avg_dest(h=3)={ka3:.3f}")
    print(f"  P(h>3)={100*ph_gt3:.3f}%")
    print(f"  Drift: {ka1:.3f} -> {ka2:.3f} -> {ka3:.3f}  "
          f"(diffs: {ka2-ka1:+.3f}, {ka3-ka2:+.3f})")
    print()

# =====================================================================
# PART 2: h=1 TRANSITION MATRIX ON 8-STATE k0 CHAIN
# P[k0_src][k0_dst] = fraction of h=1 outputs from k0_src landing in k0_dst
# =====================================================================
print("\n=== h=1 TRANSITION MATRIX (8-STATE k0 CHAIN) ===\n")

# Build transition matrix from h1_dest data
# P[k0_src][k0_dst] = (count of h1 destinations with k0(r')=k0_dst) / total_h1 from k0_src

P_h1 = np.zeros((len(K0_TYPES), len(K0_TYPES)))
k0_idx = {k0: i for i, k0 in enumerate(K0_TYPES)}

for k0 in K0_TYPES:
    h1 = type_results[k0]['h1']
    total = sum(h1.values())
    if total == 0:
        continue
    for r_dst, cnt in h1.items():
        k0_dst = k0_of[r_dst]
        P_h1[k0_idx[k0]][k0_idx[k0_dst]] += cnt / total

print("P_h1 matrix (rows=src k0, cols=dst k0, values=fraction of h=1 transitions):")
print(f"      " + "  ".join(f"k0={k0:d}" for k0 in K0_TYPES))
for i, k0_src in enumerate(K0_TYPES):
    row = "  ".join(f"{P_h1[i,j]:.3f}" for j in range(len(K0_TYPES)))
    print(f"k0={k0_src}: {row}")

# =====================================================================
# PART 3: STATIONARY DISTRIBUTION OF h=1 k0 CHAIN
# (where weights = fraction of h=1 transitions going to each k0 type)
# =====================================================================
print("\n=== STATIONARY DISTRIBUTION (h=1 k0 CHAIN) ===\n")

# Power iteration to find stationary distribution
pi = np.ones(len(K0_TYPES)) / len(K0_TYPES)
for _ in range(10000):
    pi = pi @ P_h1
    pi /= pi.sum()

print("Stationary distribution π (fraction of h=1 booster arrivals at each k0):")
for i, k0 in enumerate(K0_TYPES):
    print(f"  k0={k0} (r∈{k0_groups[k0]}): π={pi[i]:.6f}  ({100*pi[i]:.3f}%)")

avg_k0_stationary = sum(pi[i] * K0_TYPES[i] for i in range(len(K0_TYPES)))
print(f"\n  Avg k0 under stationary distribution: {avg_k0_stationary:.4f}")
print(f"  Uniform avg k0 (15 equal weights): {sum(k0_of[r] for r in BList)/15:.4f}")
print(f"  BSet k0 mean: {sum(k0_of[r] for r in BList)/15:.4f}")

# =====================================================================
# PART 4: FULL UNCONDITIONAL DESTINATION DISTRIBUTION FROM EACH k0 TYPE
# Combine h=1, h=2, h=3 (exact) + h>3 (approximated by h=1 stationary dist)
# =====================================================================
print("\n=== APPROXIMATE UNCONDITIONAL DESTINATION DISTRIBUTION BY k0 ===\n")
print("Combining P(h=1)*k1_avg + P(h=2)*k2_avg + P(h=3)*k3_avg + P(h>3)*k_stat\n")

for k0 in K0_TYPES:
    d = type_results[k0]
    ph1, ph2, ph3, ph_gt3 = d['ph1'], d['ph2'], d['ph3'], d['ph_gt3']
    ka1, ka2, ka3 = d['ka1'], d['ka2'], d['ka3']
    ka_stat = avg_k0_stationary  # for h>3, use stationary distribution k0 avg

    # Weighted average of k_avg_dest
    total_weight = ph1 + ph2 + ph3 + ph_gt3
    if total_weight > 0:
        ka_uncond = (ph1*ka1 + ph2*ka2 + ph3*ka3 + ph_gt3*ka_stat) / total_weight
    else:
        ka_uncond = float('nan')

    print(f"k0={k0}: unconditional k_avg_dest={ka_uncond:.3f}  "
          f"(h1:{ph1:.3f}*{ka1:.2f} + h2:{ph2:.3f}*{ka2:.2f} + "
          f"h3:{ph3:.3f}*{ka3:.2f} + h>3:{ph_gt3:.3f}*{ka_stat:.2f})")

# =====================================================================
# PART 5: VERIFY k-DESTINATION DRIFT ACROSS h=1,2,3
# =====================================================================
print("\n=== k-DESTINATION DRIFT SUMMARY (k_avg_dest by hop length) ===\n")
print(f"{'k0':>5}  {'h=1':>8}  {'h=2':>8}  {'h=3':>8}  {'diff(1->2)':>12}  {'diff(2->3)':>12}")
print("-"*70)
for k0 in K0_TYPES:
    d = type_results[k0]
    ka1, ka2, ka3 = d['ka1'], d['ka2'], d['ka3']
    d12 = ka2 - ka1
    d23 = ka3 - ka2
    print(f"k0={k0}:  {ka1:.4f}  {ka2:.4f}  {ka3:.4f}  {d12:+.4f}     {d23:+.4f}")

print()
print("Conjecture to verify: k_avg_dest(h) is DECREASING in h (for all k0 types).")
all_decreasing_12 = all(type_results[k0]['ka2'] < type_results[k0]['ka1'] for k0 in K0_TYPES)
all_decreasing_23 = all(type_results[k0]['ka3'] < type_results[k0]['ka2']
                        for k0 in K0_TYPES if sum(type_results[k0]['h3'].values()) > 0)
print(f"k_avg_dest(h=2) < k_avg_dest(h=1) for ALL k0: {all_decreasing_12}")
print(f"k_avg_dest(h=3) < k_avg_dest(h=2) for all k0 with h=3 data: {all_decreasing_23}")
