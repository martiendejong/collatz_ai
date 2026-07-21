"""
127_bset_scaling.py
====================
Does the BSet A/B partition (K>=5 vs K<=4) persist at larger moduli?

At mod-256: BSet has 15 elements. K>=5 = GroupA, K<=4 = GroupB.
At mod-512: The BSet analogue has more elements. Does the K threshold change?
At mod-1024: Even more BSet elements.

Key questions:
1. How many BSet elements at each modulus?
2. Do BSet elements still cluster around specific K values?
3. Does the bipartite structure (negative eigenvalue) persist?
4. What is the K threshold that defines the partition at each modulus?

Definition: BSet_N = set of residues r (mod 2^N) that appear as outputs
of the mod-2^N Markov chain and that are "gateway" residues (have
self-loops or concentrate transitions). More precisely, here we use:
BSet at mod-M = residues mod M that the chain's stationary distribution
assigns significant weight to relative to a uniform comparison.
"""
import numpy as np
import math
from collections import defaultdict

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("=" * 70)
print("PART 1: K DISTRIBUTION AT EACH MODULUS")
print("=" * 70)
print()
print("For mod-M chain, what K values are represented in odd residues?")
print("K=v2(r+1) is determined by r alone. For r in [1,M-1] odd:")
print()

for MOD in [128, 256, 512, 1024]:
    odd_res = list(range(1, MOD, 2))
    K_dist = defaultdict(int)
    for r in odd_res:
        K_dist[v2(r+1)] += 1
    print(f"mod-{MOD}: {len(odd_res)} odd residues")
    for k in sorted(K_dist)[:10]:
        cnt = K_dist[k]
        frac = cnt/len(odd_res)
        print(f"  K={k}: {cnt:5d} residues ({frac:.4f}), theory 1/2^k={1/2**k:.4f}")
    kmax = max(K_dist.keys())
    print(f"  Max K = {kmax} = log2(MOD)")
    print()

print("=" * 70)
print("PART 2: BSet ANALOGUES AT LARGER MODULI")
print("=" * 70)
print()
print("Compute stationary distribution, find residues with anomalously high weight.")
print()

NS = 256  # transitions per residue

for MOD in [256, 512]:
    odd_res = list(range(1, MOD, 2))
    N = len(odd_res)
    idx = {r: i for i, r in enumerate(odd_res)}

    # Build transition matrix
    P = np.zeros((N, N), dtype=np.float64)
    for i, r in enumerate(odd_res):
        K0 = v2(r+1)
        counts = np.zeros(N)
        valid = 0
        for k in range(NS):
            n = r + MOD*k
            if v2(n+1) != K0: continue
            n_out, _, _ = macro_step(n)
            r_out = n_out % MOD
            if r_out % 2 == 1 and r_out in idx:
                counts[idx[r_out]] += 1; valid += 1
        if valid > 0: P[i] = counts / valid

    # Stationary distribution
    pi = np.ones(N)/N
    for _ in range(1000):
        pi_new = pi @ P
        if np.max(np.abs(pi_new - pi)) < 1e-10: break
        pi = pi_new

    # Find high-weight residues (> threshold * uniform)
    threshold = 3.0
    uniform = 1/N
    high_weight = [(odd_res[i], pi[i], v2(odd_res[i]+1)) for i in range(N) if pi[i] > threshold * uniform]
    high_weight.sort(key=lambda x: -x[1])

    print(f"mod-{MOD}: {len(high_weight)} residues with pi > {threshold}x uniform ({threshold*uniform:.5f})")
    print(f"  {'r':>6} {'K':>4} {'pi':>10} {'pi/unif':>10}")
    for r, p, K in high_weight[:20]:
        print(f"  {r:>6} {K:>4} {p:>10.6f} {p/uniform:>10.2f}x")

    # K distribution among high-weight residues
    K_hw = defaultdict(int)
    for r, p, K in high_weight:
        K_hw[K] += 1
    print(f"  K distribution of high-weight residues: {dict(sorted(K_hw.items()))}")
    print()

print("=" * 70)
print("PART 3: K DISTRIBUTION CONDITIONAL ON BEING IN STATIONARY DISTRIBUTION")
print("=" * 70)

# The stationary distribution assigns weight pi(r) to each residue r.
# What is the effective K distribution when sampling from stationary?
# E_pi[f(K)] = sum_r pi(r) * f(K(r))

for MOD in [256, 512, 1024]:
    odd_res = list(range(1, MOD, 2))
    N = len(odd_res)

    # Approximate stationary with random orbit
    import random as _r
    _r.seed(17)
    n = _r.getrandbits(2000) | 1
    K_counts = defaultdict(int)
    for _ in range(50000):
        n_out, K_val, _ = macro_step(n)
        r = n % MOD
        K_counts[K_val] += 1  # This is K of CURRENT step, which equals K(r)
        n = n_out
        if n < 2: n = _r.getrandbits(2000) | 1

    total = sum(K_counts.values())
    print(f"mod-{MOD}: K distribution under stationary (50k steps):")
    for k in sorted(K_counts)[:8]:
        obs = K_counts[k]/total
        theory = 1/2**k
        print(f"  K={k}: {obs:.5f} (theory: {theory:.5f}, ratio: {obs/theory:.3f})")
    print()

print("=" * 70)
print("PART 4: CCT-BASED BSET DEFINITION")
print("=" * 70)

# A more principled definition of BSet at each modulus:
# BSet_N = CCT elements that appear in the STATIONARY distribution
# with high weight, i.e., residues that are "central" to the chain.
#
# Alternatively, BSet is defined by: r in BSet iff r appears as output
# of mod-M chain with probability significantly higher than 1/M.
#
# Let's compute: for each residue r, how often does the orbit MOD return to r?
# High return frequency = higher stationary weight = BSet candidate.

print("Checking: which residues have v2(r+1) = N-1 or N-2 (near max K)?")
print("These are the 'gateway' residues at each modulus.")
print()
for MOD in [256, 512, 1024]:
    N = int(math.log2(MOD))
    near_max = [(r, v2(r+1)) for r in range(1, MOD, 2) if v2(r+1) >= N-2]
    print(f"mod-{MOD} (N={N}): residues with K>={N-2} (near max K=N):")
    for r, K in sorted(near_max):
        print(f"  r={r:6d}, K={K}")
    print(f"  Total: {len(near_max)} residues (expected: 3)")
    print()
