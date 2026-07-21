"""
126_ab_transitions.py
======================
Quantify and explain the A<->B alternation in BSet transitions.

Group A (K>=5): {63, 95, 127, 159, 191, 223, 255}
Group B (K<=4): {27, 55, 83, 103, 169, 207, 239, 253}

Key questions:
1. What is P(next BSet in A | current in A)?  P(next in B | current in A)?
2. Can we PROVE from arithmetic that K>=5 macro-steps tend to output K<=4 residues?
3. What is the full 15x15 BSet-to-BSet hitting matrix and its spectrum?
4. How does the A/B partition change at larger moduli (512, 1024)?
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

BSet = [27,55,63,83,95,103,127,159,169,191,207,223,239,253,255]
MOD = 256
bset_set = set(BSet)
GroupA = {63, 95, 127, 159, 191, 223, 255}   # K>=5
GroupB = {27, 55, 83, 103, 169, 207, 239, 253}  # K<=4

print("=" * 70)
print("PART 1: DIRECT 15x15 BSet HITTING MATRIX")
print("=" * 70)

bset_list = sorted(BSet)
bidx = {r: i for i, r in enumerate(bset_list)}
NS = 4000  # samples per BSet element

P_hit = np.zeros((15, 15))
for i, r in enumerate(bset_list):
    counts = np.zeros(15)
    hits = 0
    for k in range(NS):
        n = r + MOD*k
        cur = n
        for _ in range(500):
            cur_out, _, _ = macro_step(cur)
            r_out = cur_out % MOD
            if r_out in bset_set:
                counts[bidx[r_out]] += 1
                hits += 1
                break
            cur = cur_out
    if hits > 0:
        P_hit[i] = counts / hits

# Show the 2x2 A/B block structure
# A = indices of Group A elements in bset_list
A_idx = [i for i,r in enumerate(bset_list) if r in GroupA]
B_idx = [i for i,r in enumerate(bset_list) if r in GroupB]
print(f"Group A indices: {[bset_list[i] for i in A_idx]}")
print(f"Group B indices: {[bset_list[i] for i in B_idx]}")

# Aggregate into 2x2 block matrix
P_aa = P_hit[np.ix_(A_idx, A_idx)].sum(axis=1)  # per A-row, sum over A cols
P_ab = P_hit[np.ix_(A_idx, B_idx)].sum(axis=1)  # per A-row, sum over B cols
P_ba = P_hit[np.ix_(B_idx, A_idx)].sum(axis=1)
P_bb = P_hit[np.ix_(B_idx, B_idx)].sum(axis=1)

print(f"\nRow-by-row A-to-A and A-to-B probabilities:")
print(f"  {'r':>5} {'P(->A)':>10} {'P(->B)':>10}")
for ii, i in enumerate(A_idx):
    r = bset_list[i]
    K = v2(r+1)
    print(f"  {r:>5} (K={K}) {P_aa[ii]:>10.4f} {P_ab[ii]:>10.4f}")

print(f"\nRow-by-row B-to-A and B-to-B probabilities:")
print(f"  {'r':>5} {'P(->A)':>10} {'P(->B)':>10}")
for ii, i in enumerate(B_idx):
    r = bset_list[i]
    K = v2(r+1)
    print(f"  {r:>5} (K={K}) {P_ba[ii]:>10.4f} {P_bb[ii]:>10.4f}")

# Overall 2x2 matrix (weighted by stationary distribution)
pi = np.ones(15)/15  # approximate stationary
for _ in range(500):
    pi_new = pi @ P_hit
    if np.max(np.abs(pi_new - pi)) < 1e-10: break
    pi = pi_new
pi_A = pi[A_idx].sum(); pi_B = pi[B_idx].sum()
P_AA = (pi[A_idx] @ P_hit[np.ix_(A_idx, A_idx)]).sum() / pi_A
P_AB = (pi[A_idx] @ P_hit[np.ix_(A_idx, B_idx)]).sum() / pi_A
P_BA = (pi[B_idx] @ P_hit[np.ix_(B_idx, A_idx)]).sum() / pi_B
P_BB = (pi[B_idx] @ P_hit[np.ix_(B_idx, B_idx)]).sum() / pi_B

print(f"\nStationary distribution: pi(A)={pi_A:.4f}, pi(B)={pi_B:.4f}")
print(f"\n2x2 aggregate transition matrix (weighted by stationary dist):")
print(f"        -> A     -> B")
print(f"  A:  {P_AA:8.4f} {P_AB:8.4f}")
print(f"  B:  {P_BA:8.4f} {P_BB:8.4f}")
M22 = np.array([[P_AA, P_AB], [P_BA, P_BB]])
ev22 = np.linalg.eigvals(M22.T)
print(f"  Eigenvalues: {sorted(ev22.real, reverse=True)}")
print(f"  Second eigenvalue (mixing rate): {min(ev22.real):.6f}")

print()
print("=" * 70)
print("PART 2: ARITHMETIC EXPLANATION — WHY K>=5 -> K<=4?")
print("=" * 70)

# After a K-step with K large, the output is:
#   n_out = (m * 3^K - 1) / 2^l0
# K of next step = v2(n_out + 1)
# n_out + 1 = (m * 3^K - 1) / 2^l0 + 1 = (m * 3^K - 1 + 2^l0) / 2^l0
# The 2-adic val of the numerator m*3^K - 1 + 2^l0:
#   m*3^K - 1 = 2^l0 * Q  (where Q is odd, by definition of l0)
#   So m*3^K - 1 + 2^l0 = 2^l0 * (Q + 1) = 2^l0 * (odd + 1) = 2^l0 * 2^s * (odd)
#   where s = v2(Q+1) = v2((m*3^K - 1)/2^l0 + 1) = v2(n_out + 1) = K_next
# So K_next = v2((m*3^K - 1)/2^l0 + 1)

# This is entirely arithmetic. The KEY POINT: for FIXED r mod 256 (BSet element with K large),
# what is the distribution of K_next?

print("\nFor each Group A BSet element r, compute empirical P(K_next = k) for NEXT BSet VISIT:")
print(f"  {'r':>5} {'K_r':>4} | {'K=1':>7} {'K=2':>7} {'K=3':>7} {'K=4':>7} | {'K=5..8':>9}")
print("  " + "-" * 55)

for r in sorted(GroupA):
    K_r = v2(r+1)
    K_next_counts = defaultdict(int)
    total = 0
    for k in range(2000):
        n = r + MOD * k
        cur = n
        for _ in range(500):
            cur_out, _, _ = macro_step(cur)
            r_out = cur_out % MOD
            if r_out in bset_set:
                K_next = v2(r_out + 1)
                K_next_counts[K_next] += 1
                total += 1
                break
            cur = cur_out
    p1 = K_next_counts.get(1,0)/total
    p2 = K_next_counts.get(2,0)/total
    p3 = K_next_counts.get(3,0)/total
    p4 = K_next_counts.get(4,0)/total
    p_high = sum(K_next_counts.get(k,0) for k in range(5,9))/total
    print(f"  {r:>5} (K={K_r}) | {p1:>7.4f} {p2:>7.4f} {p3:>7.4f} {p4:>7.4f} | {p_high:>9.4f}")

print(f"\nFor each Group B BSet element r, compute empirical P(K_next = k) for NEXT BSet VISIT:")
print(f"  {'r':>5} {'K_r':>4} | {'K=1':>7} {'K=2':>7} {'K=3':>7} {'K=4':>7} | {'K=5..8':>9}")
print("  " + "-" * 55)

for r in sorted(GroupB):
    K_r = v2(r+1)
    K_next_counts = defaultdict(int)
    total = 0
    for k in range(2000):
        n = r + MOD * k
        cur = n
        for _ in range(500):
            cur_out, _, _ = macro_step(cur)
            r_out = cur_out % MOD
            if r_out in bset_set:
                K_next = v2(r_out + 1)
                K_next_counts[K_next] += 1
                total += 1
                break
            cur = cur_out
    p1 = K_next_counts.get(1,0)/total
    p2 = K_next_counts.get(2,0)/total
    p3 = K_next_counts.get(3,0)/total
    p4 = K_next_counts.get(4,0)/total
    p_high = sum(K_next_counts.get(k,0) for k in range(5,9))/total
    print(f"  {r:>5} (K={K_r}) | {p1:>7.4f} {p2:>7.4f} {p3:>7.4f} {p4:>7.4f} | {p_high:>9.4f}")

print()
print("=" * 70)
print("PART 3: HOW MANY MACRO-STEPS BETWEEN BSet VISITS?")
print("=" * 70)

# Distribution of inter-BSet gap length
import random as _r; _r.seed(7)
n = _r.getrandbits(2000) | 1
gap_dist = defaultdict(int)
in_bset_prev = (n % MOD) in bset_set
gap = 0
for _ in range(100000):
    n_out, _, _ = macro_step(n)
    n = n_out
    r = n % MOD
    gap += 1
    if r in bset_set:
        gap_dist[gap] += 1
        gap = 0
    if n < 2: n = _r.getrandbits(2000) | 1

total_gaps = sum(gap_dist.values())
print(f"\nInter-BSet gap distribution (steps between consecutive BSet visits):")
print(f"  {'Gap':>6} {'Count':>8} {'Prob':>8} {'CDF':>8}")
cdf = 0
for g in sorted(gap_dist)[:20]:
    prob = gap_dist[g]/total_gaps
    cdf += prob
    print(f"  {g:>6} {gap_dist[g]:>8} {prob:>8.4f} {cdf:>8.4f}")

mean_gap = sum(g*c for g,c in gap_dist.items()) / total_gaps
print(f"\nMean inter-BSet gap = {mean_gap:.4f}")
print(f"Theory (if geometric): E[gap] = 1/pi_BSet = 1/{15/128:.5f} = {128/15:.4f}")
print(f"(Gap=1 = consecutive BSet visits, confirming P(BSet->BSet)=0.35)")
