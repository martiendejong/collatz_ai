"""
125_bset_kstruct.py
====================
How does the BSet structure affect K-distributions?

The BSet = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255} mod 256
consists of exactly the 15 residues where the macro-step behavior is "forced"
in the sense that K = v2(n+1) is fixed by the residue class mod 256.

Questions:
1. What is the K value for each BSet element? (It should be exactly K=7 or K=8
   for all of them, since BSet is defined as having n ≡ -1 mod 128 or mod 256.)
2. How does the transition structure from BSet differ from non-BSet?
3. Does exiting BSet always produce l0 values that match the geometric distribution?
4. What is the variance of l0 from BSet transitions? Does it differ from the
   l0 variance from non-BSet?
5. Connection: The spectral gap is 0.606 for BSet chain (|lambda_2|=0.394).
   Does the K=7 or K=8 forced value at BSet cause the oscillatory mode?
"""
import numpy as np
from collections import defaultdict
import math

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

print("=" * 70)
print("PART 1: K VALUE FOR EACH BSET ELEMENT (and l0 distribution)")
print("=" * 70)
print()

# For each BSet element, determine K (fixed by residue class)
# and compute l0 distribution from many lifts
NS = 2000  # samples per BSet element
print(f"{'Residue':>8} {'K':>4} {'E[l0]':>8} {'Var[l0]':>8} {'P(l0=1)':>8} {'P(l0=2)':>8} {'P(l0=3)':>8}")
print("-" * 65)

bset_K_vals = {}
bset_l0_stats = {}
for r in sorted(BSet):
    K = v2(r + 1)
    l0_list = []
    for k in range(NS):
        n = r + MOD * k
        _, K_actual, l0 = macro_step(n)
        if K_actual == K:  # confirm K matches
            l0_list.append(l0)
    l0_arr = np.array(l0_list)
    bset_K_vals[r] = K
    bset_l0_stats[r] = l0_arr
    p1 = (l0_arr == 1).mean()
    p2 = (l0_arr == 2).mean()
    p3 = (l0_arr == 3).mean()
    print(f"{r:>8} {K:>4} {l0_arr.mean():>8.4f} {l0_arr.var():>8.4f} {p1:>8.4f} {p2:>8.4f} {p3:>8.4f}")

print()
print(f"Theory: E[l0]=2, Var[l0]=2, P(l0=k)=1/2^k")
print(f"K values in BSet: {sorted(set(bset_K_vals.values()))}")

# For non-BSet residues
print()
print("=" * 70)
print("PART 2: K DISTRIBUTION IN NON-BSET RESIDUES")
print("=" * 70)

non_bset = [r for r in range(1, MOD, 2) if r not in BSet]
NS2 = 200
K_dist_nonbset = defaultdict(int)
l0_dist_nonbset = defaultdict(int)
for r in non_bset[:32]:  # sample first 32 non-BSet residues
    for k in range(NS2):
        n = r + MOD * k
        _, K_val, l0 = macro_step(n)
        K_dist_nonbset[K_val] += 1
        l0_dist_nonbset[l0] += 1

total_nb = sum(K_dist_nonbset.values())
print(f"\nK distribution from non-BSet residues ({len(non_bset[:32])} residues, {NS2} lifts each):")
print(f"{'K':>4} {'Observed':>10} {'Theory 1/2^K':>14}")
for k in sorted(K_dist_nonbset)[:10]:
    obs = K_dist_nonbset[k] / total_nb
    theory = 1/2**k
    print(f"{k:>4} {obs:>10.5f} {theory:>14.5f}")

print()
print("l0 distribution from non-BSet:")
total_nb2 = sum(l0_dist_nonbset.values())
for l in sorted(l0_dist_nonbset)[:8]:
    obs = l0_dist_nonbset[l] / total_nb2
    theory = 1/2**l
    print(f"l0={l}: {obs:.5f} (theory: {theory:.5f})")

print()
print("=" * 70)
print("PART 3: LYAPUNOV CONTRIBUTIONS OF EACH BSET ELEMENT")
print("=" * 70)

# The ACTUAL K and l0 values for each BSet element (from Part 1)
# Compute Lyapunov for each element using the ACTUAL E[l0] from measurements
print()
print(f"{'Residue':>8} {'K':>4} {'E[l0]':>8} {'Lyapunov':>12} {'Sign':>6}")
print("-" * 45)
for r in sorted(BSet):
    K = bset_K_vals[r]
    l0_arr = bset_l0_stats[r]
    el0 = l0_arr.mean()
    lp = K*math.log(3) - (K + el0)*math.log(2)
    sign = "+" if lp > 0 else "-"
    print(f"{r:>8} {K:>4} {el0:>8.4f} {lp:>12.6f} {sign:>6}")

print()
print("Cross-over: K*log3 - (K+E[l0])*log2 = 0 when E[l0] = K*(log3/log2-1) = K*0.585")
print("For K=1: need E[l0] = 0.585 to break even")
print("For K=2: need E[l0] = 1.170")
print("For K=3: need E[l0] = 1.754")
print("For K=4: need E[l0] = 2.339")
print("General: Lyapunov > 0 iff E[l0] < K*(log3-log2)/log2 = K*0.585")

print()
print("l0 is FIXED (Var=0) for BSet elements with K<=5:")
for r in sorted(BSet):
    K = bset_K_vals[r]
    l0_arr = bset_l0_stats[r]
    if l0_arr.var() < 0.01:
        l0_fixed = int(l0_arr[0])
        lp = K*math.log(3) - (K + l0_fixed)*math.log(2)
        print(f"  r={r:>3}: K={K}, l0={l0_fixed} (fixed), Lyapunov={lp:+.4f}")

print()
print("=" * 70)
print("PART 4: AUTOCORRELATION WITHIN BSet-CONSTRAINED SEGMENTS")
print("=" * 70)

# When an orbit enters BSet, does the NEXT step always leave BSet?
# Or can two consecutive steps both be BSet?
bset_set = set(BSet)
NS_ORBIT = 10000

# Run orbit from a large number; restart from random points when hitting 1
import random as _rng
_rng.seed(99)
n = _rng.getrandbits(2000) | 1  # 2000-bit random odd number
bset_sequence = []
for _ in range(NS_ORBIT):
    r = n % MOD
    in_bset = 1 if r in bset_set else 0
    bset_sequence.append(in_bset)
    n_out, _, _ = macro_step(n)
    n = n_out
    if n == 1: n = _rng.getrandbits(2000) | 1  # restart with fresh large number

bset_seq = np.array(bset_sequence, dtype=float)
print(f"\nBSet occupancy rate: {bset_seq.mean():.5f}")
print(f"Theory: 15/128 = {15/128:.5f} (BSet has 15 elements out of 128 odd residues mod 256)")
print(f"Actually: 15 out of {MOD//2} odd residues mod {MOD}")
print(f"Theory: 15/{MOD//2} = {15/(MOD//2):.5f}")

print("\nBSet membership autocorrelation (lag 1 = consecutive step is also BSet?):")
print(f"{'Lag':>5} {'ACF':>10}")
bset_mean = bset_seq.mean()
bset_var = bset_seq.var()
for lag in range(1, 12):
    acf = np.corrcoef(bset_seq[:-lag], bset_seq[lag:])[0,1]
    print(f"{lag:>5} {acf:>10.6f}")

# Specifically: P(next in BSet | current in BSet)?
transitions_bb = 0  # BSet -> BSet
transitions_bn = 0  # BSet -> non-BSet
transitions_nb = 0  # non-BSet -> BSet
transitions_nn = 0  # non-BSet -> non-BSet
for i in range(len(bset_sequence)-1):
    a, b = bset_sequence[i], bset_sequence[i+1]
    if a == 1 and b == 1: transitions_bb += 1
    elif a == 1 and b == 0: transitions_bn += 1
    elif a == 0 and b == 1: transitions_nb += 1
    else: transitions_nn += 1

total_from_b = transitions_bb + transitions_bn
total_from_n = transitions_nb + transitions_nn
print(f"\nTransition matrix BSet vs non-BSet:")
print(f"  P(BSet->BSet) = {transitions_bb/total_from_b:.5f}")
print(f"  P(BSet->non-BSet) = {transitions_bn/total_from_b:.5f}")
print(f"  P(non-BSet->BSet) = {transitions_nb/total_from_n:.5f}")
print(f"  P(non-BSet->non-BSet) = {transitions_nn/total_from_n:.5f}")

# Spectral analysis of 2x2 transition matrix
T22 = np.array([
    [transitions_bb/total_from_b, transitions_bn/total_from_b],
    [transitions_nb/total_from_n, transitions_nn/total_from_n]
])
eigvals = np.linalg.eigvals(T22.T)
print(f"\nEigenvalues of 2x2 [BSet,nonBSet] transition matrix:")
for ev in sorted(eigvals, key=abs, reverse=True):
    print(f"  lambda = {ev.real:.6f}  (dominant = mixing rate of BSet indicator)")

print()
print("=" * 70)
print("PART 5: NEGATIVE EIGENVECTOR OF P_BSet — OSCILLATION PARTITION")
print("=" * 70)

# Rebuild the 15x15 P_BSet matrix empirically
NS5 = 3000
P_bset = np.zeros((15, 15))
bset_list = sorted(BSet)
bidx = {r: i for i, r in enumerate(bset_list)}
for i, r in enumerate(bset_list):
    counts = np.zeros(15)
    hits = 0
    for k in range(NS5):
        n = r + MOD*k
        cur = n
        for _ in range(200):  # walk until hitting BSet again
            cur_out, _, _ = macro_step(cur)
            r_out = cur_out % MOD
            if r_out in bset_set:
                counts[bidx[r_out]] += 1
                hits += 1
                break
            cur = cur_out
    if hits > 0:
        P_bset[i] = counts / hits

# Eigenvalue decomposition
eigvals, eigvecs = np.linalg.eig(P_bset.T)
# Sort by real part, then find the most negative
sorted_idx = np.argsort(eigvals.real)[::-1]
print("\nTop eigenvalues of P_BSet:")
for idx in sorted_idx[:6]:
    ev = eigvals[idx]
    print(f"  lambda = {ev.real:+.6f} + {ev.imag:+.6f}i  |l|={abs(ev):.6f}")

# Find the dominant negative real eigenvalue
neg_idx = None
for idx in sorted_idx:
    ev = eigvals[idx]
    if ev.real < -0.1 and abs(ev.imag) < 0.01:
        neg_idx = idx
        break

if neg_idx is not None:
    neg_ev = eigvals[neg_idx]
    neg_vec = eigvecs[:, neg_idx].real
    print(f"\nDominant negative eigenvalue: {neg_ev.real:.6f}")
    print(f"\nNegative eigenvector (reveals oscillation partition):")
    print(f"  {'Residue':>8} {'K':>4} {'l0':>6} {'Lyapunov':>12} {'EigVec':>10} {'Group':>8}")
    for i, r in enumerate(bset_list):
        K = bset_K_vals[r]
        l0_arr = bset_l0_stats[r]
        el0 = l0_arr.mean()
        lp = K*math.log(3) - (K + el0)*math.log(2)
        v = neg_vec[i]
        group = "A(+)" if v > 0 else "B(-)"
        print(f"  {r:>8} {K:>4} {el0:>6.2f} {lp:>12.4f} {v:>10.4f} {group:>8}")

    groupA = [bset_list[i] for i in range(15) if neg_vec[i] > 0]
    groupB = [bset_list[i] for i in range(15) if neg_vec[i] <= 0]
    print(f"\nGroup A (positive): {groupA}")
    print(f"Group B (negative): {groupB}")
    print(f"\nMixing mechanism: from A, the chain tends to go to B, and vice versa.")
    print(f"The -0.394 oscillation means each BSet-visit alternates A<->B.")
else:
    print("No dominant negative eigenvalue found; eigenvalues:", eigvals)

print()
print("=" * 70)
print("PART 6: BSet LYAPUNOV BALANCE")
print("=" * 70)
print()
print("The Lyapunov exponent from the full orbit is -0.575.")
print("BSet contributes ~11.7% of steps. Non-BSet contributes ~88.3% of steps.")
print()
print("BSet Lyapunov contributions:")
bset_lyapunovs = []
for r in sorted(BSet):
    K = bset_K_vals[r]
    l0_arr = bset_l0_stats[r]
    el0 = l0_arr.mean()
    lp = K*math.log(3) - (K + el0)*math.log(2)
    bset_lyapunovs.append(lp)
mean_bset_lp = np.mean(bset_lyapunovs)
print(f"  Mean Lyapunov over BSet elements: {mean_bset_lp:.4f}")
print(f"  Total Lyapunov = 0.117 * {mean_bset_lp:.4f} + 0.883 * L_nonBSet = -0.575")
L_nonBSet = (-0.575 - 0.117 * mean_bset_lp) / 0.883
print(f"  => L_nonBSet = {L_nonBSet:.4f}")
print(f"  Theoretical L_nonBSet (E[K]=2, E[l0]=2): {2*math.log(3) - 4*math.log(2):.4f}")

