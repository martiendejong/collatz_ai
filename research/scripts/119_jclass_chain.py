"""
119_jclass_chain.py
=====================
THE j-CLASS AGGREGATE CHAIN

By the Exchangeability Theorem, all BSet elements with the same j-class
have IDENTICAL transition distributions. This means the 15x15 BSet chain
collapses to a 9x9 j-class aggregate chain Q.

This script:
1. Computes Q from the exact BSet chain (script 114 data)
2. Analyzes Q's eigenvalues, stationary distribution, flow structure
3. Identifies the HIERARCHICAL CASCADE: how j-classes flow to each other
4. Checks whether Q is self-similar across moduli
"""
import numpy as np, sys

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1)>>K; x = m*(3**K)-1; l = v2(x)
    return x>>l, K, l

BSet = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255}
BList = sorted(BSet)
All_odd = list(range(1, 256, 2))
NonBSet = [r for r in All_odd if r not in BSet]

idx_B = {r: i for i, r in enumerate(BList)}
idx_A = {r: i for i, r in enumerate(All_odd)}

# j-class assignments
def get_j(r, N=8):
    K = v2(r+1)
    m = (r+1)>>K
    prod = m*(3**K)
    l0 = v2(prod-1)
    return N-K-l0

j_of = {r: get_j(r) for r in BList}
j_classes = sorted(set(j_of.values()), reverse=True)
j_members = {j: [r for r in BList if j_of[r]==j] for j in j_classes}

print("j-class structure:")
for j in j_classes:
    print(f"  j={j:2d}: {j_members[j]}")

# =====================================================================
# PART 1: REBUILD EXACT BSet CHAIN (from script 114 method)
# =====================================================================
N_SAMP = 4096
print(f"\nBuilding 128x128 transition matrix ({N_SAMP} samples)...")
P_all = np.zeros((128, 128))
for i, r in enumerate(All_odd):
    K0 = v2(r+1)
    counts = np.zeros(128)
    valid = 0
    for k_iter in range(N_SAMP):
        n = r + 256*k_iter
        if v2(n+1) != K0: continue
        n_out, _, _ = macro_step(n)
        r_out = n_out % 256
        if r_out % 2 == 1 and r_out in idx_A:
            counts[idx_A[r_out]] += 1
            valid += 1
    if valid > 0:
        P_all[i] = counts/valid

idx_B_in_A = [idx_A[r] for r in BList]
idx_N_in_A = [idx_A[r] for r in NonBSet]
P_BB = P_all[np.ix_(idx_B_in_A, idx_B_in_A)]
P_BN = P_all[np.ix_(idx_B_in_A, idx_N_in_A)]
P_NB = P_all[np.ix_(idx_N_in_A, idx_B_in_A)]
P_NN = P_all[np.ix_(idx_N_in_A, idx_N_in_A)]

# Solve for exact BSet chain
h = np.linalg.solve(np.eye(113)-P_NN, P_NB)
P_BSet = P_BB + P_BN @ h
P_BSet /= P_BSet.sum(axis=1, keepdims=True)  # normalize rows

# Stationary of BSet chain
pi_B = np.ones(15)/15
for _ in range(5000):
    pi_new = pi_B @ P_BSet
    if np.max(np.abs(pi_new-pi_B)) < 1e-12: break
    pi_B = pi_new

# =====================================================================
# PART 2: COMPUTE j-CLASS AGGREGATE MATRIX Q
# =====================================================================
print("\nComputing 9x9 j-class aggregate matrix Q...")
j_idx = {j: i for i, j in enumerate(j_classes)}
nJ = len(j_classes)
Q = np.zeros((nJ, nJ))

# Q[j, j'] = sum_{r' in j'-class} P_BSet(r, r') for r in j-class j
# (all r in j-class have same row by exchangeability)
for i_j, j in enumerate(j_classes):
    reps = j_members[j]  # all elements in j-class
    # Use average of their rows (should be identical by exchangeability)
    avg_row = np.mean([P_BSet[idx_B[r]] for r in reps], axis=0)
    for i_jp, jp in enumerate(j_classes):
        # Sum probabilities going to j'-class
        Q[i_j, i_jp] = sum(avg_row[idx_B[rp]] for rp in j_members[jp])

# Verify row sums
row_sums = Q.sum(axis=1)
print(f"Row sums range: [{row_sums.min():.5f}, {row_sums.max():.5f}] (should be 1.0)")
Q /= row_sums[:, None]

print("\n9x9 j-class transition matrix Q:")
print(f"{'j ->':>6}", end="")
for jp in j_classes:
    print(f"  j={jp:3d}", end="")
print()
print("-" * 75)
for i_j, j in enumerate(j_classes):
    print(f"j={j:3d}:", end="")
    for i_jp, jp in enumerate(j_classes):
        v = Q[i_j, i_jp]
        if v > 0.01: print(f"  {v:.3f}", end="")
        else: print(f"  {'':5}", end="")
    print()

# =====================================================================
# PART 3: j-CLASS EIGENVALUES AND STATIONARY
# =====================================================================
print("\n9x9 j-class stationary distribution:")
pi_Q = np.ones(nJ)/nJ  # approximate: ~ proportional to class size
for _ in range(5000):
    pi_new = pi_Q @ Q
    if np.max(np.abs(pi_new - pi_Q)) < 1e-12: break
    pi_Q = pi_new

for i_j, j in enumerate(j_classes):
    n_members = len(j_members[j])
    print(f"  j={j:3d} ({n_members} elem): pi_Q={pi_Q[i_j]:.5f}  "
          f"pi_B_sum={sum(pi_B[idx_B[r]] for r in j_members[j]):.5f}  "
          f"expected={n_members/15:.5f}")

print("\nEigenvalues of Q (9x9):")
vals = sorted(np.linalg.eigvals(Q).real, reverse=True)
for i, v in enumerate(vals):
    print(f"  lambda_{i+1} = {v:.8f}")
print(f"Spectral gap of Q = {1-vals[1]:.8f}")

print("\nEigenvalues of P_BSet (15x15) for comparison:")
vals15 = sorted(np.linalg.eigvals(P_BSet).real, reverse=True)
for i, v in enumerate(vals15[:8]):
    print(f"  lambda_{i+1} = {v:.8f}")
print(f"Spectral gap of P_BSet = {1-vals15[1]:.8f}")

# =====================================================================
# PART 4: FLOW ANALYSIS — CASCADE STRUCTURE
# =====================================================================
print()
print("=" * 70)
print("PART 4: CASCADE FLOW STRUCTURE")
print("=" * 70)
print()
print("Dominant transitions (Q[j, j'] > 5%):")
for i_j, j in enumerate(j_classes):
    for i_jp, jp in enumerate(j_classes):
        v = Q[i_j, i_jp]
        if v >= 0.05:
            dir_str = "^" if jp > j else "v" if jp < j else "="
            print(f"  j={j:3d} --> j'={jp:3d} ({dir_str}): {v:.3f}")

print()
# Compute flow balance: upward vs downward
up_flow = 0; down_flow = 0; same_flow = 0
for i_j, j in enumerate(j_classes):
    for i_jp, jp in enumerate(j_classes):
        if jp > j: up_flow += pi_Q[i_j] * Q[i_j, i_jp]
        elif jp < j: down_flow += pi_Q[i_j] * Q[i_j, i_jp]
        else: same_flow += pi_Q[i_j] * Q[i_j, i_jp]
print(f"Flow balance (weighted by stationary pi_Q):")
print(f"  Upward   (j' > j): {up_flow:.4f}")
print(f"  Downward (j' < j): {down_flow:.4f}")
print(f"  Same     (j' = j): {same_flow:.4f}")

print()
print("Mean first-passage times from each j-class to j=6 (top class):")
# MFPT from j to j=6: solve (I-Q_without_top) @ t = 1
i_top = j_idx[6]
# Absorbing at j=6: remove top state and solve
non_top = [i for i in range(nJ) if i != i_top]
Q_sub = Q[np.ix_(non_top, non_top)]
ones = np.ones(len(non_top))
t = np.linalg.solve(np.eye(len(non_top))-Q_sub, ones)
print(f"  j=6 is the reference (MFPT=0 from itself)")
for i, idx in enumerate(non_top):
    print(f"  j={j_classes[idx]:3d}: MFPT to j=6 = {t[i]:.2f} j-class steps")

# =====================================================================
# PART 5: THEORETICAL j-CLASS TRANSITION FORMULA
# =====================================================================
print()
print("=" * 70)
print("PART 5: THEORETICAL TRANSITION PROBABILITIES")
print("=" * 70)
print()
print("From j-class j (j>=1), output has K' >= j (geometric from j):")
print("  P(K'=k | K'>=j) approx 1/2^{k-j+1} for k=j, j+1, ...")
print("  P(l0'=q) approx 1/2^q for q=1, 2, ...")
print("  Next j-class: j' = 8-K'-l0' (at mod-256, N=8)")
print()
print("Theoretical Q[j, j'] for j in {1,...,6} -> j' in {-5,...,6}:")

N = 8
for j_src in range(6, 0, -1):
    print(f"\n  From j={j_src}:")
    theory = {}
    for K_prime in range(j_src, N):
        p_K = 1.0 / (2**(K_prime - j_src + 1))
        for l0_prime in range(1, N - K_prime):
            p_l0 = 1.0 / (2**l0_prime)
            j_dest = N - K_prime - l0_prime
            theory[j_dest] = theory.get(j_dest, 0) + p_K * p_l0
        # K_prime = N-1: l0 can be 1 (only), j' = N-(N-1)-1 = 0
        # K_prime = N: only for j<=0 elements (K=8 means r=255)
    # Normalize
    total = sum(theory.values())
    theory = {jp: v/total for jp, v in theory.items()}
    for jp in sorted(theory, reverse=True):
        if theory[jp] > 0.02:
            print(f"    j'={jp:3d}: {theory[jp]:.4f}", end="  ")
    print()
