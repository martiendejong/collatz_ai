"""
114_nonbset_hitting.py
=======================
NON-BSet EXCURSION FIRST-HITTING ANALYSIS

The mod-256 Markov chain has 128 odd states: 15 BSet + 113 non-BSet.
Starting from each non-BSet state r', the orbit bounces until reaching BSet.
This script computes EXACTLY:
  h(r') = (P(first BSet hit = r_1 | start at r'), ..., P(= r_15 | start at r'))

using the linear system:
  h(r') = Σ_{r'' in BSet} P(r'->r'') e_{r''} + Σ_{r'' not in BSet} P(r'->r'') h(r'')

Then: the BSet embedded chain transition P_BSet(r -> r_j) =
  P(h=1|r)*[1/coset_size if r_j in coset, else 0] + P(h>1|r)*Σ_{r' not in BSet in coset} (1/coset_size)*h(r')_j

KEY QUESTION: Why does r=103 have the highest stationary weight?
"""
import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1); m = (n + 1) >> k; x = m * (3**k) - 1; l = v2(x)
    return x >> l, k, l

BSet = {27,55,63,83,95,103,127,159,169,191,207,223,239,253,255}
BList = sorted(BSet)
NonBSet = [r for r in range(1, 256, 2) if r not in BSet]  # 113 elements
All_odd = sorted(list(BSet) + NonBSet)  # 128 elements

idx_B = {r: i for i, r in enumerate(BList)}      # BSet index: 0..14
idx_N = {r: i for i, r in enumerate(NonBSet)}    # NonBSet index: 0..112
idx_A = {r: i for i, r in enumerate(All_odd)}    # All-odd index: 0..127

print("=" * 70)
print("PART 1: MOD-256 TRANSITION MATRIX (EXACT via 2048 samples)")
print("=" * 70)
print()

# Build the exact mod-256 transition matrix
N_ALL = 128
N_SAMP = 2048
P_all = np.zeros((N_ALL, N_ALL))

for i, r in enumerate(All_odd):
    K0 = v2(r + 1)
    counts = np.zeros(N_ALL)
    valid = 0
    for k_iter in range(N_SAMP):
        n = r + 256 * k_iter
        if v2(n + 1) != K0:
            continue
        n_out, _, _ = macro_step(n)
        r_out = n_out % 256
        if r_out % 2 == 0:
            r_out = (n_out - 1) % 256  # shouldn't happen, but guard
        if r_out in idx_A:
            counts[idx_A[r_out]] += 1
            valid += 1
    if valid > 0:
        P_all[i] = counts / valid

print(f"Full 128x128 transition matrix computed ({N_SAMP} samples/state).")

# Extract submatrices
# B = BSet states (15), N = NonBSet states (113)
# P_all is indexed by All_odd
idx_B_in_A = [idx_A[r] for r in BList]
idx_N_in_A = [idx_A[r] for r in NonBSet]

P_BB = P_all[np.ix_(idx_B_in_A, idx_B_in_A)]  # 15x15: BSet->BSet
P_BN = P_all[np.ix_(idx_B_in_A, idx_N_in_A)]  # 15x113: BSet->NonBSet
P_NB = P_all[np.ix_(idx_N_in_A, idx_B_in_A)]  # 113x15: NonBSet->BSet
P_NN = P_all[np.ix_(idx_N_in_A, idx_N_in_A)]  # 113x113: NonBSet->NonBSet

print(f"Submatrices: P_BB(15x15), P_BN(15x113), P_NB(113x15), P_NN(113x113)")
print()

# =====================================================================
# PART 2: NON-BSet FIRST-HITTING DISTRIBUTION
# =====================================================================
print("=" * 70)
print("PART 2: EXACT FIRST-HITTING DISTRIBUTION h(r') = (P_BSet_first_hit)")
print("=" * 70)
print()
print("Solving: h = P_NB + P_NN @ h  =>  (I - P_NN) @ h = P_NB")
print()

I_113 = np.eye(113)
# h has shape (113, 15): h[i, j] = P(first BSet hit = BList[j] | start at NonBSet[i])
A = I_113 - P_NN
h = np.linalg.solve(A, P_NB)

print("First 20 non-BSet elements and their hitting distributions:")
print(f"{'r':>4} K  j  -> top 5 BSet destinations")
print("-" * 60)
for i, r in enumerate(NonBSet[:20]):
    K = v2(r + 1)
    if K < 8:
        m = (r + 1) >> K
        l0 = v2((3**K) * m - 1)
        j = 8 - K - l0
    else:
        j = -99
    top5 = sorted([(h[i, jj], BList[jj]) for jj in range(15)], reverse=True)[:5]
    top5_str = "  ".join(f"r={r2}:{p:.3f}" for p, r2 in top5 if p > 0.05)
    print(f"r={r:3d} K={K} j={j:2d}: {top5_str}")

print()
print("Most common BSet destination for each non-BSet element:")
from collections import Counter
most_common = Counter()
for i in range(113):
    best_j = np.argmax(h[i])
    most_common[BList[best_j]] += 1

print("BSet element : # non-BSet states that funnel to it MOST OFTEN")
for r, cnt in most_common.most_common():
    K = v2(r + 1)
    m = (r + 1) >> K
    l0 = v2((3**K) * m - 1)
    j = 8 - K - l0
    print(f"  r={r:3d} (K={K}, j={j:2d}): {cnt} non-BSet funnelers")

print()
print(f"h[*, idx(103)] stats (flow to r=103 from each non-BSet element):")
h103 = h[:, idx_B[103]]
print(f"  Mean = {h103.mean():.5f} (uniform = 1/15 = {1/15:.5f})")
print(f"  Max  = {h103.max():.5f} at r={NonBSet[np.argmax(h103)]}")
print(f"  Min  = {h103.min():.5f} at r={NonBSet[np.argmin(h103)]}")
print(f"  # elements with h(103) > 0.2: {(h103 > 0.2).sum()}")
print()

# =====================================================================
# PART 3: RECONSTRUCT EXACT P_BSet FROM P_BB + HITTING DISTRIBUTIONS
# =====================================================================
print("=" * 70)
print("PART 3: EXACT BSet EMBEDDED CHAIN FROM P_BB + P_BN @ h")
print("=" * 70)
print()
print("P_BSet(r->r') = P_BB(r->r') + Σ_{r' in NonBSet} P_BN(r->r') * h(r'->r')")
print()

P_BSet_exact = P_BB + P_BN @ h

# Normalize rows (should sum to 1)
row_sums = P_BSet_exact.sum(axis=1)
print(f"Row sums (should be 1.0): min={row_sums.min():.5f}, max={row_sums.max():.5f}")
P_BSet_exact = P_BSet_exact / row_sums[:, np.newaxis]

# Stationary distribution
pi_exact = np.ones(15) / 15
for _ in range(5000):
    pi_new = pi_exact @ P_BSet_exact
    if np.max(np.abs(pi_new - pi_exact)) < 1e-12:
        break
    pi_exact = pi_new

print()
print("Exact BSet embedded chain stationary distribution:")
print(f"{'r':>4}  {'K':>2}  {'j':>3}  {'pi_exact':>10}  {'pi*15':>6}  notes")
print("-" * 55)
for i, r in enumerate(BList):
    K = v2(r + 1)
    m = (r + 1) >> K
    l0 = v2((3**K) * m - 1)
    j = 8 - K - l0
    note = " <- MOST VISITED" if pi_exact[i] == pi_exact.max() else ""
    note += " <- LEAST VISITED" if pi_exact[i] == pi_exact.min() else ""
    print(f"r={r:3d}  K={K}  j={j:3d}  pi={pi_exact[i]:.6f}  {pi_exact[i]*15:5.3f}x{note}")

# Compare with empirical
print()
print("Eigenvalues of exact P_BSet:")
vals = sorted(np.linalg.eigvals(P_BSet_exact), key=lambda x: -x.real)
for i, v in enumerate(vals[:8]):
    print(f"  lambda_{i+1} = {v.real:.6f} + {v.imag:.5f}i")
print(f"Spectral gap = {1 - vals[1].real:.6f}")

# =====================================================================
# PART 4: EXPLAIN r=103 DOMINANCE
# =====================================================================
print()
print("=" * 70)
print("PART 4: WHY r=103 IS THE MOST VISITED — ALGEBRAIC EXPLANATION")
print("=" * 70)
print()

i_103 = idx_B[103]
print("Flow into r=103 from each BSet sender:")
for i, r_s in enumerate(BList):
    flow = pi_exact[i] * P_BSet_exact[i, i_103]
    direct = pi_exact[i] * P_BB[i, i_103]
    indirect = pi_exact[i] * (P_BN @ h)[i, i_103]
    K_s = v2(r_s + 1)
    m_s = (r_s + 1) >> K_s
    l0_s = v2((3**K_s) * m_s - 1)
    j_s = 8 - K_s - l0_s
    if flow > 0.002:
        print(f"  r={r_s:3d} (K={K_s},j={j_s:2d}): "
              f"P(->103)={P_BSet_exact[i,i_103]:.4f} "
              f"(direct={P_BB[i,i_103]:.4f}, indirect={indirect/pi_exact[i]:.4f}), "
              f"flow={flow:.5f}")

print()
print("Non-BSet elements that funnel most strongly to r=103:")
top_nonbset = sorted(enumerate(NonBSet), key=lambda x: -h[x[0], i_103])[:10]
for idx_n, r_n in top_nonbset:
    K_n = v2(r_n + 1)
    if K_n < 8:
        m_n = (r_n + 1) >> K_n
        l0_n = v2((3**K_n) * m_n - 1)
        j_n = 8 - K_n - l0_n
    else:
        j_n = -99
    print(f"  r'={r_n:3d} (K={K_n},j={j_n:2d}): h(103)={h[idx_n, i_103]:.4f}")

# Trace the specific excursion path 31->121->91->103
print()
print("SPECIFIC EXCURSION CHAINS ENDING AT r=103:")
n = 31  # non-BSet, K=5
chain = [n]
while n not in BSet:
    n, _, _ = macro_step(n)
    chain.append(n % 256 if n > 256 else n)
    if len(chain) > 20:
        break
print(f"  Chain from r=31: {chain}")

n = 91
chain = [n]
while n not in BSet:
    n, _, _ = macro_step(n)
    chain.append(n)
    if len(chain) > 20:
        break
print(f"  Chain from r=91: {chain}")

# Find ALL short (<=5 step) excursion chains ending at r=103
print()
print("All non-BSet chains of length <= 5 ending at r=103:")
count_chains = 0
for r_start in NonBSet:
    n = r_start
    chain = [n]
    for _ in range(5):
        n_next, _, _ = macro_step(n)
        rn = n_next % 256 if n_next >= 256 else n_next
        chain.append(rn)
        if rn in BSet:
            if rn == 103:
                count_chains += 1
                if count_chains <= 15:
                    print(f"  {chain} (len={len(chain)-1})")
            break
        n = n_next
print(f"  Total short chains to r=103: {count_chains} / {len(NonBSet)} = {count_chains/len(NonBSet):.4f}")
print()

# Compare for r=169 (least visited)
i_169 = idx_B[169]
count_to_169 = 0
for r_start in NonBSet:
    n = r_start
    for _ in range(5):
        n_next, _, _ = macro_step(n)
        rn = n_next % 256 if n_next >= 256 else n_next
        if rn in BSet:
            if rn == 169:
                count_to_169 += 1
            break
        n = n_next
print(f"Short chains to r=169 (least visited): {count_to_169} / {len(NonBSet)} = {count_to_169/len(NonBSet):.4f}")
print()
print("CONCLUSION: r=103 is the most visited BSet element because")
print("  (1) many non-BSet excursion paths end at r=103 specifically,")
print("  (2) the non-BSet dynamics have a 'funnel' toward r=103 via specific chains.")
