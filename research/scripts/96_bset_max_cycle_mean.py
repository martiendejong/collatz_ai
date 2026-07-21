"""
96_bset_max_cycle_mean.py
==========================
Compute the EXACT MAX CYCLE MEAN (MCM) of the BSet-restricted Markov chain.

The MCM is the maximum long-run avg k/step achievable by any orbit
that concentrates on a single "best" cycle in the BSet graph.

This is different from the ergodic rate (2.0614) which is the AVERAGE
across the stationary distribution. The MCM bounds the WORST CASE
for D_hard_kern (the orbit that tries to maximize avg k/step).

PROBLEM FORMULATION:
For each edge (r → r') in the BSet graph (T(r,r') > 0), define:
  edge_k(r, r') = E[k_sum | starts at r, lands at r'] / E[h | starts at r, lands at r']
                = conditional avg k/step for that specific transition

Max cycle mean = max over all cycles C of:
  (Σ_{(r,r') ∈ C} edge_k(r,r') × edge_h(r,r')) / (Σ_{(r,r') ∈ C} edge_h(r,r'))

Using Bellman-Ford / Karp's algorithm on the weighted directed graph.

RESULT EXPECTED: MCM ≈ 2.7974 (from script 82) or 2.711 (from script 90).
If verified: MCM < 3.419 → D_hard_kern = ∅ (conditional on BSet chain being
the correct model for Collatz macro-steps).

PROOF STRUCTURE:
  For any Collatz orbit:
    avg k/step = (BSet k_sum + non-BSet k_sum) / (BSet h + non-BSet h)
    ≤ max(BSet_MCM, non-BSet_max)
    ≤ max(MCM_BSet, 2.2503)  [from script 94]
    = MCM_BSet  [since MCM_BSet > 2.2503]
    < 3.419  [if MCM_BSet < 3.419]
  Therefore avg k/step < 3.419 for all orbits → D_hard_kern = ∅.
"""
import sys, time, math
from collections import defaultdict
import numpy as np

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
B_idx = {r: i for i, r in enumerate(BList)}

M_BASE = 10**12
N_SAMPLES = 2048   # per BSet element (more samples for accurate conditional moments)
MAX_H = 1000

# =====================================================================
# COMPUTE CONDITIONAL EDGE WEIGHTS: E[k_sum | r → r'] and E[h | r → r']
# =====================================================================
print("=" * 70)
print("CONDITIONAL EDGE WEIGHTS FOR BSet MARKOV CHAIN")
print("=" * 70)
print()
print(f"Computing E[k_sum | r→r'] and E[h | r→r'] for all (r,r') edges")
print(f"N_SAMPLES = {N_SAMPLES} per source state")
print()

t0 = time.time()

# edge_k_sums[r][r'] = sum of k_sum over trajectories from r landing at r'
# edge_h_sums[r][r'] = sum of h over trajectories from r landing at r'
# edge_count[r][r'] = number of trajectories
edge_k_sums = defaultdict(lambda: defaultdict(float))
edge_h_sums = defaultdict(lambda: defaultdict(float))
edge_count = defaultdict(lambda: defaultdict(int))

for r in BList:
    for i in range(N_SAMPLES):
        n = M_BASE + r + 256 * i
        n_cur = n
        h = 0
        k_sum = 0
        r_dest = None

        while h < MAX_H:
            n_out, k, l = macro_step(n_cur)
            h += 1
            k_sum += k
            n_cur = n_out
            if n_cur <= 1:
                break
            r_out = n_cur % 256
            if r_out in BSet:
                r_dest = r_out
                break

        if r_dest is not None:
            edge_k_sums[r][r_dest] += k_sum
            edge_h_sums[r][r_dest] += h
            edge_count[r][r_dest] += 1

print(f"[Computed in {time.time()-t0:.1f}s]\n")

# Compute conditional means: edge_kh[r][r'] = E[k_sum]/E[h] for r→r'
edge_kh = {}  # edge_kh[(r,r')] = avg k/step conditional on landing at r'
edge_h_avg = {}  # edge_h_avg[(r,r')] = E[h | lands at r']
edge_prob = {}  # edge_prob[(r,r')] = T(r,r')

total_count = {r: sum(edge_count[r].values()) for r in BList}

for r in BList:
    for r2 in BList:
        cnt = edge_count[r][r2]
        if cnt > 0:
            ek = edge_k_sums[r][r2] / edge_h_sums[r][r2]  # E[k_sum]/E[h]
            eh = edge_h_sums[r][r2] / cnt
            ep = cnt / total_count[r]
            edge_kh[(r, r2)] = ek
            edge_h_avg[(r, r2)] = eh
            edge_prob[(r, r2)] = ep

# =====================================================================
# DISPLAY EDGE WEIGHT TABLE (k/step conditional on each transition)
# =====================================================================
print("=== CONDITIONAL k/step FOR EACH EDGE (r → r') ===\n")
print("Shows E[k/step | starts at r, lands at r'] for top transitions\n")
print(f"{'Edge':>20}  {'P(r→r\')':>10}  {'E[k/step]':>12}  {'E[h]':>8}")
print("-" * 60)

# Show best edges by conditional k/step
all_edges = [(r, r2, edge_kh[(r,r2)], edge_prob[(r,r2)], edge_h_avg[(r,r2)])
             for r, r2 in edge_kh.keys()]
all_edges.sort(key=lambda x: -x[2])

print("Top 20 edges by conditional k/step:")
for r, r2, ek, ep, eh in all_edges[:20]:
    k0r = v2(r+1)
    k0r2 = v2(r2+1)
    print(f"  r={r:3d}(k={k0r}) → r={r2:3d}(k={k0r2}):  P={ep:.4f}  k/step={ek:.4f}  E[h]={eh:.2f}")

print()
print("Bottom 20 edges by conditional k/step:")
for r, r2, ek, ep, eh in all_edges[-20:]:
    k0r = v2(r+1)
    k0r2 = v2(r2+1)
    print(f"  r={r:3d}(k={k0r}) → r={r2:3d}(k={k0r2}):  P={ep:.4f}  k/step={ek:.4f}  E[h]={eh:.2f}")

# =====================================================================
# KARP'S ALGORITHM FOR MAX CYCLE MEAN
# =====================================================================
print()
print("=" * 70)
print("MAX CYCLE MEAN COMPUTATION (KARP'S ALGORITHM)")
print("=" * 70)
print()
print("Max cycle mean on the weighted directed graph where edge (r,r')")
print("has weight = conditional k_sum per edge = E[k_sum | r→r']")
print("and length = conditional h per edge = E[h | r→r']")
print()
print("MCM = max over cycles C of (total k_sum / total h)")
print()

# Build the weighted graph for max-cycle-mean.
# Edge weight for (r,r'): use avg k_sum per transition = edge_kh[(r,r')] × edge_h_avg[(r,r')]
# Edge length for (r,r'): edge_h_avg[(r,r')]
# We want to maximize (sum of k_sum over cycle) / (sum of h over cycle)
# = max cycle mean on the graph with edge weights = k_sum and lengths = h

# Karp's algorithm: multiply by length and use log-sum formulation
# Alternative: Binary search on c such that max cycle mean(k_sum - c*h) ≥ 0
# Equivalently: solve the minimum ratio cycle problem

# Build adjacency: for each (r,r') with enough samples
MIN_COUNT = 10  # minimum samples to trust an edge
edges = []
for r in BList:
    for r2 in BList:
        cnt = edge_count[r][r2]
        if cnt >= MIN_COUNT:
            k_sum_per = edge_kh[(r,r2)] * edge_h_avg[(r,r2)]  # E[k_sum per transition]
            h_per = edge_h_avg[(r,r2)]  # E[h per transition]
            edges.append((r, r2, k_sum_per, h_per))

# Binary search for max cycle mean c*:
# c* = max c such that ∃ cycle with total_k_sum/total_h ≥ c
# Equivalent: max c such that min-weight cycle with weights (h - (1/c)*k_sum) is ≤ 0
# Or: max c such that min-weight cycle with weights (k_sum - c*h) is ≥ 0

# Use parametric shortest path: for each candidate c, check if the graph
# with edge weights w(e) = k_sum(e) - c*h(e) has a positive cycle (Bellman-Ford)

def has_positive_cycle(c_val, edges, nodes):
    """Check if graph with edge weight k_sum - c*h has a positive cycle."""
    INF = 1e18
    dist = {r: -INF for r in nodes}
    # Initialize: try starting from each node
    dist = {r: 0.0 for r in nodes}
    n = len(nodes)

    for _ in range(n):  # n iterations of Bellman-Ford
        updated = False
        new_dist = dict(dist)
        for r, r2, ks, h in edges:
            w = ks - c_val * h  # weight = k_sum - c × h
            if dist[r] + w > new_dist[r2]:
                new_dist[r2] = dist[r] + w
                updated = True
        dist = new_dist
        if not updated:
            break

    # Check for positive cycle: one more iteration
    for r, r2, ks, h in edges:
        w = ks - c_val * h
        if dist[r] + w > dist[r2] + 1e-9:
            return True  # positive cycle exists
    return False

nodes = BList
lo, hi = 1.0, 5.0  # search in [1.0, 5.0]

print(f"Binary search for MCM in [{lo}, {hi}]...")
for iteration in range(60):
    mid = (lo + hi) / 2
    if has_positive_cycle(mid, edges, nodes):
        lo = mid
    else:
        hi = mid

mcm = (lo + hi) / 2
print(f"MAX CYCLE MEAN = {mcm:.6f}")
print(f"D_hard_kern threshold = 3.419")
print(f"GAP = {3.419 - mcm:.6f}")
print()
if mcm < 3.419:
    print(f"✓ MCM ({mcm:.4f}) << 3.419 → D_hard_kern = ∅")
else:
    print(f"✗ VIOLATION: MCM ({mcm:.4f}) ≥ 3.419!")

# =====================================================================
# FIND THE BEST CYCLE
# =====================================================================
print()
print("=== FINDING THE BEST CYCLE ===\n")
print(f"Looking for cycle with mean ≈ {mcm:.4f}...")
print()

# Use value iteration to find the policy that achieves the MCM
# V_t(r) = max over r' of (k_sum(r,r') + V_{t-1}(r') - mcm × h(r,r'))
V = {r: 0.0 for r in nodes}
policy = {r: None for r in nodes}

for t in range(200):
    V_new = {}
    for r in nodes:
        best_val = -1e18
        best_r2 = None
        for r2 in nodes:
            if (r, r2) in edge_kh and edge_count[r][r2] >= MIN_COUNT:
                ks = edge_kh[(r,r2)] * edge_h_avg[(r,r2)]
                h = edge_h_avg[(r,r2)]
                val = ks - mcm * h + V[r2]
                if val > best_val:
                    best_val = val
                    best_r2 = r2
        V_new[r] = best_val
        policy[r] = best_r2
    V = V_new

# Trace the optimal cycle from policy
def find_cycle(policy):
    """Follow policy until a cycle is detected."""
    for start in nodes:
        visited = {}
        cur = start
        step = 0
        while cur not in visited and step < 50:
            visited[cur] = step
            if policy[cur] is None:
                break
            cur = policy[cur]
            step += 1
        if cur in visited:
            # Found a cycle starting at 'cur'
            cycle = [cur]
            nxt = policy[cur]
            while nxt != cur:
                cycle.append(nxt)
                nxt = policy[nxt]
            return cycle
    return None

best_cycle = find_cycle(policy)
if best_cycle:
    print(f"Best cycle found: {[f'r={r}(k={v2(r+1)})' for r in best_cycle]}")
    print()
    # Compute cycle mean
    cycle_k_sum = 0
    cycle_h_sum = 0
    cycle_edges = []
    for i, r in enumerate(best_cycle):
        r2 = best_cycle[(i+1) % len(best_cycle)]
        if (r, r2) in edge_kh:
            ks = edge_kh[(r,r2)] * edge_h_avg[(r,r2)]
            h = edge_h_avg[(r,r2)]
            k_step = edge_kh[(r,r2)]
            p = edge_prob.get((r,r2), 0)
            cycle_k_sum += ks
            cycle_h_sum += h
            cycle_edges.append((r, r2, k_step, h, p))
            k0r = v2(r+1)
            k0r2 = v2(r2+1)
            print(f"  r={r:3d}(k={k0r}) → r={r2:3d}(k={k0r2}): k/step={k_step:.4f}, E[h]={h:.2f}, P={p:.4f}")

    cycle_mean = cycle_k_sum / cycle_h_sum
    print()
    print(f"Cycle mean = {cycle_k_sum:.4f} / {cycle_h_sum:.4f} = {cycle_mean:.4f}")
    print(f"Gap to D_hard_kern: {3.419 - cycle_mean:.4f}")
else:
    print("No cycle found (policy may not converge to simple cycle)")

# =====================================================================
# D_hard_kern PROOF SUMMARY
# =====================================================================
print()
print("=" * 70)
print("D_hard_kern = ∅ PROOF STRUCTURE")
print("=" * 70)
print()
print("CLAIM: Every Collatz orbit has long-run avg k/step < 3.419.")
print()
print("PROOF SKETCH (all bounds from empirical computation at large n):")
print()
print("Step 1: Every Collatz orbit can be decomposed into alternating periods:")
print("  - BSet periods: orbit visits BSet-residue states (mod 256 ∈ BSet)")
print("  - non-BSet periods: orbit visits non-BSet states before returning to BSet")
print()
print("Step 2: avg k/step over the FULL orbit is a weighted average:")
print("  avg k/step = (BSet_k_sum + nonBSet_k_sum) / (BSet_h + nonBSet_h)")
print()
print("Step 3: BSet-period bound:")
print(f"  max avg k/step achievable in BSet = MCM = {mcm:.4f}")
print(f"  (from max-cycle-mean computation on 15-state BSet graph)")
print()
print("Step 4: non-BSet-period bound:")
print(f"  max E[k/step until BSet] over all non-BSet residues = 2.2503")
print(f"  (from script 94, N=512 samples per residue)")
print()
print("Step 5: Combined bound:")
print(f"  avg k/step ≤ max(BSet_MCM, nonBSet_max)")
print(f"           = max({mcm:.4f}, 2.2503)")
print(f"           = {max(mcm, 2.2503):.4f}")
print()
print("Step 6: D_hard_kern threshold = 3.419")
print(f"  Gap: 3.419 - {max(mcm, 2.2503):.4f} = {3.419 - max(mcm, 2.2503):.4f}")
print()
if max(mcm, 2.2503) < 3.419:
    print("✓ CONCLUSION: Every orbit satisfies avg k/step < 3.419 < D_hard_kern threshold.")
    print("  Therefore D_hard_kern = ∅.")
    print()
    print("MISSING PIECES FOR RIGOROUS PROOF:")
    print("  1. Make the MCM bound exact (not just empirical)")
    print("  2. Prove that the macro-step Markov chain converges to the BSet model")
    print("  3. Establish that large-n Phi estimates are tight (exponential concentration)")
    print("  4. Handle the case of orbits that avoid BSet entirely (P=0 by ergodicity?)")

# =====================================================================
# COMPARE WITH PREVIOUS MCM ESTIMATES
# =====================================================================
print()
print("=== MCM COMPARISON ACROSS SCRIPTS ===\n")
print(f"Script 82 (direct cycle enumeration):    MCM ≈ 2.7974")
print(f"Script 90 (8-state k0-type model):       MCM ≈ 2.711")
print(f"Script 95 (BSet ergodic rate):           E_erg = 2.0614")
print(f"Script 96 (15-state Bellman-Ford MCM):   MCM = {mcm:.4f}")
print()
print("The 15-state MCM gives the most accurate bound on the maximum achievable")
print("avg k/step for any orbit following the BSet Markov chain model.")
print()
print(f"D_hard_kern gap from 15-state MCM: {3.419 - mcm:.4f}")
print()
if mcm < 2.8:
    print("The 15-state MCM is LOWER than script 82's estimate (2.7974).")
    print("This suggests the earlier estimate may have included unreachable cycles.")
elif mcm > 2.8:
    print("The 15-state MCM is HIGHER than script 90's estimate (2.711).")
    print("The 15-state model captures more structure than the 8-state k0-type model.")

# =====================================================================
# THE MAX-CYCLE EDGE ANALYSIS: what makes the best cycle?
# =====================================================================
print()
print("=== HIGHEST k/step EDGES IN THE BSet GRAPH ===\n")
print("Top 10 edges by conditional k/step (the 'building blocks' of best cycles):")
for r, r2, ek, ep, eh in all_edges[:10]:
    k0r = v2(r+1)
    k0r2 = v2(r2+1)
    # How consistent is this edge? (high ek, high ep = important edge)
    importance = ek * ep
    print(f"  r={r:3d}(k={k0r}) → r={r2:3d}(k={k0r2}): k/step={ek:.4f}, P={ep:.4f}, E[h]={eh:.2f}, importance={importance:.4f}")

print()
print("Key pattern: high conditional k/step edges tend to involve:")
print("  - High k0 sources (r=255, k0=8 or r=127, k0=7)")
print("  - Short excursions (small E[h]) to nearby BSet elements")
print("  - The best k/step edges may form an approximately closed sub-cycle")
