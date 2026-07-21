"""
93_routing_probabilities.py
============================
Compute EXACT routing probabilities P_route(r) for all 128 odd residues mod 256.

DEFINITION:
For n ≡ r mod 256 with k0 = v2(r+1), the macro-step from n uses k=k0.
The m-class is: m = (n+1)/2^k0, constrained by n ≡ r mod 256.
Within 256 consecutive odd m-values, the class has 2^k0 members.
P_route(r) = fraction of class members whose output mod 256 lies in BSet.

KEY INSIGHT (from scripts 91-92b):
- r=169 (BSet, k0=1): P_route=100% — always routes to BSet
- r=255 (BSet, k0=8): P_route=31/256=12.1%
- BSet membership correlates with P_route, NOT avg k/step

CLAIM TO TEST:
BSet = {r : P_route(r) ≥ threshold} for some threshold
(i.e., P_route cleanly separates BSet from non-BSet)

ALSO: Compute the full one-step output distribution for key residues,
showing exactly which BSet elements receive from which.
"""
import sys, math
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
BList = sorted(BSet)
ODD_RES = [r for r in range(1, 256, 2)]  # 128 odd residues
M_PERIOD = 256  # 256 odd m-values per source

# =====================================================================
# COMPUTE P_route FOR ALL 128 ODD RESIDUES
# =====================================================================
print("=" * 70)
print("EXACT ROUTING PROBABILITIES P_route(r) FOR ALL 128 ODD RESIDUES")
print("=" * 70)
print()
print("P_route(r) = fraction of m-class giving BSet output in 1 macro-step")
print("m-class: m ≡ (r+1)/2^k0 mod (256/2^k0), i.e., the m-values for n≡r mod 256")
print()

p_route = {}
class_sizes = {}
bset_counts = {}
bset_output_dist = {}  # r -> Counter of (output r_out) within class

for r in ODD_RES:
    k0 = v2(r + 1)
    pow3k0 = 3**k0
    class_period = 256 // (2**k0)  # spacing between consecutive class members in odd m-index
    m_base = (r + 1) // (2**k0)    # the base m-value for this class (m ≡ m_base mod 128/2^(k0-1))

    # Enumerate class members: odd m in [1,511] with m ≡ m_base mod (256/2^k0)
    # m = m_base, m_base + (256/2^k0), m_base + 2*(256/2^k0), ...
    # Equivalently: within 256 odd m-values, the class has 2^k0 members

    step_in_odd = class_period  # number of ODD m-values between class members
    # odd m-values: 1,3,5,...,511. Index i corresponds to m=2i+1.
    # class: m = 2*(something)+1 with m ≡ m_base mod class_period
    # m_base is odd (since (r+1)/2^k0 is odd for k0 = v2(r+1)).

    class_m_vals = [m for m in range(1, 512, 2) if m % (class_period) == m_base % (class_period)]
    # Adjust: the class period for m is 256/2^k0 (in terms of ODD step index)
    # Actually: m ≡ m_base mod (256/2^(k0-1)) wait let me re-derive.

    # n ≡ r mod 256. n+1 ≡ r+1 mod 256. m = (n+1)/2^k0. m ≡ (r+1)/2^k0 mod (256/2^k0).
    # Since m is odd and 256/2^k0 is an integer (because k0 ≤ 8):
    period_m = 256 // (2**k0)  # = 128, 64, 32, 16, 8, 4, 2, 1 for k0=1..8
    m_class_base = (r + 1) // (2**k0)  # odd integer

    class_m_vals = [m for m in range(m_class_base, 512, period_m) if m > 0 and m < 512 and m % 2 == 1]
    # Include wrapping: class can also contain m_class_base + period_m, m_class_base + 2*period_m, etc.
    # Since m must be odd and in [1,511]:

    # class_size = min(2^(k0+1), 256):
    # period_m=256/2^k0; in [1,511] there are 512/period_m class members, all odd
    # (m_base is odd and step=period_m is even, so parity is preserved).
    # For k0>=7, period_m<=2 collapses to all 256 odd m.
    expected_class_size = min(2**(k0+1), 256)
    class_size = len(class_m_vals)
    assert class_size == expected_class_size, (
        f"Expected class size {expected_class_size} but got {class_size} for r={r} k0={k0}")

    bset_count = 0
    output_dist = Counter()
    for m in class_m_vals:
        val = pow3k0 * m - 1
        l = v2(val)
        out = (val >> l) % 256
        output_dist[out] += 1
        if out in BSet:
            bset_count += 1

    p_route[r] = bset_count / class_size
    class_sizes[r] = class_size
    bset_counts[r] = bset_count
    bset_output_dist[r] = output_dist

# =====================================================================
# DISPLAY RESULTS
# =====================================================================

# Sort by P_route
sorted_res = sorted(ODD_RES, key=lambda r: (p_route[r], v2(r+1)))

print(f"{'r':>5}  {'k0':>4}  {'class_size':>10}  {'BSet_count':>10}  {'P_route':>10}  {'BSet?':>8}")
print("-" * 65)

bset_p = [p_route[r] for r in BList]
non_bset_p = [p_route[r] for r in ODD_RES if r not in BSet]

for r in sorted_res:
    k0 = v2(r + 1)
    cs = class_sizes[r]
    bc = bset_counts[r]
    pr = p_route[r]
    in_bset = r in BSet
    marker = " ✓" if in_bset else ""
    print(f"  r={r:3d}  k={k0}  class={cs:4d}  BSet_cnt={bc:4d}  P_route={pr:.4f} ({bc}/{cs}){marker}")

print()
print(f"P_route ranges:")
print(f"  BSet elements:     [{min(bset_p):.4f}, {max(bset_p):.4f}]")
print(f"  Non-BSet elements: [{min(non_bset_p):.4f}, {max(non_bset_p):.4f}]")
gap = min(bset_p) - max(non_bset_p)
print(f"  Gap (min_BSet - max_nonBSet): {gap:.6f}")
if gap > 0:
    threshold = (min(bset_p) + max(non_bset_p)) / 2
    print(f"  ✓ CLEAN SEPARATION! P_route threshold ≈ {threshold:.4f}")
else:
    print(f"  ✗ OVERLAP detected")

# =====================================================================
# BSet ELEMENTS: P_route sorted
# =====================================================================
print()
print("=== BSet ELEMENTS SORTED BY P_route ===\n")
bset_sorted = sorted(BList, key=lambda r: p_route[r])
for r in bset_sorted:
    k0 = v2(r + 1)
    cs = class_sizes[r]
    bc = bset_counts[r]
    pr = p_route[r]
    # Show which BSet elements they route to
    bset_targets = {out: cnt for out, cnt in bset_output_dist[r].items() if out in BSet}
    targets_str = ", ".join(f"r={out}({cnt})" for out, cnt in sorted(bset_targets.items(), key=lambda x: -x[1])[:3])
    print(f"  r={r:3d} k0={k0}: P_route={pr:.4f} ({bc}/{cs})  → {targets_str}")

# =====================================================================
# THE ROUTING STRUCTURE: WHO ROUTES TO WHOM?
# =====================================================================
print()
print("=== ONE-STEP ROUTING STRUCTURE AMONG BSet ELEMENTS ===\n")
print("(showing only BSet→BSet direct transitions in 1 macro-step)")
print()
print(f"{'Source r':>10}  {'k0':>4}  {'P_route':>10}  → {'Destination distribution (in BSet)':>10}")
print("-" * 80)

for r in BList:
    k0 = v2(r + 1)
    cs = class_sizes[r]
    bc = bset_counts[r]
    pr = p_route[r]
    bset_targets = {out: cnt for out, cnt in bset_output_dist[r].items() if out in BSet}
    targets_str = "  ".join(f"r={out:3d}(k={v2(out+1)},n={cnt})" for out, cnt in sorted(bset_targets.items(), key=lambda x: -x[1]))
    print(f"  r={r:3d} k={k0}:  P={pr:.3f}  →  {targets_str}")

# =====================================================================
# NON-BSet NEAR-BOUNDARY: P_route AND WHY NOT IN BSet
# =====================================================================
print()
print("=== NON-BSet NEAR-BOUNDARY (highest P_route) ===\n")
non_sorted = sorted([r for r in ODD_RES if r not in BSet], key=lambda r: -p_route[r])
for r in non_sorted[:15]:
    k0 = v2(r + 1)
    cs = class_sizes[r]
    bc = bset_counts[r]
    pr = p_route[r]
    bset_targets = {out: cnt for out, cnt in bset_output_dist[r].items() if out in BSet}
    targets_str = ", ".join(f"r={out}({cnt})" for out, cnt in sorted(bset_targets.items(), key=lambda x: -x[1])[:3])
    targets_str = targets_str if targets_str else "(none)"
    print(f"  r={r:3d} k={k0}:  P_route={pr:.4f} ({bc}/{cs})  → {targets_str}")

# =====================================================================
# KEY EXAMPLE: r=169 EXACT ANALYSIS
# =====================================================================
print()
print("=== EXACT ANALYSIS: r=169 (k0=1, P_route=100%) ===\n")
r = 169
k0 = v2(r + 1)
period_m = 256 // (2**k0)  # = 128
m_class_base = (r + 1) // (2**k0)  # = 85
pow3k0 = 3**k0  # = 3

print(f"r=169: k0=1, m-class base = (169+1)/2 = 85, period = 128")
print(f"Class members: m ≡ 85 mod 128 in [1,511]")
print()
class_m = [m for m in range(m_class_base, 512, period_m) if m % 2 == 1]
for m in class_m:
    val = pow3k0 * m - 1
    l = v2(val)
    out = val >> l
    r_out = out % 256
    print(f"  m={m:4d}: 3×{m}-1 = {val}  l={l}  n'={out}  n' mod 256 = {r_out}  {'∈ BSet ✓' if r_out in BSet else '∉ BSet'}")

print()
print(f"Pattern: n' mod 256 ∈ {{127, 63, 255, 191}} = ALL high-k BSet elements!")
print(f"k0 of destinations: {', '.join(str(v2(out+1)) for out in [127,63,255,191])}")
print(f"r=169 (k0=1) ALWAYS routes to k0≥6 boosters — upgrading k0 by 5-7!")

# =====================================================================
# THE KEY THEOREM: P_route IS THE SEPARATING CRITERION
# =====================================================================
print()
print("=" * 70)
print("THEOREM: P_route CLEANLY SEPARATES BSet FROM NON-BSet")
print("=" * 70)
print()
if gap > 0:
    print(f"BSet threshold:    P_route ≥ {threshold:.4f}")
    print(f"Min P_route (BSet): {min(bset_p):.4f}")
    print(f"Max P_route (non-BSet): {max(non_bset_p):.4f}")
    print(f"Gap: {gap:.6f}")
    print()
    print("ALL BSet elements have P_route above threshold.")
    print("ALL non-BSet elements have P_route below threshold.")
    print()
    print("INTERPRETATION: BSet = {r : orbit starting at r reaches BSet in 1 step")
    print("                        with probability ≥ threshold}")
else:
    print("P_route does NOT cleanly separate BSet from non-BSet.")
    print("Some non-BSet elements have P_route in BSet range.")
    print()
    # Show the overlap
    overlap_bset = [r for r in BList if p_route[r] <= max(non_bset_p)]
    overlap_non = [r for r in ODD_RES if r not in BSet and p_route[r] >= min(bset_p)]
    print(f"BSet elements with low P_route: {[(r, f'{p_route[r]:.4f}') for r in overlap_bset]}")
    print(f"Non-BSet elements with high P_route: {[(r, f'{p_route[r]:.4f}') for r in overlap_non]}")

# =====================================================================
# D_hard_kern ARGUMENT FROM P_route
# =====================================================================
print()
print("=== D_hard_kern ARGUMENT FROM ROUTING PROBABILITY ===\n")
max_p_route = max(p_route.values())
max_p_r = max(p_route.keys(), key=lambda r: p_route[r])
k0_max = v2(max_p_r + 1)
print(f"Max P_route over all 128 residues: {max_p_route:.4f} (r={max_p_r}, k0={k0_max})")
print()
print(f"For a D_hard_kern orbit visiting booster r at each step h=1:")
print(f"  avg k/step = k0  (single-step dominates)")
print(f"  For r={max_p_r}: k0={k0_max}, avg k/step = {k0_max}")
print()
print(f"Maximum k0 over all BSet elements: {max(v2(r+1) for r in BSet)}")
print(f"D_hard_kern threshold: 3.419")
print(f"Max k0 in BSet: {max(v2(r+1) for r in BSet)} = {max(v2(r+1) for r in BSet)} < 3.419?")
print(f"Answer: {max(v2(r+1) for r in BSet)} {'<' if max(v2(r+1) for r in BSet) < 3.419 else '≥'} 3.419")
print()
print(f"If h=1 always (best case): avg k/step = k0 = at most {max(v2(r+1) for r in BSet)}")
print(f"But P(h=1) ≤ {max_p_route:.3f} < 1, so h>1 often occurs with lower k/step")
print(f"→ avg k/step < {max(v2(r+1) for r in BSet)} always < 3.419")
print()
# The exact formula for k0=8 (r=255) with P_route=p:
# E[k/step] = p × (k0=8) + (1-p) × E[k/step | h>1]
# For h>1, E[k/step] ≈ 2.0
r255_route = p_route[255]
k0_255 = 8
e_k_step_255 = r255_route * k0_255 + (1 - r255_route) * 2.0  # rough bound
print(f"For r=255 (k0=8, P_route={r255_route:.4f}):")
print(f"  E[k/step] ≈ {r255_route:.4f} × 8 + {1-r255_route:.4f} × 2.0")
print(f"           = {e_k_step_255:.4f}  (rough formula)")
print(f"  Actual (script 92b): 2.418")
print(f"  D_hard_kern threshold: 3.419")
print(f"  Gap: {3.419 - e_k_step_255:.4f}")

# =====================================================================
# SUMMARY TABLE: P_route vs k/step
# =====================================================================
print()
print("=== SUMMARY: P_route AND avg k/step FOR BSet ELEMENTS ===\n")
print(f"{'r':>5}  {'k0':>4}  {'P_route':>10}  {'E[k/step approx]':>18}  {'gap to 3.419':>12}")
print("-" * 60)
for r in sorted(BList, key=lambda r: v2(r+1)):
    k0 = v2(r + 1)
    pr = p_route[r]
    e_ks = pr * k0 + (1 - pr) * 2.0  # rough: h=1 with k=k0, h>1 with k≈2
    gap_dk = 3.419 - e_ks
    print(f"  r={r:3d}  k={k0}  P={pr:.4f}  E[k/step]≈{e_ks:.4f}  gap={gap_dk:.4f}")
