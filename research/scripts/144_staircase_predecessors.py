"""
144_staircase_predecessors.py
==============================
Map the predecessor tree of the phantom staircase.

The canonical terminal staircase (Obs 287):
47->121->91->103->175->445->167->283->319->911->577->433->325->61->23->5->1

For each element p in the staircase, find ALL predecessors q (values where
macro_step(q) = p) up to some bound. This reveals:
1. Which values "feed into" each staircase node (the funnel tributaries)
2. Why some N=9 phantom elements (95, 253, 399) are rarely visited
3. The exact branching structure of the staircase attractor

Also: verify that side branches of the N=9 phantom cycle (253->95->91)
explain the entry-level discrepancy.
"""
import math
from collections import defaultdict, Counter

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

def predecessors_of(target, max_n=100000):
    """Find all odd n <= max_n such that macro_step(n) = target."""
    result = []
    # We need m*3^K - 1 = 2^l0 * target.
    # For each K and l0, compute m = (2^l0 * target + 1) / 3^K.
    for K in range(1, 20):
        for l0 in range(1, 40):
            numerator = (1 << l0) * target + 1
            denom = 3**K
            if numerator % denom != 0:
                continue
            m = numerator // denom
            if m <= 0 or m % 2 == 0:  # m must be odd positive
                continue
            n = m * (1 << K) - 1
            if n > max_n or n <= 0:
                continue
            # Verify
            n_out, K_act, l0_act = macro_step(n)
            if n_out == target and K_act == K and l0_act == l0:
                result.append((n, K, l0))
    return sorted(result)

# Canonical staircase (from the dominant terminal path)
STAIRCASE = [47, 121, 91, 103, 175, 445, 167, 283, 319, 911, 577, 433, 325, 61, 23, 5, 1]
# Note: these are in FORWARD order (47 is T-16, 1 is T-0)

STAIRCASE_SET = set(STAIRCASE)

# All phantom elements
PHANTOMS = {
    7:  {47, 91, 103, 121},
    8:  {71, 91, 103, 121, 175, 189},
    9:  {91, 95, 103, 167, 175, 253, 283, 319, 399, 445},
    10: {703, 937},
}
ALL_PHANTOM = set().union(*PHANTOMS.values())

print("=" * 70)
print("PART 1: PREDECESSORS OF EACH STAIRCASE ELEMENT")
print("=" * 70)
print()
print("For each element p in the canonical staircase, find predecessors q")
print("(values where macro_step(q) = p), categorized as:")
print("  [canon] = canonical staircase predecessor (from above)")
print("  [phantom] = other phantom element")
print("  [external] = neither staircase nor phantom")
print()

# The canonical predecessor of each staircase element
CANONICAL_PRED = {}
for i in range(len(STAIRCASE)-1):
    CANONICAL_PRED[STAIRCASE[i+1]] = STAIRCASE[i]

# Track: for T-8 to T-16 elements, what are the tributaries?
for target in STAIRCASE[:-1]:  # all except n=1
    preds = predecessors_of(target, max_n=200000)
    canon_preds = [q for q, K, l0 in preds if q in STAIRCASE_SET]
    phantom_preds = [q for q, K, l0 in preds if q in ALL_PHANTOM and q not in STAIRCASE_SET]
    external_preds = [q for q, K, l0 in preds if q not in ALL_PHANTOM and q not in STAIRCASE_SET]

    t_k = STAIRCASE.index(target)
    t_label = f"T-{len(STAIRCASE)-1-t_k}"

    print(f"p={target} ({t_label}): {len(preds)} predecessors <= 200000")
    if canon_preds:
        print(f"  [canon]:    {canon_preds}")
    if phantom_preds:
        print(f"  [phantom]:  {phantom_preds}")
    if external_preds:
        print(f"  [external]: {external_preds[:10]}{'...' if len(external_preds)>10 else ''} ({len(external_preds)} total)")
    print()

print()
print("=" * 70)
print("PART 2: SIDE BRANCHES OF THE N=9 PHANTOM CYCLE")
print("=" * 70)
print()
print("The N=9 phantom has 3 elements NOT in the canonical path: {95, 253, 399}")
print("Tracing their real orbits to see where they merge with the canonical path.")
print()

side_elements = [95, 253, 399]
for p in side_elements:
    orbit = [p]
    n = p
    while n > 1 and len(orbit) < 50:
        n, _, _ = macro_step(n)
        orbit.append(n)
        if n in STAIRCASE_SET:
            break
    merge_val = orbit[-1] if orbit[-1] in STAIRCASE_SET else None
    merge_step = len(orbit) - 1 if merge_val else None

    preds = predecessors_of(p, max_n=200000)
    print(f"n={p} (N=9 phantom, NOT in canonical path):")
    print(f"  Real orbit: {orbit}")
    if merge_val:
        t_k_merge = STAIRCASE.index(merge_val)
        t_label = f"T-{len(STAIRCASE)-1-t_k_merge}"
        print(f"  Merges into canonical path at n={merge_val} ({t_label}) after {merge_step} steps")
    print(f"  Predecessors (<= 200000): {preds[:5]}{'...' if len(preds)>5 else ''}")
    print()

print()
print("=" * 70)
print("PART 3: SIDE BRANCHES OF THE N=7/8 PHANTOM CYCLES")
print("=" * 70)
print()

side_n78 = [71, 189]  # N=8 elements not in canonical path
for p in side_n78:
    orbit = [p]
    n = p
    while n > 1 and len(orbit) < 50:
        n, _, _ = macro_step(n)
        orbit.append(n)
        if n in STAIRCASE_SET:
            break
    merge_val = orbit[-1] if orbit[-1] in STAIRCASE_SET else None
    merge_step = len(orbit) - 1 if merge_val else None

    preds = predecessors_of(p, max_n=200000)
    print(f"n={p} (N=8 phantom, NOT in canonical path):")
    print(f"  Real orbit: {orbit}")
    if merge_val:
        t_k_merge = STAIRCASE.index(merge_val)
        t_label = f"T-{len(STAIRCASE)-1-t_k_merge}"
        print(f"  Merges into canonical path at n={merge_val} ({t_label}) after {merge_step} steps")
    print(f"  Predecessors (<= 200000): {preds[:5]}{'...' if len(preds)>5 else ''}")
    print()

print()
print("=" * 70)
print("PART 4: WHY n=399 HAS ZERO PASSAGE RATE")
print("=" * 70)
print()
print("n=399 is in the N=9 phantom but has 0% passage rate.")
print("Finding ALL predecessors with no bound:")
print()

# Find predecessors of 399 with very high K
for K in range(1, 30):
    for l0 in range(1, 50):
        numerator = (1 << l0) * 399 + 1
        denom = 3**K
        if numerator % denom != 0:
            continue
        m = numerator // denom
        if m <= 0 or m % 2 == 0:
            continue
        n = m * (1 << K) - 1
        # Verify
        n_out, K_act, l0_act = macro_step(n)
        if n_out == 399:
            print(f"  Predecessor: n={n} (K={K}, l0={l0}, m={m})")
            if n.bit_length() <= 60:
                # Trace forward a few steps
                path = [n]
                x = n
                for _ in range(5):
                    x, _, _ = macro_step(x)
                    path.append(x)
                    if x in {1, 91, 103, 175, 319}:
                        break
                print(f"    orbit: {path}")
            break
    else:
        continue
    break

# Actually count predecessors with no bound
pred_sizes = {}
for target in ALL_PHANTOM:
    preds = predecessors_of(target, max_n=10**15)
    pred_sizes[target] = len(preds)

print()
print("Number of predecessors (<= 10^15) for each phantom element:")
print(f"{'n':>8} {'in_canon':>10} {'pred_count':>12} {'phantom_N':>12}")
print("-" * 50)
for p in sorted(ALL_PHANTOM):
    in_can = "yes" if p in STAIRCASE_SET else "no"
    Ns = [N for N, elems in PHANTOMS.items() if p in elems]
    lvl = ','.join(str(N) for N in sorted(Ns))
    print(f"{p:>8} {in_can:>10} {pred_sizes[p]:>12} {lvl:>12}")

print()
print("=" * 70)
print("PART 5: COUNT OF SUCCESSORS (FAN-OUT) OF STAIRCASE ELEMENTS")
print("=" * 70)
print()
print("How many predecessors does each staircase element have at LARGE n?")
print("We compute the EXPECTED NUMBER of predecessors of p in range [2^{b-1}, 2^b].")
print()
print("Formula: n = m * 2^K - 1, m * 3^K = 2^l0 * p + 1.")
print("For fixed K, l0: exactly one value of m if (2^l0*p+1) divisible by 3^K.")
print("Count of valid (K, l0) pairs with m in [2^{b-K-1}, 2^{b-K}) and m odd.")
print()

for p in [47, 121, 91, 103, 175, 445, 167, 283, 319]:
    # Count valid (K, l0) pairs giving n in b-bit range
    count_preds_b = {}
    for b in [20, 50, 100, 200]:
        count = 0
        for K in range(1, b):
            for l0 in range(1, b):
                numerator = (1 << l0) * p + 1
                denom = 3**K
                if numerator % denom != 0:
                    continue
                m = numerator // denom
                if m % 2 == 0:
                    continue  # m must be odd
                n = m * (1 << K) - 1
                if (1 << (b-1)) <= n < (1 << b):
                    count += 1
        count_preds_b[b] = count
    print(f"p={p}: preds in b-bit range: " +
          ", ".join(f"b={b}:{count_preds_b[b]}" for b in [20, 50, 100, 200]))

print()
print("Expected: count of predecessors grows with b (more b-bit numbers feed into each staircase element).")
print("This explains why passage rates are high even for large starting numbers.")

print()
print("=" * 70)
print("PART 6: STAIRCASE SUMMARY - FULL STRUCTURE")
print("=" * 70)
print()
print("The phantom staircase as an attractor tree:")
print()
print("  Incoming tributaries (side branches merging into canonical path)")
print()

# Build the full tree: for each staircase node, list its predecessors
max_search = 10000
for i, target in enumerate(STAIRCASE[:9]):  # T-16 to T-8
    t_label = f"T-{len(STAIRCASE)-1-i}"
    preds = predecessors_of(target, max_n=max_search)
    canon_pred = CANONICAL_PRED.get(target)
    non_canon = [(q, K, l0) for q, K, l0 in preds if q != canon_pred]
    phantom_non_canon = [(q, K, l0) for q, K, l0 in non_canon if q in ALL_PHANTOM]
    external_non_canon = [(q, K, l0) for q, K, l0 in non_canon if q not in ALL_PHANTOM]
    print(f"{t_label}: n={target}")
    if canon_pred:
        print(f"  <- [canon] {canon_pred}")
    for q, K, l0 in phantom_non_canon:
        print(f"  <- [phantom] {q} (K={K},l0={l0})")
    for q, K, l0 in external_non_canon[:5]:
        print(f"  <- [external] {q} (K={K},l0={l0})")
    if len(external_non_canon) > 5:
        print(f"  <- [external] ...and {len(external_non_canon)-5} more")
    print()
