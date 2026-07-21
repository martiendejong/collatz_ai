"""
132_collatz_tree.py
====================
Structure of the Collatz tree: predecessor counts and tree depth.

Inverse of the macro-step:
  If macro_step(n) = n', then n = 2^K * m - 1 where:
    m * 3^K = 2^{l0} * n' + 1
    K = v2(n + 1), l0 = v2(m * 3^K - 1)

For each n', the predecessors are parameterized by (K, l0) pairs where:
  1. 3^K | (2^{l0} * n' + 1)
  2. m = (2^{l0} * n' + 1) / 3^K is a positive odd integer
  3. K = v2(m + 1) matches... wait, K = v2(2^K * m) - ... let me reconsider.

Actually: n = 2^K * m - 1 means n+1 = 2^K * m, so v2(n+1) = K (since m is odd). ✓

Key question: how many predecessors does n' have? How does this vary?
"""
import math
from collections import defaultdict
import random as _r
import numpy as np

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

def predecessors(n_prime, max_K=20, max_l0=30):
    """Find all predecessors n s.t. macro_step(n) = n'."""
    preds = []
    for l0 in range(1, max_l0+1):
        target = (n_prime << l0) + 1  # = 2^l0 * n' + 1
        # Find all K s.t. 3^K | target
        # Also need: m = target / 3^K is odd
        for K in range(1, max_K+1):
            pk = 3**K
            if target % pk != 0: continue
            m = target // pk
            if m <= 0 or m % 2 == 0: continue  # m must be positive odd
            n_pred = (m << K) - 1  # = 2^K * m - 1
            # Verify:
            v_check = v2(n_pred + 1)
            if v_check != K: continue  # K must match v2(n_pred+1)
            # Double-check macro_step
            n_out, K_out, l0_out = macro_step(n_pred)
            if n_out == n_prime:
                preds.append((n_pred, K, l0))
    return preds

print("=" * 70)
print("PART 1: PREDECESSORS OF SMALL ODD NUMBERS")
print("=" * 70)
print()

for n_prime in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27]:
    preds = predecessors(n_prime, max_K=15, max_l0=20)
    if preds:
        print(f"n'={n_prime:>4}: {len(preds)} predecessors")
        for n_p, K, l0 in sorted(preds)[:8]:
            print(f"           n={n_p:>8} (K={K}, l0={l0})")
        if len(preds) > 8: print(f"           ... ({len(preds)-8} more)")
    else:
        print(f"n'={n_prime:>4}: 0 predecessors (NO entry from outside the tree!)")

print()
print("=" * 70)
print("PART 2: PREDECESSOR COUNT DISTRIBUTION")
print("=" * 70)
print()

# Count predecessors for odd numbers up to some bound
pred_counts = defaultdict(int)
pred_count_vals = {}
for n_prime in range(1, 500, 2):  # odd numbers 1..499
    preds = predecessors(n_prime, max_K=12, max_l0=15)
    cnt = len(preds)
    pred_counts[cnt] += 1
    pred_count_vals[n_prime] = cnt

print(f"Predecessor count distribution (n'=1..499 odd, max_K=12, max_l0=15):")
print(f"{'Count':>8} {'#(n) with this count':>22} {'Fraction':>10}")
for cnt in sorted(pred_counts):
    frac = pred_counts[cnt] / 250
    print(f"{cnt:>8} {pred_counts[cnt]:>22} {frac:>10.4f}")

# Find the most-connected numbers
top_n = sorted(pred_count_vals.items(), key=lambda x: -x[1])[:10]
print(f"\nTop 10 most-connected odd numbers (by predecessor count):")
for n_prime, cnt in top_n:
    print(f"  n'={n_prime:>5}: {cnt} predecessors")

print()
print("=" * 70)
print("PART 3: TREE STRUCTURE — HOW MANY LEVELS UP FROM 1?")
print("=" * 70)
print()
print("BFS from n=1: how many odd numbers are at each depth?")
print("Depth 0: {1}")
print("Depth 1: predecessors of 1 (excluding 1)")
print("Depth 2: predecessors of depth-1 nodes")
print()

from collections import deque
visited = {1}
current_level = [1]
for depth in range(1, 7):
    next_level = []
    for n_prime in current_level:
        preds = predecessors(n_prime, max_K=15, max_l0=20)
        for n_p, K, l0 in preds:
            if n_p not in visited:
                visited.add(n_p)
                next_level.append(n_p)
    print(f"Depth {depth}: {len(next_level)} new nodes | examples: {sorted(next_level)[:10]}")
    current_level = next_level
    if len(next_level) == 0: break

print()
print("=" * 70)
print("PART 4: ARITHMETIC STRUCTURE OF PREDECESSORS")
print("=" * 70)
print()
print("For each (K, l0), which n' have a valid predecessor?")
print("Condition: 3^K | (2^{l0} * n' + 1)")
print("=> 2^{l0} * n' = -1 mod 3^K")
print("=> n' = -2^{-l0} mod 3^K")
print("=> n' belongs to a specific residue class mod 3^K")
print()

# For small K, what is the residue class mod 3^K?
print(f"{'K':>4} {'l0':>4} {'Residue n mod 3^K':>22} {'3^K':>8}")
for K in range(1, 6):
    pk = 3**K
    for l0 in range(1, min(K+3, 8)):
        # n' ≡ -pow(2, l0, pk)^{-1} * 1 mod pk? Wait:
        # 2^{l0} * n' + 1 ≡ 0 mod pk
        # n' ≡ -1 * modinv(2^l0, pk) mod pk
        from sympy import mod_inverse
        try:
            inv = mod_inverse(pow(2, l0, pk), pk)
            res = (-inv) % pk
            print(f"{K:>4} {l0:>4} {res:>22d} (mod {pk})")
        except:
            print(f"K={K}, l0={l0}: no inverse")

print()
print("=" * 70)
print("PART 5: DENSITY OF NUMBERS WITH NO PREDECESSORS")
print("=" * 70)
print()
print("Is every odd number reachable from the Collatz tree root?")
print("(If yes, this would be equivalent to Collatz conjecture)")
print()
# Numbers with 0 predecessors: n' such that there's no n with macro_step(n)=n'
# (within our search window)
no_preds = [n for n, cnt in pred_count_vals.items() if cnt == 0]
print(f"Odd numbers with 0 predecessors (n'=1..499, max_K=12, max_l0=15):")
print(f"  Count: {len(no_preds)}")
print(f"  Examples: {sorted(no_preds)[:30]}")
print()
print("Note: These may have predecessors with larger K or l0 (beyond our search).")
print("For n'=3: 3*3^K must divide 2^l0*3+1 = 2^l0*3+1. For K=1: 2^l0*3+1 must")
print("  be divisible by 3. But 2^l0*3+1 = 1 mod 3 always (never 0). So n'=3 has no predecessor!")
print("  Wait: let's check n'=3 predecessors more carefully...")
p3 = predecessors(3, max_K=20, max_l0=30)
print(f"  n'=3: {len(p3)} predecessors = {[(n,K,l0) for n,K,l0 in p3[:5]]}")
print()
print("Checking: 2^l0 * 3 + 1 mod 3 = 0+1 = 1 ≠ 0 for all l0. So K=1 fails.")
print("2^l0 * 3 + 1 mod 9 = 2^l0*3+1. For l0=1: 7. For l0=2: 13=4. For l0=3: 25=7.")
print("For l0=4: 49=4. Pattern: 7,4,7,4,... none is 0 mod 9. So K=2 also fails!")
print("In fact: 2^l0*3+1 mod 3^K. Since 3|(2^l0*3+1) iff 3|1 which is false.")
print("So 3 can NEVER be the target of a macro-step! n'=3 is a 'dead end' (root-like).")
