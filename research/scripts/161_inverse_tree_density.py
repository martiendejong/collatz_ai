"""
161_inverse_tree_density.py
============================
ATTACKING WALL 2 WITH THE FIELD'S BEST UNCONDITIONAL TOOL:
the inverse-tree density method (Applegate-Lagarias I, 1995; refined to
x^0.84 by Krasikov-Lagarias 2003 via difference inequalities).

Method (tree-search variant, faithful small-scale reimplementation):
  Work on odd integers with V(n) = (3n+1)/2^e. The inverse tree from 1:
  children of m are n = (m*2^e - 1)/3 for e >= 1 with m*2^e == 1 (mod 3)
  (e even iff m == 1 mod 3, e odd iff m == 2 mod 3); such n is automatically
  odd; nodes n == 0 (mod 3) are leaves. Distinct paths give distinct n
  (each n has a unique forward orbit), so
      f(x) := #{odd n <= x reaching 1}  >=  #{enumerated nodes <= x},
  UNCONDITIONALLY. Truncating e <= E and depth <= K gives a computable
  certificate; the exponent gamma(y) = log2 count(<=2^y) / y lower-bounds
  the density exponent achievable by this truncation.

Published frontier: 0.81 (tree-search, 1995) -> 0.84 (difference
inequalities, 2003). This script reproduces the method and measures how the
certificate scales with (E, K) -- mapping the compute frontier for pushing
beyond 0.84.

BONUS: prime-D scan along the convergent line (where our prime-D reduction
applies): for A <= 300, test D = 2^ceil(A theta) - 3^A for primality.
"""
from math import log2, ceil
import numpy as np

THETA = log2(3)

def enumerate_tree(E, max_nodes=25_000_000, max_depth=64):
    """BFS of the inverse tree from 1 with halving-exponent cap e <= E.
    Returns histogram of log2(value) over all enumerated nodes."""
    level = [5]  # children of 1 with e<=E, excluding the trivial self-loop:
    # e even needed (1 mod 3): e=2 -> 1 (self), e=4 -> 5. (e=6>E for E<6)
    extra = []
    for e in range(6, E+1, 2):
        extra.append(((1 << e) - 1)//3)
    level += extra
    hist_bits = []
    total = len(level)
    depth = 1
    while level and total < max_nodes and depth < max_depth:
        hist_bits.extend(log2(m) for m in level)
        nxt = []
        for m in level:
            r = m % 3
            if r == 0:
                continue                      # leaf
            e0 = 2 if r == 1 else 1           # parity of valid e
            for e in range(e0, E+1, 2):
                n = (m*(1 << e) - 1)//3
                nxt.append(n)
        total += len(nxt)
        level = nxt
        depth += 1
    hist_bits.extend(log2(m) for m in level)
    return np.array(hist_bits), depth, total

print("="*74)
print("PART 1: INVERSE-TREE DENSITY CERTIFICATES (Applegate-Lagarias method)")
print("="*74)
print(f"{'E':>3} {'depth':>6} {'nodes':>12} {'best gamma':>11} {'at y=log2(x)':>13}")
for E in [3, 4, 5, 6, 8]:
    bits, depth, total = enumerate_tree(E)
    bits.sort()
    ymax = int(bits.max())
    ys = np.arange(min(15, max(2, ymax - 1)), ymax + 1)
    counts = np.searchsorted(bits, ys, side='right')
    gammas = np.log2(np.maximum(counts, 1)) / ys
    k = int(np.argmax(gammas))
    print(f"{E:>3} {depth:>6} {total:>12,} {gammas[k]:>11.4f} {ys[k]:>13}")
print()
print("Published frontier: 0.81 (tree-search 1995), 0.84 (Krasikov-Lagarias 2003).")
print("The certificate gamma rises with E and depth: the scaling map above shows")
print("what raw enumeration buys; passing 0.84 requires the difference-inequality")
print("machinery at higher congruence levels (mod 3^k), not enumeration alone.")
print()

print("="*74)
print("PART 2: PRIME-D SCAN -- where the prime-D reduction applies")
print("="*74)
def is_probable_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    import random
    for _ in range(24):
        a = random.randrange(2, n-1)
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else:
            return False
    return True

primes = []
for A in range(2, 301):
    B = ceil(A*THETA)
    D = (1 << B) - 3**A
    if D <= 0:
        continue
    if is_probable_prime(D):
        primes.append((A, B, D.bit_length()))
print(f"signatures with D prime (A <= 300): {len(primes)} of 299")
for A, B, bits in primes:
    print(f"  (A,B)=({A},{B})  D has {bits} bits")
print()
print(f"empirical density ~ {len(primes)/299:.3f}; heuristic ~ c/ln(D) ~ 1.1/A;")
print("at each such signature the cycle problem is ONE Gauss-period value.")
