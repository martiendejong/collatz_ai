"""
169_desert_theorem.py
======================
THE PREDECESSOR DESERT THEOREM (candidate) -- from LP abstraction to real
integers.

Exact 3-adic algebra (one line each):
  shrink-map s(a) = (2a-1)/3 (the predecessor smaller than a, exists iff
  a == 2 mod 3):
    a == -1 (mod 3^k)  ==>  s(a) == -1 (mod 3^(k-1))   [tower is shrink-closed]
    a == -4 (mod 3^k)  ==>  s(a) == -3 (mod 3^(k-1)), i.e. == 0 mod 3 [DEAD]

Prediction (from the certificate funnel law, Obs 318): the number of odd
predecessors of a below a fixed budget 2^y * a is suppressed for the -4
tower relative to the -1 tower by a factor growing like (3/2)^k.

Test on ACTUAL integers: for k = 2..7, take odd representatives
a == -1 (mod 3^k) and a == -4 (mod 3^k) of comparable size, BFS the inverse
V-tree (predecessors n = (2^e a - 1)/3, e-parity forced by a mod 3), count
all odd predecessors <= 2^y * a, and compare tower ratios.
"""
import numpy as np
from math import log2

def count_predecessors(a0, y, cap_nodes=6_000_000):
    """# odd n in the inverse V-tree of a0 with n <= 2^y * a0."""
    budget = a0 << y
    level = [a0]
    total = 0
    while level and total < cap_nodes:
        nxt = []
        for m in level:
            r = m % 3
            if r == 0:
                continue
            e = 2 if r == 1 else 1
            while True:
                n = (m * (1 << e) - 1) // 3
                if n > budget:
                    break
                nxt.append(n)
                e += 2
        total += len(nxt)
        level = nxt
    return total

def reps(target, k, size_min, count=4):
    """Odd representatives == target mod 3^k, the first `count` above size_min."""
    M = 3**k
    t = target % M
    out = []
    a = t + ((size_min - t) // M + 1) * M
    while len(out) < count:
        if a % 2 == 1:
            out.append(a)
        a += M
    return out

Y = 22
SIZE = 3**8          # representatives ~ 6561+, comparable across k
print(f"budget: predecessors <= 2^{Y} * a;  representatives ~ {SIZE}+")
print(f"{'k':>3} {'N(-1 tower)':>12} {'N(-4 tower)':>12} {'ratio':>8} "
      f"{'ratio growth':>13} {'(3/2 pred.)':>11}")
prev = None
for k in range(2, 8):
    n1 = np.mean([count_predecessors(a, Y) for a in reps(-1, k, SIZE)])
    n4 = np.mean([count_predecessors(a, Y) for a in reps(-4, k, SIZE)])
    ratio = n1 / n4
    growth = ratio / prev if prev else float('nan')
    print(f"{k:>3} {n1:>12.0f} {n4:>12.0f} {ratio:>8.2f} "
          f"{growth:>13.3f} {'1.500':>11}")
    prev = ratio
print()
print("Algebra check (exact):")
for k in [5, 8, 12]:
    M = 3**k
    a1 = (-1) % M
    a4 = (-4) % M
    s1 = (2*a1 - 1) // 3 if (2*a1 - 1) % 3 == 0 else None
    s4 = (2*a4 - 1) // 3 if (2*a4 - 1) % 3 == 0 else None
    print(f"  k={k}: s(-1) == -1 mod 3^{k-1}: {s1 % 3**(k-1) == (-1) % 3**(k-1)};"
          f"  s(-4) == 0 mod 3: {s4 % 3 == 0}")
