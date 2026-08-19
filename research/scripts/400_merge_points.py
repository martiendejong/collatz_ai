# 400: MERGE POINTS — what happens when two orbits merge, and where?
# Exact local rule to verify: n and 4n+1 ALWAYS merge in one Syracuse step
# (3(4n+1)+1 = 4(3n+1): same odd part). Global question: distribution of the
# first common value of two random orbits. PREREGISTERED prediction from the
# harmonic law (P(orbit passes w) ~ C/w): P(merge value > X) ~ 1/X.
import random
random.seed(400)

def syr(n):
    m = 3*n + 1
    return m >> ((m & -m).bit_length() - 1)

# (A) exact rule
ok = 0
for _ in range(2000):
    n = random.getrandbits(64) | 1
    a = syr(n); b = syr(4*n + 1)
    ok += (a == b)
print(f"(A) regel n ~ 4n+1: {ok}/2000 mergen in EEN stap")

# (B) random pairs: first common odd value
def orbit_set(n, cap=100000):
    s = set()
    while n != 1 and len(s) < cap:
        s.add(n)
        n = syr(n)
    s.add(1)
    return s

merge_vals = []
merge_rel = []
for _ in range(400):
    a = random.getrandbits(80) | 1
    b = random.getrandbits(80) | 1
    A = orbit_set(a)
    n = b
    while n not in A:
        n = syr(n)
    merge_vals.append(n)
merge_vals.sort()
import math
med = merge_vals[len(merge_vals)//2]
print(f"\n(B) 400 willekeurige 80-bit-paren:")
print(f"    mediaan merge-waarde: {med}")
qs = [50, 75, 90, 95, 99]
for q in qs:
    print(f"    p{q}: {merge_vals[int(len(merge_vals)*q/100)-1]}")
from collections import Counter
top = Counter(merge_vals).most_common(8)
print(f"    meest voorkomende merge-punten: {top}")
# tail test: P(merge > X) ~ 1/X ?
for X in [10, 100, 1000, 10000, 100000]:
    frac = sum(1 for v in merge_vals if v > X)/len(merge_vals)
    print(f"    P(merge > {X:>6}) = {frac:.4f}   (1/X-wet ratio: {frac*X:.1f})")

# (C) 2-adisch nabije paren: mergen ze hoger?
print(f"\n(C) nabije paren (n, n + 2^k * oneven):")
for k in [8, 16, 24, 32]:
    vals = []
    for _ in range(150):
        n = random.getrandbits(80) | 1
        m = n + (1 << k)*(2*random.getrandbits(8) + 1)
        A = orbit_set(n)
        x = m
        while x not in A:
            x = syr(x)
        vals.append(x)
    vals.sort()
    print(f"    k={k:2d}: mediaan merge-waarde {vals[len(vals)//2]}, p90 {vals[int(0.9*len(vals))]}")
