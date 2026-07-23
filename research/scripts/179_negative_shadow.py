"""
179_negative_shadow.py
======================
SHADOW DECOMPOSITION OF CLIMBS (follow-up to Obs 340/341).

Hypothesis: climb phases of positive orbits are 2-adic SHADOWS of the negative
cycles of T (the inhabitants of the divergence fractal D, Obs 337):
    C1 = {-1}                                   mean e = 1     (burst climber)
    C2 = {-5, -7}                               mean e = 3/2   (the K=2 highway
                                                fixed point, Obs 340)
    C3 = {-17,-25,-37,-55,-41,-61,-91}          mean e = 11/7  (slow climber)
Shadowing n == c (mod 2^d) is preserved by T with depth loss exactly e per
step (|T'|_2 = 2^e), so a depth-d shadow predicts ~d/e_cycle steps of the
cycle's climb rate log2(3) - B/A > 0.

 S1: shadow-walk of 27 and 703: per odd iterate the deepest negative-cycle
     point and depth; compare with climb/fall phases.
 S2: statistics over record orbits vs random controls: mean max-depth during
     climb vs after peak; baseline E[max depth over 10 points] ~ log2(10)+1.
 S3: predictive test: bucket all steps by shadow depth d; report mean e on the
     NEXT step and mean growth over the next 3 odd steps. Prediction: deep
     shadow => next-step statistics of the shadowed cycle (e.g. e=1 for -1,
     alternation 1,2 for C2), i.e. growth increasing with depth.
"""
import random
from math import log2

def v2(x):
    return (x & -x).bit_length() - 1

def T(n):
    q = 3 * n + 1
    e = v2(q)
    return q >> e, e

NEG = [-1, -5, -7, -17, -25, -37, -55, -41, -61, -91]
LABEL = {-1: "C1:-1", -5: "C2:-5", -7: "C2:-7"}
for c in NEG[3:]:
    LABEL[c] = f"C3:{c}"

def shadow(n):
    best_d, best_c = -1, None
    for c in NEG:
        d = v2(n - c)
        if d > best_d:
            best_d, best_c = d, c
    return best_d, best_c

def orbit(n):
    seq = [n]
    while n != 1:
        n, e = T(n)
        seq.append(n)
    return seq

# ---------------- S1: shadow walks ------------------------------------------
for n0 in (27, 703):
    seq = orbit(n0)
    peak = max(seq)
    ipeak = seq.index(peak)
    print(f"S1: shadow-walk n0={n0} (odd-peak {peak} at step {ipeak}):")
    for i, n in enumerate(seq[:-1]):
        d, c = shadow(n)
        nxt, e = T(n)
        phase = "CLIMB" if i < ipeak else "fall"
        mark = " <== deep" if d >= 5 else ""
        if i <= ipeak + 3 or d >= 5:
            print(f"    t={i:<3} n={n:<10} e={e} shadow={LABEL[c]:<7} depth={d:<2}"
                  f" {phase}{mark}")
    print()

# ---------------- S2: depth statistics --------------------------------------
def phase_depths(n0):
    seq = orbit(n0)
    peak = max(seq)
    ipeak = seq.index(peak)
    climb = [shadow(n)[0] for n in seq[:ipeak]]
    fall = [shadow(n)[0] for n in seq[ipeak:-1]]
    return climb, fall

print("S2: mean shadow depth, climb vs fall:")
print("    n0        E[d|climb]  E[d|fall]   max(d,climb)")
for n0 in (27, 703, 26623, 626331, 837799):
    climb, fall = phase_depths(n0)
    print(f"    {n0:<9} {sum(climb)/len(climb):8.2f}  {sum(fall)/len(fall):8.2f}"
          f"     {max(climb)}")
random.seed(11)
alld = []
for _ in range(500):
    n = random.getrandbits(20) | 1
    for m in orbit(n)[:-1]:
        alld.append(shadow(m)[0])
print(f"    random controls: E[d] = {sum(alld)/len(alld):.2f}  "
      f"(baseline max-of-10 ~ {log2(10)+1:.1f})")

# ---------------- S3: predictive test ----------------------------------------
print()
print("S3: shadow depth d  ->  next-step e and 3-step growth (all record orbits"
      " + 2000 random):")
from collections import defaultdict
buck_e = defaultdict(list)
buck_g = defaultdict(list)
pool = [27, 703, 26623, 626331, 837799] + \
       [random.getrandbits(24) | 1 for _ in range(2000)]
for n0 in pool:
    seq = orbit(n0)
    for i in range(len(seq) - 4):
        d, c = shadow(seq[i])
        _, e = T(seq[i])
        g = log2(seq[i + 3] / seq[i])
        db = min(d, 9)
        buck_e[db].append(e)
        buck_g[db].append(g)
print("    d    #steps     E[e_next]   E[growth 3 odd-steps] (bits)")
for d in sorted(buck_e):
    ne = buck_e[d]
    ng = buck_g[d]
    print(f"    {d:<4} {len(ne):<10,} {sum(ne)/len(ne):8.3f}   {sum(ng)/len(ng):10.3f}")
print("    (Haar: E[e]=2.0, growth=3*(log2 3 - 2) = -1.245 bits)")
