"""
182_tail_typicality.py
======================
RUNG t14: TAIL TYPICALITY -- do the zero tails behave Haar-typically at the
record extremes?

Key numerology: a record excursion E at seed scale K needs t ~ E/0.2516 odd
steps, consuming S ~ (4/3)t ~ 5.3E halving-bits >> K: records live far BEYOND
the free bit-window of the seed. If the deterministic zero-tail continuation
follows Haar statistics (tail typicality), the free-window LD theory still
applies and predicts E*(K) = K - c*log K (Gumbel log correction).
If seeds could encode luck, E*(K)/K would exceed 1; if tails were adversarial,
records would fall well short.

 R1: record scan n < 2^23, fit E*(K) = a*K - b*log2(K): prediction a = 1.
 R2: climb e-histogram of top records pooled vs the Gibbs-tilted law
     Geom(3/4): P(e=j) = (3/4)(1/4)^{j-1} -- the conditional-limit theorem
     made testable (rigorous under Haar via the exact-i.i.d. camouflage
     theorem, Obs 341; here tested on actual integer records = tail regime).
"""
from math import log2
from collections import Counter

def v2(x):
    return (x & -x).bit_length() - 1

def T(n):
    q = 3 * n + 1
    e = v2(q)
    return q >> e, e

LIM = 1 << 23
best = 0.0
records = []
for n in range(3, LIM, 2):
    m, peak = n, n
    while m >= n and m != 1:
        m, _ = T(m)
        if m > peak:
            peak = m
    if log2(peak / n) > best - 2:
        m2, p2 = n, n
        while m2 != 1:
            m2, _ = T(m2)
            p2 = max(p2, m2)
        peak = p2
    E = log2(peak / n)
    if E > best:
        best = E
        records.append((n, peak, E))

print("R1: excursion records n < 2^23")
print("    n          K=log2 n   E          E/K       E/(K - 1.15*log2 K)")
rows = []
for n, peak, E in records:
    K = log2(n)
    if K < 4:
        continue
    rows.append((n, K, E))
    print(f"    {n:<10} {K:8.2f}  {E:8.2f}   {E/K:7.3f}   "
          f"{E/(K - 1.15*log2(K)):7.3f}")
# least squares fit E = a*K + b*log2(K) on the record points (K >= 8)
pts = [(K, E) for _, K, E in rows if K >= 8]
import statistics
# fit via normal equations for [K, logK, 1]
X = [(K, log2(K), 1.0) for K, _ in pts]
Y = [E for _, E in pts]
# solve 3x3
def solve3(X, Y):
    A = [[sum(x[i] * x[j] for x in X) for j in range(3)] for i in range(3)]
    b = [sum(x[i] * y for x, y in zip(X, Y)) for i in range(3)]
    # gaussian elimination
    for i in range(3):
        p = A[i][i]
        for j in range(i + 1, 3):
            f = A[j][i] / p
            for k2 in range(3):
                A[j][k2] -= f * A[i][k2]
            b[j] -= f * b[i]
    x = [0.0] * 3
    for i in (2, 1, 0):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, 3))) / A[i][i]
    return x
a, bb, c = solve3(X, Y)
print(f"    fit: E*(K) = {a:.3f}*K + {bb:.3f}*log2(K) + {c:.2f}   "
      f"(tail-typicality prediction: slope a = 1)")

print()
print("R2: pooled climb e-histogram of top records vs Gibbs law Geom(3/4):")
ehist = Counter()
tot = 0
for n, peak, E in records:
    if E < 8:
        continue
    m = n
    seq = [n]
    while m != 1:
        m, _ = T(m)
        seq.append(m)
    ip = seq.index(max(seq))
    m = n
    for i in range(ip):
        m, e = T(m)
        ehist[min(e, 6)] += 1
        tot += 1
print("    e    measured   Geom(3/4)")
for e in range(1, 7):
    pred = (3 / 4) * (1 / 4) ** (e - 1) if e < 6 else (1 / 4) ** 5
    print(f"    {e}    {ehist[e]/tot:8.4f}   {pred:8.4f}")
print(f"    ({tot} pooled climb steps from records with E >= 8; "
      f"mean e = {sum(e*c for e,c in ehist.items())/tot:.3f} vs 4/3 = 1.333)")
