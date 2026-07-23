"""
178_d7_dimension.py
===================
The channel-7 highway D7 (Obs 340): 2-adic points whose macro-chain stays in
channel 7 (l=1, K'>=2) with K >= 2 forever -- perpetual climbers.

 M1: dimension of D7 by survival counting: among odd residues mod 2^W, count
     those whose first j macro-steps are all channel-7 with K>=2 (each step
     determined by the residue while consumed bits S_j = sum(K_i + 1) <= W).
     dim ~ log2(#survivors * classes) / ... measured as (W - cost_j)/W with
     cost_j = W - log2 N_j at fixed consumed depth.
 M2: the D7-shadow integers: actual positive integers n < 2^W with the longest
     channel-7/K>=2 prefixes. Prediction (Obs 337/340): these are the seeds of
     record excursions.
 M3: analytic check: per macro-step Haar cost of (c=7 AND K>=2):
     P = P(K>=2)*P(c=7) = 1/8, consumed bits E[K+1 | K>=2] = 4
     -> naive dim(D7) = 1 - 3/4 = ... measure the true value.
"""
from math import log2
from collections import Counter

def v2(x):
    return (x & -x).bit_length() - 1

W = 22
LIMIT = 1 << W

def d7_prefix(n):
    """Length of the channel-7/K>=2 macro-prefix of n, and bits consumed,
    stopping when the next step is no longer determined by n mod 2^W
    (consumed + K + 1 > W) or the chain leaves the highway."""
    steps, consumed = 0, 0
    while True:
        K = v2(n + 1)
        a = (n + 1) >> K
        if consumed + K + 1 > W:
            return steps, consumed, True    # censored (ran out of window)
        x = 3**K * a
        c = x & 7
        if not (K >= 2 and c == 7):
            return steps, consumed, False
        n = (x - 1) >> 1                     # l=1 on the highway
        steps += 1
        consumed += K + 1

surv = Counter()          # j -> number of odd n with prefix >= j (uncensored counting)
best = []                 # (steps, n) records
maxsteps = 0
for n in range(3, LIMIT, 2):
    s, cons, censored = d7_prefix(n)
    for j in range(1, s + 1):
        surv[j] += 1
    if s > maxsteps and not censored:
        maxsteps = s
        best.append((s, n))
    elif s == maxsteps and len(best) < 40 and s > 0:
        best.append((s, n))

total = LIMIT // 2
print(f"W = {W}: {total:,} odd residues")
print("M1: survival of channel-7/K>=2 prefixes:")
print("    j   N_j        frac        cost-bits   cost/step")
prev = total
for j in range(1, maxsteps + 1):
    N = surv[j]
    frac = N / total
    cost = -log2(frac)
    print(f"    {j:<3} {N:<10,} {frac:.6f}    {cost:7.3f}     {cost/j:.3f}")
print(f"    (Haar prediction M3: P(step) = 1/8 -> 3.000 bits/step;")
print(f"     naive dim(D7) = 1 - 3/E[K+1|K>=2] = 1 - 3/4 = 0.25)")

print()
print("M2: D7-shadow integers (longest highway prefixes, n < 2^22):")
top = sorted(best, reverse=True)[:15]
for s, n in top:
    # show their macro chain start and their eventual excursion
    m, peak = n, n
    while m != 1:
        q = 3 * m + 1
        e = v2(q)
        m = q >> e
        peak = max(peak, m)
    print(f"    n={n:<9} prefix={s}  odd-peak={peak:<12}"
          f" ratio=2^{log2(peak/n):.1f}")
