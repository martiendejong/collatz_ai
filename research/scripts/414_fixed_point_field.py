"""Obs 621 part 2: the fixed-point FIELD vs the measured visit field.

Local features explain R2 = 0.36 of log(v*P) (script 413). Test the
chain-accumulation hypothesis: compute the memoized fixed-point field
P_fp(v) (recursion P = sum of children, boundary 1/u beyond VCUT) for
every odd non-spring v < 3^8 and correlate log P_fp with measured
log P over the same 1837 values. High correlation = the correction
layer is fully chain-computable (deep but not mysterious); the gap to
1 = boundary noise + measurement noise + true depth.
"""
import math
import sys
import numpy as np
from collections import Counter
sys.setrecursionlimit(100000)

MMAX = 64
VCUT = 1_000_000
memo = {}

def P(v):
    if v in memo:
        return memo[v]
    s = 0.0
    for m in range(1, MMAX + 1):
        t = (1 << m) * v - 1
        if t % 3 == 0:
            u = t // 3
            if u % 2 == 1 and u % 3 != 0 and u > 1:
                s += P(u) if u < VCUT else 1.0 / u
    memo[v] = s
    return s

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

NSTARTS = 300_000
hits = Counter()
for n0 in range(1_000_001, 1_000_001 + 2 * NSTARTS, 2):
    n = n0
    while n != 1:
        if n < 6561:
            hits[n] += 1
        m = 3 * n + 1
        n = m >> v2(m)

vs = sorted(v for v in hits if v % 3 != 0 and hits[v] >= 150 and v >= 5)
meas = np.array([math.log(hits[v] / NSTARTS) for v in vs])
fp = np.array([math.log(P(v)) for v in vs])
r = np.corrcoef(meas, fp)[0, 1]
a, b = np.polyfit(fp, meas, 1)
resid = meas - (a * fp + b)
print(f"{len(vs)} values; corr(log P_meas, log P_fp) = {r:.4f}")
print(f"slope = {a:.3f} (1 = exact proportionality); residual sd = {resid.std():.3f}")
print(f"(compare: local-feature R2 was 0.358; this R2 = {r * r:.3f})")
print(f"largest residuals (v, meas/fp ratio):")
for i in np.argsort(np.abs(resid))[-5:][::-1]:
    print(f"  v={vs[i]:>5}: measured {math.exp(meas[i]):.5f}  fp {math.exp(fp[i]):.5f}")
