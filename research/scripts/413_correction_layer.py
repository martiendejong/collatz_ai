"""Obs 621: regression anatomy of the multifractal correction layer.

The bridge verdict (Obs 617): visit measure = mean-field (1/v) times a
pointwise correction field invisible to K-L. Obs 620 found its levers.
Here: measure P(hit v) for every odd non-spring v < 3^8 and regress
log(v * P) on the structural features now in hand:
  - direction: v mod 3 (2 -> descending m=1 rung, 1 -> ascending m=2)
  - first alive rung index m1 (spring pattern on low rungs)
  - springness of the second rung
  - ternary trailing-ones count (in-degree capacity, stairway law)
  - log2 v  (residual size trend)
  - mod-9 roulette class dummies
Report incremental R^2; the unexplained remainder is the true depth of
the correction layer at this scale.
"""
import math
import numpy as np
from collections import Counter

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def tern_trailing_ones(v):
    t = 0
    while v % 3 == 1:
        v //= 3
        t += 1
    return t

def first_rungs(v):
    """(m1, alive1, m2, alive2) for the two lowest valid rungs."""
    out = []
    m = 1
    while len(out) < 2:
        t = (1 << m) * v - 1
        if t % 3 == 0 and (t // 3) % 2 == 1:
            u = t // 3
            out.append((m, u % 3 != 0))
        m += 1
    return out[0][0], out[0][1], out[1][0], out[1][1]

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
y = np.array([math.log(v * hits[v] / NSTARTS) for v in vs])
print(f"{len(vs)} values; log(v*P): mean {y.mean():.3f}, sd {y.std():.3f}")

def r2(X):
    X = np.column_stack([np.ones(len(vs))] + X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return 1 - resid.var() / y.var()

f_dir = [np.array([1.0 if v % 3 == 2 else 0.0 for v in vs])]
fr = [first_rungs(v) for v in vs]
f_m1 = f_dir + [np.array([float(a[0]) for a in fr])]
f_al = f_m1 + [np.array([float(a[1]) for a in fr]),
               np.array([float(a[3]) for a in fr])]
f_t1 = f_al + [np.array([float(min(tern_trailing_ones(v), 5)) for v in vs])]
f_lv = f_t1 + [np.array([math.log2(v) for v in vs])]
mod9 = f_lv + [np.array([1.0 if v % 9 == r else 0.0 for v in vs])
               for r in (1, 2, 4, 5, 7)]
steps = [("direction (v mod 3)", f_dir), ("+ first rung m1", f_m1),
         ("+ rung springness", f_al), ("+ ternary trailing ones", f_t1),
         ("+ log2 v", f_lv), ("+ mod-9 dummies", mod9)]
prev = 0.0
for name, X in steps:
    r = r2(X)
    print(f"  {name:<26} R2 = {r:.3f}  (+{r - prev:.3f})")
    prev = r
