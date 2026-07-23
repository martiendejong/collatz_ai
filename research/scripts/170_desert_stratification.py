"""
170_desert_stratification.py
=============================
THE DECISIVE TEST: does desert depth explain the entire certificate field?

Theory consolidated (Obs 319 + paper algebra):
  s(a) = (2a-1)/3 on {a == 2 mod 3}:
   - unique fixed point a = -1  (2a-1 = 3a <=> a = -1)
   - 3-adically expanding: |s(a)-(-1)|_3 = 3 |a-(-1)|_3
   - death condition: s(a) == 0 mod 3 <=> a == 5 mod 9
   - depth law: P(d >= j) = (2/3)^j exactly (verify below)
   - fertile set = s-invariant Cantor set, dim = log_3 2 = 0.6309...

Predictions to test on the k=13 certificate:
  P1: empirical P(d >= j) = (2/3)^j exactly (combinatorial identity).
  P2: log2 v stratifies by d(m): mean log2 v increases ~ log2(3/2) = 0.585
      bits per fertility level (this DERIVES the funnel law H3).
  P3: d(m) explains most of the variance of log2 v (high R^2) -- if so, the
      certificate field is essentially a function of desert depth, and the
      limit theorem reduces to a 1-D strata model with (2/3)^j weights.
"""
import numpy as np
from math import log2

k = 13
K3 = 3**k
v = np.load(f"certificate_k{k}.npy").astype(np.float64)
v /= v.max()
N = 3**(k-1)
i = np.arange(N, dtype=np.int64)
m = (3*i + 2) % K3          # class representatives mod 3^k

# desert depth: iterate s on residues; s well-defined mod 3^(level-1)
depth = np.zeros(N, dtype=np.int32)
alive = np.ones(N, dtype=bool)
a = m.copy()
level = k
maxd = k - 2
for j in range(1, maxd + 1):
    # death test: a == 5 mod 9  (needs level >= 2)
    dies = alive & (a % 9 == 5)
    depth[dies] = j
    alive &= ~dies
    if not alive.any() or level < 3:
        break
    a = np.where(alive, (2*a - 1) // 3 % 3**(level-1), a)
    level -= 1
depth[alive] = maxd + 1      # survived all resolvable levels

print("P1: depth distribution vs (2/3)^j")
print(f"{'j':>3} {'P(d>=j) emp':>12} {'(2/3)^j':>10}")
for j in range(1, 9):
    emp = float((depth >= j).mean())
    print(f"{j:>3} {emp:>12.6f} {(2/3)**j:>10.6f}")
print()

print("P2/P3: stratification of log2 v by desert depth")
lv = np.log2(v)
print(f"{'d':>3} {'#classes':>9} {'mean log2 v':>12} {'std':>6} {'delta':>7}")
prev = None
means, weights = [], []
for j in range(1, maxd + 2):
    sel = depth == j
    n = int(sel.sum())
    if n == 0:
        continue
    mu = float(lv[sel].mean()); sd = float(lv[sel].std())
    d = mu - prev if prev is not None else float('nan')
    print(f"{j:>3} {n:>9} {mu:>12.3f} {sd:>6.3f} {d:>7.3f}")
    prev = mu
    means.append(mu); weights.append(n)
print(f"  predicted delta per level: log2(3/2) = {log2(1.5):.3f}")
# R^2 of log2 v explained by depth strata
grand = lv.mean()
ss_tot = float(((lv - grand)**2).sum())
ss_between = sum(w * (mu - grand)**2 for mu, w in zip(means, weights))
print(f"\nP3: variance of log2 v explained by depth strata: R^2 = "
      f"{ss_between/ss_tot:.3f}")
