"""
167_bottleneck_analysis.py
===========================
WHAT REMAINS TO BE SOLVED, PART A: anatomy of the certificate vectors.

The Perron/certificate vectors c^m (k = 12..17) encode exactly what the
Krasikov hierarchy can and cannot certify: classes with SMALL c^m are the
bottleneck -- the method guarantees the fewest integers there, and they
control the rate 1 - gamma(k).

Hypotheses to test:
 H1 (D2-starvation): bottleneck classes are those whose 4m-orbit spends the
    most steps in the D2 branch (m == 5 mod 9), which contributes only
    lambda^-2 with NO second term. Prediction: bottom classes have a
    systematically higher D2-frequency along their 4m-orbit.
 H2 (self-similarity): the bottleneck set is nested across k (the level-k+1
    bottom classes project onto the level-k bottom classes) -- i.e. the
    bottleneck converges to a 3-adic closed set (a subshift), which would
    explain the power-law rate 1-gamma(k) ~ k^-0.85 and open the road to
    PROVING the limit theorem.
 H3 (spread): the ratio max c / min c grows with k (the LP stretches ever
    further to cover the bottleneck).
"""
import numpy as np
from math import log2

def load(k):
    v = np.load(f"certificate_k{k}.npy").astype(np.float64)
    return v / v.max()

def orbit_D2_fraction(i0, k, steps=60):
    """Fraction of D2-steps (i == 1 mod 3) along the 4m-orbit from compact i0."""
    N = 3**(k-1)
    i = i0
    d2 = 0
    for _ in range(steps):
        if i % 3 == 1:
            d2 += 1
        i = (4*i + 2) % N
    return d2 / steps

ks = [12, 13, 14, 15, 16, 17]
print("H3 -- spread of the certificate vectors:")
print(f"{'k':>3} {'min c':>12} {'log2(max/min)':>14}")
data = {}
for k in ks:
    try:
        v = load(k)
    except FileNotFoundError:
        continue
    data[k] = v
    print(f"{k:>3} {v.min():>12.3e} {log2(1.0/v.min()):>14.2f}")
print()

print("H1 -- D2-starvation at the bottleneck (k=15):")
k = 15
v = data[k]
N = 3**(k-1)
order = np.argsort(v)
bottom = order[:2000]
top = order[-2000:]
rng = np.random.default_rng(0)
rand = rng.integers(0, N, 2000)
for name, idx in [("bottom-2000", bottom), ("random-2000", rand), ("top-2000", top)]:
    fr = np.mean([orbit_D2_fraction(int(i), k) for i in idx])
    print(f"  {name:>12}: mean D2-fraction along 4m-orbit = {fr:.4f}")
print("  (uniform expectation = 1/3; higher at bottom supports H1)")
print()

print("H2 -- nesting of the bottleneck across k:")
for k in [12, 13, 14, 15, 16]:
    if k not in data or k+1 not in data:
        continue
    va, vb = data[k], data[k+1]
    Na, Nb = 3**(k-1), 3**k
    fa = 0.02          # bottom 2% at each level
    bota = set(np.argsort(va)[:int(fa*Na)].tolist())
    botb = np.argsort(vb)[:int(fa*Nb)]
    # project level-(k+1) compact indices to level-k: i mod 3^(k-1)
    proj = set((botb % Na).tolist())
    inter = len(bota & proj)
    # random baseline: projecting a random 2% set covers <= 2%*3 classes
    print(f"  k={k}->k={k+1}: |bottom_k cap proj(bottom_k+1)| / |bottom_k| "
          f"= {inter/len(bota):.3f}  (random baseline ~ {min(1.0, 3*fa):.3f})")
print()

print("Distribution shape (k=15): quantiles of log2 c")
q = np.quantile(np.log2(data[15]), [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
print("  q[0.1%,1%,10%,50%,90%,99%,99.9%] =", np.round(q, 2))
