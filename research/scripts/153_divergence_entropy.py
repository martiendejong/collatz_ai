"""
153_divergence_entropy.py
==========================
THE DIVERGENCE CASE: rigorous constraints + the unified entropy gap.

Setup (V-map on odd n): n_{i+1} = (3 n_i + 1) / 2^{e_i},  e_i = v2(3 n_i + 1).
After m odd steps with b_m = e_1 + ... + e_m halvings:

    log2 n_m  =  log2 n_0 + m*theta - b_m + sum_i log2(1 + 1/(3 n_i)),
    theta = log2(3).

D1 (rigorous): every element of a divergent orbit exceeds the verification
    limit N0 (otherwise it converges), so the correction term is < m*delta,
    delta = log2(1 + 1/(3 N0)) ~ 2^-69 -- the SAME delta as the cycle case.

D2 (rigorous): divergence (n_m >= n_0 for all m) forces the parity vector
    into the thin slab
        m <= b_m <= m*(theta + delta)   for all m.
    A generic 2-adic point has b_m/m -> 2; divergence demands b_m/m <= 1.585.

D3 (rigorous, Terras-style with sharp constant): the first m exponents
    (e_1..e_m) determine n_0 mod 2^{b_m}; each pattern occupies exactly one
    residue class of density 2^{-b_m}. Hence the density of starting values
    still in the slab after m odd steps is at most

        S_m = sum_{b <= m*theta'} binom(b-1, m-1) * 2^{-b},   theta' = theta+delta.

    Large deviations: with b = beta*m,  log2(binom(beta m, m) 2^{-beta m})/m
    -> g(beta) = beta*H(1/beta) - beta, maximized on [1, theta] at beta=theta:

        g(theta) = theta*H(1/theta) - theta  =  -(1 - H(1/theta)) * theta.

    THE SAME CONSTANT AS THE CYCLE GAP, converted per odd step:
    cycle gap 1-H(1/theta) bits/halving  x  theta halvings/odd-step.

This script: (1) exact rate computation; (2) exact finite-m slab densities
S_m (rigorous upper bounds on the density of not-yet-dropped starters);
(3) Monte-Carlo confirmation of the decay exponent on 128-bit numbers.
"""
from math import log2, comb, log
import random

THETA = log2(3)

def H(x):
    return -x*log2(x) - (1-x)*log2(1-x)

print("="*74)
print("PART 1: THE UNIFIED ENTROPY GAP")
print("="*74)
cycle_gap_per_halving = 1 - H(1/THETA)
div_gap_per_oddstep = THETA - THETA*H(1/THETA)
print(f"cycle gap:      1 - H(1/theta)          = {cycle_gap_per_halving:.6f} bits/halving")
print(f"divergence gap: theta*(1 - H(1/theta))  = {div_gap_per_oddstep:.6f} bits/odd-step")
print(f"identity check: {cycle_gap_per_halving:.6f} * {THETA:.6f} = "
      f"{cycle_gap_per_halving*THETA:.6f}  == divergence gap: "
      f"{abs(cycle_gap_per_halving*THETA - div_gap_per_oddstep) < 1e-12}")
print()
print("g(beta) = beta*H(1/beta) - beta on [1, theta] (must be increasing, max at theta):")
for beta in [1.0001, 1.2, 1.4, THETA]:
    g = beta*H(1/beta) - beta
    print(f"  g({beta:.4f}) = {g:+.6f}")
print()

print("="*74)
print("PART 2: EXACT SLAB DENSITIES S_m (rigorous upper bounds)")
print("="*74)
print("S_m = density of odd n whose orbit has b_j <= j*theta for all j<=m is")
print("bounded by the single-time bound sum_{b<=m*theta} C(b-1,m-1) 2^-b:")
print(f"{'m':>6} {'log2(S_m bound)':>16} {'rate/m':>10}  (-> {-div_gap_per_oddstep:.4f})")
from math import lgamma
def l2comb(n, k):
    return (lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)) / log(2)
for m in [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
    bmax = int(m*THETA)
    terms = [l2comb(b-1, m-1) - b for b in range(m, bmax+1)]
    mx = max(terms)
    l2 = mx + log2(sum(2.0**(t - mx) for t in terms))
    print(f"{m:>6} {l2:>16.2f} {l2/m:>10.5f}")
print()

print("="*74)
print("PART 3: MONTE-CARLO — drop-time tail on 128-bit odd integers")
print("="*74)
print("m*(n) = number of odd steps until the orbit first drops below n.")
print("Prediction: P(m* > m) decays ~ 2^(-0.0793 m) (polynomial corrections).")
random.seed(2026)
N_SAMPLES = 400_000
MAXM = 120
surv = [0]*(MAXM+1)
for _ in range(N_SAMPLES):
    n0 = random.getrandbits(128) | 1
    n = n0
    for m in range(1, MAXM+1):
        x = 3*n + 1
        e = (x & -x).bit_length() - 1
        n = x >> e
        if n < n0:
            break
        surv[m] += 1

print(f"{'m':>5} {'P(m*>m)':>12} {'log2 P':>9} {'slope (last 10)':>16}")
import math
prev = None
for m in [5, 10, 20, 30, 40, 50, 60, 80, 100]:
    p = surv[m]/N_SAMPLES
    if p == 0:
        print(f"{m:>5} {'0':>12}")
        continue
    l2 = log2(p)
    print(f"{m:>5} {p:>12.6f} {l2:>9.2f}")
# fit slope on m in [30, 70]
import numpy as np
ms = [m for m in range(30, 71) if surv[m] > 50]
ys = [log2(surv[m]/N_SAMPLES) for m in ms]
slope = np.polyfit(ms, ys, 1)[0]
print(f"\nfitted decay exponent on m in [30,70]: {slope:.4f} bits/odd-step")
print(f"theoretical (Cramer, boundary-free):   {-div_gap_per_oddstep:.4f} bits/odd-step")
print("""
(The empirical slope should be slightly steeper than the single-time Cramer
rate because survival requires b_j <= j*theta for ALL j <= m, a ballot-type
constraint; agreement of leading order is the check.)""")
