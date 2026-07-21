"""
134_padic_orbit.py
===================
p-adic structure along Collatz orbits for primes p = 5, 7, 11, 13.

Questions:
1. What is the stationary distribution of v_p(n+1) for p = 5, 7, 11?
2. Is it different from the p-adic valuation of a random odd number?
3. Is K = v2(n+1) independent of v_p(n+1) for all primes p?
4. For p=3, we have the exact formula: v3(2^L-1) = 1+v3(L/2) for L even.
   What is the analogue for p=5, 7, 11?
   - For p=5: ord_5(2) = 4. 5 | 2^L - 1 iff 4 | L.
     v5(2^L-1) = v5(2^4-1) + v5(L/4) = 1 + v5(L/4) for 4|L, else 0.
   - For p=7: ord_7(2) = 3. 7 | 2^L - 1 iff 3 | L.
     v7(2^L-1) = 1 + v7(L/3) for 3|L, else 0.
   - For p=11: ord_11(2) = 10. 11 | 2^L-1 iff 10|L.
     v11(2^L-1) = 1 + v11(L/10) for 10|L, else 0.

5. Transition rule generalization:
   v_p(n_out+1) = min(K, v_p(2^{l0}-1)) where v_p(2^{l0}-1) = 1+v_p(l0/d)
   for d = ord_p(2) and d | l0, else 0.
"""
import numpy as np
from collections import defaultdict
import random as _r
import math

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def vp(x, p):
    if x == 0: return 999
    c = 0
    while x % p == 0: x //= p; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

PRIMES = [5, 7, 11, 13]

print("=" * 70)
print("PART 1: IDENTITY v_p(2^L - 1) FOR VARIOUS PRIMES p")
print("=" * 70)
print()
print("Verify: v_p(2^L - 1) = 1 + v_p(L/d) when d | L, else 0")
print("where d = ord_p(2) = multiplicative order of 2 mod p.")
print()

from sympy import n_order, mod_inverse

for p in PRIMES:
    d = n_order(2, p)
    print(f"p={p}, ord_{p}(2) = {d}:")
    print(f"  L  {'v_p(2^L-1)':>12} {'1+v_p(L/d)':>12} {'d|L?':>6} {'Match':>6}")
    for L in [d, 2*d, 3*d, 4*d, d*3, d*5, d*9, 2*d+1, d*15, d*25]:
        actual = vp(2**L - 1, p)
        if d == 0 or L % d != 0:
            formula = 0
        else:
            formula = 1 + vp(L // d, p)
        divides = "yes" if L % d == 0 else "no"
        match = "OK" if actual == formula else "FAIL"
        print(f"  {L:>3}  {actual:>12}  {formula:>12}  {divides:>6}  {match:>6}")
    print()

print()
print("=" * 70)
print("PART 2: EMPIRICAL DISTRIBUTION OF v_p(n+1) ALONG ORBITS")
print("=" * 70)
print()

_r.seed(42)
n = _r.getrandbits(3000) | 1
n_steps = 200000
K_vals = []
vp_vals = {p: [] for p in PRIMES}
l0_vals = []

for _ in range(n_steps):
    K = v2(n+1)
    K_vals.append(K)
    l0_vals.append(0)  # placeholder, will compute in macro_step
    for p in PRIMES:
        vp_vals[p].append(vp(n+1, p))
    n_out, _, l0 = macro_step(n)
    l0_vals[-1] = l0
    n = n_out
    if n < 2:
        n = _r.getrandbits(3000) | 1

K_arr = np.array(K_vals)
l0_arr = np.array(l0_vals)
vp_arrays = {p: np.array(vp_vals[p]) for p in PRIMES}

print("Distribution of v_p(n+1) along Collatz orbits:")
print(f"{'p':>4} {'d':>4} {'j':>4} {'P(vp=j)':>12} {'Theory Geom':>14} {'Ratio':>8}")

for p in PRIMES:
    d = n_order(2, p)
    J_arr = vp_arrays[p]
    p_uniform = (p-1)/p  # P(vp=0) for random number = (p-1)/p
    print(f"\np={p}, ord_{p}(2)={d}:")
    total_check = 0
    for j in range(6):
        prob = (J_arr == j).mean()
        total_check += prob
        # Theory for random odd number: P(vp(n+1)=j) = (p-1)/p * (1/p)^j = (p-1)/p^{j+1}
        # (since n+1 = 2^K * m, and vp(n+1) = vp(m) since gcd(2,p)=1 for p odd prime)
        # Among odd numbers: vp(m) has same distribution as for random integers
        # (since restricting to odd doesn't change p-adic structure for p=5,7,11,13)
        theory_random = ((p-1) / p) * (1/p)**j
        ratio = prob / theory_random if theory_random > 0 else float('nan')
        print(f"  {p:>4} {d:>4} {j:>4} {prob:>12.6f} {theory_random:>14.6f} {ratio:>8.4f}")
    print(f"  Total prob checked: {total_check:.6f}")

print()
print("=" * 70)
print("PART 3: IS v_p(n+1) INDEPENDENT OF K = v2(n+1)?")
print("=" * 70)
print()

for p in PRIMES:
    J_arr = vp_arrays[p]
    corr = np.corrcoef(K_arr, J_arr)[0, 1]
    print(f"p={p}: Pearson corr(K, v_p(n+1)) = {corr:.6f}")

print()
print("Conditional P(v_p(n+1)=0 | K=k) for p=5,7:")
for p in [5, 7]:
    J_arr = vp_arrays[p]
    print(f"\np={p}:")
    for k in range(1, 6):
        mask = K_arr == k
        if mask.sum() > 100:
            prob = (J_arr[mask] == 0).mean()
            print(f"  K={k}: P(vp=0)={prob:.4f} (vs overall {(J_arr==0).mean():.4f})")

print()
print("=" * 70)
print("PART 4: v_p(n+1) TRANSITION DYNAMICS")
print("=" * 70)
print()
print("Transition rule: v_p(n_out+1) governed by v_p(2^{l0}-1)")
print("  = 0 if d nmid l0  (where d = ord_p(2))")
print("  = 1 + v_p(l0/d) if d | l0")
print()
print("So: d | l0 with probability P(d|l0) for l0 ~ Geom(1/2)")
print()

for p in PRIMES:
    d = n_order(2, p)
    # P(d | l0) where l0 ~ Geom(1/2)
    p_d_divides = sum(1/2**k for k in range(d, 300, d))
    print(f"p={p}: d={d}, P(d | l0) = {p_d_divides:.6f} = 1/(2^d - 1) = {1/(2**d-1):.6f}")

print()
print("Empirical P(v_p(n_out+1) = 0):")
print("Theory: P(d nmid l0) = 1 - 1/(2^d - 1)")

for p in PRIMES:
    d = n_order(2, p)
    # Compute transitions empirically
    J_out_vals = []
    _r.seed(111)
    n = _r.getrandbits(2000) | 1
    for _ in range(100000):
        n_out, K, l0 = macro_step(n)
        J_out_vals.append(vp(n_out+1, p))
        n = n_out
        if n < 2:
            n = _r.getrandbits(2000) | 1
    J_out_arr = np.array(J_out_vals)

    p_j0_emp = (J_out_arr == 0).mean()
    p_j0_theory = 1 - 1/(2**d - 1)
    print(f"  p={p}: d={d}: P(J_out=0) = {p_j0_emp:.6f}, theory = {p_j0_theory:.6f}")

print()
print("=" * 70)
print("PART 5: v_p(n+1) AUTOCORRELATION")
print("=" * 70)
print()

for p in PRIMES:
    J_arr = vp_arrays[p][:50000]
    for lag in [1, 2]:
        acf = np.corrcoef(J_arr[:-lag], J_arr[lag:])[0, 1]
        print(f"p={p}, lag-{lag} ACF: {acf:.6f}")
    print()

print()
print("=" * 70)
print("PART 6: THEORY — STATIONARY DISTRIBUTION OF v_p(n+1)")
print("=" * 70)
print()
print("Prediction from the transition rule:")
print("  J_out = 0 with prob 1 - 1/(2^d - 1) [= when d nmid l0]")
print("  J_out >= 1 with prob 1/(2^d - 1)      [= when d | l0]")
print()
print("  Conditional on d | l0:")
print("  l0/d ~ Geom(1 - 1/2^d) on {1,2,3,...} [from Geom(1/2) conditioned on d|l0]")
print("  v_p(l0/d) follows a distribution depending on l0/d")
print()

for p in PRIMES:
    d = n_order(2, p)
    # P(d | l0) = sum_{k=1}^{inf} (1/2)^{dk} = (1/2^d) / (1 - 1/2^d) = 1/(2^d - 1)
    p_d_divides = 1/(2**d - 1)
    # Given d | l0: l0 = d*k where k ~ Geom(1 - 1/2^d), i.e., P(k=j) = (1/2^d)^j * (1-1/2^d) / (1/(2^d-1))
    # = (1 - 1/2^d) * (2^d - 1) * (1/2^d)^j  ... hmm let me redo
    # l0 ~ Geom(1/2): P(l0 = k) = 1/2^k for k=1,2,3,...
    # Given d | l0: l0 = d, 2d, 3d, ...
    # P(l0 = jd | d | l0) = (1/2^{jd}) / p_d_divides = (1/2^{jd}) * (2^d - 1)
    # This is Geometric with parameter q = 1/2^d: P(l0/d = j) = (1-q) * q^{j-1} ... wait:
    # Sum_{j=1}^{inf} (1/2^{jd}) = (1/2^d)/(1-1/2^d) = 1/(2^d-1) = p_d_divides. Good.
    # P(l0/d = j | d|l0) = (1/2^{jd}) / p_d_divides = (1/2^{jd}) * (2^d-1)
    # This is Geom on {1,2,...} with "success prob" 1-1/2^d... let me check:
    # = (1 - 1/2^d) * (1/2^d)^{j-1} * (1/2^d) ... no.
    # Actually (1/2^{jd}) * (2^d-1) = (2^d-1)/2^{jd}.
    # P(l0/d=1) = (2^d-1)/2^d = 1 - 1/2^d
    # P(l0/d=2) = (2^d-1)/2^{2d}
    # P(l0/d=j) = (2^d-1)/2^{jd} = (1-1/2^d) * (1/2^d)^{j-1} * (1/2^d) ...
    # = (1/2^d)^j * (2^d-1) which sums to (2^d-1)/2^d / (1-1/2^d) = ...
    # Let r = 1/2^d. P(j) = r^j*(1-r)/r = (1-r)*r^{j-1}? No:
    # Sum P(j) = (2^d-1) * Sum r^j = (1/r - 1) * r/(1-r) = (1-r)/r * r/(1-r) = 1. Good.
    # So l0/d ~ Geom(1/2^d) on {1,2,...} means P(l0/d=j) = (1-1/2^d)*(1/2^d)^{j-1}.
    # But (1/2^{jd})*(2^d-1) = (1-1/2^d)*(1/2^d)^{j-1} * (1/2^d) ... No:
    # (1-1/2^d) * (1/2^d)^{j-1} with j=1: (1-1/2^d). But we want (2^d-1)/2^d = 1-1/2^d for j=1. ✓
    # j=2: (2^d-1)/4^d = (1-1/2^d)/2^d. And (1-1/2^d)*(1/2^d)^1 = (1-1/2^d)/2^d. ✓
    # So l0/d | (d|l0) ~ Geom(1/2^d) on {1,...}: P(l0/d=j) = (1-1/2^d)*(1/2^d)^{j-1}.

    # P(v_p(l0/d) = r | d|l0):
    # = Sum_{j: v_p(j)=r} (1-1/2^d)*(1/2^d)^{j-1}
    p_vp0 = sum((1-1/2**d)*(1/2**d)**(j-1) for j in range(1, 200) if vp(j, p) == 0)
    p_vp1 = sum((1-1/2**d)*(1/2**d)**(j-1) for j in range(1, 200) if vp(j, p) == 1)
    p_vp2 = sum((1-1/2**d)*(1/2**d)**(j-1) for j in range(1, 200) if vp(j, p) == 2)

    print(f"p={p}, d={d}:")
    print(f"  P(d | l0) = 1/(2^d-1) = {p_d_divides:.6f}")
    print(f"  Given d|l0: P(v_p(l0/d)=0) = {p_vp0:.6f}")
    print(f"              P(v_p(l0/d)=1) = {p_vp1:.6f}")
    print(f"              P(v_p(l0/d)=2) = {p_vp2:.6f}")

    # Now: P(J_out = 0) = P(d nmid l0) = 1 - p_d_divides
    # P(J_out >= 1) = p_d_divides
    # P(J_out = j | J_out>=1) = P(min(K, 1+v_p(l0/d)) = j)
    # This is the same structure as for p=3 but with d instead of 2.
    p_J0 = 1 - p_d_divides
    print(f"  Predicted P(J_out=0) = {p_J0:.6f}")
    J_arr = vp_arrays[p]
    print(f"  Empirical P(J_out=0) = {(J_arr==0).mean():.6f}")
    print()
