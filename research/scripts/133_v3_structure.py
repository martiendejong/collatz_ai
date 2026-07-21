"""
133_v3_structure.py
====================
3-adic structure along Collatz orbits: distribution of v3(n+1).

Key question: what is the joint distribution of (v2(n+1), v3(n+1)) = (K, J)
along Collatz orbits? Is K independent of J? What governs J transitions?

Background:
- K = v2(n+1) is the 2-adic valuation (the macro-step index)
- J = v3(n+1) is the 3-adic valuation
- For the inverse map: n has a predecessor from (K, l0) iff 3^K | (2^l0 * n' + 1)
  This depends on l0 mod (2 * 3^{K-1}), which links l0 to v3(2^l0 - 1).

From theory:
  macro_step(n) = (m * 3^K - 1) / 2^l0
  n_out + 1 = (m * 3^K + 2^l0 - 1) / 2^l0

  v3(n_out + 1) = v3(m * 3^K + 2^l0 - 1)

  Key: v3(2^l0 - 1) = 0 if l0 odd; = 1 + v3(l0/2) if l0 even.
  (Uses: ord_3(2) = 2, and lifting-the-exponent formula.)

  If l0 odd: v3(n_out + 1) = 0 (since 3^K | (m * 3^K) but 3 nmid (2^l0-1) when l0 odd)
  If l0 even: v3(n_out + 1) = min(K, 1 + v3(l0/2)) [approximately]
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

def v3(x):
    if x == 0: return 999
    c = 0
    while x % 3 == 0: x //= 3; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("=" * 70)
print("PART 1: DISTRIBUTION OF v3(n+1) ALONG ORBITS")
print("=" * 70)
print()

_r.seed(42)
n = _r.getrandbits(3000) | 1
K_vals = []; J_vals = []; l0_vals = []
n_steps = 200000

for _ in range(n_steps):
    K = v2(n+1); J = v3(n+1)
    K_vals.append(K); J_vals.append(J)
    n_out, _, l0 = macro_step(n)
    l0_vals.append(l0)
    n = n_out
    if n < 2:
        n = _r.getrandbits(3000) | 1

K_arr = np.array(K_vals); J_arr = np.array(J_vals); l0_arr = np.array(l0_vals)

print("Distribution of J = v3(n+1):")
print(f"  {'J':>4} {'Count':>10} {'Prob':>10} {'Theory?':>12}")
for j in range(8):
    cnt = (J_arr == j).sum()
    prob = cnt / n_steps
    # Theory: J=0 with prob 2/3? Let's see
    print(f"  {j:>4} {cnt:>10} {prob:>10.6f}")

print()
print("Mean J:", J_arr.mean(), "  Var J:", J_arr.var())
print("Compare K: mean K=", K_arr.mean(), "  Var K=", K_arr.var())

print()
print("=" * 70)
print("PART 2: JOINT DISTRIBUTION OF (K, J)")
print("=" * 70)
print()

joint = defaultdict(int)
for K, J in zip(K_arr, J_arr):
    joint[(K, J)] += 1

print("Joint P(K,J):")
print(f"  {'K':>4} {'J':>4} {'P(K,J)':>10} {'P(K)':>8} {'P(J|K)':>10}")
for K in range(1, 7):
    for J in range(0, 5):
        cnt = joint.get((K, J), 0)
        cnt_K = (K_arr == K).sum()
        if cnt_K > 0:
            print(f"  {K:>4} {J:>4} {cnt/n_steps:>10.6f} {cnt_K/n_steps:>8.4f} {cnt/cnt_K:>10.6f}")

# Test independence: P(K,J) = P(K)*P(J)?
print()
print("Testing K-J independence: ratio P(K,J) / [P(K)*P(J)]")
P_K = {k: (K_arr==k).sum()/n_steps for k in range(1, 8)}
P_J = {j: (J_arr==j).sum()/n_steps for j in range(0, 6)}
for K in range(1, 5):
    for J in range(0, 4):
        p_kj = joint.get((K, J), 0) / n_steps
        p_k = P_K.get(K, 0)
        p_j = P_J.get(J, 0)
        if p_k * p_j > 0:
            ratio = p_kj / (p_k * p_j)
            print(f"  K={K}, J={J}: ratio={ratio:.4f}")

print()
print("=" * 70)
print("PART 3: v3(n+1) TRANSITION DYNAMICS")
print("=" * 70)
print()
print("Empirical: J_out = v3(n_out+1) as function of J_in and l0")
print()

# Record (J_in, l0_parity, K, J_out) transitions
J_in_vals = []; l0_par_vals = []; K_in_vals = []; J_out_vals = []

_r.seed(123)
n = _r.getrandbits(2000) | 1
for _ in range(100000):
    J_in = v3(n+1)
    K = v2(n+1)
    n_out, _, l0 = macro_step(n)
    J_out = v3(n_out+1)
    J_in_vals.append(J_in); l0_par_vals.append(l0 % 2)
    K_in_vals.append(K); J_out_vals.append(J_out)
    n = n_out
    if n < 2:
        n = _r.getrandbits(2000) | 1

J_in_arr = np.array(J_in_vals)
l0_par_arr = np.array(l0_par_vals)
K_in_arr = np.array(K_in_vals)
J_out_arr = np.array(J_out_vals)

print("P(J_out | l0 parity):")
for par in [0, 1]:
    mask = l0_par_arr == par
    J_outs = J_out_arr[mask]
    par_name = "l0 odd" if par == 1 else "l0 even"
    print(f"\n  {par_name} (n={mask.sum()}):")
    for j in range(5):
        prob = (J_outs == j).mean()
        print(f"    J_out={j}: {prob:.6f}")

print()
print("Theory prediction:")
print("  l0 odd  (prob 2/3): J_out = 0 always (predicted)")
print("  l0 even (prob 1/3): J_out = min(K, 1+v3(l0/2)) (predicted)")

print()
print("P(J_out | l0 even, K) — verifying J_out = min(K, 1+v3(l0/2)):")
# l0 even cases only
mask_even = l0_par_arr == 0
K_even = K_in_arr[mask_even]
J_out_even = J_out_arr[mask_even]
l0_even_vals = np.array(l0_vals[:100000])[mask_even]  # need l0 values
# Actually let's collect directly

J_out_theory = []
J_out_actual = []
_r.seed(456)
n = _r.getrandbits(2000) | 1
for _ in range(200000):
    K = v2(n+1); m = (n+1) >> K
    n_out, _, l0 = macro_step(n)
    if l0 % 2 == 0:  # l0 even
        v3_l0_2 = v3(l0 // 2)
        theory_J_out = min(K, 1 + v3_l0_2)
        actual_J_out = v3(n_out + 1)
        J_out_theory.append(theory_J_out)
        J_out_actual.append(actual_J_out)
    n = n_out
    if n < 2:
        n = _r.getrandbits(2000) | 1

J_out_theory = np.array(J_out_theory)
J_out_actual = np.array(J_out_actual)
match_rate = (J_out_theory == J_out_actual).mean()
print(f"\nTheory accuracy for l0-even cases: {match_rate:.6f} ({len(J_out_theory)} samples)")
print(f"Mean theory J_out: {J_out_theory.mean():.4f}")
print(f"Mean actual J_out: {J_out_actual.mean():.4f}")
print()
# Show where theory fails
wrong_mask = J_out_theory != J_out_actual
if wrong_mask.sum() > 0:
    print(f"Theory fails in {wrong_mask.sum()} cases ({100*wrong_mask.mean():.2f}%)")
    print("When J_out_theory != J_out_actual:")
    print("  theory values:", np.unique(J_out_theory[wrong_mask], return_counts=True))
    print("  actual values:", np.unique(J_out_actual[wrong_mask], return_counts=True))
else:
    print("Theory is EXACT: v3(n_out+1) = min(K, 1+v3(l0/2)) for all l0-even cases. CONFIRMED.")

print()
print("=" * 70)
print("PART 4: THEORETICAL DISTRIBUTION OF J = v3(n+1)")
print("=" * 70)
print()
print("From the transition rule:")
print("  l0 odd (prob 2/3): J_out = 0")
print("  l0 even (prob 1/3): J_out = min(K, 1+v3(l0/2))")
print()
print("Distribution of v3(l0/2) for l0 ~ Geom(1/2) conditioned on l0 even:")
print("  l0 ~ Geom(1/2): P(l0=k) = 1/2^k, P(l0 even) = 1/3")
print("  Given l0 even: P(l0=2j) = (1/2^{2j}) / (1/3) = 3/4^j for j=1,2,3,...")
print("  So l0/2 ~ distribution with P(l0/2=j) = 3/4^j")
print()

# P(v3(l0/2) = r) for l0 ~ Geom(1/2) conditioned on l0 even
# P(l0/2 = j) = 3/4^j for j=1,2,3,...
# P(v3(l0/2) = r) = P(3^r | l0/2, 3^{r+1} nmid l0/2)
# = sum_{j: v3(j)=r} 3/4^j

from functools import lru_cache

def P_v3_l0half_cond(r, max_j=200):
    """P(v3(l0/2) = r | l0 even)"""
    total = sum(3.0/4**j for j in range(1, max_j+1) if v3(j) == r)
    return total

print("P(v3(l0/2) = r | l0 even):")
for r in range(6):
    p = P_v3_l0half_cond(r)
    print(f"  r={r}: {p:.6f}")

print()
print("Theoretical P(J_out = j):")
print("  j=0: 2/3 (l0 odd) + 1/3 * P(min(K,1+v3(l0/2))=0) = 2/3 + 1/3 * 0 = 2/3")
print("       (since K>=1 always and 1+v3(l0/2) >= 1, min(K,...) >= 1)")
# So J=0 always when l0 odd (prob 2/3)
# J=1 when l0 even AND min(K, 1+v3(l0/2)) = 1
#         = P(l0 even) * P(K=1 or v3(l0/2)=0 | l0 even)
# Actually min(K, 1+v3(l0/2)) = 1 iff K=1 OR (K>=2 AND v3(l0/2)=0)
#         = 1 iff v3(l0/2)=0 OR K=1
# Hmm, min(K, 1+v3(l0/2)) = 1 iff min(K,j+1) = 1 where j=v3(l0/2)
# = 1 iff K=1 (and then min(1,j+1)=1) OR K>=2 and j+1=1 (j=0, impossible)
# Wait: min(K, 1+j) = 1 iff K=1 or 1+j=1 (j=0, but j=v3(l0/2)>=0, so j=0 gives 1+j=1; and K>=1 means min(K,1)=1 iff K=1)
# Actually if j=0: min(K, 1) = 1 for all K>=1. So j=0 always gives J_out=1.
# If j>=1: min(K, j+1) = 1 iff K=1.

p_l0_even = 1/3
p_j0 = sum(3.0/4**j for j in range(1, 200) if v3(j) == 0) / (1/3)  # P(v3(l0/2)=0 | l0 even)
# Ah wait I need to renormalize. Let me redo:
v3_probs = {}
for r in range(6):
    v3_probs[r] = P_v3_l0half_cond(r)  # Already conditional

# P(J=0) = P(l0 odd) = 2/3
p_J0 = 2/3
# P(J>=1) = P(l0 even) = 1/3
# P(J_out = j | l0 even) = P(min(K, 1+v3(l0/2)) = j)
# = P(v3(l0/2)=j-1, K>=j) + P(v3(l0/2)>=j, K=j)
# (Both conditions give min=j)
p_K = {k: 1/2**k for k in range(1, 20)}  # P(K=k) = 1/2^k

print("\nComputing theoretical P(J_out=j) from transition rule:")
p_J = {0: 2/3}
# P(J_out=j | l0 even) for j=1,2,3,...
for j in range(1, 6):
    # min(K, 1+v3(l0/2)) = j iff:
    # Case 1: v3(l0/2) = j-1 AND K >= j (so 1+v3(l0/2)=j <= K)
    # Case 2: v3(l0/2) >= j AND K = j (so K < 1+v3(l0/2), min = K = j)
    prob_case1 = v3_probs.get(j-1, 0) * sum(1/2**k for k in range(j, 20))  # P(v3=j-1)*P(K>=j)
    prob_case2 = sum(v3_probs.get(r,0) for r in range(j, 10)) * (1/2**j)   # P(v3>=j)*P(K=j)
    p_J_given_even_j = prob_case1 + prob_case2
    p_J[j] = p_l0_even * p_J_given_even_j

for j in sorted(p_J):
    empirical = (J_arr == j).mean()
    print(f"  J={j}: theory={p_J.get(j,0):.6f}, empirical={empirical:.6f}")

print()
print("=" * 70)
print("PART 5: CORRELATION BETWEEN J=v3(n+1) AND K=v2(n+1)")
print("=" * 70)
print()

corr = np.corrcoef(K_arr[:100000], J_arr[:100000])[0,1]
print(f"Pearson correlation K vs J: {corr:.6f}")
print()

# Are K and J independent?
# Theory says: K depends on binary representation, J depends on ternary representation
# The Collatz map mixes them via n+1 -> m*3^K
print("Conditional distribution P(J=0 | K=k):")
for k in range(1, 7):
    mask = K_arr == k
    if mask.sum() > 100:
        p_j0_given_k = (J_arr[mask] == 0).mean()
        p_j1_given_k = (J_arr[mask] == 1).mean()
        p_j2_given_k = (J_arr[mask] >= 2).mean()
        print(f"  K={k}: P(J=0)={p_j0_given_k:.4f}, P(J=1)={p_j1_given_k:.4f}, P(J>=2)={p_j2_given_k:.4f}")

print()
print("=" * 70)
print("PART 6: v3(n+1) AUTOCORRELATION ALONG ORBITS")
print("=" * 70)
print()

# Is v3(n+1) autocorrelated along orbits?
J_small = J_arr[:50000]
for lag in [1, 2, 3, 4, 5]:
    acf = np.corrcoef(J_small[:-lag], J_small[lag:])[0,1]
    print(f"  Lag-{lag} ACF of J=v3(n+1): {acf:.6f}")

print()
print("Compare: K autocorrelation (should be ~0):")
K_small = K_arr[:50000]
for lag in [1, 2]:
    acf = np.corrcoef(K_small[:-lag], K_small[lag:])[0,1]
    print(f"  Lag-{lag} ACF of K=v2(n+1): {acf:.6f}")

print()
print("=" * 70)
print("PART 7: v3(l0) DISTRIBUTION AND v3(n+1) DYNAMICS")
print("=" * 70)
print()

# In each step, v3(n_out+1) is determined by:
# - parity of l0
# - if l0 even: min(K, 1+v3(l0/2))
# Let's verify the formula v3(2^L - 1) = 1 + v3(L/2) for L even

print("Verification: v3(2^L - 1) = 1 + v3(L/2) for L even:")
print(f"  {'L':>6} {'v3(2^L-1)':>12} {'1+v3(L/2)':>12} {'Match':>8}")
for L in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 30, 36]:
    actual = v3(2**L - 1)
    formula = 1 + v3(L // 2)
    match = "OK" if actual == formula else "FAIL"
    print(f"  {L:>6} {actual:>12} {formula:>12} {match:>8}")

print()
print("This identity is exact: v3(2^L - 1) = 1 + v3(L/2) for all even L.")
print("Proof: ord_3(2) = 2, so by LTE: v3(2^L - 1) = v3(2^2 - 1) + v3(L/2)")
print("  for L = 2*(L/2) even, = 1 + v3(L/2). QED.")
