"""
141_phantom_k5.py
=================
Extend phantom spectrum to K=5 (D3=179) and verify.
Also: investigate the information decay rate in Collatz orbits.
"""
import math
from sympy import n_order

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

print("=" * 70)
print("PART 1: K=5, l0=1 PHANTOM SPECTRUM (D3=179)")
print("=" * 70)
print()

K, l0 = 5, 1
D3 = 3**K - 2**(K+l0)
print(f"K={K}, l0={l0}: D3 = 3^{K} - 2^{K+l0} = {3**K} - {2**(K+l0)} = {D3}")
print(f"Validity check: 2^{{K+l0}} = {2**(K+l0)} < D3 = {D3}? {2**(K+l0) < D3}")
print()

D = int(n_order(2, D3))
print(f"ord_{{D3}}(2) = ord_{{{D3}}}(2) = {D}")
print(f"Phantom N values: N = {D}*j - 1 for j=1,2,3,...")
print(f"First few: N = {', '.join(str(D*j-1) for j in range(1, 8))}")
print()

# Find valid phantoms
valid_phantoms = []
for j in range(1, 20):
    N = D * j - 1
    numerator = (1 << l0) * ((1 << N) - 1) + 1
    if numerator % D3 != 0: continue
    m = numerator // D3
    if m % 2 == 0: continue  # m must be odd
    n = m * (1 << K) - 1
    if n >= (1 << N): continue  # must be valid representative

    # Verify l0
    x = m * (3**K) - 1
    actual_l0 = v2(x)
    if actual_l0 != l0: continue

    # Verify directly
    n_out, K_act, l0_act = macro_step(n)
    is_phantom = (n_out == n + (1 << N))
    valid_phantoms.append((N, m, n, n_out, is_phantom))
    print(f"j={j}: N={N}, m={m}, n={n}")
    print(f"  macro_step({n}) = {n_out}")
    print(f"  n + 2^{N} = {n + (1<<N)}")
    print(f"  Phantom confirmed: {is_phantom}")
    print()
    if len(valid_phantoms) >= 4:
        break

if not valid_phantoms:
    print("No valid phantoms found for this type within j=1..20.")

print()
print("=" * 70)
print("PART 2: COMPLETE PHANTOM SPECTRUM SURVEY")
print("=" * 70)
print()
print("All (K,l0) types with D3 > 2^{K+l0} and first phantom N <= 100:")
print()
print(f"{'(K,l0)':>12} {'D3':>8} {'ord':>8} {'N_1':>8} {'N_2':>8} {'N_3':>8}")
print("-" * 60)

spectrum = []
for K in range(1, 15):
    for l0 in range(1, 8):
        D3 = 3**K - 2**(K+l0)
        if D3 <= 1 or 2**(K+l0) >= D3:
            continue
        try:
            D = int(n_order(2, D3))
        except Exception:
            continue

        # Find first few valid phantoms
        phantoms_j = []
        for j in range(1, 50):
            N = D * j - 1
            if N > 200: break
            numerator = (1 << l0) * ((1 << N) - 1) + 1
            if numerator % D3 != 0: continue
            m = numerator // D3
            if m % 2 == 0: continue
            n = m * (1 << K) - 1
            if n >= (1 << N): continue
            x = m * (3**K) - 1
            if v2(x) != l0: continue
            phantoms_j.append(N)

        if phantoms_j and phantoms_j[0] <= 100:
            spectrum.append((phantoms_j[0], K, l0, D3, D, phantoms_j[:3]))

spectrum.sort()
for N1, K, l0, D3, D, Ns in spectrum:
    Ns_str = ', '.join(str(x) for x in Ns)
    print(f"  ({K},{l0}):   {D3:>8}  {D:>8}  {Ns_str}")

print()
print(f"Total types with first phantom N <= 100: {len(spectrum)}")

print()
print("=" * 70)
print("PART 3: INFORMATION DECAY IN COLLATZ ORBITS")
print("=" * 70)
print()
print("How quickly does orbit n_T 'forget' its starting value n_0?")
print()
print("Method: Sample many pairs (n_0, n_0+2^k) for large k (differ in 1 bit at position k).")
print("Compute correlation of final orbit positions after T steps.")
print()
import random as _r
_r.seed(42)

b = 100  # bit length of starting numbers

print(f"Starting numbers: {b}-bit random odd integers.")
print(f"Perturbation: flip bit k of n_0 to get n_0' = n_0 XOR 2^k.")
print()

for k in [0, 5, 10, 20, 50, 95]:
    # For each perturbation position k, compute how many macro-steps
    # until the orbit of n_0 and n_0 XOR 2^k diverge (n_T != n_T')
    n_pairs = 200
    divergence_steps = []

    for _ in range(n_pairs):
        n0 = _r.getrandbits(b) | 1
        n0_prime = n0 ^ (1 << k)
        if n0_prime % 2 == 0: n0_prime ^= 1  # keep odd

        n_a = n0; n_b = n0_prime
        t = 0
        max_t = 5000

        while n_a != n_b and t < max_t:
            n_a, _, _ = macro_step(n_a)
            n_b, _, _ = macro_step(n_b)
            t += 1
            if n_a < 2 or n_b < 2:
                break

        divergence_steps.append(t if t < max_t else max_t)

    mean_div = sum(divergence_steps) / len(divergence_steps)
    max_div = max(divergence_steps)
    # Note: "divergence" here means when they FIRST MEET (collide), not when they diverge!
    # After collision, both orbits are identical. Let's re-measure collision time.
    print(f"k={k:>3} (bit {k} flipped): mean collision time = {mean_div:.1f}, max = {max_div}")

print()
print("Expected: collision time ~ (b - k) * 1.2 / (something)")
print("When k=0 (LSB flipped): two odd numbers differ in bit 0 -- but both are odd, so bit 0 is always 1. Hmm.")
print()

print("Actually: measuring how quickly two orbits COLLIDE (visit same value).")
print("This happens when both orbits first visit the same value n_t = n_t'.")
print()
print("Collision time should scale with how DIFFERENT n_0 and n_0' are.")
print("Difference |n_0 - n_0'| = 2^k. By coupling: collision time ~ k/4 steps? Or ~k * 1.2?")

print()
print("=" * 70)
print("PART 4: FIRST COLLISION DISTRIBUTION")
print("=" * 70)
print()
_r.seed(123)
b = 200

print(f"Testing collision time for {b}-bit starting numbers (200 pairs each k)")
print()
print(f"{'k (bit flipped)':>18} {'|diff|':>12} {'mean T_coll':>14} {'T_coll/k':>12}")
print("-" * 60)

for k in [1, 5, 10, 20, 50, 100, 150, 190]:
    n_pairs = 200
    collision_times = []

    for _ in range(n_pairs):
        n0 = _r.getrandbits(b) | 1
        # Perturb bit k
        n0_prime = n0 ^ (1 << k)
        if n0_prime % 2 == 0: n0_prime ^= 1

        n_a = n0; n_b = n0_prime
        t = 0; collided = False
        max_t = 5000

        while t < max_t:
            n_a, _, _ = macro_step(n_a)
            n_b, _, _ = macro_step(n_b)
            t += 1
            if n_a == n_b:
                collided = True
                break
            if n_a < 2 and n_b < 2:
                collided = True
                break

        collision_times.append(t)

    mean_T = sum(collision_times) / len(collision_times)
    print(f"{k:>18} {2**k:>12.2e} {mean_T:>14.1f} {mean_T/k:>12.4f}")

print()
print("If T_coll ~ C * k: ratio T_coll/k should be constant = C.")
print("Theory from 2-adic expansion rate: C = 1/(4?) since the 2-adic valuation")
print("of the difference grows by ~4 per step (Obs 274 formula).")
