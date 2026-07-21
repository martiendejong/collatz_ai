"""
130_orbit_coupling.py
======================
Orbit coupling: how quickly do two Collatz orbits merge?

Two starting points n and n + 2^L (differing by 2^L in 2-adic distance)
eventually merge if they produce the same orbit segment.

Key question: What is the expected COUPLING TIME T(n, n+2^L)?
That is, the first step t such that n_t = m_t (both orbits coincide).

This gives a direct measure of "2-adic mixing" that complements the
spectral gap analysis.

Also: ORBIT COUPLING THEOREM
If n ≡ m mod 2^N, then after ONE macro-step, n_out and m_out satisfy:
  n_out ≡ m_out mod 2^{N-K} (approximately)
So the 2-adic distance DECREASES by a factor of about 2^K per step.
The coupling time should be ~N / E[K] = N/2 steps.
"""
import numpy as np
from collections import defaultdict
import random as _r

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l = v2(x)
    return x >> l, K, l

print("=" * 70)
print("PART 1: ORBIT MERGING — HOW QUICKLY DO NEARBY ORBITS COALESCE?")
print("=" * 70)
print()

# Start n1 and n2 = n1 + 2^L, run both orbits until they merge
# or for max_steps steps.

def coupling_time(n1, n2, max_steps=10000):
    """Return (merge_time, did_merge) for orbits starting at n1 and n2."""
    for t in range(1, max_steps+1):
        n1_out, _, _ = macro_step(n1)
        n2_out, _, _ = macro_step(n2)
        n1 = n1_out; n2 = n2_out
        if n1 == n2:
            return t, True
        if n1 < 2 or n2 < 2:
            # If one reached 1 before merging, they merge at 1
            if n1 == n2:
                return t, True
            # Continue until both reach 1
    return max_steps, False

print("Testing coupling times for pairs (n, n + 2^L):")
print(f"{'L (gap)':>8} {'N_trials':>9} {'Mean T':>8} {'Max T':>8} {'P(merge)':>10}")
print("-" * 50)

_r.seed(123)
for L in [1, 2, 4, 8, 16, 32, 64]:
    times = []
    N_trials = 200
    for _ in range(N_trials):
        n1 = _r.getrandbits(200) | 1  # 200-bit random odd
        n2 = n1 + (1 << L)
        # Ensure n2 is also odd (n1 + 2^L: if L>=1, n1 odd, n2 = odd + even = odd) ✓
        t, merged = coupling_time(n1, n2, max_steps=5000)
        times.append(t)
    times = np.array(times)
    p_merge = (times < 5000).mean()
    print(f"{L:>8} {N_trials:>9} {times.mean():>8.2f} {times.max():>8} {p_merge:>10.3f}")

print()
print("=" * 70)
print("PART 2: 2-ADIC DISTANCE CONTRACTION PER STEP")
print("=" * 70)
print()
print("If n1 = n2 mod 2^N, after one macro-step do n1_out = n2_out mod 2^{N-K}?")
print("Theory: macro_step(n) = (m*3^K - 1) >> l0. If n1 = n2 mod 2^N:")
print("  n1 + 1 = a * 2^K1,  n2 + 1 = b * 2^K2 (family, K)")
print("  If K1 = K2 = K (same leading zeros): m1 = m2 mod 2^{N-K}")
print("  Then m1*3^K - 1 = m2*3^K - 1 mod 2^{N-K}")
print("  l0_1 may equal l0_2 or differ...")
print()

# Verify empirically: start n1, n2 with n1 ≡ n2 mod 2^N
# After 1 step, measure v2(n1_out - n2_out) as a function of N

print("2-adic distance v2(n1_out - n2_out) when v2(n1 - n2) = N:")
print(f"{'N':>4} {'E[v2(out diff)]':>18} {'Std':>8} {'Theory N-E[K]':>16}")
print("-" * 50)

_r.seed(456)
for N in [4, 8, 12, 16, 20]:
    diffs_out = []
    for _ in range(1000):
        n1 = _r.getrandbits(200) | 1
        # Create n2 ≡ n1 mod 2^N but n2 ≠ n1 (differ at bit N)
        n2 = n1 + (1 << N)  # n2 - n1 = 2^N exactly
        n1_out, K1, _ = macro_step(n1)
        n2_out, K2, _ = macro_step(n2)
        if n1_out != n2_out:
            diff = abs(n1_out - n2_out)
            d = v2(diff)
            diffs_out.append(d)
    if diffs_out:
        arr = np.array(diffs_out)
        print(f"{N:>4} {arr.mean():>18.4f} {arr.std():>8.4f} {N-2:>16.0f}")
    else:
        print(f"{N:>4} (all merged in 1 step)")

print()
print("=" * 70)
print("PART 3: COUPLING TIME vs BIT LENGTH")
print("=" * 70)
print()
print("For fixed gap L=1 (n, n+2), does coupling time grow with bit length b?")
print("Theory: T_couple ~ b/E[K] = b/2 (each step reduces 2-adic distance by E[K]=2)")

_r.seed(789)
for b in [20, 50, 100, 200, 500]:
    times = []
    N_trials = 100
    for _ in range(N_trials):
        n1 = _r.getrandbits(b) | 1
        n2 = n1 + 2  # differ by 2 = 2^1, so 2-adic distance is 1
        t, merged = coupling_time(n1, n2, max_steps=10000)
        times.append(t)
    times = np.array(times)
    print(f"b={b:>4}: mean T={times.mean():.1f}, std={times.std():.1f}, "
          f"max T={times.max()}, theory b/2={b/2:.0f}")

print()
print("=" * 70)
print("PART 4: WHEN DO ORBITS MERGE? — THE CONVERGENCE MECHANISM")
print("=" * 70)
print()
print("When n1 and n2 merge, what is the structure just before merging?")
print("Look at the K and l0 values in the last few steps before merge.")

_r.seed(321)
last_steps_before_merge = []
K_merge = defaultdict(int)  # K value at the merge step

for _ in range(500):
    n1 = _r.getrandbits(100) | 1
    n2 = n1 + 2  # start 2 apart
    history1 = [(n1, None, None)]
    history2 = [(n2, None, None)]
    for t in range(5000):
        n1_out, K1, l1 = macro_step(n1)
        n2_out, K2, l2 = macro_step(n2)
        history1.append((n1_out, K1, l1))
        history2.append((n2_out, K2, l2))
        n1 = n1_out; n2 = n2_out
        if n1 == n2:
            # Just merged! Look at last 3 steps
            last_steps_before_merge.append(t+1)
            K_merge[K1] += 1
            break

if last_steps_before_merge:
    arr = np.array(last_steps_before_merge)
    print(f"\nMerging time statistics (n2=n1+2, 100-bit starts):")
    print(f"  Mean T: {arr.mean():.2f}")
    print(f"  Std T:  {arr.std():.2f}")
    print(f"  Min/Max: {arr.min()}/{arr.max()}")
    print(f"\nK value at merge step:")
    for K in sorted(K_merge)[:10]:
        print(f"  K={K}: {K_merge[K]} times ({100*K_merge[K]/len(last_steps_before_merge):.1f}%)")
    # Theory: merging happens when n1 ≡ n2 mod 2^L for some L, then one more step makes them equal
    # The merge step typically has large l0 that bridges the difference
    print(f"\nMerge typically happens when K and l0 conspire to bring both orbits to same value.")

print()
print("=" * 70)
print("PART 5: COUPLING RATE EXPONENT")
print("=" * 70)
print()

# Fit T_couple = a * b^alpha
import math
b_vals = [20, 50, 100, 200, 500]
T_vals = []
_r.seed(777)
for b in b_vals:
    times = []; N_trials = 200
    for _ in range(N_trials):
        n1 = _r.getrandbits(b) | 1; n2 = n1 + 2
        t, _ = coupling_time(n1, n2, max_steps=20000)
        times.append(t)
    T_vals.append(np.mean(times))
    print(f"b={b:>4}: E[T]={np.mean(times):.2f}")

log_b = np.log(b_vals)
log_T = np.log(T_vals)
alpha, intercept = np.polyfit(log_b, log_T, 1)
print(f"\nFit: T ~ b^alpha, alpha = {alpha:.4f}")
print(f"Predicted T for b=1000: {math.exp(intercept) * 1000**alpha:.1f}")
print(f"Theory (T~b/2): alpha=1.00")
