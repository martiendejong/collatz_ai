"""
142_last_mile.py
================
Distribution of values in the final steps of Collatz orbits.

The orbit is a random walk ending at 1. What values does it visit
in the last 10, 20, 50 steps before reaching 1?

Key questions:
1. Is the distribution of n_{T-k} (k steps before end) concentrated?
2. Which small values are "gateway" values that almost all orbits pass through?
3. How many distinct "terminal paths" (sequences of last k values) exist?
4. Is there a "funnel" structure where all large orbits converge to a
   small set of values before reaching 1?
"""
import random as _r
import math
from collections import Counter, defaultdict

def v2(x):
    if x == 0: return 999
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

def collatz_orbit(n0):
    """Full orbit from n0 to 1, returns list of visited values."""
    orbit = [n0]
    n = n0
    while n > 1:
        n, _, _ = macro_step(n)
        orbit.append(n)
    return orbit

print("=" * 70)
print("PART 1: FREQUENCY OF LAST-k VALUES")
print("=" * 70)
print()
print("For 10000 random 200-bit starting numbers, record the last k steps.")
print()

_r.seed(42)
b = 200
n_trials = 5000

# Record last-k values for k = 1, 2, 5, 10, 20
last_k_counters = {k: Counter() for k in [1, 2, 5, 10, 20, 50]}

for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    for k in last_k_counters:
        if len(orbit) >= k:
            val = orbit[-(k+1)]  # value k+1 steps from end (since orbit ends at 1 = orbit[-1])
            last_k_counters[k][val] += 1

print(f"Based on {n_trials} random {b}-bit starting numbers:")
print()
for k in [1, 2, 5, 10, 20, 50]:
    counter = last_k_counters[k]
    total = sum(counter.values())
    top = counter.most_common(5)
    print(f"n_(T-{k}): {len(counter)} distinct values. Top 5:")
    for val, cnt in top:
        print(f"  n={val:>10}: {cnt:>6} times ({100*cnt/total:.2f}%)")
    print()

print()
print("=" * 70)
print("PART 2: FUNNEL VALUES — WHAT FRACTION OF ALL ORBITS PASS THROUGH n?")
print("=" * 70)
print()
print("Count how many of {1,...,1000} are visited by 10000 random 1000-bit orbits.")
print()

b = 1000; n_trials = 5000
# Count passage frequency through small values
passage_counts = Counter()

_r.seed(123)
for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    n = n0
    visited_small = set()
    while n > 1:
        if n <= 10000:
            visited_small.add(n)
        n, _, _ = macro_step(n)
    for v in visited_small:
        passage_counts[v] += 1

# Show the most frequently visited small values (all orbits pass through these)
top_small = [(v, passage_counts.get(v, 0)) for v in range(1, 201, 2) if passage_counts.get(v, 0) > 0]
top_small.sort(key=lambda x: -x[1])

print(f"Most frequently visited odd values n <= 200 (by {n_trials} random {b}-bit orbits):")
print(f"{'n':>8} {'count':>8} {'%':>8} {'log2(n)':>10}")
print("-" * 45)
for v, cnt in top_small[:30]:
    print(f"{v:>8} {cnt:>8} {100*cnt/n_trials:>8.2f}% {math.log2(v):>10.2f}")

print()
print("Values with 100% passage rate (universal funnels):")
universal = [v for v, cnt in top_small if cnt == n_trials]
print(f"  {universal}")

print()
print("=" * 70)
print("PART 3: DISTINCT TERMINAL PATHS (last k macro-steps)")
print("=" * 70)
print()

_r.seed(456)
b = 200
n_trials = 10000

terminal_paths = {k: Counter() for k in [5, 10, 15, 20]}

for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    for k in terminal_paths:
        if len(orbit) > k:
            path = tuple(orbit[-(k+1):-1])  # last k values before final 1
            terminal_paths[k][path] += 1

print(f"Based on {n_trials} random {b}-bit starting numbers:")
print()
for k, counter in sorted(terminal_paths.items()):
    n_distinct = len(counter)
    top3 = counter.most_common(3)
    print(f"Distinct terminal paths (last {k} steps): {n_distinct}")
    print(f"  Top 3 paths and their freq:")
    for path, cnt in top3:
        print(f"    {path[:5]}... : {cnt} times ({100*cnt/n_trials:.2f}%)")
    print()

print()
print("=" * 70)
print("PART 4: LAST-k BIT LENGTH DISTRIBUTION")
print("=" * 70)
print()
print("At step T-k, how large is n? (bit length distribution)")
print()

_r.seed(789)
b = 200
n_trials = 5000

for k in [0, 1, 2, 5, 10, 20, 50]:
    bits_at_k = []
    for _ in range(n_trials):
        n0 = _r.getrandbits(b) | 1
        orbit = collatz_orbit(n0)
        if len(orbit) > k:
            val = orbit[-(k+1)] if k > 0 else orbit[-1]
            bits_at_k.append(val.bit_length())

    if bits_at_k:
        import statistics
        mean_bits = statistics.mean(bits_at_k)
        std_bits = statistics.stdev(bits_at_k)
        max_bits = max(bits_at_k)
        print(f"T-{k:>2}: mean bits = {mean_bits:.2f}, std = {std_bits:.2f}, max = {max_bits}")

print()
print("Expected from random walk with drift -0.575:")
print("At T-k, the orbit is k steps from 1, so mean log2(n) ~ k * 0.693 / 0.575 ~ k * 1.2")
print("(since each step increases bit length by 0.693/0.575 ~ 1.2 in expectation if we reverse)")
for k in [0, 1, 2, 5, 10, 20, 50]:
    print(f"  T-{k:>2}: expected mean bits ~ {k * 1.205:.1f}")

print()
print("=" * 70)
print("PART 5: GATEWAY VALUE CENSUS")
print("=" * 70)
print()
print("Which small odd n have macro_step(n) > n (growth steps)?")
print("These are 'expansion nodes' where the orbit temporarily grows.")
print()

expansion_nodes = []
for n in range(1, 500, 2):
    n_out, K, l0 = macro_step(n)
    if n_out > n:
        expansion_nodes.append((n, n_out, K, l0, n_out/n))

print(f"{'n':>6} {'n_out':>8} {'K':>4} {'l0':>4} {'ratio':>8}")
print("-" * 40)
for n, n_out, K, l0, ratio in expansion_nodes[:20]:
    print(f"{n:>6} {n_out:>8} {K:>4} {l0:>4} {ratio:>8.3f}")
print(f"... (total {len(expansion_nodes)} expansion nodes in [1, 499])")
