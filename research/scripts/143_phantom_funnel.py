"""
143_phantom_funnel.py
=====================
Rigorously verify the phantom-funnel correspondence:
  - Phantom cycle elements at level N are the dominant gateway values
    at specific T-k depths in large Collatz orbits.

Key hypothesis: the phantom cycles form a "staircase" of attractor channels.
Orbits visiting N=9 phantom elements are ~10 steps from 1.
Orbits visiting N=7/8 phantom elements are ~13-16 steps from 1.
Orbits visiting N=10 phantom elements (703, 937) are even further.

Questions:
1. For each known phantom element, what is its mean T-k "funnel depth"?
2. Do higher-N phantoms correspond to larger T-k (further from 1)?
3. What fraction of orbits are absorbed into each phantom transit?
4. Are there phantom elements NOT in the funnel, or funnel nodes NOT in any phantom?
5. What is the real orbit length of the N=20 phantom fixed point (684783)?
6. Connection to N=10 phantom: where does {703, 937} sit in the funnel?
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

def collatz_orbit(n0, max_steps=100000):
    orbit = [n0]
    n = n0
    for _ in range(max_steps):
        if n == 1: break
        n, _, _ = macro_step(n)
        orbit.append(n)
    return orbit

# Known phantom elements from scripts 137-139
PHANTOMS = {
    7:  [47, 91, 103, 121],
    8:  [71, 91, 103, 121, 175, 189],
    9:  [91, 95, 103, 167, 175, 253, 283, 319, 399, 445],
    10: [703, 937],
}

ALL_PHANTOM_ELEMENTS = sorted(set(e for elems in PHANTOMS.values() for e in elems))

# Additional gateway candidates from script 142
EXTRA_GATEWAY = [5, 13, 23, 61, 433, 325, 577, 911]

print("=" * 70)
print("PART 1: FUNNEL DEPTH OF EACH PHANTOM ELEMENT")
print("=" * 70)
print()
print("For each phantom element p, run 5000 large random orbits.")
print("For each orbit passing through p, record the step T-k at which it visited p.")
print()

_r.seed(42)
b = 500
n_trials = 5000

# For each value of interest, record which T-k step it was visited at
# T-k means: orbit[-1] = 1, orbit[-2] = T-1 value, orbit[-(k+1)] = T-k value
# So if p appears at index i in orbit (0-indexed), T-k = len(orbit) - 1 - i

funnel_depths = defaultdict(list)  # element -> list of (T-k steps before 1)
passage_counts = Counter()

for trial in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    # Build index of where each value appears
    T = len(orbit) - 1  # total steps (orbit[T] = 1)
    visited_in_orbit = {}
    for i, v in enumerate(orbit):
        if v in visited_in_orbit:
            continue  # only record first visit
        if v in set(ALL_PHANTOM_ELEMENTS) or v in set(EXTRA_GATEWAY):
            visited_in_orbit[v] = i
    for v, i in visited_in_orbit.items():
        k_steps_from_end = T - i  # T-k where k = T-i
        funnel_depths[v].append(k_steps_from_end)
        passage_counts[v] += 1

print(f"Based on {n_trials} random {b}-bit starting numbers:")
print()
print(f"{'Element':>10} {'Phantom N':>12} {'Passage%':>10} {'MeanT-k':>10} {'StdT-k':>8} {'MedianT-k':>10}")
print("-" * 70)

import statistics

all_probe_values = sorted(set(ALL_PHANTOM_ELEMENTS + EXTRA_GATEWAY))

for v in all_probe_values:
    depths = funnel_depths.get(v, [])
    phantom_levels = [N for N, elems in PHANTOMS.items() if v in elems]
    level_str = ','.join(str(N) for N in sorted(phantom_levels)) if phantom_levels else 'none'
    pct = 100 * passage_counts[v] / n_trials
    if depths:
        mean_d = statistics.mean(depths)
        std_d = statistics.stdev(depths) if len(depths) > 1 else 0
        med_d = statistics.median(depths)
        print(f"{v:>10} {level_str:>12} {pct:>10.2f}% {mean_d:>10.2f} {std_d:>8.2f} {med_d:>10.1f}")
    else:
        print(f"{v:>10} {level_str:>12} {pct:>10.2f}%  (not visited)")

print()
print("=" * 70)
print("PART 2: PHANTOM N vs MEAN FUNNEL DEPTH (staircase test)")
print("=" * 70)
print()
print("Average funnel depth (mean T-k) grouped by phantom level N:")
print()

for N in sorted(PHANTOMS.keys()):
    elems = PHANTOMS[N]
    depths_for_N = []
    for v in elems:
        depths_for_N.extend(funnel_depths.get(v, []))
    if depths_for_N:
        mean_d = statistics.mean(depths_for_N)
        print(f"N={N}: elements={elems}, mean T-k = {mean_d:.2f}")
    else:
        print(f"N={N}: elements={elems}, (not visited)")

print()
print("If staircase hypothesis holds: mean T-k should INCREASE with N.")

print()
print("=" * 70)
print("PART 3: THE DOMINANT TERMINAL PATH RECONSTRUCTION")
print("=" * 70)
print()
print("Reconstruct the full dominant terminal path from n=1 backwards.")
print("At each step T-k, report the dominant value and whether it is a phantom element.")
print()

# Use 10000 orbits for this
_r.seed(123)
b = 500
n_trials = 10000
last_k_counters = {}
for k in range(1, 26):
    last_k_counters[k] = Counter()

for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1
    for k in range(1, 26):
        if T >= k:
            val = orbit[T - k]
            last_k_counters[k][val] += 1

print(f"{'T-k':>6} {'Top value':>12} {'Freq%':>8} {'Phantom N':>12} {'Cumul path dominance'}")
print("-" * 70)

dominant_path = []
for k in range(1, 26):
    counter = last_k_counters[k]
    total = sum(counter.values())
    if not total:
        continue
    top_val, top_cnt = counter.most_common(1)[0]
    pct = 100 * top_cnt / total
    phantom_levels = [N for N, elems in PHANTOMS.items() if top_val in elems]
    level_str = ','.join(str(N) for N in sorted(phantom_levels)) if phantom_levels else '---'
    dominant_path.append((k, top_val, pct, level_str))
    print(f"{k:>6} {top_val:>12} {pct:>8.2f}% {level_str:>12}")

print()
print("Phantom elements appear at T-k:")
for k, val, pct, lvl in dominant_path:
    if lvl != '---':
        print(f"  T-{k}: n={val} (phantom N={lvl}, passage={pct:.1f}%)")

print()
print("=" * 70)
print("PART 4: REAL ORBIT OF N=10 PHANTOM ELEMENTS {703, 937}")
print("=" * 70)
print()
print("How many steps does each N=10 phantom element take to reach 1?")
print("At what T-k step do large orbits visit these?")
print()

for p in [703, 937]:
    orbit = collatz_orbit(p)
    T = len(orbit) - 1
    print(f"n={p}: orbit length = {T} steps")
    max_val = max(orbit)
    print(f"  Max value: {max_val} ({max_val.bit_length()} bits)")
    print(f"  Orbit (first 20 values): {orbit[:20]}")
    print(f"  Orbit enters phantom path at: {[v for v in orbit if v in set(ALL_PHANTOM_ELEMENTS)]}")
    print()

print()
print("=" * 70)
print("PART 5: REAL ORBIT OF N=20 PHANTOM FIXED POINT (n=684783)")
print("=" * 70)
print()
p = 684783
orbit = collatz_orbit(p)
T = len(orbit) - 1
print(f"n={p}: orbit length = {T} steps")
print(f"Max value: {max(orbit)} ({max(orbit).bit_length()} bits)")
print(f"Orbit (all {T} steps):")
for i, v in enumerate(orbit):
    phantom_levels = [N for N, elems in PHANTOMS.items() if v in elems]
    lvl = f" <-- phantom N={','.join(str(N) for N in phantom_levels)}" if phantom_levels else ""
    print(f"  T-{T-i:>3} (step {i:>3}): n={v:>12}{lvl}")

print()
print("=" * 70)
print("PART 6: PHANTOM ELEMENT ABSORBER ANALYSIS")
print("=" * 70)
print()
print("How much of the 'funnel flow' passes through each phantom element?")
print("In what fraction of orbits is each phantom element the FIRST phantom-level element visited?")
print()

_r.seed(456)
b = 500
n_trials = 5000

all_phantom_set = set(ALL_PHANTOM_ELEMENTS)
first_phantom_hit = Counter()  # first phantom element encountered (from large n toward 1)

for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    for v in orbit:
        if v in all_phantom_set:
            first_phantom_hit[v] += 1
            break  # first phantom hit only

print(f"First phantom element hit (out of {n_trials} orbits):")
print(f"{'Element':>10} {'Count':>8} {'%':>8} {'Phantom N':>12}")
print("-" * 45)
for v, cnt in sorted(first_phantom_hit.items(), key=lambda x: -x[1]):
    pct = 100 * cnt / n_trials
    phantom_levels = [N for N, elems in PHANTOMS.items() if v in elems]
    lvl = ','.join(str(N) for N in sorted(phantom_levels))
    print(f"{v:>10} {cnt:>8} {pct:>8.2f}% {lvl:>12}")

print()
total_absorbed = sum(first_phantom_hit.values())
print(f"Total orbits hitting ANY phantom element: {total_absorbed}/{n_trials} = {100*total_absorbed/n_trials:.1f}%")

print()
print("=" * 70)
print("PART 7: EXPANSION NODES -- WHAT IS SPECIAL ABOUT n=127?")
print("=" * 70)
print()
print("n=127 = 2^7-1: K=v2(128)=7, m=1, x=3^7-1=2186, l0=v2(2186)=1, n_out=1093")
print("Ratio 1093/127 = 8.606 -- the highest ratio for small n")
print()
print("General expansion: n=2^K-1 for various K:")
print(f"{'K':>4} {'n=2^K-1':>10} {'n_out':>12} {'ratio':>8}")
print("-" * 40)
for K in range(1, 15):
    n = (1 << K) - 1
    n_out, K_act, l0 = macro_step(n)
    ratio = n_out / n if n > 0 else 0
    print(f"{K:>4} {n:>10} {n_out:>12} {ratio:>8.4f}")

print()
print("Pattern: n=2^K-1 has m=1 (since n+1=2^K), so n_out = (3^K-1)/2^l0.")
print("l0 = v2(3^K-1). The expansion ratio is n_out/n = (3^K-1) / (2^l0 * (2^K-1)).")
print()
print("For which K is l0=1? Need 3^K-1 = 2 mod 4, i.e., 3^K = 3 mod 4.")
print("3^K mod 4: 3,1,3,1,... -> K odd -> l0 = v2(3^K-1) = 1 when K is odd")
print("           K even -> 3^K-1 = 0 mod 4 -> l0 >= 2")
print()
print("Ratio for K odd (l0=1): (3^K-1) / (2*(2^K-1)) ~ 3^K/2^{K+1} = (3/2)^K / 2")
print("This grows exponentially! K=7: 3^7/2^8 = 2187/256 = 8.55 ~ ratio 8.606. [confirmed]")
print()

# Table of l0 values for n=2^K-1
print("l0 = v2(3^K - 1) for K = 1..20:")
for K in range(1, 21):
    x = (3**K) - 1
    l0 = v2(x)
    ratio = x / (l0 > 0 and (1 << l0) * ((1 << K) - 1) or 1)
    print(f"  K={K:>2}: 3^K-1={x:>10}, l0={l0}, ratio ~ {(3**K - 1) / ((1<<l0) * ((1<<K) - 1)):.4f}")
