"""
145_two_channels.py
====================
Two-channel structure of Collatz orbits.

The dominant T-2 values split orbits into two major channels:
  - "23-channel": T-2 = n=23 (macro_step chain: ...->61->23->5->1)
  - "13-channel": T-2 = n=13 (macro_step chain: ...->7/11->13->5->1)

Key question: does the 13-channel have its OWN phantom staircase at T-8 to T-16,
or is it structurally different from the 23-channel?

If yes: there should be a second phantom staircase corresponding to the 13-channel,
formed by DIFFERENT phantom cycle elements (or possibly the same ones).
If no: the 13-channel orbits are diffuse and don't concentrate at phantom elements.
"""
import random as _r
from collections import Counter, defaultdict
import math

def v2(x):
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

# All known phantom elements
PHANTOMS = {
    7:  {47, 91, 103, 121},
    8:  {71, 91, 103, 121, 175, 189},
    9:  {91, 95, 103, 167, 175, 253, 283, 319, 399, 445},
    10: {703, 937},
}
ALL_PHANTOM = set().union(*PHANTOMS.values())

print("=" * 70)
print("PART 1: SEPARATING THE TWO CHANNELS")
print("=" * 70)
print()
print("Classification of 500-bit orbits by T-2 value:")
print("  23-channel: T-2 = 23 (via ...->61->23->5->1)")
print("  13-channel: T-2 = 13 (via ...->7/11->13->5->1)")
print()

_r.seed(42)
b = 500
n_trials = 10000

channel_23_counters = {k: Counter() for k in range(1, 26)}
channel_13_counters = {k: Counter() for k in range(1, 26)}
other_counters = {k: Counter() for k in range(1, 26)}
channel_counts = Counter()

for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1

    if T < 2:
        continue

    t2_val = orbit[T - 2]
    if t2_val == 23:
        channel = "23"
    elif t2_val == 13:
        channel = "13"
    else:
        channel = "other"
    channel_counts[channel] += 1

    for k in range(1, 26):
        if T >= k:
            val = orbit[T - k]
            if channel == "23":
                channel_23_counters[k][val] += 1
            elif channel == "13":
                channel_13_counters[k][val] += 1
            else:
                other_counters[k][val] += 1

total = sum(channel_counts.values())
print(f"Channel distribution ({n_trials} trials):")
for ch, cnt in sorted(channel_counts.items()):
    print(f"  {ch}-channel: {cnt} ({100*cnt/total:.1f}%)")
print()

print()
print("=" * 70)
print("PART 2: DOMINANT T-k VALUES BY CHANNEL")
print("=" * 70)
print()
print(f"{'T-k':>6} {'23-ch top':>14} {'23%':>8} {'Phantom?':>10} | {'13-ch top':>14} {'13%':>8} {'Phantom?':>10}")
print("-" * 80)

for k in range(1, 26):
    c23 = channel_23_counters[k]
    c13 = channel_13_counters[k]
    n23 = sum(c23.values())
    n13 = sum(c13.values())

    if n23 > 0:
        top23, cnt23 = c23.most_common(1)[0]
        pct23 = 100 * cnt23 / n23
        ph23_N = [N for N, elems in PHANTOMS.items() if top23 in elems]
        ph23_str = ','.join(str(N) for N in sorted(ph23_N)) if ph23_N else "---"
    else:
        top23, pct23, ph23_str = "N/A", 0, "---"

    if n13 > 0:
        top13, cnt13 = c13.most_common(1)[0]
        pct13 = 100 * cnt13 / n13
        ph13_N = [N for N, elems in PHANTOMS.items() if top13 in elems]
        ph13_str = ','.join(str(N) for N in sorted(ph13_N)) if ph13_N else "---"
    else:
        top13, pct13, ph13_str = "N/A", 0, "---"

    print(f"{k:>6} {str(top23):>14} {pct23:>8.2f}% {ph23_str:>10} | {str(top13):>14} {pct13:>8.2f}% {ph13_str:>10}")

print()
print("=" * 70)
print("PART 3: FULL DOMINANT T-k PATH FOR EACH CHANNEL")
print("=" * 70)
print()
print("Reconstructing the dominant terminal path for each channel:")
print()

print("23-CHANNEL dominant path:")
for k in range(1, 20):
    c23 = channel_23_counters[k]
    if not c23: continue
    n23 = sum(c23.values())
    top23, cnt23 = c23.most_common(1)[0]
    pct23 = 100 * cnt23 / n23
    ph23_N = [N for N, elems in PHANTOMS.items() if top23 in elems]
    ph23_str = f" [phantom N={','.join(str(N) for N in sorted(ph23_N))}]" if ph23_N else ""
    print(f"  T-{k:>2}: n={top23:>8} ({pct23:.1f}%){ph23_str}")

print()
print("13-CHANNEL dominant path:")
for k in range(1, 20):
    c13 = channel_13_counters[k]
    if not c13: continue
    n13 = sum(c13.values())
    top13, cnt13 = c13.most_common(1)[0]
    pct13 = 100 * cnt13 / n13
    ph13_N = [N for N, elems in PHANTOMS.items() if top13 in elems]
    ph13_str = f" [phantom N={','.join(str(N) for N in sorted(ph13_N))}]" if ph13_N else ""
    print(f"  T-{k:>2}: n={top13:>8} ({pct13:.1f}%){ph13_str}")

print()
print("=" * 70)
print("PART 4: PHANTOM ELEMENT PASSAGE RATES BY CHANNEL")
print("=" * 70)
print()
print("For each phantom staircase element, what fraction of EACH CHANNEL passes through it?")
print()
print(f"{'Element':>8} {'All%':>8} {'23-ch%':>10} {'13-ch%':>10} {'Phantom N':>12}")
print("-" * 55)

# Count passage rates within each channel
PROBE_VALS = sorted(ALL_PHANTOM) + [5, 13, 23, 61, 325, 433, 577, 911]

_r.seed(123)
b = 500
n_trials2 = 5000
ch23_passage = Counter()
ch13_passage = Counter()
all_passage = Counter()
ch23_total = 0
ch13_total = 0
all_total = 0

for _ in range(n_trials2):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1
    if T < 2: continue

    t2_val = orbit[T - 2]
    channel = "23" if t2_val == 23 else ("13" if t2_val == 13 else "other")
    all_total += 1
    if channel == "23": ch23_total += 1
    elif channel == "13": ch13_total += 1

    visited = set(orbit)
    for v in PROBE_VALS:
        if v in visited:
            all_passage[v] += 1
            if channel == "23": ch23_passage[v] += 1
            elif channel == "13": ch13_passage[v] += 1

for v in sorted(PROBE_VALS):
    ph_N = [N for N, elems in PHANTOMS.items() if v in elems]
    lvl = ','.join(str(N) for N in sorted(ph_N)) if ph_N else "---"
    pct_all = 100 * all_passage[v] / all_total
    pct_23 = 100 * ch23_passage[v] / ch23_total if ch23_total else 0
    pct_13 = 100 * ch13_passage[v] / ch13_total if ch13_total else 0
    print(f"{v:>8} {pct_all:>8.1f}% {pct_23:>10.1f}% {pct_13:>10.1f}% {lvl:>12}")

print()
print(f"Channel totals: 23-ch = {ch23_total} ({100*ch23_total/all_total:.1f}%), 13-ch = {ch13_total} ({100*ch13_total/all_total:.1f}%)")

print()
print("=" * 70)
print("PART 5: THE 13-CHANNEL PREDECESSORS")
print("=" * 70)
print()
print("The 13-channel uses paths ending ...->7->13->5->1 or ...->11->13->5->1.")
print("What are the full dominant paths for the 13-channel?")
print()

# Trace the backward tree from n=13
# macro_step(7) = 13? Check: v2(8)=3, m=1, x=27-1=26, l0=1, n_out=13. Yes!
# macro_step(11) = 13? Check: v2(12)=2, m=3, x=27-1=26, l0=1, n_out=13. Yes!
# macro_step(17) = 13? Check: v2(18)=1, m=9, x=27-1=26, l0=1, n_out=13. Yes!

def predecessors_of(target, max_n=100000):
    result = []
    for K in range(1, 20):
        for l0 in range(1, 40):
            numerator = (1 << l0) * target + 1
            denom = 3**K
            if numerator % denom != 0: continue
            m = numerator // denom
            if m <= 0 or m % 2 == 0: continue
            n = m * (1 << K) - 1
            if n > max_n or n <= 0: continue
            n_out, K_act, l0_act = macro_step(n)
            if n_out == target:
                result.append((n, K, l0))
    return sorted(result)

print("Predecessors of n=13 (T-3 values for 13-channel):")
p13 = predecessors_of(13, 200000)
for n, K, l0 in p13[:15]:
    print(f"  n={n} (K={K}, l0={l0}, mod3={n%3})")

print()
print("The 13-channel at T-3 level: which values are most common?")
c13_k3 = channel_13_counters[3]
n13 = sum(c13_k3.values())
print(f"Top T-3 values in 13-channel ({n13} orbits):")
for val, cnt in c13_k3.most_common(10):
    print(f"  n={val}: {cnt} ({100*cnt/n13:.1f}%), mod3={val%3}")

print()
print("Does the 13-channel pass through phantom staircase elements?")
# Check if 13-channel orbits go through {47,91,103,121,...}
_r.seed(789)
phantom_in_13ch = Counter()
n_13ch = 0
for _ in range(5000):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1
    if T < 2: continue
    if orbit[T-2] != 13: continue
    n_13ch += 1
    for v in orbit:
        if v in ALL_PHANTOM:
            phantom_in_13ch[v] += 1

print(f"(Based on {n_13ch} 13-channel orbits)")
print(f"Phantom elements visited in 13-channel orbits:")
for v, cnt in sorted(phantom_in_13ch.items(), key=lambda x: -x[1])[:15]:
    ph_N = [N for N, elems in PHANTOMS.items() if v in elems]
    lvl = ','.join(str(N) for N in sorted(ph_N))
    print(f"  n={v} (phantom N={lvl}): {cnt}/{n_13ch} = {100*cnt/n_13ch:.1f}%")
