"""
147_other_channel.py
=====================
Investigate the "other-channel" (T-2 != 13 or 23): the 9.2% of large Collatz
orbits that reach n=1 via a different penultimate path.

From Obs 286: T-1=5 for 94% of orbits (via n=5->1).
The remaining 6% go through:
  n=85 (2.24%), n=151 (1.96%), n=227 (1.42%), n=341 (0.22%), etc.

Key questions:
1. What are the T-2 values in the other-channel, and do they form organized sub-channels?
2. Do other-channel orbits have ANY phantom structure (N=7-10)?
3. What is the T-k structure of the other-channel for k=3 to 25?
4. Can we algebraically characterize the predecessors of 1 and their sub-channel structure?
"""
import random as _r
from collections import Counter
import math

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

def collatz_orbit(n0, max_steps=200000):
    orbit = [n0]
    n = n0
    for _ in range(max_steps):
        if n == 1: break
        n, _, _ = macro_step(n)
        orbit.append(n)
    return orbit

PHANTOMS = {
    7:  {47, 91, 103, 121},
    8:  {71, 91, 103, 121, 175, 189},
    9:  {91, 95, 103, 167, 175, 253, 283, 319, 399, 445},
    10: {703, 937},
}
ALL_PHANTOM = set().union(*PHANTOMS.values())

print("=" * 70)
print("PART 1: STRUCTURE OF THE OTHER-CHANNEL (T-2 != 13 or 23)")
print("=" * 70)
print()
print("Predecessors of n=1: values p where macro_step(p) = 1")
print("These are: p = (m * 3^K - 1) / 2^l0 with m*3^K = 2^l0 + 1 (i.e., m*3^K = 2^l0 + 1)")
print()

# Find all p with macro_step(p) = 1 (predecessors of 1), up to 10^6
preds_of_1 = []
for K in range(1, 25):
    for l0 in range(1, 80):
        val = (1 << l0) + 1  # 2^l0 + 1, needed for (m*3^K - 1)/2^l0 = 1 means m*3^K = 2^l0 + 1
        if val % (3**K) != 0:
            continue
        m = val // (3**K)
        if m <= 0 or m % 2 == 0:
            continue
        n = m * (1 << K) - 1
        if n > 10**6:
            continue
        # Verify
        n_out, K_act, l0_act = macro_step(n)
        if n_out == 1:
            preds_of_1.append((n, K, l0))

preds_of_1.sort()
print(f"Predecessors of 1 up to 10^6 ({len(preds_of_1)} found):")
for n, K, l0 in preds_of_1[:15]:
    print(f"  n={n:>10} (K={K}, l0={l0}), n mod 3 = {n%3}")
    # Verify T-2 value in real orbit from large numbers
    # T-2 of n = n itself (since n -> 1 directly)

print()
print("Pattern: predecessors of 1 are n where m*3^K = 2^{l0} + 1.")
print("The smallest predecessor with each K:")
for K in range(1, 8):
    preds_K = [(n, l0) for n, K_, l0 in preds_of_1 if K_==K]
    if preds_K:
        n_min, l0_min = preds_K[0]
        print(f"  K={K}: smallest n = {n_min} (l0={l0_min})")

print()
print("=" * 70)
print("PART 2: OTHER-CHANNEL T-k DISTRIBUTION")
print("=" * 70)
print()

_r.seed(42)
b = 500
n_trials = 20000

# Classify orbits by T-2 value
t2_counter = Counter()
other_channel_tk = {k: Counter() for k in range(1, 30)}
channel_23_size = 0
channel_13_size = 0
other_size = 0

other_t2_subchannels = {}  # t2_val -> {k: Counter}

for _ in range(n_trials):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1
    if T < 2: continue

    t2_val = orbit[T - 2]
    t2_counter[t2_val] += 1

    if t2_val == 23:
        channel_23_size += 1
    elif t2_val == 13:
        channel_13_size += 1
    else:
        other_size += 1
        for k in range(1, 30):
            if T >= k:
                other_channel_tk[k][orbit[T-k]] += 1
        # Track by subchannel (T-2 value)
        if t2_val not in other_t2_subchannels:
            other_t2_subchannels[t2_val] = {k: Counter() for k in range(1, 30)}
        for k in range(1, 30):
            if T >= k:
                other_t2_subchannels[t2_val][k][orbit[T-k]] += 1

total = sum(t2_counter.values())
print(f"T-2 distribution ({n_trials} trials):")
for t2_val, cnt in sorted(t2_counter.items(), key=lambda x: -x[1])[:15]:
    n_out, K, l0 = macro_step(t2_val)
    print(f"  T-2={t2_val:>8}: {cnt:>6} ({100*cnt/total:.2f}%), macro_step->{n_out} (K={K},l0={l0})")
print(f"  Total 23-ch: {channel_23_size} ({100*channel_23_size/total:.1f}%)")
print(f"  Total 13-ch: {channel_13_size} ({100*channel_13_size/total:.1f}%)")
print(f"  Total other: {other_size} ({100*other_size/total:.1f}%)")

print()
print("=" * 70)
print("PART 3: OTHER-CHANNEL T-k DOMINANT VALUES (phantoms?)")
print("=" * 70)
print()
print(f"{'T-k':>6} {'Dominant n':>14} {'Freq%':>8} {'Phantom?':>12}")
print("-" * 50)

for k in range(1, 26):
    c = other_channel_tk[k]
    total_k = sum(c.values())
    if not total_k: continue
    top_val, top_cnt = c.most_common(1)[0]
    pct = 100 * top_cnt / total_k
    ph_N = [N for N, elems in PHANTOMS.items() if top_val in elems]
    ph_str = ','.join(str(N) for N in sorted(ph_N)) if ph_N else "---"
    print(f"{k:>6} {top_val:>14} {pct:>8.2f}% {ph_str:>12}")

print()
print("=" * 70)
print("PART 4: PHANTOM CONTENT OF OTHER-CHANNEL ORBITS")
print("=" * 70)
print()
print("Do other-channel orbits pass through any phantom elements?")
print()

_r.seed(123)
b = 500
n_trials2 = 5000
other_phantom_counts = Counter()
other_total = 0

for _ in range(n_trials2):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1
    if T < 2: continue
    t2_val = orbit[T-2]
    if t2_val == 13 or t2_val == 23: continue
    other_total += 1
    for v in orbit:
        if v in ALL_PHANTOM:
            other_phantom_counts[v] += 1

print(f"Based on {other_total} other-channel orbits:")
print()
if other_phantom_counts:
    print(f"{'Element':>8} {'Count':>8} {'%':>8} {'Phantom N':>12}")
    for v, cnt in sorted(other_phantom_counts.items(), key=lambda x: -x[1])[:15]:
        ph_N = [N for N, elems in PHANTOMS.items() if v in elems]
        lvl = ','.join(str(N) for N in sorted(ph_N))
        print(f"{v:>8} {cnt:>8} {100*cnt/other_total:>8.2f}% {lvl:>12}")
else:
    print("No phantom elements found in other-channel orbits!")

print()
print("=" * 70)
print("PART 5: ALGEBRAIC CHARACTERIZATION OF PREDECESSORS OF 1")
print("=" * 70)
print()
print("All T-1 values (predecessors of 1) share the property: m*3^K - 1 = 2^l0.")
print("This means m*3^K = 2^l0 + 1.")
print()
print("Key constraints:")
print("1. m must be odd positive")
print("2. K >= 1 (since K = v2(n+1) >= 1 for odd n)")
print("3. l0 >= 1 (since l0 = v2(m*3^K - 1) >= 1 -- the output IS 1, so any l0 works)")
print()

# By modular arithmetic: m * 3^K = 2^l0 + 1
# mod 3: 0 = 2^l0 + 1 (mod 3) since 3^K = 0 mod 3
# So 2^l0 = -1 = 2 (mod 3), which means l0 is ODD.
print("From mod 3: 2^l0 + 1 = 0 (mod 3) requires 2^l0 = 2 (mod 3), so l0 is ODD.")
print()
print("From mod 9 (for K >= 2): 0 = m*9*3^{K-2} - 1 ... complex analysis")
print()
print("Values with l0 = 1 (x=2^1+1=3, m*3^K=3): K=1, m=1, n=1*2-1=1 (trivial)")
print("Values with l0 = 3 (x=2^3+1=9=3^2): K=2, m=1, n=1*4-1=3 (standard predecessor)")
print("                                      or K=1, m=3, n=3*2-1=5 [n=5: K=1,m=3: confirmed]")

# Double check n=5: K=v2(6)=1, m=3, x=9-1=8=2^3, l0=3, n_out=1. Yes!
print()
print("  n=5: K=1, l0=3, m=3 -> 3*3 - 1 = 8 = 2^3. [confirmed]")
print()

# What are l0=3 predecessors?
# m*3^K = 2^3 + 1 = 9 = 3^2
# K=1: m=3 -> n=5 ✓
# K=2: m=1 -> n=3 (trivial orbit 3->1)
# K=3 or more: 3^K | 9 -> K<=2 only
print("l0=3 (target=9=3^2): K=1 gives n=5; K=2 gives n=3.")
print()

# l0=7 (target=128+1=129=3*43): K=1, m=43, n=85. K>1: 3^K|129? 3|129=3*43, 9|129? 129/9=14.33 no.
print("l0=7 (target=129=3*43): K=1, m=43, n=85 [ok]; K>=2: 9 does not divide 129.")
print()

# l0=9 (target=512+1=513=3^3*19): K=1, m=171, n=341; K=2, m=57, n=227; K=3, m=19, n=151.
print("l0=9 (target=513=3^3*19): K=1->n=341; K=2->n=227; K=3->n=151; K>=4: 81|513? 513/81=6.33 no.")
print()

# This explains the "other-channel" T-2 values beautifully!
# n=5: K=1, l0=3 (most common)
# n=85: K=1, l0=7 (next: 2.24%)
# n=151: K=3, l0=9 (1.96%)
# n=227: K=2, l0=9 (1.42%)
# n=341: K=1, l0=9 (0.22%)

print("Complete T-1 value structure (by l0):")
print()
print(f"  l0 = 3:  2^3+1 = 9 = 3^2.    T-1 values: n=5 (K=1,m=3), n=3 (K=2,m=1 [trivial])")
print(f"  l0 = 5:  2^5+1 = 33 = 3*11.  T-1 values: n=21 (K=1,m=11). [But 21->1 directly?]")

# Check n=21: K=v2(22)=1, m=11, x=33-1=32=2^5, l0=5, n_out=1. Yes!
print(f"           Actually n=21: K=1, l0=5, 21->1 [confirmed] (rarely visited: T-2=21 is rare)")
print(f"  l0 = 7:  2^7+1 = 129 = 3*43. T-1 values: n=85 (K=1,m=43).")
print(f"  l0 = 9:  2^9+1 = 513 = 3^3*19. T-1 values: n=341 (K=1), 227 (K=2), 151 (K=3).")
print(f"  l0 = 11: 2^11+1 = 2049 = 3*683. T-1 values: n=1365 (K=1,m=683).")
print(f"  l0 = 13: 2^13+1 = 8193 = 3*2731. T-1 values: n=5461 (K=1,m=2731).")

print()
print("Pattern: T-1 = n values are {(2^{l0}+1) / 3^K * 2^K - 1} for l0 odd and 3^K | 2^{l0}+1.")
print("The l0 values are the ODD numbers where 3 | 2^{l0}+1.")
print()
print(f"2^l0 + 1 mod 3: {[((1<<l0)+1) % 3 for l0 in range(1,20)]}")
print(f"l0 values where 3 | 2^l0+1: {[l0 for l0 in range(1,30) if ((1<<l0)+1)%3==0]}")
print()
print("These are l0 = 1, 3, 5, 7, 9, 11, ... (all odd l0!)")
print("And within each odd l0, the exact factorization of 2^l0+1 by powers of 3")
print("determines how many T-1 values exist at that l0 level.")

print()
print("PASSAGE RATES: proportion of b-bit orbits with each T-1 value:")
print("  n=5 (l0=3): 94% -- dominant because l0=3 is the smallest odd l0 giving a nontrivial path")
print("  n=21 (l0=5): very rare (not in our T-2 data at >0.1%)")
print("  n=85 (l0=7): 2.24%")
print("  n=151,227,341 (l0=9): 1.96%+1.42%+0.22%=3.60%")
print("  n=1365 (l0=11): very rare")
print()
print("Why l0=5 (n=21) is rare but l0=7 (n=85) is 2.24%?")
print("Because n=21 and n=85 have different predecessor densities.")
print("n=21: small number (4 bits), needs 5-bit output to reach 1.")
print("n=85: 7-bit number, needs 8-bit output (128) to reach 1.")
print()

# Find predecessors of n=21 vs n=85
def predecessors_of(target, max_n=1000000):
    result = []
    for K in range(1, 25):
        for l0 in range(1, 50):
            num = (1<<l0)*target + 1
            denom = 3**K
            if num % denom != 0: continue
            m = num // denom
            if m%2==0: continue
            n = m*(1<<K)-1
            if 0 < n <= max_n:
                n_out,_,_ = macro_step(n)
                if n_out == target:
                    result.append((n,K,l0))
    return sorted(result)

for target in [5, 21, 85, 151, 227, 341]:
    preds = predecessors_of(target, max_n=100000)
    t2_ch = "23" if any(macro_step(p)[0]==23 for p,_,__ in preds) else "?"
    print(f"  n={target}: {len(preds)} predecessors up to 100K: {[p for p,_,__ in preds[:5]]}")
print()
print("n=21 has fewer predecessors than n=85 -> lower passage rate (T-2=21 is very rare).")
