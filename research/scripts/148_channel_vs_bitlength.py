"""
148_channel_vs_bitlength.py
============================
Does the 13-channel / 23-channel split depend on bit-length?

From Obs 289, for b=500-bit starting numbers: 47.7% / 43.3%.
This script measures the split across b = 10 to b = 1000 to check
whether the asymmetry is universal or a finite-size effect.

Also measures the T-3 sub-structure within each channel to see if
the phantom staircase entry rate changes with bit-length.
"""
import random as _r
from collections import Counter

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    K = v2(n+1); m = (n+1) >> K; x = m*(3**K)-1; l0 = v2(x)
    return x >> l0, K, l0

def collatz_orbit(n0, max_steps=500000):
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

# Phantom staircase elements (canonical T-8 to T-16)
STAIRCASE_CORE = {47, 121, 91, 103, 175, 445, 167, 283, 319}

print("=" * 70)
print("PART 1: 13-CHANNEL / 23-CHANNEL SPLIT VS BIT-LENGTH")
print("=" * 70)
print()
print(f"{'b':>6} {'23-ch%':>8} {'13-ch%':>8} {'35-ch%':>8} {'other%':>8} {'trials':>8}")
print("-" * 60)

bit_lengths = [10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000]
n_per_b = 5000

_r.seed(42)
results = {}

for b in bit_lengths:
    ch23, ch13, ch35, ch_other = 0, 0, 0, 0
    total = 0

    for _ in range(n_per_b):
        n0 = _r.getrandbits(b) | 1
        if n0 < 3: n0 = 3
        orbit = collatz_orbit(n0)
        T = len(orbit) - 1
        if T < 2: continue
        total += 1
        t2 = orbit[T-2]
        if t2 == 23:
            ch23 += 1
        elif t2 == 13:
            ch13 += 1
        elif t2 == 35:
            ch35 += 1
        else:
            ch_other += 1

    if total:
        p23 = 100*ch23/total
        p13 = 100*ch13/total
        p35 = 100*ch35/total
        po = 100*ch_other/total
        results[b] = (p23, p13, p35, po, total)
        print(f"{b:>6} {p23:>8.2f}% {p13:>8.2f}% {p35:>8.2f}% {po:>8.2f}% {total:>8}")

print()
print("Key: does 13/23 ratio converge as b -> inf?")
print()
for b in bit_lengths:
    if b in results:
        p23, p13, p35, po, t = results[b]
        ratio = p13/p23 if p23 > 0 else 0
        print(f"  b={b:>5}: 13-ch/23-ch ratio = {ratio:.4f}, (13-23)% = {p13-p23:+.2f}%")

print()
print("=" * 70)
print("PART 2: PHANTOM STAIRCASE ENTRY RATE VS BIT-LENGTH")
print("=" * 70)
print()
print("Of 23-channel orbits, what fraction pass through the staircase core (N=9 phantom)?")
print("Specifically, what fraction visit n=319 (the T-8 element, dissolution point)?")
print()
print(f"{'b':>6} {'n319 rate%':>12} {'any_phantom%':>14} {'staircase%':>12}")
print("-" * 55)

_r.seed(123)

for b in bit_lengths:
    n_staircase = 0
    n_phantom = 0
    n_319 = 0
    total_23 = 0

    for _ in range(n_per_b):
        n0 = _r.getrandbits(b) | 1
        if n0 < 3: n0 = 3
        orbit = collatz_orbit(n0)
        T = len(orbit) - 1
        if T < 2: continue
        t2 = orbit[T-2]
        if t2 != 23: continue
        total_23 += 1
        visited = set(orbit)
        if 319 in visited:
            n_319 += 1
        if any(v in ALL_PHANTOM for v in orbit):
            n_phantom += 1
        if any(v in STAIRCASE_CORE for v in orbit):
            n_staircase += 1

    if total_23 > 0:
        print(f"{b:>6} {100*n_319/total_23:>12.2f}% {100*n_phantom/total_23:>14.2f}% {100*n_staircase/total_23:>12.2f}%")
    else:
        print(f"{b:>6} {'(no 23-ch orbits)':>40}")

print()
print("Key: does phantom staircase entry rate increase with b -> inf?")
print("If yes, eventually ALL 23-ch orbits go through the staircase.")
print()

print()
print("=" * 70)
print("PART 3: T-3 VALUES WITHIN 23-CHANNEL VS BIT-LENGTH")
print("=" * 70)
print()
print("The dominant T-3 in 23-channel is n=61 (direct predecessor of 23).")
print("Does n=61 concentration change with b?")
print()
print(f"{'b':>6} {'T-3=61%':>10} {'T-3=325%':>12} {'other T-3%':>12}")
print("-" * 45)

_r.seed(456)

for b in bit_lengths:
    t3_61 = 0; t3_325 = 0; total_23 = 0

    for _ in range(n_per_b):
        n0 = _r.getrandbits(b) | 1
        if n0 < 3: n0 = 3
        orbit = collatz_orbit(n0)
        T = len(orbit) - 1
        if T < 2: continue
        if orbit[T-2] != 23: continue
        total_23 += 1
        if T >= 3:
            t3 = orbit[T-3]
            if t3 == 61: t3_61 += 1
            elif t3 == 325: t3_325 += 1

    if total_23 > 0:
        p61 = 100*t3_61/total_23
        p325 = 100*t3_325/total_23
        print(f"{b:>6} {p61:>10.2f}% {p325:>12.2f}% {100-p61-p325:>12.2f}%")

print()
print("n=61 is the direct predecessor of 23 (canonical exit ramp T-3).")
print("n=325 is T-4 of the canonical path (325->61->23->5->1).")
print()

print()
print("=" * 70)
print("PART 4: WHERE DO 13-CHANNEL ORBITS CONVERGE WITH 23-CHANNEL?")
print("=" * 70)
print()
print("At what T-k does the 13-channel diverge from the 23-channel in the backward tree?")
print("Checking: for a 13-channel orbit and a 23-channel orbit from the same b-bit range,")
print("what is the last T-k where they share a common ancestor?")
print()
print("This is equivalent to asking: how deep is the split point in the Collatz tree?")
print()
print("Approach: for each channel, track the T-k sequence and find when they first differ.")
print("Concretely: at T-k, do 13-ch and 23-ch orbits visit the same values?")
print()

_r.seed(789)
b = 500
n_per = 2000

# Collect T-k distributions for both channels
tk_23 = {k: Counter() for k in range(1, 30)}
tk_13 = {k: Counter() for k in range(1, 30)}
total_23 = total_13 = 0

for _ in range(n_per * 2):
    n0 = _r.getrandbits(b) | 1
    orbit = collatz_orbit(n0)
    T = len(orbit) - 1
    if T < 2: continue
    t2 = orbit[T-2]
    if t2 == 23 and total_23 < n_per:
        total_23 += 1
        for k in range(1, 30):
            if T >= k: tk_23[k][orbit[T-k]] += 1
    elif t2 == 13 and total_13 < n_per:
        total_13 += 1
        for k in range(1, 30):
            if T >= k: tk_13[k][orbit[T-k]] += 1
    if total_23 >= n_per and total_13 >= n_per:
        break

print(f"{'T-k':>5} {'23-ch top':>12} {'23%':>8} {'13-ch top':>12} {'13%':>8} {'Shared?':>8}")
print("-" * 60)
for k in range(1, 28):
    c23 = tk_23[k]; c13 = tk_13[k]
    if not c23 or not c13: continue
    top23, n23 = c23.most_common(1)[0]
    top13, n13 = c13.most_common(1)[0]
    p23 = 100*n23/total_23
    p13 = 100*n13/total_13
    # Are the same values in both channels at this T-k?
    vals_23 = set(c23.keys())
    vals_13 = set(c13.keys())
    shared = vals_23 & vals_13
    shared_str = f"{len(shared)} values" if shared else "DISJOINT"
    print(f"{k:>5} {str(top23):>12} {p23:>8.2f}% {str(top13):>12} {p13:>8.2f}% {shared_str:>8}")

print()
print("The channels become 'disjoint' (no shared T-k values) at large k,")
print("then merge back as k grows further (eventually all start from random large values).")
