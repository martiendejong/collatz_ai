"""R359-366: CHAMPION GENEALOGY. Hypothesis: laggard records propagate along
backward edges -- next record is often a preimage (2n or (4n-1)/3 or (2n-1)/3)
of a previous record; genuinely NEW seeds are rare. Decompose A006877 (measured
to 2^24) into backward-chains vs seeds; also measure orbit-merge structure:
where does each seed's orbit merge into the previous seed's orbit?"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

N = 1 << 24
d = np.full(N, -1, dtype=np.int32)
d[1] = 0
for n in range(2, N):
    if d[n] >= 0: continue
    path = []; m = n
    while m >= N or d[m] < 0:
        path.append(m); m = m // 2 if m % 2 == 0 else (3 * m + 1) // 2
    base = d[m]
    for x in reversed(path):
        base += 1
        if x < N: d[x] = base

recs = []
best = -1
for n in range(2, N):
    if d[n] > best: best = d[n]; recs.append(n)

# preimage relations
def is_preimage(child, parent):
    # child -> parent in ONE T-step?
    if child == 2 * parent: return "2n"
    if child % 2 == 1 and (3 * child + 1) // 2 == parent: return "up"
    return None

def backward_reach(child, parent, max_steps=6):
    # is child reachable from parent by <= max_steps backward T-steps?
    frontier = {parent}
    for s in range(1, max_steps + 1):
        nxt = set()
        for v in frontier:
            nxt.add(2 * v)
            if (2 * v - 1) % 3 == 0 and ((2 * v - 1) // 3) % 2 == 1:
                nxt.add((2 * v - 1) // 3)
        if child in nxt: return s
        frontier = {x for x in nxt if x < (1 << 26)}
    return None

print("record | relation to previous record | seed?")
seeds = []
prev = None
for n in recs:
    if n < 20:
        prev = n; continue
    rel = None
    if prev:
        s = backward_reach(n, prev, 8)
        if s: rel = f"preimage-chain, {s} backward steps from {prev}"
    if rel is None:
        seeds.append(n)
        rel = "*** NEW SEED ***"
    print(f"  {n:>10} (d={d[n]:>3}): {rel}")
    prev = n

print(f"\nseeds: {seeds}")
print(f"{len(recs)-2} records, {len(seeds)} genuine seeds -> "
      f"{100*(1 - len(seeds)/(len(recs)-2)):.0f}% of records are backward-propagation")

# orbit merge structure of seeds: where does seed_i's orbit merge into seed_{i-1}'s?
def orbit_set(n, cap=100000):
    o = []; m = n
    while m != 1:
        o.append(m); m = m // 2 if m % 2 == 0 else (3 * m + 1) // 2
    return o
print("\nseed orbit merges (first common value with previous seed's orbit, and % of orbit shared):")
prev_orb = None
for s in seeds:
    orb = orbit_set(s)
    if prev_orb:
        ps = set(prev_orb)
        merge = next((i for i, v in enumerate(orb) if v in ps), None)
        if merge is not None:
            shared = len(orb) - merge
            print(f"  seed {s:>10}: merges with previous seed at value {orb[merge]:>8} "
                  f"after {merge} own steps; shares {shared}/{len(orb)} = {shared/len(orb):.0%} of its orbit")
    prev_orb = orb
