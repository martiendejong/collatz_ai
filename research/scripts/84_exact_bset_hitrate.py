"""
84_exact_bset_hitrate.py
=========================
Exact P(h=1) for each booster: fraction of departures that land in BSet in 1 step.

The output residue mod 256 from booster r (with k0=v2(r+1)) depends on m mod 256
(or 512 for the full v2-analysis). The exact P(h=1) requires 256 odd m values
(one full period of the output distribution mod 256 in the high-precision sense).

For r with k0 = v2(r+1):
  n+1 = 2^k0 * m, m odd.
  output = (3^k0 * m - 1) / 2^{v2(3^k0 * m - 1)}
  P(h=1) = #{m in [1,3,...,511]: output mod 256 in BSet} / 256

Also: show which BSet destinations are reachable from each booster in 1 step.
"""
import sys
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n+1)
    m = (n+1) >> k
    x = m * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
BList = sorted(BSet)

print("=== EXACT P(h=1) FOR ALL BOOSTERS (256 odd-m period) ===\n")

# For each booster r, use 256 odd m values: m = 1, 3, ..., 511
# This gives the full period of output mod 256
M_PERIOD = 256  # 256 odd m values = one full output period

results = {}

for r in BList:
    k0 = v2(r + 1)
    pow3k = 3**k0

    bset_hits = Counter()  # destination -> count
    non_bset_outputs = Counter()  # non-BSet output mod 256 -> count
    all_outputs = Counter()

    for m_idx in range(M_PERIOD):
        m = 2 * m_idx + 1   # m = 1, 3, ..., 511

        # Compute macro-step output for n = 2^k0 * m - 1
        val = pow3k * m - 1
        l = v2(val)
        output = val >> l
        r_out = output % 256

        all_outputs[r_out] += 1
        if r_out in BSet:
            bset_hits[r_out] += 1
        else:
            non_bset_outputs[r_out] += 1

    bset_total = sum(bset_hits.values())
    p_h1 = bset_total / M_PERIOD
    results[r] = {'k0': k0, 'p_h1': p_h1, 'hits': bset_hits, 'non_bset': non_bset_outputs}

    print(f"r={r:3d} k={k0}: P(h=1)={bset_total}/{M_PERIOD} = {100*p_h1:.3f}%")
    if bset_hits:
        dests = ", ".join(f"r'={r2}(k={v2(r2+1)}):{cnt}" for r2, cnt in sorted(bset_hits.items()))
        print(f"  Destinations: {dests}")
    print()

print()
print("=== SUMMARY TABLE ===\n")
print(f"{'r':>5} {'k':>4} {'P(h=1)':>10}  {'Destinations (BSet r prime, count)':}")
print("-"*70)
for r in BList:
    d = results[r]
    k0 = d['k0']
    p = d['p_h1']
    hits = d['hits']
    total_hits = sum(hits.values())
    dest_str = " ".join(f"{r2}x{c}" for r2, c in sorted(hits.items(), key=lambda x: -x[1]))
    print(f"r={r:3d} k={k0}: {total_hits:3d}/{M_PERIOD} = {100*p:.3f}%  -> {dest_str}")

print()
print("=== SORTED BY P(h=1) ===")
sorted_r = sorted(BList, key=lambda r: -results[r]['p_h1'])
for r in sorted_r:
    d = results[r]
    print(f"  r={r:3d} k={d['k0']}: P(h=1)={100*d['p_h1']:.3f}%")

print()
print("=== THEORETICAL IMPLICATION ===")
print("For D_hard_kern: need ~40% consecutive-booster rate from r=255.")
print("Actual P(h=1) values from each booster:")
max_p = max(d['p_h1'] for d in results.values())
max_r = max(BList, key=lambda r: results[r]['p_h1'])
print(f"  Maximum P(h=1): {100*max_p:.3f}% (from r={max_r})")
print(f"  Required: ~40%")
print(f"  Gap factor: {0.40/max_p:.2f}x (max P(h=1) is {0.40/max_p:.2f}x too small)")

# Bonus: for each booster, list ALL BSet elements reachable in 1 step and their k values
print()
print("=== BOOSTER-TO-BOOSTER 1-STEP REACHABILITY ===")
for r in BList:
    d = results[r]
    k0 = d['k0']
    if d['hits']:
        k_weighted = sum(v2(r2+1)*c for r2,c in d['hits'].items()) / sum(d['hits'].values())
        print(f"r={r:3d}(k={k0}) -> {sorted(d['hits'].keys())} : k_avg_dest={k_weighted:.2f}")
    else:
        print(f"r={r:3d}(k={k0}) -> [] (no direct BSet hits)")
