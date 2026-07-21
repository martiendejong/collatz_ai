"""
89_large_n_h_dist_r255.py
==========================
Large-n h distribution for r=255.

For LARGE m values (m ~ 10^12), compute the exact h distribution over 256 odd-m values.
This avoids the convergence artifact of the small-n computation (script 88).

Key comparison:
  Small-n (script 88): m=1..511, 31.6% convergence, E[h|BSet hit]=5.83
  Large-n (script 82, empirical N=5000): avg_h = 9.2
  Large-n (this script): ???

If spectral gap = 0.99 correctly predicts geometric h distribution in large-n:
  E[h] ≈ 1/P_stat = 1/0.109 = 9.17
  P(h=j) ≈ 0.109 * (1-0.109)^(j-1) for j >= 2
  P(h=1) = 31/256 (exact, arithmetic)

Also: get the exact large-n 255->127 transition stats.
"""
import sys, time
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step(n):
    k = v2(n + 1)
    m = (n + 1) >> k
    x = m * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
BList = sorted(BSet)

# Large-n starting points: n = 256 * m - 1 for m = M_BASE + 1, M_BASE + 3, ..., M_BASE + 511
# (256 consecutive odd m-values starting from a large base)
M_BASE = 10**12  # 1 trillion -- ensures n >> 10^14, negligible convergence

M_PERIOD = 256  # 256 odd m-values per period
MAX_H = 500  # max hop steps before giving up

print(f"=== LARGE-n h DISTRIBUTION FOR r=255 ===")
print(f"m in [{M_BASE+1}, {M_BASE+511}] (256 odd values)")
print(f"n = 256m - 1 in [{256*(M_BASE+1)-1}, {256*(M_BASE+511)-1}]")
print(f"(n ~ 2.56 * 10^14, negligible convergence expected)\n")

t0 = time.time()

h_vals = []
k_sums = []
dests = []
converged_count = 0
timeout_count = 0

for m_idx in range(M_PERIOD):
    m = M_BASE + 2 * m_idx + 1  # m = M_BASE+1, M_BASE+3, ..., M_BASE+511
    n = 256 * m - 1

    # Verify k=8 (since v2(256m) = 8 for odd m)
    k_check = v2(n + 1)
    if k_check != 8:
        print(f"WARNING: k={k_check} for m={m}, expected 8")
        continue

    # Trace to BSet hit
    n_cur = n
    h = 0
    k_sum = 0
    hit = None

    while h < MAX_H:
        n_cur, k, l = macro_step(n_cur)
        h += 1
        k_sum += k
        if n_cur <= 1:
            converged_count += 1
            hit = 'converged'
            break
        if n_cur % 256 in BSet:
            hit = n_cur % 256
            break

    if hit is None:
        timeout_count += 1
        hit = 'timeout'

    if isinstance(hit, int):
        h_vals.append(h)
        k_sums.append(k_sum)
        dests.append(hit)

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")
print(f"BSet hits: {len(h_vals)}/256, converged: {converged_count}, timeout: {timeout_count}\n")

# h distribution
h_counter = Counter(h_vals)
total = M_PERIOD  # 256

print("=== HOP LENGTH DISTRIBUTION (large-n) ===\n")
print(f"{'h':>5}  {'count':>6}  {'P(h=j)':>10}  {'cumP':>10}")
print("-"*50)
cum = 0
for h in sorted(h_counter.keys())[:30]:
    cnt = h_counter[h]
    p = cnt / total
    cum += cnt
    print(f"  h={h:3d}:  {cnt:4d}  {100*p:8.4f}%   cumP={100*cum/total:7.4f}%")

if len(h_counter) > 30:
    print(f"  ... ({len(h_counter)-30} more h values)")

# Key statistics
E_h = sum(h_vals) / total
E_h_bset = sum(h_vals) / len(h_vals) if h_vals else 0
E_k_step = sum(k_sums) / sum(h_vals) if h_vals else 0
print(f"\nE[h] (over all 256 starting points) = {E_h:.4f}")
print(f"E[h | BSet hit] = {E_h_bset:.4f}")
print(f"E[k/step] = {E_k_step:.4f}")

# Compare with small-n and geometric model
p_h1 = h_counter.get(1, 0) / total
p_h2 = h_counter.get(2, 0) / total
print(f"\nP(h=1) = {h_counter.get(1,0)}/{total} = {100*p_h1:.4f}%  (small-n exact: 31/256 = 12.11%)")
print(f"P(h=2) = {h_counter.get(2,0)}/{total} = {100*p_h2:.4f}%  (small-n exact: 21/256 = 8.20%)")
print(f"\nGeometric model prediction (p=10.9%):")
import math
p_geom = 0.109
for j in range(1, 8):
    p_pred = p_geom * (1-p_geom)**(j-1)
    actual = h_counter.get(j, 0) / total
    print(f"  P(h={j}) predicted = {100*p_pred:.4f}%  actual = {100*actual:.4f}%")

# Destination distribution
dest_counter = Counter(dests)
print(f"\n=== DESTINATION DISTRIBUTION (large-n, all h) ===\n")
print(f"{'dest r':>8}  {'k0':>4}  {'count':>6}  {'fraction':>10}")
for r, cnt in sorted(dest_counter.items(), key=lambda x: -x[1]):
    k0 = v2(r+1)
    print(f"  r={r:3d}  k={k0}: {cnt:5d}  {100*cnt/len(h_vals):.3f}%")

# 255->127 specific
paths_to_127 = [(h_vals[i], k_sums[i]) for i in range(len(h_vals)) if dests[i] == 127]
print(f"\n=== 255->127 TRANSITION (large-n, {M_PERIOD} starting points) ===")
print(f"Paths reaching r=127 first: {len(paths_to_127)}/{total} = {100*len(paths_to_127)/total:.3f}%")
if paths_to_127:
    E_h_127 = sum(h for h,k in paths_to_127) / len(paths_to_127)
    E_k_127 = sum(k for h,k in paths_to_127) / sum(h for h,k in paths_to_127)
    print(f"E[h(255->127)] = {E_h_127:.4f}")
    print(f"E[k/step(255->127)] = {E_k_127:.4f}")
    h_counter_127 = Counter(h for h,k in paths_to_127)
    print(f"h distribution:")
    for h in sorted(h_counter_127.keys()):
        k_here = sum(k for hh,k in paths_to_127 if hh==h)
        print(f"  h={h}: {h_counter_127[h]} paths  k_total={k_here}  avg_k/step={k_here/h_counter_127[h]/h:.3f}")

# Check: P(h=1) should be exactly 31/256 regardless of n
print(f"\n=== ARITHMETIC VERIFICATION ===")
print(f"P(h=1) = {h_counter.get(1,0)}/256: {'31/256 CONFIRMED' if h_counter.get(1,0) == 31 else 'UNEXPECTED!'}")
print(f"(Expected 31/256 from arithmetic — mod-256 structure is n-invariant)")
