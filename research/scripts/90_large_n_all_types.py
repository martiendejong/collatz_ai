"""
90_large_n_all_types.py
========================
Large-n h distribution and transition statistics for ALL 8 booster types.
Using one representative booster per k0 type, with m ~ 10^12.

Gives the complete large-n 8-state transition matrix with:
- E[h(k0 -> k0')] per edge
- E[k_sum(k0 -> k0')] per edge
- max cycle mean from this matrix

Also: confirms the Mersenne pattern in BSet (r=63, 127, 255 are 2^k-1).
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict

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
k0_of = {r: v2(r+1) for r in BList}

# k0 -> representative booster
k0_groups = defaultdict(list)
for r in BList:
    k0_groups[k0_of[r]].append(r)
K0_TYPES = sorted(k0_groups.keys())

rep = {k0: k0_groups[k0][0] for k0 in K0_TYPES}

# Large-n starting points
M_BASE = 10**12
M_PERIOD = 256
MAX_H = 500

print("=== LARGE-n TRANSITION MATRIX FOR ALL 8 BOOSTER TYPES ===\n")
print(f"m in [M_BASE+1, M_BASE+511] (256 odd values per type)")
print(f"M_BASE = {M_BASE} (n ~ 2.56e14)\n")

t0 = time.time()

# edge_data[k0_src][k0_dst] = (count, total_steps, total_k)
edge_data = {k0: defaultdict(lambda: [0, 0, 0]) for k0 in K0_TYPES}
# Also track per-booster destination (not just k0 type)
dest_by_type = {k0: Counter() for k0 in K0_TYPES}

type_stats = {}  # k0 -> {'ph1': float, 'avg_h': float, 'avg_k': float, ...}

for k0 in K0_TYPES:
    r = rep[k0]
    step = 2**(k0 + 1)  # spacing between n values with k=k0 exactly

    # n = step * m - 1 for odd m, gives k=k0 exactly
    # For m ~ M_BASE/step (large), n = step * m - 1 ~ M_BASE

    h_vals = []
    k_sums = []
    bset_hits = 0
    converged_count = 0

    for m_idx in range(M_PERIOD):
        m = M_BASE + 2 * m_idx + 1  # 256 consecutive odd m-values from M_BASE+1
        n = step * m - 1

        # Verify k
        k_check = v2(n + 1)
        if k_check != k0:
            # This can happen for very small m, but shouldn't for m >> 1
            # Try a different base
            n = (n // step + 1) * step + (r % step)
            if n < M_BASE:
                n += step * ((-n + M_BASE) // step + 1)
            k_check = v2(n + 1)

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
            r_out = n_cur % 256
            if r_out in BSet:
                hit = r_out
                break

        if isinstance(hit, int):
            h_vals.append(h)
            k_sums.append(k_sum)
            k0_dst = k0_of[hit]
            edge_data[k0][k0_dst][0] += 1
            edge_data[k0][k0_dst][1] += h
            edge_data[k0][k0_dst][2] += k_sum
            dest_by_type[k0][hit] += 1
            bset_hits += 1

    n_hits = len(h_vals)
    if n_hits > 0:
        ph1 = Counter(h_vals).get(1, 0) / M_PERIOD
        avg_h = sum(h_vals) / n_hits
        avg_k = sum(k_sums) / sum(h_vals)
        type_stats[k0] = {'ph1': ph1, 'avg_h': avg_h, 'avg_k': avg_k, 'n_hits': n_hits,
                           'n_converged': converged_count}
        print(f"k0={k0} (r={r}): n_hits={n_hits}/256  avg_h={avg_h:.3f}  avg_k/step={avg_k:.4f}  "
              f"P(h=1)={100*ph1:.3f}%  converged={converged_count}")

print(f"\n[Done in {time.time()-t0:.1f}s]\n")

# =====================================================================
# 8-STATE TRANSITION MATRIX FROM LARGE-n DATA
# =====================================================================
print("\n=== 8-STATE k0 TRANSITION MATRIX (E[h] per edge) ===\n")
print(f"Rows = source k0, Cols = dest k0")
print(f"Value = E[h | k0_src -> k0_dst] (avg steps for paths from k0_src ending at k0_dst)\n")

print("       " + "  ".join(f"k0={k0}" for k0 in K0_TYPES))
for k0_src in K0_TYPES:
    row = []
    for k0_dst in K0_TYPES:
        cnt, tot_s, tot_k = edge_data[k0_src][k0_dst]
        if cnt >= 1:
            row.append(f"{tot_s/cnt:6.2f}")
        else:
            row.append(f"  --  ")
    print(f"k0={k0_src}:  " + "  ".join(row))

print("\n=== 8-STATE k0 TRANSITION MATRIX (E[k/step] per edge) ===\n")
print(f"Value = E[k/step | k0_src -> k0_dst]\n")
print("       " + "  ".join(f"k0={k0}" for k0 in K0_TYPES))
for k0_src in K0_TYPES:
    row = []
    for k0_dst in K0_TYPES:
        cnt, tot_s, tot_k = edge_data[k0_src][k0_dst]
        if cnt >= 1 and tot_s > 0:
            row.append(f"{tot_k/tot_s:6.4f}")
        else:
            row.append(f"  --   ")
    print(f"k0={k0_src}:  " + "  ".join(row))

# =====================================================================
# MAX CYCLE MEAN (KARP'S ALGORITHM)
# =====================================================================
print("\n=== MAX CYCLE MEAN FROM LARGE-n TRANSITION DATA ===\n")

def has_positive_cycle(lam, K0_TYPES, edge_data):
    # Bellman-Ford: check if max-weight cycle exists
    dist = {k0: 0.0 for k0 in K0_TYPES}
    for _ in range(len(K0_TYPES) + 1):
        updated = False
        new_dist = dict(dist)
        for k0_src in K0_TYPES:
            for k0_dst in K0_TYPES:
                cnt, tot_s, tot_k = edge_data[k0_src][k0_dst]
                if cnt < 1: continue
                k_avg = tot_k / cnt
                s_avg = tot_s / cnt
                w = k_avg - lam * s_avg
                cand = dist[k0_src] + w
                if cand > new_dist.get(k0_dst, -1e18):
                    new_dist[k0_dst] = cand
                    updated = True
        dist = new_dist
        if not updated: break
    # Check for positive-weight cycle
    for k0_src in K0_TYPES:
        for k0_dst in K0_TYPES:
            cnt, tot_s, tot_k = edge_data[k0_src][k0_dst]
            if cnt < 1: continue
            k_avg = tot_k / cnt
            s_avg = tot_s / cnt
            w = k_avg - lam * s_avg
            if dist[k0_src] + w > dist.get(k0_dst, -1e18) + 1e-9:
                return True
    return False

lo, hi = 0.5, 8.0
for _ in range(80):
    mid = (lo + hi) / 2
    if has_positive_cycle(mid, K0_TYPES, edge_data):
        lo = mid
    else:
        hi = mid

max_cycle_mean_large_n = (lo + hi) / 2
DK_THRESHOLD = 3.419
print(f"Max cycle mean (large-n data): lambda* = {max_cycle_mean_large_n:.6f}")
print(f"D_hard_kern threshold: {DK_THRESHOLD}")
print(f"Gap: {DK_THRESHOLD - max_cycle_mean_large_n:.6f}")
print(f"{'lambda* < threshold (D_hard_kern excluded)' if max_cycle_mean_large_n < DK_THRESHOLD else 'CAUTION: lambda* >= threshold!'}")

# Best 2-cycles
print("\n=== BEST 2-CYCLES FROM LARGE-n DATA ===\n")
best_2cyc = []
for k0_a in K0_TYPES:
    for k0_b in K0_TYPES:
        c1, s1, k1 = edge_data[k0_a][k0_b]
        c2, s2, k2 = edge_data[k0_b][k0_a]
        if c1 >= 1 and c2 >= 1 and s1 > 0 and s2 > 0:
            lam = (k1/c1 + k2/c2) / (s1/c1 + s2/c2)
            h1 = s1/c1
            h2 = s2/c2
            best_2cyc.append((lam, k0_a, k0_b, c1, c2, h1, h2))

best_2cyc.sort(reverse=True)
for lam, k0_a, k0_b, c1, c2, h1, h2 in best_2cyc[:10]:
    print(f"  k0={k0_a}<->k0={k0_b}: lambda={lam:.4f}  h_a={h1:.2f}  h_b={h2:.2f}  n={c1}/{c2}")

# =====================================================================
# MERSENNE PATTERN IN BSET
# =====================================================================
print("\n=== MERSENNE PATTERN IN BSET ===\n")
print("BSet elements of the form 2^k - 1 (Mersenne numbers):")
for r in BList:
    k0 = k0_of[r]
    if r + 1 == 2**k0:  # r = 2^k0 - 1 (pure Mersenne)
        print(f"  r={r:3d} = 2^{k0}-1 (k0={k0})  IN BSet ✓  [MERSENNE]")
    elif (r+1) == 2**(k0+1) - 1 or r == 2**k0 - 1:
        pass
    else:
        r1 = r+1
        print(f"  r={r:3d}: r+1={r1} = 2^{k0} * {r1>>k0}  (k0={k0})")

print()
print("All 2^k-1 values in 1..255 and BSet membership:")
for k in range(1, 9):
    r = 2**k - 1
    k0 = v2(r+1)
    in_bset = r in BSet
    print(f"  2^{k}-1 = {r:3d}: k0={k0}  {'IN BSet ✓' if in_bset else 'NOT in BSet'}")
