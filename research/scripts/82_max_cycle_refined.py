"""
82_max_cycle_refined.py
========================
Refined max-cycle-mean computation with N=5000 samples per booster.
Also: full per-edge conditional stats (not just marginals).
"""
import math, sys, time
from collections import defaultdict

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
DK_THRESHOLD = 3.419

def trace_to_next_booster(n_start, max_steps=5000):
    n = n_start
    steps = 0; k_sum = 0
    while True:
        n, k, l = macro_step(n)
        steps += 1; k_sum += k
        if n <= 1: return None, steps, k_sum
        if n % 256 in BSet: return n % 256, steps, k_sum
        if steps >= max_steps: return None, steps, k_sum

BASE = 512 * 2_000_000
N_PER_BOOSTER = 5000

print(f"N={N_PER_BOOSTER} per booster, BASE={BASE:,}")
t0 = time.time()

# edge_data[r][r2] = [count, total_steps, total_k]
edge_data = {r: defaultdict(lambda: [0, 0, 0]) for r in BList}

for r in BList:
    k0 = v2(r + 1)
    step = 2**(k0 + 1)
    n_base = (BASE // step) * step + (r % step)
    if n_base < BASE: n_base += step
    _, k_check, _ = macro_step(n_base)
    assert k_check == k0, f"k mismatch r={r}"

    for j in range(N_PER_BOOSTER):
        n = n_base + j * step
        r2, steps, k_sum = trace_to_next_booster(n)
        if r2 is not None:
            ed = edge_data[r][r2]
            ed[0] += 1; ed[1] += steps; ed[2] += k_sum

print(f"Simulation done [{time.time()-t0:.0f}s]")

# Per-booster summary
print("\n=== PER-BOOSTER UNCONDITIONAL STATS ===")
for r in BList:
    k0 = v2(r + 1)
    tc = sum(v[0] for v in edge_data[r].values())
    ts = sum(v[1] for v in edge_data[r].values())
    tk = sum(v[2] for v in edge_data[r].values())
    if tc > 0:
        print(f"r={r:3d} k={k0}: avg_steps={ts/tc:.3f}  avg_k/step={tk/ts:.4f}  reach={tc}/{N_PER_BOOSTER}")

# Self-loop analysis
print("\n=== SELF-LOOP ANALYSIS (r -> r) ===")
for r in BList:
    if r in edge_data[r]:
        cnt, tot_s, tot_k = edge_data[r][r]
        if cnt > 0:
            ak = tot_k / tot_s
            print(f"r={r:3d}: self-loop count={cnt}/{N_PER_BOOSTER} ({100*cnt/N_PER_BOOSTER:.1f}%)  "
                  f"avg_h={tot_s/cnt:.2f}  avg_k/step={ak:.4f}  cycle_lambda={ak:.4f}")

# Max cycle mean via binary search
def has_positive_cycle(lam, BList, edge_data):
    INF = float('inf')
    dist = {v: 0.0 for v in BList}
    for _ in range(len(BList) + 1):
        updated = False
        new_dist = dict(dist)
        for r in BList:
            for r2, (cnt, tot_s, tot_k) in edge_data[r].items():
                if cnt < 2: continue
                k_avg = tot_k / cnt
                s_avg = tot_s / cnt
                w = k_avg - lam * s_avg
                cand = dist[r] + w
                if r2 in new_dist and cand > new_dist[r2]:
                    new_dist[r2] = cand
                    updated = True
        dist = new_dist
        if not updated: break
    for r in BList:
        for r2, (cnt, tot_s, tot_k) in edge_data[r].items():
            if cnt < 2: continue
            k_avg = tot_k / cnt
            s_avg = tot_s / cnt
            w = k_avg - lam * s_avg
            if r2 in dist and dist[r] + w > dist[r2] + 1e-10:
                return True
    return False

lo, hi = 0.5, 8.0
for _ in range(80):
    mid = (lo + hi) / 2
    if has_positive_cycle(mid, BList, edge_data):
        lo = mid
    else:
        hi = mid

max_cycle_mean = (lo + hi) / 2
print(f"\n=== MAX CYCLE MEAN (N={N_PER_BOOSTER}) ===")
print(f"lambda* = {max_cycle_mean:.6f}")
print(f"D_hard_kern threshold = {DK_THRESHOLD:.4f}")
print(f"Gap = {DK_THRESHOLD - max_cycle_mean:.6f}")
print(f"Max cycle mean {'<' if max_cycle_mean < DK_THRESHOLD else '>='} threshold")

# Best 2-cycles
print("\n=== BEST 2-CYCLES (r->r'->r) ===")
best_2cyc = []
for r in BList:
    for r2 in BList:
        if r2 in edge_data[r] and r in edge_data[r2]:
            c1, s1, k1 = edge_data[r][r2]; c2, s2, k2 = edge_data[r2][r]
            if c1 >= 2 and c2 >= 2:
                lam = (k1/c1 + k2/c2) / (s1/c1 + s2/c2)
                best_2cyc.append((lam, r, r2, c1, c2, s1/c1, s2/c2))

best_2cyc.sort(reverse=True)
for lam, r, r2, c1, c2, h1, h2 in best_2cyc[:10]:
    print(f"  {r}->{r2}->{r}: lambda={lam:.4f}  h1={h1:.1f}  h2={h2:.1f}  "
          f"n1={c1}  n2={c2}")

print(f"\n[Total {time.time()-t0:.0f}s]")
