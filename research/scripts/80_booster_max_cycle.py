"""
80_booster_max_cycle.py
========================
Compute the booster-to-booster transition matrix and the maximum cycle mean
(max achievable avg k in any infinite orbit that follows booster chains).

Key question: can any infinite orbit achieve avg k >= 3.419 (D_hard_kern threshold)?
If max_cycle_mean < 3.419, then D_hard_kern = empty.
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
LOG23 = math.log2(3) - 1  # 0.58496
DK_THRESHOLD = 3.419

def trace_to_next_booster(n_start, max_steps=2000):
    """From n_start (a booster-class number), walk until hitting the next booster.
    Returns (next_r, steps, k_sum) or (None, steps, k_sum) if converged/timed out."""
    n = n_start
    steps = 0
    k_sum = 0
    while True:
        n, k, l = macro_step(n)
        steps += 1
        k_sum += k
        if n <= 1:
            return None, steps, k_sum
        r = n % 256
        if r in BSet:
            return r, steps, k_sum
        if steps >= max_steps:
            return None, steps, k_sum

BASE = 256 * 2_000_000   # large enough to avoid small-n artifacts
N_PER_BOOSTER = 1000     # samples per booster

print(f"=== BOOSTER MAX CYCLE ANALYSIS ===")
print(f"Samples per booster: {N_PER_BOOSTER}, BASE={BASE:,}\n")

t0 = time.time()

# edge_data[r][r2] = (count, total_steps, total_k)
edge_data = {r: defaultdict(lambda: [0, 0, 0]) for r in BList}
all_hop_avg_k = []   # avg k per step for EVERY individual hop
max_hop_k = 0.0
max_hop_info = None

for r in BList:
    k0 = v2(r + 1)
    step = 2**(k0 + 1)   # stride; n ≡ r mod step ensures exactly k=k0

    # Align to r mod step at or above BASE
    n_base = (BASE // step) * step + (r % step)
    if n_base < BASE:
        n_base += step

    # Verify k
    _, k_check, _ = macro_step(n_base)
    if k_check != k0:
        # Try with offset
        n_base += step
        _, k_check, _ = macro_step(n_base)
        assert k_check == k0, f"Cannot align r={r}: got k={k_check}"

    converged = 0
    timeout = 0

    for j in range(N_PER_BOOSTER):
        n = n_base + j * step
        r2, steps, k_sum = trace_to_next_booster(n)
        if r2 is not None:
            ed = edge_data[r][r2]
            ed[0] += 1
            ed[1] += steps
            ed[2] += k_sum
            hop_avg_k = k_sum / steps
            all_hop_avg_k.append(hop_avg_k)
            if hop_avg_k > max_hop_k:
                max_hop_k = hop_avg_k
                max_hop_info = (r, r2, steps, k_sum, n)
        elif n <= 1:
            converged += 1
        else:
            timeout += 1

    total_cnt = sum(v[0] for v in edge_data[r].values())
    total_steps = sum(v[1] for v in edge_data[r].values())
    total_k = sum(v[2] for v in edge_data[r].values())

    if total_cnt > 0:
        avg_steps = total_steps / total_cnt
        avg_k_step = total_k / total_steps
        print(f"r={r:3d} k={k0}: avg_steps={avg_steps:.2f}  avg_k/step={avg_k_step:.4f}  "
              f"reach_rate={total_cnt}/{N_PER_BOOSTER}  conv={converged}  timeout={timeout}")

print(f"\n[elapsed {time.time()-t0:.0f}s]")

# ==========================================
# MAX CYCLE MEAN (Karp's algorithm variant)
# ==========================================
# Build graph with edge weights = (k_sum_avg, steps_avg)
# Want: max over all cycles of (sum k) / (sum steps)

# Build edge weight matrix
# W[r][r2] = (mean_k_sum, mean_steps) for edge r->r2
W = {}
for r in BList:
    W[r] = {}
    for r2, (cnt, tot_steps, tot_k) in edge_data[r].items():
        if cnt >= 3:  # require at least 3 observations
            W[r][r2] = (tot_k / cnt, tot_steps / cnt)

# Karp's algorithm for max cycle mean:
# For each source s, compute V[k][v] = max total k-sum over paths of length exactly k steps
# reaching v from s.
# max_cycle_mean = max_v max_k (V[N][v] - V[k][v]) / (N - k)
# where N = number of nodes.

# Since graph is sparse, use policy iteration (simpler for small graphs):
# Value iteration to find max avg k cycle.

nodes = [r for r in BList if W.get(r)]
N = len(nodes)
node_idx = {r: i for i, r in enumerate(nodes)}
INF = float('inf')

# V[k][i] = best total k-sum over paths of exactly k hops from any start, ending at node i
# We also need total steps for the ratio.
# Use (k_sum, steps) pairs and maximize k_sum / steps.

# Parametric approach: binary search on lambda*
# For given lam, check if max cycle (k_sum - lam*steps) > 0.
# This is positive iff there's a cycle with avg_k > lam.

def has_positive_cycle(lam, nodes, W, node_idx):
    """Bellman-Ford: check if there's a positive-weight cycle
    using edge weight = k_sum - lam * steps."""
    n = len(nodes)
    dist = {v: 0.0 for v in nodes}  # start at 0 (reachable from any node)

    for _ in range(n + 1):
        updated = False
        new_dist = dict(dist)
        for r in nodes:
            if r not in W:
                continue
            for r2, (k_avg, s_avg) in W[r].items():
                if r2 in dist:
                    w = k_avg - lam * s_avg
                    cand = dist[r] + w
                    if cand > new_dist.get(r2, -INF):
                        new_dist[r2] = cand
                        updated = True
        dist = new_dist
        if not updated:
            break

    # Check for positive cycle: run one more relaxation
    for r in nodes:
        if r not in W:
            continue
        for r2, (k_avg, s_avg) in W[r].items():
            if r2 in dist:
                w = k_avg - lam * s_avg
                if dist[r] + w > dist.get(r2, -INF) + 1e-10:
                    return True
    return False

# Binary search for max cycle mean
lo, hi = 0.5, 10.0
for _ in range(60):
    mid = (lo + hi) / 2
    if has_positive_cycle(mid, nodes, W, node_idx):
        lo = mid
    else:
        hi = mid

max_cycle_mean = (lo + hi) / 2

print(f"\n=== MAX CYCLE MEAN ===")
print(f"Max achievable avg k in any booster chain: {max_cycle_mean:.4f}")
print(f"D_hard_kern threshold:                     {DK_THRESHOLD:.4f}")
print(f"Gap to threshold:                          {DK_THRESHOLD - max_cycle_mean:.4f}")
if max_cycle_mean < DK_THRESHOLD:
    print(f"  --> MAX CYCLE MEAN < THRESHOLD: D_hard_kern CANNOT be maintained via booster chains!")
else:
    print(f"  --> WARNING: max cycle mean meets threshold")

print(f"\n=== MAX SINGLE HOP avg k ===")
print(f"Maximum avg k/step over any single booster hop: {max_hop_k:.4f}")
if max_hop_info:
    r, r2, steps, k_sum, n = max_hop_info
    print(f"  From r={r} -> r'={r2}, steps={steps}, k_sum={k_sum}, start n={n}")

# Distribution of hop avg k values
import statistics
if all_hop_avg_k:
    print(f"\n=== HOP avg k DISTRIBUTION (N={len(all_hop_avg_k)} hops) ===")
    print(f"  Mean:    {statistics.mean(all_hop_avg_k):.4f}")
    print(f"  Median:  {statistics.median(all_hop_avg_k):.4f}")
    print(f"  Stdev:   {statistics.stdev(all_hop_avg_k):.4f}")
    sorted_k = sorted(all_hop_avg_k, reverse=True)
    print(f"  Top 5:   {[round(x,3) for x in sorted_k[:5]]}")
    above_dk = sum(1 for x in all_hop_avg_k if x >= DK_THRESHOLD)
    print(f"  Hops with avg_k >= {DK_THRESHOLD}: {above_dk}/{len(all_hop_avg_k)}")

# Best cycles (enumerate short cycles explicitly)
print(f"\n=== SHORT CYCLE AVG K (explicit enumeration) ===")
best_cycle = None
best_cycle_ak = 0.0

for r_start in nodes:
    # Self-loops (1-step cycles)
    if r_start in W.get(r_start, {}):
        k_sum, steps = W[r_start][r_start]
        ak = k_sum / steps
        if ak > best_cycle_ak:
            best_cycle_ak = ak
            best_cycle = [r_start, r_start]
        print(f"  Cycle {r_start}->{r_start}: avg_k={ak:.4f} ({steps:.1f} steps)")

    # 2-cycles
    for r2 in W.get(r_start, {}):
        if r_start in W.get(r2, {}):
            k1, s1 = W[r_start][r2]
            k2, s2 = W[r2][r_start]
            ak = (k1 + k2) / (s1 + s2)
            if ak > best_cycle_ak:
                best_cycle_ak = ak
                best_cycle = [r_start, r2, r_start]
            if ak > max_cycle_mean - 0.1:  # only print near-best
                print(f"  Cycle {r_start}->{r2}->{r_start}: avg_k={ak:.4f} ({s1+s2:.1f} steps)")

print(f"\nBest explicit cycle: {best_cycle} with avg_k={best_cycle_ak:.4f}")
print(f"\n[Total elapsed: {time.time()-t0:.0f}s]")
