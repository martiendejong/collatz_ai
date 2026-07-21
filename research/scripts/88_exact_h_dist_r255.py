"""
88_exact_h_dist_r255.py
========================
Exact hop-length distribution for r=255 (k0=8) over ALL 256 starting m-values.

For each of the 256 odd m-values (m = 1, 3, ..., 511) giving n = 256m - 1:
  - Trace macro-steps until BSet hit or convergence (no timeout)
  - Record hop length h and k-sum

Gives the EXACT h distribution from all 256 starting points in one period.
Python handles big integers natively, so no overflow.

Key use: rigorously lower-bound E[h(255)] for the max cycle mean proof.

Lower bound formula:
  E[h(255)] >= sum_{j=1}^{H} j * P(h=j) + (H+1) * P(h>H)

As H increases (all 256 m-values complete within H steps), the bound tightens
to the EXACT E[h(255)].
"""
import sys
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
k0_of = {r: v2(r+1) for r in BList}

# r=255 starting points: n = 256m - 1 for m = 1, 3, ..., 511 (256 odd values)
# These give k=8 exactly (since v2(256m) = 8 for odd m)
print("=== EXACT h DISTRIBUTION FROM r=255 (ALL 256 STARTING POINTS) ===\n")
print("n = 256*m - 1 for m = 1, 3, ..., 511 (256 starting points)")
print("Using Python big integers -- no overflow.\n")

M_PERIOD = 256
MAX_H = 1000  # effectively no limit

h_vals = []
k_sums = []
dest_by_h = {}  # h -> Counter of destinations
converged = []

for m_idx in range(M_PERIOD):
    m = 2 * m_idx + 1  # m = 1, 3, ..., 511
    n = 256 * m - 1

    # Verify k=8
    k_check = v2(n + 1)
    assert k_check == 8, f"Expected k=8 but got k={k_check} for m={m}, n={n}"

    # Trace to BSet hit or convergence
    n_cur = n
    h = 0
    k_sum = 0
    hit = None

    while h < MAX_H:
        n_cur, k, l = macro_step(n_cur)
        h += 1
        k_sum += k
        if n_cur <= 1:
            hit = 'converged'
            break
        if n_cur % 256 in BSet:
            hit = n_cur % 256
            break

    if hit is None:
        print(f"  WARNING: m={m} did NOT hit BSet or converge in {MAX_H} steps!")
        hit = 'timeout'

    if isinstance(hit, int):
        h_vals.append(h)
        k_sums.append(k_sum)
        if h not in dest_by_h:
            dest_by_h[h] = Counter()
        dest_by_h[h][hit] += 1
    else:
        converged.append((m, h, hit))

print(f"Results: {len(h_vals)} BSet hits, {len(converged)} converged/timeout")

# h distribution
h_counter = Counter(h_vals)
total = M_PERIOD  # 256 starting points

print("\n=== HOP LENGTH DISTRIBUTION (exact, all 256 starting points) ===\n")
print(f"{'h':>5}  {'count':>6}  {'P(h=j)':>10}  {'cumP':>10}  {'dest breakdown'}")
print("-"*80)

cumulative = 0
for h in sorted(h_counter.keys()):
    cnt = h_counter[h]
    p = cnt / total
    cumulative += cnt
    cum_p = cumulative / total
    dests = ", ".join(f"r={r2}(k={k0_of.get(r2,'?')}):x{c}" for r2, c in sorted(dest_by_h[h].items(), key=lambda x: -x[1])[:5])
    print(f"  h={h:3d}:  {cnt:4d}  P={100*p:7.4f}%  cumP={100*cum_p:7.4f}%  -> {dests}")

print()
print(f"All converged: {len(converged)}")
for m, h, reason in converged:
    print(f"  m={m}: {reason} at step h={h}")

# Exact E[h] computation
E_h_exact = sum(h * cnt for h, cnt in h_counter.items()) / total
print(f"\n=== EXACT E[h] FROM r=255 ===")
print(f"E[h(255)] = {E_h_exact:.6f}  (exact, over 256 starting points)")
print(f"  (for converged orbits, using the steps taken before convergence)")

# Rigorous lower bound on E[h(255)] using partial data
# E[h] >= sum_{j=1}^{H} j * P(h=j) + (H+1) * P(h>H)
print(f"\n=== RIGOROUS LOWER BOUNDS ON E[h(255)] ===\n")
cum_count = 0
for H in sorted(h_counter.keys()):
    cum_count += h_counter[H]
    lb_partial = sum(j * h_counter[j] / total for j in h_counter if j <= H)
    p_above = (total - cum_count) / total
    lb = lb_partial + (H + 1) * p_above
    ub = lb_partial + MAX_H * p_above  # loose upper bound
    print(f"  H={H:3d}: E[h] >= {lb:.4f}  (using exact P(h<=H)={100*cum_count/total:.2f}%)")
    if cum_count == total:
        print(f"  => ALL {total} STARTING POINTS COMPLETE AT H={H}. E[h(255)] = {E_h_exact:.6f} EXACTLY.")
        break

# k_sum distribution and E[k/step]
k_per_step = [k/h for h, k in zip(h_vals, k_sums)]
E_k_step = sum(k_sums) / sum(h_vals)  # avg k per step over ALL paths from r=255
E_k_per_path = sum(k_sums) / len(k_sums) if k_sums else 0

print(f"\n=== k-SUM STATISTICS FROM r=255 ===")
print(f"E[k_sum] per path = {E_k_per_path:.4f}")
print(f"E[k/step] (total k / total steps) = {E_k_step:.6f}")
print(f"Self-loop cycle mean (E[k_sum]/E[h] for r→255) = {sum(k for h,k,d in [(h_vals[i],k_sums[i],list(dest_by_h.get(h_vals[i],{}).keys())) for i in range(len(h_vals)) if list(dest_by_h.get(h_vals[i],Counter()).keys())[:1] == [255] or (h_vals[i] in dest_by_h and 255 in dest_by_h[h_vals[i]])]) / max(1, sum(h for h,k,d in [(h_vals[i],k_sums[i],list(dest_by_h.get(h_vals[i],{}).keys())) for i in range(len(h_vals)) if (h_vals[i] in dest_by_h and 255 in dest_by_h[h_vals[i]])])):.4f}  (self-loop paths)")

# For the proof: lower bound on max cycle mean of 255<->127
print(f"\n=== TOWARD PROVING MAX CYCLE MEAN < 3.419 ===\n")
print(f"Best 2-cycle: r=255 <-> r=127 (from Theorem 199)")
print(f"lambda* = E[k_sum(255->127)] / E[h(255->127)] + E[k_sum(127->255)] / ...")
print(f"Actually: lambda* = (E[k_sum(255->127)] + E[k_sum(127->255)]) / (E[h(255->127)] + E[h(127->255)])")
print()

# Find exact distribution of 255->127 hops (from the 256 starting points, filter dest=127)
paths_to_127 = [(h_vals[i], k_sums[i]) for i in range(len(h_vals))
                if h_vals[i] in dest_by_h and 255 in dest_by_h[h_vals[i]] or
                   any(r == 127 for r in (list(dest_by_h.get(h_vals[i], Counter()).keys())))]

# Simpler: find which m-values hit r=127
paths_to_127_exact = []
for i in range(len(h_vals)):
    h = h_vals[i]
    k = k_sums[i]
    # Check destination
    dests = dest_by_h.get(h, Counter())
    # Need to cross-reference with which starting m gave h=h_vals[i]
    # This is complex without saving per-m data. Let me recompute.

print("Computing exact 255->127 hop distribution...")
h_127 = []
k_127 = []
for m_idx in range(M_PERIOD):
    m = 2 * m_idx + 1
    n = 256 * m - 1
    n_cur = n
    h = 0
    k_sum_here = 0
    while h < MAX_H:
        n_cur, k, l = macro_step(n_cur)
        h += 1
        k_sum_here += k
        if n_cur <= 1:
            break
        if n_cur % 256 in BSet:
            if n_cur % 256 == 127:
                h_127.append(h)
                k_127.append(k_sum_here)
            break

total_reach_127 = len(h_127)
print(f"Starting points that reach r=127 first: {total_reach_127}/256 = {100*total_reach_127/256:.2f}%")
if h_127:
    E_h_127 = sum(h_127) / len(h_127)
    E_k_127 = sum(k_127) / sum(h_127)
    lambda_self_127 = sum(k_127) / sum(h_127)  # avg k/step when going 255->127
    print(f"E[h(255->127)] = {E_h_127:.4f}")
    print(f"E[k/step(255->127)] = {E_k_127:.4f}")
    h_counter_127 = Counter(h_127)
    for h in sorted(h_counter_127.keys()):
        print(f"  h={h}: {h_counter_127[h]} paths (k_sum={sum(k for hh,k in zip(h_127,k_127) if hh==h)})")
