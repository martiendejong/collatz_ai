"""
83_exact_h_distribution.py
===========================
Exact one-period h-distribution for each booster r.

For each r in BSet with k0 = v2(r+1):
  - One period = 2^k0 starting m values (all odd m in [1, 2^(k0+1)-1])
  - For each m, n = 2^k0 * m - 1 gives exactly k=k0
  - Trace until BSet hit or convergence, record h

This gives an exact (non-probabilistic) sample of all arithmetic orbits
within one period of the starting structure. The distribution is periodic
in m and represents the large-n regime exactly.
"""
import math, sys, time
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
DK_THRESHOLD = 3.419

def one_period_orbits(r, max_h=200):
    """
    Enumerate all orbits from one full period of starting values at r.
    Returns list of (h, k_sum) for each orbit, where h = steps to next BSet hit.
    """
    k0 = v2(r + 1)
    period = 2**k0   # number of odd m values in one period
    step = 2**(k0 + 1)  # n spacing for exactly k=k0

    # m values: 1, 3, 5, ..., 2^(k0+1) - 1 = period odd values
    # n = 2^k0 * m - 1 for m odd in [1, step-1]
    results = []
    k_prefix = k0  # the macro-step AT the booster r has k=k0 (or larger — see below)

    for m_idx in range(period):
        m = 2 * m_idx + 1   # m = 1, 3, 5, ..., 2*period-1
        n = (2**k0) * m - 1
        assert n % 2 == 1, f"n should be odd: n={n}, m={m}, k0={k0}"
        # Verify n % (step) == r % step
        # n + 1 = 2^k0 * m, and v2(n+1) = k0 since m is odd. Good.

        # Trace to next BSet hit
        n_cur = n
        h = 0
        k_sum = 0
        hit = None
        while h < max_h:
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
            hit = 'timeout'
        results.append((h, k_sum, hit, m))

    return results, k0, period

# Main analysis
print("=== EXACT ONE-PERIOD h-DISTRIBUTION FOR ALL BOOSTERS ===\n")
t0 = time.time()

summary = {}

for r in sorted(BSet):
    results, k0, period = one_period_orbits(r)

    # Filter: only hits (not converged/timeout)
    hits = [(h, k_sum, dest, m) for h, k_sum, dest, m in results if isinstance(dest, int)]
    timeouts = sum(1 for _, _, d, _ in results if d == 'timeout')
    converged = sum(1 for _, _, d, _ in results if d == 'converged')

    h_values = [h for h, _, _, _ in hits]
    k_values = [k_sum for _, k_sum, _, _ in hits]

    h_counter = Counter(h_values)
    total = period  # total starting points in one period

    # Cumulative distribution
    cumulative = {}
    running = 0
    for h in range(1, max(h_values)+1 if h_values else 1):
        running += h_counter.get(h, 0)
        cumulative[h] = running

    if hits:
        avg_h = sum(h_values) / len(hits)
        avg_k_step = sum(k_values) / sum(h_values)
        max_h = max(h_values)
        min_h = min(h_values)
    else:
        avg_h = float('nan')
        avg_k_step = float('nan')
        max_h = 0
        min_h = 0

    summary[r] = {
        'k0': k0, 'period': period, 'hits': len(hits), 'converged': converged,
        'timeout': timeouts, 'avg_h': avg_h, 'avg_k_step': avg_k_step,
        'h_counter': h_counter, 'max_h': max_h, 'min_h': min_h
    }

    p1 = h_counter.get(1, 0) / total
    p2 = h_counter.get(2, 0) / total
    p3 = h_counter.get(3, 0) / total
    p_leq3 = (h_counter.get(1,0)+h_counter.get(2,0)+h_counter.get(3,0)) / total

    print(f"r={r:3d} k={k0} period={period:3d}: avg_h={avg_h:.3f}  avg_k/s={avg_k_step:.4f}  "
          f"P(h=1)={p1:.4f}  P(h=2)={p2:.4f}  P(h=3)={p3:.4f}  P(h<=3)={p_leq3:.4f}")

    # Detailed h distribution
    max_show = min(15, max_h+1)
    dist_str = "  h: " + "  ".join(f"{h}:{h_counter.get(h,0):3d}/{total}"
                                    for h in range(1, max_show+1) if h_counter.get(h,0)>0)
    print(dist_str)

    if converged > 0 or timeouts > 0:
        print(f"  [converged={converged}, timeout={timeouts}]")
    print()

print("=== E[h] LOWER BOUNDS (from one-period exact distribution) ===\n")
print("Lower bound: E[h] >= sum_j j * P(h=j) for j=1..H, plus (H+1)*(1 - P(h<=H))")
print()

for r in sorted(BSet):
    s = summary[r]
    total = s['period']
    hc = s['h_counter']
    k0 = s['k0']

    max_h_obs = max(hc.keys()) if hc else 0

    # Sum j * P(h=j) for all observed j
    lower_bound = sum(j * hc[j] / total for j in hc)
    # Add (max_h+1) * P(h > max_h) as lower bound contribution
    p_above = (total - sum(hc.values())) / total  # converged/timeout fraction
    lower_bound += (max_h_obs + 1) * p_above

    print(f"r={r:3d} k={k0}: E[h] >= {lower_bound:.4f}  (from period={total} exact orbits)")

print(f"\n[Elapsed: {time.time()-t0:.2f}s]")

# === MAX CYCLE MEAN FROM EXACT PERIOD DATA ===
print("\n=== EXACT MAX CYCLE MEAN FROM ONE-PERIOD DATA ===\n")
print("Using avg_h and avg_k/step from exact one-period computation:")
print()

# For each r, self-loop cycle mean = avg_k/step (since cycle_lambda = total_k/total_steps)
best_cycle = None
best_lam = 0.0

for r in sorted(BSet):
    s = summary[r]
    if s['hits'] > 0:
        lam = s['avg_k_step']  # this IS the cycle lambda for self-loop
        print(f"r={r:3d}: self-loop lambda = {lam:.4f}  avg_h = {s['avg_h']:.3f}")
        if lam > best_lam:
            best_lam = lam
            best_cycle = r

print(f"\nBest self-loop: r={best_cycle} with lambda={best_lam:.4f}")
print(f"D_hard_kern threshold: {DK_THRESHOLD}")
print(f"Gap: {DK_THRESHOLD - best_lam:.4f}")

# Also check 2-cycles using exact h data... harder without cross-booster data
# But we can note: the exact one-period avg_h gives the same result as the large-n simulation
