"""
85_exact_h2_distribution.py
============================
Exact P(h=2) for each booster: fraction of departures that first miss BSet
(h>1) then land in BSet on the second macro-step.

Also: joint distribution P(h=1, dest=r') and P(h=2, dest=r'') for each source.

Key question: is P(h=2) also near-uniform across boosters, like P(h=1)?
Is P(h=1) + P(h=2) also approximately 2 * 15/128?

For exact computation:
  h=1 case: output = (3^k0 * m - 1)/2^l, count m-values where output%256 in BSet
  h=2 case: for each m where output%256 NOT in BSet, apply one more macro-step
            and check if THAT output%256 in BSet

The period for the h=2 distribution depends on the period of (n mod 512) after
the first macro-step. Since we use 256 odd m values (one period for output mod 256),
the h=2 computation requires tracing each non-BSet output one more step.
"""
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def v2(x):
    c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def macro_step_val(val_in):
    """Given 3^k0 * m - 1 (already computed), do ONE macro-step.
    Input: an odd positive integer n.
    Output: (n', k, l) where n' is the next odd value.
    """
    k = v2(val_in + 1)
    m2 = (val_in + 1) >> k
    x = m2 * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

def macro_step(n):
    k = v2(n + 1)
    m = (n + 1) >> k
    x = m * (3**k) - 1
    l = v2(x)
    return x >> l, k, l

BSet = {27, 55, 63, 83, 95, 103, 127, 159, 169, 191, 207, 223, 239, 253, 255}
BList = sorted(BSet)
M_PERIOD = 256  # 256 odd m values = one full period for output mod 256

print("=== EXACT P(h=1) AND P(h=2) FOR ALL BOOSTERS ===\n")
print("Using 256 odd m-values (one output period mod 256).\n")

results = {}

for r in BList:
    k0 = v2(r + 1)
    pow3k = 3**k0

    h1_hits = Counter()   # h=1 destination -> count
    h2_hits = Counter()   # h=2 destination -> count
    h_gt2 = 0             # count of m-values with h > 2

    for m_idx in range(M_PERIOD):
        m = 2 * m_idx + 1   # m = 1, 3, ..., 511

        # First macro-step from booster r
        val1 = pow3k * m - 1
        l1 = v2(val1)
        out1 = val1 >> l1
        r1 = out1 % 256

        if r1 in BSet:
            h1_hits[r1] += 1
        elif out1 <= 1:
            pass  # converged, not counted
        else:
            # Second macro-step
            k2 = v2(out1 + 1)
            m2 = (out1 + 1) >> k2
            val2 = m2 * (3**k2) - 1
            l2 = v2(val2)
            out2 = val2 >> l2
            r2 = out2 % 256

            if r2 in BSet:
                h2_hits[r2] += 1
            else:
                h_gt2 += 1

    h1_total = sum(h1_hits.values())
    h2_total = sum(h2_hits.values())
    p_h1 = h1_total / M_PERIOD
    p_h2 = h2_total / M_PERIOD
    p_h12 = (h1_total + h2_total) / M_PERIOD

    results[r] = {
        'k0': k0, 'h1': h1_hits, 'h2': h2_hits,
        'p_h1': p_h1, 'p_h2': p_h2, 'h_gt2': h_gt2
    }

    print(f"r={r:3d} k={k0}: P(h=1)={h1_total}/{M_PERIOD}={100*p_h1:.3f}%  "
          f"P(h=2)={h2_total}/{M_PERIOD}={100*p_h2:.3f}%  "
          f"P(h<=2)={100*p_h12:.3f}%  "
          f"h>2:{h_gt2}")

print()
print("=== THEORETICAL PREDICTIONS ===")
print(f"P(h=1) theory: 15/128 = {100*15/128:.3f}%")
# For P(h=2): the non-BSet fraction is ~113/128. Of those outputs (113 residues),
# each then takes one more macro-step, landing uniformly over 128 odd residues.
# P(h=2) ≈ (113/128) × (15/128) ≈ 10.33%
# P(h<=2) ≈ 15/128 + (113/128)*(15/128) = (15/128)*(1 + 113/128)
theory_h1 = 15/128
theory_h2 = (113/128) * (15/128)
theory_h12 = theory_h1 + theory_h2
print(f"P(h=2) theory: (113/128)*(15/128) = {100*theory_h2:.3f}%")
print(f"P(h<=2) theory: {100*theory_h12:.3f}%")
print(f"P(h<=T) theory (geometric): 1 - (113/128)^T * (1-15/128)")

print()
print("=== SUMMARY: SORTED BY P(h=2) ===")
sorted_r = sorted(BList, key=lambda r: -results[r]['p_h2'])
for r in sorted_r:
    d = results[r]
    print(f"  r={r:3d} k={d['k0']}: P(h=1)={100*d['p_h1']:.3f}%  P(h=2)={100*d['p_h2']:.3f}%  "
          f"P(h<=2)={100*(d['p_h1']+d['p_h2']):.3f}%")

print()
print("=== JOINT DISTRIBUTION: P(h=2, dest=r'') FROM EACH BOOSTER ===")
print("(only top destinations shown)\n")
for r in BList:
    d = results[r]
    k0 = d['k0']
    h2 = d['h2']
    if sum(h2.values()) > 0:
        top = sorted(h2.items(), key=lambda x: -x[1])[:5]
        top_str = ", ".join(f"r''={r2}(k={v2(r2+1)}):{c}" for r2,c in top)
        k_avg = sum(v2(r2+1)*c for r2,c in h2.items()) / sum(h2.values())
        print(f"r={r:3d}(k={k0}): h2 dests: {top_str} | avg k_dest2={k_avg:.2f}")

print()
print("=== COMPARISON: P(h=1) vs P(h=2) DESTINATION k-AVERAGES ===")
print(f"{'r':>5} {'k':>4}  {'k_avg_h1':>10}  {'k_avg_h2':>10}  {'diff':>8}")
for r in BList:
    d = results[r]
    k0 = d['k0']
    h1 = d['h1']
    h2 = d['h2']
    if sum(h1.values()) > 0:
        ka1 = sum(v2(r2+1)*c for r2,c in h1.items()) / sum(h1.values())
    else:
        ka1 = float('nan')
    if sum(h2.values()) > 0:
        ka2 = sum(v2(r2+1)*c for r2,c in h2.items()) / sum(h2.values())
    else:
        ka2 = float('nan')
    diff = ka2 - ka1 if ka1 == ka1 and ka2 == ka2 else float('nan')
    print(f"r={r:3d} k={k0}:  k_avg_h1={ka1:.3f}   k_avg_h2={ka2:.3f}   diff={diff:+.3f}")

print()
print("=== CUMULATIVE P(h<=T) FROM EXACT FIRST TWO STEPS ===")
print("Using exact h=1 and h=2 counts, plus geometric model for h>2.\n")

# For each booster, compute P(h<=1), P(h<=2), and geometric extrapolation
# P(h=3) ≈ P(h>2) * P(h=1) (geometric approximation for h>2)
for r in BList:
    d = results[r]
    ph1 = d['p_h1']
    ph2 = d['p_h2']
    ph_gt2 = d['h_gt2'] / M_PERIOD
    # Geometric model: P(h=j) ≈ P(h>j-1) * r_BSet for j >= 2
    # r_BSet ≈ 15/128 (probability of hitting BSet from any residue)
    # P(h=2) exact; P(h=j) ≈ ph_gt2 * (113/128)^(j-2) * (15/128) for j>=3
    # E[h] = 1*ph1 + 2*ph2 + sum_j>=3 j * ph_gt2 * (113/128)^(j-2) * (15/128)
    # E[h | h>=3] ≈ 3 + (113/128)/(15/128) = 3 + 113/15 = 3 + 7.53 = 10.53
    q = 113/128  # probability of missing BSet
    p = 15/128   # probability of hitting BSet
    e_h_given_gte3 = 3 + q/p  # expected h conditional on h>=3 (geometric)
    e_h = 1*ph1 + 2*ph2 + ph_gt2 * e_h_given_gte3
    print(f"r={r:3d}: P(h=1)={100*ph1:.2f}%  P(h=2)={100*ph2:.2f}%  "
          f"P(h>=3)={100*ph_gt2:.2f}%  E[h]_approx={e_h:.3f}")
