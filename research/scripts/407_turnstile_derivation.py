"""Obs 616: turnstile weights from first principles (S5 candidates 5+6).

Measure, over all odd starts n < N:
  (1) the GATEWAY distribution: the alternator A_y = (4^y-1)/3 that an
      orbit visits as its last odd value before 1;
  (2) hit probabilities P(orbit visits v) for a small target set
      (alternator ladder + rungs of gate 5), to test the visit-law
      prediction P(hit v) ~ c/v.

First-principles model (geometric backward mint + visit measure 1/v):
  d(A_y) proportional to 1/A_y over non-spring alternators
  (springs y = 0 mod 3 are divisible by 3: no odd predecessors, their
  asymptotic traffic is exactly 0). Prediction:
  d(5) = (1/5) / (1/5 + 1/85 + 1/341 + 1/5461 + 1/21845 + ...) = 93.06%.
Rung version for gate 5 (direct predecessors 3, 13, 53, 213, 853, ...;
springs 3 and 213 blocked): share via rung r proportional to 1/r.
"""
import sys
sys.setrecursionlimit(10000)

N = 4_000_001          # odd starts n < N
TARGETS = [5, 13, 53, 853, 3413, 85, 341, 5461]
TBIT = {v: 1 << i for i, v in enumerate(TARGETS)}

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

# memo[n] = (gateway, hitmask) for odd n < N
memo = {1: (1, 0)}

def process(n0):
    path = []
    n = n0
    while True:
        if n < N and n in memo:
            g, hm = memo[n]
            break
        path.append(n)
        m = 3 * n + 1
        nxt = m >> v2(m)
        if nxt == 1:
            g, hm = n, TBIT.get(n, 0)   # n is the gateway
            memo[n] = (g, hm) if n < N else (g, hm)
            path.pop()
            break
        n = nxt
    for v in reversed(path):
        hm |= TBIT.get(v, 0)
        if v < N:
            memo[v] = (g, hm)
    return memo.get(n0, (g, hm))

from collections import Counter
gate_count = Counter()
hit_count = Counter()
total = 0
for n in range(3, N, 2):
    g, hm = process(n)
    gate_count[g] += 1
    for v, b in TBIT.items():
        if hm & b or n == v:
            hit_count[v] += 1
    total += 1

print(f"odd starts: {total}")
print("\n(1) GATEWAY distribution (top 8):")
alt = [(4 ** y - 1) // 3 for y in range(2, 12)]
for g, c in gate_count.most_common(8):
    y = next((i + 2 for i, a in enumerate(alt) if a == g), "?")
    print(f"  gate {g:>8} (y={y}): {c / total * 100:7.3f}%")
pred = {}
weights = [(y, (4 ** y - 1) // 3) for y in range(2, 14) if y % 3 != 0]
Z = sum(1 / a for _, a in weights)
print("  1/v-model prediction vs measured:")
for y, a in weights[:5]:
    meas = gate_count.get(a, 0) / total * 100
    print(f"    gate {a:>6}: predicted {1 / a / Z * 100:7.3f}%  measured {meas:7.3f}%")

print("\n(2) hit probabilities, testing P(hit v) ~ c/v  (show v * P):")
for v in TARGETS:
    p = hit_count[v] / total
    print(f"  v={v:>6}: P(hit)={p:8.5f}   v*P = {v * p:7.4f}")

print("\n(3) rungs of gate 5 (share of 5-traffic per direct predecessor):")
t5 = hit_count[5]
for r in (13, 53, 853, 3413):
    print(f"  via {r:>5}: {hit_count[r] / t5 * 100:6.2f}% of 5-traffic "
          f"(1/v-model: {(1 / r) / (1 / 13 + 1 / 53 + 1 / 853 + 1 / 3413) * 100:6.2f}%)")
