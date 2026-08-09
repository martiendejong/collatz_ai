# 356: Martien's intuition "the fractional 1/n logically never returns to n"
# made exact, in two layers.
# LAYER 1 (the subtlety): along the orbit the fractions TELESCOPE perfectly:
#   prod (3 n_i + 1) = 2^N prod n_i  is automatic for any cyclic orbit — each
#   n_i cancels between consecutive steps. So "denominators can't line up" is
#   NOT the obstruction. The whole constraint sits at the endpoint closure:
#   for a v-pattern (v_1..v_K), the unique rational fixed point is
#     n_1 = c / (2^N - 3^K),  c = sum_{i=0}^{K-1} 3^(K-1-i) * 2^(v_1+..+v_i)
#   and a cycle exists iff that rational is an odd integer (sign per side).
# LAYER 2 (the quantification): each composition gives one candidate rational;
#   heuristic P(integer) ~ 1/|2^N - 3^K|. Expected number of K-cycles:
#     E_K = sum_N C(N-1, K-1) / (2^N - 3^K)   (positive side 2^N > 3^K)
#   Calibration: E_1 should be ~1 (the trivial cycle exists!).
from math import comb, log2
from fractions import Fraction
from itertools import product as iproduct

def fixed_point(vs):
    K = len(vs)
    N = sum(vs)
    c = 0
    pref = 0
    for i in range(K):
        c += 3**(K - 1 - i) * 2**pref
        pref += vs[i]
    den = 2**N - 3**K
    return Fraction(c, den) if den != 0 else None

# LAYER 1 demo: recover ALL four known cycles from the endpoint formula
print("endpoint fixed points for small v-patterns (K<=7, N<=11):")
found = []
for K in range(1, 8):
    for N in range(1, 12):
        for vs in iproduct(range(1, N + 1), repeat=K):
            if sum(vs) != N:
                continue
            fp = fixed_point(vs)
            if fp is not None and fp.denominator == 1 and fp.numerator % 2 != 0 and fp != 0:
                found.append((K, N, vs, int(fp)))
seen = set()
for K, N, vs, n in found:
    if n in seen:
        continue
    seen.add(n)
    print(f"  K={K}, N={N}, pattern {vs}: n = {n}")
print()

# LAYER 2: expected number of positive K-cycles
print("expected number of positive K-cycles (heuristic E_K), positive side:")
tot = 0.0
tail = 0.0
for K in range(1, 41):
    EK = 0.0
    Nmin = int(K * log2(3)) + 1
    for N in range(Nmin, Nmin + 40):
        den = 2**N - 3**K
        EK += comb(N - 1, K - 1) / den
    tot += EK
    if K >= 2:
        tail += EK
    if K <= 12 or K % 5 == 0:
        print(f"  K={K:>2}: E_K = {EK:.6f}")
print(f"  sum K=1..40:  {tot:.4f}   (K=1 term ~1 = the trivial cycle, correctly predicted)")
print(f"  sum K=2..40 (nontrivial): {tail:.4f}")
# decay rate: E_K ~ 2^(-(1-H(1/log2 3)) N) -> per-K factor
import math
H = lambda p: -p * math.log2(p) - (1 - p) * math.log2(1 - p)
rate = (1 - H(1 / log2(3))) * log2(3)
print(f"  asymptotic decay: E_K ~ 2^(-{rate:.4f} K) per extra odd step")
