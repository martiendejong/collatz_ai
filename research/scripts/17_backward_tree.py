"""E27/G1: backward-tree minimum growth (Applegate-Lagarias style) — foundation for the
Krasikov-Lagarias 0.84 program.

Backward Collatz: m -> 2m (always), m -> (m-1)/3 (iff m even, m = 1 mod 3; child is odd).
State: (parity p, ternary residue r mod 3^j) with budget j (each /3-step consumes one trit).
L(p, r, j, d) = minimal number of depth-d leaves (worst case over unknown deeper trits).
Nodes r = 0 mod 3 are chains (only doubling) -> L = 1. Budget 0 -> assume no branch (worst case).
Exponent lower bound gamma(d) = log2(min-leaves over admissible roots)/d.
"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8")
sys.setrecursionlimit(100000)
from functools import lru_cache

P3 = [3 ** i for i in range(14)]

def run(K, D):
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def L(p, r, j, d):
        if d == 0: return 1
        if j > 0 and r % 3 == 0: return 1          # 3 | m: chain forever
        # child 1: 2m — parity 0, ternary 2r mod 3^j, budget j
        best = L(0, (2 * r) % P3[j], j, d - 1) if j > 0 else 1
        if j == 0:
            return 1  # no ternary info: worst case no branching ever
        # child 2 exists iff m even (p==0) and r = 1 mod 3
        if p == 0 and r % 3 == 1:
            r2 = ((r - 1) // 3) % P3[j - 1]
            best = best + L(1, r2, j - 1, d - 1)
        return best
    worst = None
    import math
    for r in range(P3[K]):
        if r % 3 == 0: continue
        for p in (0, 1):
            v = L(p, r, K, D)
            if worst is None or v < worst: worst = v
    import math
    return worst, math.log2(worst) / D

out = {}
for K, D in [(4, 12), (5, 16), (6, 20), (7, 24), (8, 28), (8, 32)]:
    worst, gamma = run(K, D)
    out[f"K={K},D={D}"] = dict(min_leaves=worst, gamma=round(gamma, 4))
    print(f"K={K} D={D}: min leaves {worst}, gamma >= {gamma:.4f}", flush=True)
print(json.dumps(out, indent=1))
