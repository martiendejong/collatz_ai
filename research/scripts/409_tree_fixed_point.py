"""Obs 617 part 2: derive the turnstile prefactors as a tree fixed point.

Exact identity: orbits hit v through exactly one direct predecessor
(forward orbits are injective; predecessor branches are disjoint), so
    P(hit v) = sum over odd non-spring predecessors u of P(hit u),
with u = (2^m v - 1)/3 over valid rungs; springs (u = 0 mod 3)
contribute 0 asymptotically (no predecessors, single-point mass).

Scheme: recurse only into children u < VCUT (this subtree is finite:
a backward path cannot repeat values, so it must eventually exceed
VCUT); for u >= VCUT apply the boundary ansatz P = 1/u (the mean visit
law). P is LINEAR in the boundary values, so target RATIOS are
invariant under the boundary's overall scale; only its 1/u shape
matters. Sensitivity probe: also run boundary shapes u^-0.9 and u^-1.1.
"""
import sys
sys.setrecursionlimit(20000)

MMAX = 64

def solve(vcut, expo):
    def P(v):
        s = 0.0
        for m in range(1, MMAX + 1):
            t = (1 << m) * v - 1
            if t % 3 == 0:
                u = t // 3
                if u % 2 == 1 and u % 3 != 0 and u > 1:
                    s += P(u) if u < vcut else u ** (-expo)
        return s
    return {v: P(v) for v in (5, 13, 53, 85, 341, 853, 3413, 5461, 21845)}

if __name__ == "__main__":
    meas = {"341/85": 1.582, "21845/5461": 5.7, "53/13": 0.9586,
            "5abs": 0.93793}
    print(f"{'vcut':>9} {'shape':>6} | P(341)/P(85)  P(21845)/P(5461)  P(53)/P(13)")
    for vcut in (10 ** 4, 10 ** 5, 10 ** 6, 4 * 10 ** 6):
        for expo in (1.0,):
            r = solve(vcut, expo)
            print(f"{vcut:>9} {expo:>6.2f} | {r[341] / r[85]:12.4f}  "
                  f"{r[21845] / r[5461]:16.3f}  {r[53] / r[13]:11.4f}")
    for expo in (0.9, 1.1):
        r = solve(10 ** 6, expo)
        print(f"{10**6:>9} {expo:>6.2f} | {r[341] / r[85]:12.4f}  "
              f"{r[21845] / r[5461]:16.3f}  {r[53] / r[13]:11.4f}")
    print(f"\nmeasured (Obs 616): 341/85 = 1.582, 21845/5461 = 5.7, 53/13 = 0.959")
