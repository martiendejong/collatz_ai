"""Obs 620: mechanism of the pair law + fixed-point robustness.

(1) Rung decomposition: split the fixed-point mass P(v) by first
    backward rung, for the pair 85 vs 341 (and 5461 vs 21845). The
    hypothesis from ternary structure: the pair asymmetry is made by
    SPRING PLACEMENT on the low rungs of the dominant chains (a
    mod-9 roulette property), e.g. 85's main child 113 loses its own
    m=1 child to the spring 75, while 341's main child 227 keeps 151.
(2) Dominant-chain trace: follow the lowest alive rung 12 deep for
    both pair members; print each node mod 9 and springness pattern.
(3) Boundary-noise Monte Carlo (honest version of "convergence"):
    multiply every boundary value 1/u by iid lognormal noise
    exp(sigma*Z), sigma = 1.0 (an order of magnitude per node!) and
    recompute the ratios 20x. If the ratios are stable, the fixed
    point is statistically robust against the multifractal boundary
    spread even though no uniform pointwise bracket exists.
"""
import random

MMAX = 64
VCUT = 1_000_000

def children(v):
    out = []
    for m in range(1, MMAX + 1):
        t = (1 << m) * v - 1
        if t % 3 == 0:
            u = t // 3
            if u % 2 == 1 and u > 1:
                out.append((m, u))
    return out

def P(v, bnd):
    s = 0.0
    for _, u in children(v):
        if u % 3 == 0:
            continue
        s += P(u, bnd) if u < VCUT else bnd(u)
    return s

if __name__ == "__main__":
    plain = lambda u: 1.0 / u
    print("(1) rung decomposition (share of P(v) per first rung):")
    for v in (85, 341, 5461, 21845):
        tot = P(v, plain)
        parts = []
        for m, u in children(v)[:6]:
            share = 0.0 if u % 3 == 0 else (P(u, plain) if u < VCUT else 1 / u) / tot
            parts.append(f"m={m}:u={u}{'(SPRING)' if u % 3 == 0 else ''} {share*100:.1f}%")
        print(f"  v={v:>6}: " + "  ".join(parts))

    print("\n(2) dominant chain (lowest alive rung), node mod 9:")
    for v in (85, 341):
        chain, cur = [], v
        for _ in range(12):
            alive = [(m, u) for m, u in children(cur) if u % 3 != 0]
            dead_low = [(m, u) for m, u in children(cur)[:2] if u % 3 == 0]
            tag = "!" if dead_low and dead_low[0][0] < alive[0][0] else ""
            cur = alive[0][1]
            chain.append(f"{cur}{tag}(m9={cur % 9})")
        print(f"  {v}: " + " -> ".join(chain[:8]))

    print("\n(3) boundary-noise MC (sigma=1.0 lognormal per boundary node):")
    rng = random.Random(11)
    ratios = {"341/85": [], "53/13": [], "21845/5461": []}
    for rep in range(20):
        cache = {}
        def noisy(u):
            if u not in cache:
                cache[u] = (1.0 / u) * pow(2.718281828, rng.gauss(0, 1.0))
            return cache[u]
        vals = {v: P(v, noisy) for v in (5, 13, 53, 85, 341, 5461, 21845)}
        ratios["341/85"].append(vals[341] / vals[85])
        ratios["53/13"].append(vals[53] / vals[13])
        ratios["21845/5461"].append(vals[21845] / vals[5461])
    for k, xs in ratios.items():
        mu = sum(xs) / len(xs)
        sd = (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5
        print(f"  {k:>12}: {mu:.4f} +- {sd:.4f}  (plain: see Obs 617)")
