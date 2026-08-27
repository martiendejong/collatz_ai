"""Obs 619: champion anatomy vs the tilted-coin model (agenda item 5).

The S4 census (Obs 613) cleared the survivor ENSEMBLE. Here: the
individual record orbits (stopping-time champions). Under the coin
model, an orbit conditioned to climb at rate delta bits/step has
halving counts w that are iid TILTED geometric: P(w) ~ 2^-w z^w with
z matched so that E[w] = log2(3) - delta. Test per champion, over the
climb window (start -> maximum):
  - match z to the observed mean w (one parameter, the drift);
  - compare VARIANCE and P(w=1) with the tilted model (shape test);
  - serial 2-gram MI vs permutation floor (independence test).
Any stable surplus = structure the coin model misses; agreement = the
champions are lucky, not clever (fair-coin runs, cf. R541).
"""
import math
import random
from collections import Counter

CHAMPIONS = [27, 703, 10087, 35655, 270271, 1126015, 56924955, 63728127,
             217740015, 1200991791, 1827397567, 2788008987]

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def climb_ws(n0):
    """w-sequence from start to the orbit maximum (climb window)."""
    n, ws, hs = n0, [], []
    while n != 1:
        m = 3 * n + 1
        w = v2(m)
        n = m >> w
        ws.append(w)
        hs.append(n)
    peak = max(range(len(hs)), key=lambda i: hs[i])
    return ws[: peak + 1]

def tilted_moments(z, wcap=60):
    ps = [(2.0 ** -w) * (z ** w) for w in range(1, wcap)]
    Z = sum(ps)
    ps = [p / Z for p in ps]
    mean = sum((w + 1) * p for w, p in enumerate(ps))
    var = sum((w + 1) ** 2 * p for w, p in enumerate(ps)) - mean ** 2
    return mean, var, ps[0]

def match_z(target_mean):
    lo, hi = 1e-6, 1.999
    for _ in range(60):
        mid = (lo + hi) / 2
        if tilted_moments(mid)[0] < target_mean:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def mi_pairs(seq, rng, cap=5):
    pairs = [(min(a, cap), min(b, cap)) for a, b in zip(seq, seq[1:])]
    def mi(ps):
        n = len(ps)
        cj, c1, c2 = Counter(ps), Counter(a for a, _ in ps), Counter(b for _, b in ps)
        return sum(c / n * math.log2(c * n / (c1[a] * c2[b])) for (a, b), c in cj.items())
    real = mi(pairs)
    fl = []
    for _ in range(20):
        sh = [b for _, b in pairs]
        rng.shuffle(sh)
        fl.append(mi(list(zip([a for a, _ in pairs], sh))))
    mu = sum(fl) / len(fl)
    sd = (sum((x - mu) ** 2 for x in fl) / len(fl)) ** 0.5 or 1e-9
    return real, mu, sd

if __name__ == "__main__":
    rng = random.Random(7)
    print(f"{'champion':>11} {'len':>4} {'w_mean':>6} | {'var':>5} {'var_mod':>7} "
          f"| {'P(w=1)':>6} {'model':>6} | {'MI':>6} {'floor':>6} {'z':>5}")
    for c in CHAMPIONS:
        ws = climb_ws(c)
        if len(ws) < 15:
            continue
        wm = sum(ws) / len(ws)
        z = match_z(wm)
        m_mean, m_var, m_p1 = tilted_moments(z)
        var = sum(w * w for w in ws) / len(ws) - wm ** 2
        p1 = ws.count(1) / len(ws)
        mi, fmu, fsd = mi_pairs(ws, rng)
        zsc = (mi - fmu) / fsd
        print(f"{c:>11} {len(ws):>4} {wm:>6.3f} | {var:>5.2f} {m_var:>7.2f} "
              f"| {p1:>6.3f} {m_p1:>6.3f} | {mi:>6.4f} {fmu:>6.4f} "
              f"{'z=%+.1f' % zsc:>7}")
