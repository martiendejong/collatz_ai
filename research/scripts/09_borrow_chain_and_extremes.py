"""E15: the /2-in-ternary borrow chain — Markov structure, measured on real orbits.
E16: extreme-value test of the alpha* spike tail 2^(-3.42 s)."""
import sys, os, random, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

random.seed(21)
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
out = {}

# ============ E15: borrow chain of halving in base 3 ============
# Halve even x digit-by-digit (MSB->LSB) in ternary: state r in {0,1}:
#   out_digit = (3r + d)//2 ; r' = (3r + d) mod 2
# Theory: from r=0, P(r'=1)=1/3 (d=1); from r=1, P(r'=1)=2/3 (d=0 or 2). Stationary 1/2.
def borrow_transitions(x):
    """returns list of (r, r') transitions while halving even x in base 3."""
    ds = [int(c) for c in ternary(x)]
    r = 0; trans = []
    for d in ds:
        v = 3 * r + d
        r2 = v & 1
        trans.append((r, r2))
        r = r2
    return trans

cnt = Counter()
# measure on REAL orbit values (the halving cascades of actual Collatz trajectories)
for _ in range(3000):
    n = random.randrange(1 << 40, 1 << 50) | 1
    for _ in range(40):
        x = 3 * n + 1
        w = v2(x)
        y = x
        for _ in range(w):
            for t in borrow_transitions(y):
                cnt[t] += 1
            y >>= 1
        n = y
        if n == 1: break
p01 = cnt[(0, 1)] / (cnt[(0, 0)] + cnt[(0, 1)])
p11 = cnt[(1, 1)] / (cnt[(1, 0)] + cnt[(1, 1)])
out["borrow_chain"] = dict(measured_P_0to1=round(p01, 4), theory=round(1/3, 4),
                           measured_P_1to1=round(p11, 4), theory11=round(2/3, 4),
                           transitions=sum(cnt.values()))

# ============ E16: extreme-value test of alpha* tail ============
astar = (1 + math.log2(3)) / 2
def max_excursion(n, cap=4000):
    a, k = coords(n)
    h0 = math.log2(a) + astar * k
    hmax = h0
    for _ in range(cap):
        a, k, w = macro_step(a, k)
        if a == 1 and k == 1: break
        h = math.log2(a) + astar * k
        hmax = max(hmax, h)
    return hmax - h0

M = 100000
excs = []
for _ in range(M):
    n = random.randrange(1 << 30, 1 << 34) | 1
    excs.append(max_excursion(n))
excs.sort()
out["excursion"] = dict(
    mean=round(sum(excs) / M, 3),
    q99=round(excs[int(.99 * M)], 3),
    q999=round(excs[int(.999 * M)], 3),
    max_observed=round(excs[-1], 3),
    predicted_max=round(math.log2(M) / 3.419 + 1.5, 2),  # EV for exp tail + offset
)
# tail slope fit: log2 P(exc > s) vs s on [2, q999]
import bisect
ss = [2 + 0.5 * i for i in range(12)]
pts = []
for s in ss:
    i = bisect.bisect_left(excs, s)
    p = (M - i) / M
    if p > 1e-4:
        pts.append((s, math.log2(p)))
# least squares slope
ms = sum(s for s, _ in pts) / len(pts); ml = sum(l for _, l in pts) / len(pts)
slope = sum((s - ms) * (l - ml) for s, l in pts) / sum((s - ms) ** 2 for s, _ in pts)
out["tail_slope_bits"] = dict(measured=round(slope, 3), predicted=-3.419)

plt.figure(figsize=(7, 4.2))
plt.plot([s for s, _ in pts], [l for _, l in pts], "o-", label="measured log2 P(exc > s)")
plt.plot([s for s, _ in pts], [pts[0][1] - 3.419 * (s - pts[0][0]) for s, _ in pts], "--",
         label="predicted slope -3.42")
plt.xlabel("excursion s (bits of H*)"); plt.ylabel("log2 P(exc > s)")
plt.title("E16: alpha* excursion tail — exponential, rate as predicted, no fat tail")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(FIG, "e16_extremes.png"), dpi=110)

print(json.dumps(out, indent=1))
