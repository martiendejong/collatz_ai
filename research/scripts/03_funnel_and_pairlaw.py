"""E4: funnel decay across the three systems. E5: pair-merge law census + alternation fuel."""
import sys, os, random, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

random.seed(3)
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
out = {}

# ============ E4: funnel decay, 3 systems, N=20000 ============
N = 20000
curves = {}
for name, f in [("3n+1", collatz_map), ("3n-1", map_3nm1), ("5n+1", map_5np1)]:
    vals = set(range(1, N + 1))
    c = []
    for step in range(1, 301):
        vals = {f(n) for n in vals}
        c.append(len(vals))
    curves[name] = c
out["funnel_at"] = {nm: {s: curves[nm][s - 1] for s in (10, 50, 100, 200, 300)} for nm in curves}

plt.figure(figsize=(7.5, 4.4))
for nm, c in curves.items():
    plt.semilogy(range(1, 301), c, label=nm)
plt.xlabel("step"); plt.ylabel("surviving distinct streams (of 20000)")
plt.title("E4: identical pair-merge funnel, three different fates")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(FIG, "e4_funnel.png"), dpi=110)

# ============ E5a: pair-merge law census ============
# true pair = consecutive (k,k+1) in family a with y=(a*3^k-1)/2 odd. Verify merge within
# (segments + slack) steps, and record family pairing offset.
def trace(n, L):
    s = [n]
    for _ in range(L):
        n = collatz_map(n)
        s.append(n)
        if n == 1: break
    return s
tot = merged = 0
offsets = {0: 0, 1: 0}
for a in range(1, 2000, 2):
    if a % 3 == 0: continue
    got = None
    for k in range(1, 10):
        y = (a * 3 ** k - 1) // 2
        if y % 2 == 1:
            if got is None:
                got = k % 2
            n1, n2 = starter(a, k), starter(a, k + 1)
            t1 = set(trace(n1, 4 * (k + 40)))
            t2 = set(trace(n2, 4 * (k + 44)))
            tot += 1
            merged += bool(t1 & t2)
    if got is not None:
        offsets[got] += 1
out["pair_law"] = dict(true_pairs=tot, merged=merged, family_offset_counts=offsets)

# ============ E5b: alternation fuel: 3 * (01)^m = (1)^2m identity + champions ============
fuel = []
for m in (2, 4, 8):
    x = (4 ** m - 1) // 3          # binary 0101...01, m ones
    fuel.append((bin(x)[2:], bin(3 * x)[2:]))
out["alternation_fuel_identity"] = fuel

# how much alternation do long-streak numbers carry vs random?
def alt_score(n):
    b = bin(n)[2:]
    return sum(1 for i in range(1, len(b)) if b[i] != b[i - 1]) / (len(b) - 1)
def streak_max(n, steps=80):
    L = ternlen(n); cur = mx = 0
    for _ in range(steps):
        n, w = syracuse(n)
        if n == 1: break
        L2 = ternlen(n)
        if L2 >= L: cur += 1; mx = max(mx, cur)
        else: cur = 0
        L = L2
    return mx
hi, lo = [], []
for _ in range(40000):
    n = random.randrange(1 << 20, 1 << 22) | 1
    s = streak_max(n)
    (hi if s >= 25 else lo).append(alt_score(n))
out["alternation_in_streakers"] = dict(
    n_high=len(hi), mean_alt_high=round(sum(hi) / len(hi), 4) if hi else None,
    mean_alt_rest=round(sum(lo) / len(lo), 4))

print(json.dumps(out, indent=1))
