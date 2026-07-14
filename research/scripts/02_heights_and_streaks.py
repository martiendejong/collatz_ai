"""E2: height-function sweep over macro-steps (with 5n+1 control).
E3: ternary-length streak analysis (the Foster-Lyapunov candidate)."""
import sys, os, random, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

random.seed(2)
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
out = {}

# ============ E2: H_alpha(a,k) = log2(a) + alpha*k over macro-steps ============
alphas = [0.0, 0.4, 0.585, 0.8, 1.0, 1.2, 1.585, 2.0]
N = 120000
samples = []
for _ in range(N):
    n = random.randrange(1 << 30, 1 << 46) | 1
    a, k = coords(n)
    a2, k2, w = macro_step(a, k)
    if a2 == 1 and k2 == 1:
        continue
    samples.append((math.log2(a), k, math.log2(a2), k2))

res = {}
for al in alphas:
    d = [(la2 + al * k2) - (la + al * k) for la, k, la2, k2 in samples]
    inc = sum(1 for x in d if x > 1e-9)
    res[al] = dict(mean=round(sum(d) / len(d), 4),
                   p_inc=round(inc / len(d), 4),
                   max_inc=round(max(d), 2))
out["height_sweep"] = res

# blocks: does H decrease over blocks of B macro-steps? (alpha = log2 3 - 1 = 0.585)
al = 0.585
blocks = {}
for B in (1, 2, 5, 10, 20):
    incs = 0; tot = 0; mx = 0
    for _ in range(4000):
        n = random.randrange(1 << 40, 1 << 60) | 1
        a, k = coords(n)
        h0 = math.log2(a) + al * k
        for _ in range(B):
            a, k, w = macro_step(a, k)
            if a == 1 and k == 1:
                break
        h1 = math.log2(a) + al * k
        tot += 1
        if h1 > h0 + 1e-9:
            incs += 1
            mx = max(mx, h1 - h0)
    blocks[B] = dict(p_inc=round(incs / tot, 4), max_inc=round(mx, 1))
out["blocks_alpha585"] = blocks

# 5n+1 control: same H, same macro structure (segments consume 2 ones per step)
ctrl = []
for _ in range(60000):
    n = random.randrange(1 << 30, 1 << 46) | 1
    m, w = syracuse_c(n, 5)
    a, k = coords(n); a2, k2 = coords(m)
    ctrl.append((math.log2(a) + 0.585 * k, math.log2(a2) + 0.585 * k2))
d = [b - a for a, b in ctrl]
out["control_5np1_alpha585"] = dict(mean=round(sum(d) / len(d), 4),
                                    p_inc=round(sum(1 for x in d if x > 0) / len(d), 4))

# ============ E3: ternary-length streaks ============
# streak = maximal run of consecutive Syracuse steps where ternlen does NOT decrease
def streaks_of(n, steps):
    s = []; cur = 0; L = ternlen(n)
    for _ in range(steps):
        n, w = syracuse(n)
        if n == 1: break
        L2 = ternlen(n)
        if L2 >= L: cur += 1
        else:
            if cur: s.append(cur)
            cur = 0
        L = L2
    if cur: s.append(cur)
    return s

allst = Counter()
for _ in range(20000):
    n = random.randrange(1 << 30, 1 << 40) | 1
    for st in streaks_of(n, 60):
        allst[st] += 1
tot = sum(allst.values())
out["streak_dist"] = {j: round(allst[j] / tot, 5) for j in range(1, 13)}
# geometric fit: P(streak >= j+1 | >= j)
rat = []
for j in range(1, 9):
    ge_j = sum(v for s, v in allst.items() if s >= j)
    ge_j1 = sum(v for s, v in allst.items() if s >= j + 1)
    rat.append(round(ge_j1 / ge_j, 3) if ge_j else None)
out["streak_continuation_ratios"] = rat

# longest-streak champions up to 2^22: what is their binary structure?
best = []
for n in range(3, 1 << 22, 2):
    s = streaks_of(n, 80)
    if s and max(s) >= 12:
        best.append((max(s), n, bin(n)[2:]))
best.sort(reverse=True)
out["streak_champions"] = [(s, n, b, "trail1s=" + str(trailing_ones(n))) for s, n, b in best[:8]]

# figure: streak distribution vs geometric
xs = list(range(1, 13))
ys = [allst[j] / tot for j in xs]
p = 1 - 1 / (sum(allst[j] * j for j in allst) / tot)  # crude
plt.figure(figsize=(7, 4.2))
plt.semilogy(xs, ys, "o-", label="observed streak length")
g0 = ys[0]
plt.semilogy(xs, [g0 * (rat[2] or .35) ** (j - 1) for j in xs], "--", label=f"geometric, ratio~{rat[2]}")
plt.xlabel("ternary-length non-decrease streak"); plt.ylabel("frequency")
plt.title("E3: ternlen streaks die geometrically (no structural cap found)")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(FIG, "e3_streaks.png"), dpi=110)

print(json.dumps(out, indent=1))
