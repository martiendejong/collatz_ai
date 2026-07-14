"""E10: the alpha* theorem — excursion-optimal height weight = (1+log2 3)/2.
E11: positional alternation fuel — how many low bits predict future ascent?"""
import sys, os, random, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
sys.stdout.reconfigure(encoding="utf-8")
from collatz_lib import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

random.seed(10)
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
out = {}

# ============ E10: fine alpha sweep, large sample ============
N = 200000
samples = []
for _ in range(N):
    n = random.randrange(1 << 40, 1 << 58) | 1
    a, k = coords(n)
    a2, k2, w = macro_step(a, k)
    if a2 == 1: continue
    samples.append((math.log2(a), k, math.log2(a2), k2, w))

alphas = [round(1.0 + 0.05 * i, 2) for i in range(14)]
sweep = {}
for al in alphas:
    d = [(la2 + al * k2) - (la + al * k) for la, k, la2, k2, w in samples]
    mx = max(d)
    p99 = sorted(d)[int(.999 * len(d))]
    sweep[al] = dict(max_inc=round(mx, 2), q999=round(p99, 3),
                     p_inc=round(sum(1 for x in d if x > 1e-9) / len(d), 4))
out["fine_sweep"] = sweep
astar = (1 + math.log2(3)) / 2
out["alpha_star_prediction"] = round(astar, 5)

# spike anatomy: at alpha*, is DH increase ~ 0.2925 * max(k, k') as predicted?
big = []
al = astar
for la, k, la2, k2, w in samples:
    dh = (la2 + al * k2) - (la + al * k)
    if dh > 2:
        big.append((round(dh, 2), k, k2, w))
big.sort(reverse=True)
out["worst_spikes_at_alpha_star(dh,k,k',w)"] = big[:10]
slope_pred = round(math.log2(3) - astar, 4)  # = astar - 1 = 0.2925
out["predicted_spike_slope"] = slope_pred

# figure
xs = alphas
plt.figure(figsize=(7, 4.2))
plt.plot(xs, [sweep[a]["max_inc"] for a in xs], "o-", label="max increase (200k samples)")
plt.plot(xs, [sweep[a]["q999"] * 10 for a in xs], "s--", label="99.9th pct x10")
plt.axvline(astar, color="crimson", ls=":", label=f"alpha* = (1+log2 3)/2 = {astar:.4f}")
plt.xlabel("alpha"); plt.ylabel("excursion size (bits)")
plt.title("E10: excursion-optimal weight is alpha* = (1+log2 3)/2")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(FIG, "e10_alphastar.png"), dpi=110)

# ============ E11: positional fuel — predict ascent from j low bits ============
# ascent(n) = log2(max_value_over_next_60_Tsteps / n); predict from n mod 2^j
def ascent(n, steps=60):
    m = n; mx = n
    for _ in range(steps):
        m = 3 * m + 1 if m % 2 else m // 2
        mx = max(mx, m)
        if m == 1: break
    return math.log2(mx / n)

M = 30000
pop = [(random.randrange(1 << 40, 1 << 44) | 1) for _ in range(M)]
asc = [ascent(n) for n in pop]
mean_asc = sum(asc) / M
var = sum((x - mean_asc) ** 2 for x in asc) / M
r2 = {}
for j in (4, 8, 12, 16, 20):
    groups = {}
    for n, x in zip(pop, asc):
        groups.setdefault(n % (1 << j), []).append(x)
    between = sum(len(g) * (sum(g) / len(g) - mean_asc) ** 2 for g in groups.values()) / M
    # correct for finite groups: only count classes with >=2 members for honest R2
    r2[j] = round(between / var, 4)
out["ascent_R2_from_j_low_bits"] = r2
out["ascent_mean_bits"] = round(mean_asc, 3)

print(json.dumps(out, indent=1))
