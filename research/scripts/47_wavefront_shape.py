"""R349-353: THE WAVEFRONT SHAPE rendered + edge laws.
The reached set B_t lives in a wedge in (t, log2 n) space:
- UPPER edge: pure doubling slides, log2 n = t exactly (slope 1)
- LOWER edge: the laggards, log2 n ~ t/c* -- c* measured, compared to the
  Lagarias-Weiss stochastic prediction for total-stopping-time records.
Figure: density heatmap of (d(n), log2 n) for all n < 2^24 + both edges +
laggard records marked. Uses same d-array as 46; graph/collatz_backward.graphml
carries the identical structure (depth attr) for interactive use."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding="utf-8")

N = 1 << 24
d = np.full(N, -1, dtype=np.int32)
d[1] = 0
for n in range(2, N):
    if d[n] >= 0: continue
    path = []; m = n
    while m >= N or d[m] < 0:
        path.append(m); m = m // 2 if m % 2 == 0 else (3 * m + 1) // 2
    base = d[m]
    for x in reversed(path):
        base += 1
        if x < N: d[x] = base

n_arr = np.arange(1, N, dtype=np.float64)
dd = d[1:].astype(np.float64)
l2 = np.log2(n_arr)

# edges per t
tmax = int(dd.max())
min_n = np.full(tmax + 1, np.inf); max_n = np.zeros(tmax + 1)
for t in range(tmax + 1):
    pass
# vectorized: for each n, update
np.minimum.at(min_n, d[1:], n_arr)
np.maximum.at(max_n, d[1:], n_arr)

# laggard records
recs = []
best = -1
for n in range(2, N):
    if d[n] > best: best = d[n]; recs.append((n, d[n]))

print("edge slopes: (t, log2 min_excluded ~ laggard edge)")
c_star = []
for t in (100, 150, 200, 250, 300, 350, 400, 441):
    if t <= tmax and np.isfinite(min_n[t]):
        # min VALUE at depth exactly t approximates the laggard frontier
        c = t / np.log2(min_n[t])
        c_star.append(c)
        print(f"  t={t:>3}: min value at depth t = {int(min_n[t]):>9,}  c* = d/log2(n) = {c:.2f}")
print("Lagarias-Weiss stochastic record prediction: c -> ~28.9 per log2 (41.68/ln 2... slow convergence)")

plt.figure(figsize=(12, 7))
sub = np.random.default_rng(349).choice(N - 1, 400_000, replace=False)
plt.hexbin(dd[sub], l2[sub], gridsize=120, cmap="viridis", bins="log", mincnt=1)
ts = np.arange(0, tmax + 1)
valid = np.isfinite(min_n) & (min_n > 1)
plt.plot(ts[valid], np.log2(min_n[valid]), "r-", lw=1.2, label="min value per depth (laggard edge)")
plt.plot([0, 24], [0, 24], "w--", lw=1.5, label="pure-doubling edge: log2 n = t")
rn = np.array([r[0] for r in recs[6:]]); rd = np.array([r[1] for r in recs[6:]])
plt.plot(rd, np.log2(rn), "wo", ms=4, mec="red", label="laggard records (A006877)")
plt.xlabel("d(n) = T-steps to reach 1  (wavefront time t)")
plt.ylabel("log2 n")
plt.title("The backward wavefront: reached set B_t for all n < 2^24\n"
          "wedge between the doubling slide (slope 1) and the laggard edge (slope ~1/18 and flattening)")
plt.colorbar(label="log10 count"); plt.legend(loc="lower right")
plt.tight_layout(); plt.savefig("figures/wavefront_shape.png", dpi=130)
print("\nfigure: figures/wavefront_shape.png")

# laggard family pattern
def tail_ones(x):
    k = 0
    while x & 1: x >>= 1; k += 1
    return k
odd_recs = [(n, dn) for n, dn in recs if n % 2 == 1 and n > 20]
ks = [tail_ones(n) for n, _ in odd_recs]
print(f"\nladder-fuel of the {len(odd_recs)} odd laggard records: mean k = {np.mean(ks):.2f} "
      f"(population mean over odd n = 2.00); max k = {max(ks)}")
big = [(n, tail_ones(n)) for n, _ in odd_recs if tail_ones(n) >= 6]
print("records with k >= 6 (heavy ladders):", big)
