"""Obs 622: LP-specialist probe of the Open Lemma on cert_k13.

The Open Lemma (drift_lowpass, open:kappa) claims uniform attenuation
kappa <= kappa_max < 1 at the K-L fixed point, given inf CV_1 > 0, via:
binding entries (argmin of each triple) are pinned by strictly positive
dual weights (complementary slackness) with concentration bounded below.

What an optimization specialist demands to see first, on the actual
optimum (k = 13, lambda = 1.818):
  (1) STRICT-COMPLEMENTARITY MARGINS: distribution of relative triple
      gaps g = (second smallest - smallest)/smallest over all min-rows.
      Degeneracy (g ~ 0) is where the binding set is unstable and the
      sensitivity argument (step 3) breaks.
  (2) DUAL CONCENTRATION: left Perron vector u of the argmin
      linearization M (row m: W0 at i4(m), W_b at the argmin member).
      The lemma needs u's pinning strength bounded below uniformly:
      report the concentration profile of u and the joint margin
      u_i * g_i (rows that are BOTH heavy and degenerate are the
      obstruction).
"""
import numpy as np

K = 13
N = 3 ** (K - 1)
M3 = 3 ** K
Mc = 3 ** (K - 1)
lam = 1.818
beta = np.log2(3)
W0 = lam ** -2
W2 = lam ** (beta - 2)
W8 = lam ** (beta - 1)

C = np.load(r"E:\projects\collatz\research\certificates\cert_k13.npy").astype(np.float64)
idx = np.arange(N, dtype=np.int64)
m = 3 * idx + 2
i4 = (((4 * m) % M3) - 2) // 3
mod9 = m % 9
rows = {"2": (mod9 == 2, 4, 2, W2), "8": (mod9 == 8, 2, 1, W8)}

tri_min_idx = np.full(N, -1, dtype=np.int64)
gap = np.full(N, np.nan)
Wb_row = np.zeros(N)
for tag, (mask, mul, r, Wb) in rows.items():
    mm = m[mask]
    t = ((mul * mm - r) // 3) % Mc
    j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3])
    vals = C[j]                      # 3 x nrows
    order = np.argsort(vals, axis=0)
    v1 = np.take_along_axis(vals, order[0:1], 0)[0]
    v2 = np.take_along_axis(vals, order[1:2], 0)[0]
    tri_min_idx[mask] = np.take_along_axis(j, order[0:1], 0)[0]
    gap[mask] = (v2 - v1) / v1
    Wb_row[mask] = Wb

g = gap[~np.isnan(gap)]
print(f"(1) strict-complementarity margins over {len(g):,} min-rows:")
for q in (1, 5, 10, 25, 50, 75, 90):
    print(f"    p{q:02d} gap = {np.percentile(g, q):.4f}")
for thr in (1e-4, 1e-3, 1e-2, 0.05):
    print(f"    fraction gap < {thr:g}: {np.mean(g < thr):.4f}")

# (2) left Perron vector of the argmin linearization (power iteration on M^T)
u = np.ones(N) / N
has_min = tri_min_idx >= 0
for it in range(200):
    nu = np.zeros(N)
    np.add.at(nu, i4, W0 * u)
    np.add.at(nu, tri_min_idx[has_min], Wb_row[has_min] * u[has_min])
    mu = nu.sum() / u.sum()
    nu /= nu.sum()
    if np.abs(nu - u).sum() < 1e-12:
        break
    u = nu
print(f"\n(2) left Perron: eigenvalue mu = {mu:.5f} (target ~ lambda-normalized 1), iters {it}")
u_scaled = u * N   # 1 = uniform
for q in (1, 5, 25, 50, 75, 95, 99):
    print(f"    p{q:02d} u/uniform = {np.percentile(u_scaled, q):.4f}")
print(f"    top-1% dual mass share: {np.sort(u)[-N // 100:].sum():.4f}")
print(f"    bottom-50% dual mass share: {np.sort(u)[: N // 2].sum():.4f}")

# joint obstruction metric: heavy AND degenerate rows
mask_r = has_min
w = u[mask_r] * N
gg = gap[mask_r]
heavy_deg = np.mean((w > 1.0) & (gg < 1e-2))
print(f"\n    joint: fraction of min-rows heavy (u>uniform) AND near-degenerate (gap<1e-2): {heavy_deg:.4f}")
c1 = np.corrcoef(np.log(w + 1e-12), np.log(gg + 1e-12))[0, 1]
print(f"    corr(log dual weight, log gap) = {c1:.3f}")
