"""Obs 626: the pathwise clipping inequality, verified on the optimum.

LEMMA (pathwise clipping, unconditional): for a sorted triple
v1 <= v2 <= v3 with gap g = v2 - v1 and any perturbation xi,
    min_i(v_i + xi_i) <= v1 + xi_1 - (xi_1 - xi_2 - g)^+ .
Proof: min <= v2 + xi_2 = v1 + xi_1 - (xi_1 - xi_2 - g); combine with
min <= v1 + xi_1. QED (two lines, no hypotheses).

Consequence: the response of the min to a fine perturbation field is
attenuated by at least E[(dxi - g)^+], where dxi = xi_1 - xi_2. The
attenuation fraction lower bound is
    lam_bound = E_mu[(dxi - g)^+] / E[xi_1-scale],
positive whenever the joint law of (dxi, g) has mass at dxi > g ·
which is exactly the anti-concentration scalar (dual-weighted shape
margin, measured flat ~0.53).

This script computes lam_bound on cert_k13/k15 with the system's OWN
fine-field statistics (dxi resampled from the empirical intra-triple
relative differences, independently of the row = part (A) imposed by
shuffling), dual-weighted, and compares with the measured
lam_clip = 0.085-0.140. Nonvacuous bound = the lemma has real content
on the actual optimum.
"""
import numpy as np

def run(K, lam, path, nsamp=2_000_000, seed=17):
    N = 3 ** (K - 1); M3 = 3 ** K; Mc = 3 ** (K - 1)
    beta = np.log2(3); W0 = lam ** -2; W2 = lam ** (beta - 2); W8 = lam ** (beta - 1)
    C = np.asarray(np.load(path, mmap_mode='r')[:]).astype(np.float64)
    idx = np.arange(N, dtype=np.int64); m = 3 * idx + 2
    i4 = (((4 * m) % M3) - 2) // 3; mod9 = m % 9
    tri = np.full(N, -1, dtype=np.int64)
    g_rel = np.full(N, np.nan); d21 = np.full(N, np.nan); Wb = np.zeros(N)
    for mask, mul, r, W in ((mod9 == 2, 4, 2, W2), (mod9 == 8, 2, 1, W8)):
        mm = m[mask]; t = ((mul * mm - r) // 3) % Mc
        j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3])
        vals = C[j]; vs = np.sort(vals, axis=0)
        g_rel[mask] = (vs[1] - vs[0]) / vs[0]
        d21[mask] = (vs[2] - vs[1]) / vs[0]
        o = np.argmin(vals, axis=0)
        tri[mask] = j[o, np.arange(j.shape[1])]
        Wb[mask] = W
        del mm, t, j, vals, vs, o
    hm = tri >= 0
    u = np.full(N, 1.0 / N)
    for it in range(150):
        nu = np.bincount(i4, weights=W0 * u, minlength=N)
        nu += np.bincount(tri[hm], weights=Wb[hm] * u[hm], minlength=N)
        nu /= nu.sum()
        if np.abs(nu - u).sum() < 1e-12:
            u = nu; break
        u = nu
    rng = np.random.default_rng(seed)
    w = u[hm] / u[hm].sum()
    rows = rng.choice(np.where(hm)[0], size=nsamp, p=w)          # dual-sampled rows
    g = g_rel[rows]
    # perturbation differences: resampled from the empirical pool of
    # intra-triple relative differences, INDEPENDENT of the row (part A)
    pool = np.concatenate([g_rel[hm], d21[hm]])
    dxi = rng.choice(pool, size=nsamp) - rng.choice(pool, size=nsamp)
    clip = np.maximum(dxi - g, 0)
    scale = np.mean(np.abs(dxi))
    lam_bound = clip.mean() / scale
    frac_active = np.mean(dxi > g)
    print(f"k={K}: E[(dxi-g)^+] = {clip.mean():.5f}, scale E|dxi| = {scale:.5f}")
    print(f"      lam_bound = {lam_bound:.4f}   P(dxi > g) = {frac_active:.4f}")
    print(f"      (measured lam_clip on this system: 0.085-0.140; bound nonvacuous if > 0)")

run(13, 1.818, r"E:\projects\collatz\research\certificates\cert_k13.npy")
run(15, 1.841, r"E:\projects\collatz\research\certificates\cert_k15.npy")
