"""Obs 623: the decisive scalar at k=17 + first P1 measurement.

(1) The repaired Open Lemma's health scalar: dual-weighted gap median
    vs k. Measured 0.0327 (k=13), 0.0255 (k=15). Saturation (Prop 23
    pattern) = lemma true-shaped; geometric decay = kappa creep.
    Three points allow an exact (c_inf, A, r) fit of median = c_inf
    + A r^k; k=19 can then confirm out-of-sample.
(2) P1 opening move: does tempering hold on the LEFT? Coarse-grain the
    dual (left Perron) mass and the certificate (right vector) mass to
    mod 27/81/243 and measure corr + fitted exponent of
    log u-mass ~ log C-mass.
"""
import numpy as np

def probe(K, lam, path, iters=120):
    N = 3 ** (K - 1); M3 = 3 ** K; Mc = 3 ** (K - 1)
    beta = np.log2(3); W0 = lam ** -2; W2 = lam ** (beta - 2); W8 = lam ** (beta - 1)
    C = np.asarray(np.load(path, mmap_mode='r')[:]).astype(np.float64)
    idx = np.arange(N, dtype=np.int64)
    m = 3 * idx + 2
    i4 = (((4 * m) % M3) - 2) // 3
    mod9 = m % 9
    tri = np.full(N, -1, dtype=np.int64)
    gap = np.full(N, np.nan)
    Wb = np.zeros(N)
    for mask, mul, r, W in ((mod9 == 2, 4, 2, W2), (mod9 == 8, 2, 1, W8)):
        mm = m[mask]
        t = ((mul * mm - r) // 3) % Mc
        j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3])
        vals = C[j]
        o = np.argsort(vals, axis=0)
        v1 = np.take_along_axis(vals, o[0:1], 0)[0]
        v2 = np.take_along_axis(vals, o[1:2], 0)[0]
        tri[mask] = np.take_along_axis(j, o[0:1], 0)[0]
        gap[mask] = (v2 - v1) / v1
        Wb[mask] = W
        del mm, t, j, vals, o, v1, v2
    hm = tri >= 0
    u = np.full(N, 1.0 / N)
    w0u = np.empty(N)
    for it in range(iters):
        nu = np.bincount(i4, weights=W0 * u, minlength=N)
        nu += np.bincount(tri[hm], weights=Wb[hm] * u[hm], minlength=N)
        nu /= nu.sum()
        if np.abs(nu - u).sum() < 1e-11:
            u = nu
            break
        u = nu
    w = u[hm]; g = gap[hm]
    order = np.argsort(g)
    cw = np.cumsum(w[order]); cw /= cw[-1]
    def wq(q): return g[order][np.searchsorted(cw, q)]
    med_w = wq(0.50)
    print(f"k={K}: dual-weighted gap p05/p25/p50 = {wq(.05):.4f}/{wq(.25):.4f}/{med_w:.4f} "
          f"(flat p50 {np.median(g):.4f}); iters {it}")
    print(f"      dual mass on gap<1e-3: {w[g < 1e-3].sum() / w.sum():.4f} "
          f"(flat {np.mean(g < 1e-3):.4f}); "
          f"corr(log u, log gap) = {np.corrcoef(np.log(w + 1e-15), np.log(g + 1e-15))[0, 1]:.3f}")
    us = np.sort(u)
    print(f"      top-1% dual share: {us[-N // 100:].sum():.4f}")
    # (2) left-right coarse comparison
    for mod in (27, 81, 243):
        r_ = m % mod
        um = np.bincount(r_, weights=u, minlength=mod)
        cm = np.bincount(r_, weights=C, minlength=mod)
        sel = (um > 0) & (cm > 0)
        lu, lc = np.log(um[sel]), np.log(cm[sel])
        cc = np.corrcoef(lu, lc)[0, 1]
        a, b = np.polyfit(lc, lu, 1)
        print(f"      left-vs-right mod {mod:>3}: corr(log) = {cc:.3f}, exponent = {a:.3f}")
    return med_w

if __name__ == "__main__":
    med17 = probe(17, 1.86168, r"E:\projects\collatz\research\certificates\cert_k17.npy")
    meds = {13: 0.0327, 15: 0.0255, 17: med17}
    # exact 3-point fit med(k) = c + A r^((k-13)/2)
    d1 = meds[15] - meds[13]; d2 = meds[17] - meds[15]
    r = d2 / d1
    A = d1 / (r - 1)
    c = meds[13] - A
    print(f"\n3-point fit: med(k) = {c:.4f} + {A:.4f} * {r:.3f}^((k-13)/2)")
    print(f"asymptote c_inf = {c:.4f}  ({'SATURATION > 0' if c > 0.005 else 'DECAY TOWARD 0'})")
    print(f"prediction k=19: {c + A * r ** 3:.4f}")
