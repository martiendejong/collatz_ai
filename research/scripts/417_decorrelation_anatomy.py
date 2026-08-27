"""Obs 624: scale-separation anatomy of Open Lemma'' part (A).

Claim to substantiate: the triple-shape field s and the dual measure u
decorrelate because they live on different digit scales · s is a
finest-trit functional, u is coarse-dominated. Measurable version on
cert_k13/k15:
  - R2 of log u explained by coarse residue dummies (mod 27, mod 243):
    should be substantial and GROW with the coarse resolution;
  - R2 of shape s explained by the same coarse dummies: should be ~0;
  - the cross-correlation then vanishes structurally (a function of
    coarse digits is uncorrelated with a coarse-independent field).
Also: which scale DOES determine s? R2 of s on the finest trits of the
row label (mod 27 of the class index t, the intra-triple coordinate).
"""
import numpy as np

def anatomy(K, lam, path):
    N = 3 ** (K - 1); M3 = 3 ** K; Mc = 3 ** (K - 1)
    beta = np.log2(3); W0 = lam ** -2; W2 = lam ** (beta - 2); W8 = lam ** (beta - 1)
    C = np.asarray(np.load(path, mmap_mode='r')[:]).astype(np.float64)
    idx = np.arange(N, dtype=np.int64); m = 3 * idx + 2
    i4 = (((4 * m) % M3) - 2) // 3; mod9 = m % 9
    tri = np.full(N, -1, dtype=np.int64); s = np.full(N, np.nan); Wb = np.zeros(N)
    for mask, mul, r, W in ((mod9 == 2, 4, 2, W2), (mod9 == 8, 2, 1, W8)):
        mm = m[mask]; t = ((mul * mm - r) // 3) % Mc
        j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3])
        vals = C[j]
        vs = np.sort(vals, axis=0)
        s[mask] = (vs[1] - vs[0]) / np.maximum(vs[2] - vs[0], 1e-300)
        o = np.argmin(vals, axis=0)
        tri[mask] = j[o, np.arange(j.shape[1])]
        Wb[mask] = W
        del mm, t, j, vals, vs, o
    hm = tri >= 0
    u = np.full(N, 1.0 / N)
    for it in range(120):
        nu = np.bincount(i4, weights=W0 * u, minlength=N)
        nu += np.bincount(tri[hm], weights=Wb[hm] * u[hm], minlength=N)
        nu /= nu.sum()
        if np.abs(nu - u).sum() < 1e-11:
            u = nu; break
        u = nu

    def r2_dummies(y, cls, ncls):
        mean_all = y.mean()
        ss_tot = ((y - mean_all) ** 2).sum()
        sums = np.bincount(cls, weights=y, minlength=ncls)
        cnts = np.bincount(cls, minlength=ncls)
        means = np.where(cnts > 0, sums / np.maximum(cnts, 1), 0)
        ss_exp = (cnts * (means - mean_all) ** 2).sum()
        return ss_exp / ss_tot

    lu = np.log(u[hm] + 1e-300)
    sh = s[hm]
    mh = m[hm]
    print(f"k={K}:")
    for mod in (27, 243, 2187):
        print(f"  coarse mod {mod:>5}: R2(log u) = {r2_dummies(lu, (mh % mod) // 1, mod):.3f}   "
              f"R2(shape s) = {r2_dummies(sh, (mh % mod) // 1, mod):.4f}")
    # fine end: top trits of the class (the intra-triple/branch coordinate)
    top = (idx[hm] * 27) // N   # top-3-trit coordinate of the class index
    print(f"  fine (top-3 trits of class): R2(log u) = {r2_dummies(lu, top, 27):.4f}   "
          f"R2(shape s) = {r2_dummies(sh, top, 27):.4f}")
    print(f"  direct corr(log u, s) = {np.corrcoef(lu, sh)[0, 1]:.4f}")

anatomy(13, 1.818, r"E:\projects\collatz\research\certificates\cert_k13.npy")
anatomy(15, 1.841, r"E:\projects\collatz\research\certificates\cert_k15.npy")
