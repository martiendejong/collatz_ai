"""Obs 628: the two normalization jobs of Prop 40, measured to the bone.

(iv-b) ENERGY-CLIPPING DECOMPOSITION (exact identity, then measured):
    min response R = xi_1 - X, X = max(0, xi1-xi2-g2, xi1-xi3-g3);
    energy loss L = E[X(2 xi_1 - X)] = E[X(xi_1 + xi_j)] + E[X g_j]
    (j = switching index), EXACT. For symmetric fields the first term
    vanishes (sum indep of difference); it is the candidate identity
    behind the measured +-0.026 antisymmetric correction. Measure both
    terms dual-weighted on cert_k13/k15 with the system's own field.

(iv-a) CLOSURE MARGINS per scale on cert_k13 (Thm 16 channels, script
    65 conventions, dual-weighted rows): channel decorrelation
    corr(T, B) per scale P, energy profile V(P), two-term fit
    V_P = a V_(P-1) + c V_(P+1), residual e(P) vs margin (c-a) V_P.
"""
import math
import numpy as np

def load(K, lam, path):
    N = 3 ** (K - 1); M3 = 3 ** K; Mc = 3 ** (K - 1)
    beta = math.log2(3); W0 = lam ** -2; W2 = lam ** (beta - 2); W8 = lam ** (beta - 1)
    C = np.asarray(np.load(path, mmap_mode='r')[:]).astype(np.float64)
    C /= C.mean()
    idx = np.arange(N, dtype=np.int64); m = 3 * idx + 2
    i4 = (((4 * m) % M3) - 2) // 3; mod9 = m % 9
    tri = np.full(N, -1, dtype=np.int64)
    g2 = np.full(N, np.nan); g3 = np.full(N, np.nan); Wb = np.zeros(N)
    for mask, mul, r, W in ((mod9 == 2, 4, 2, W2), (mod9 == 8, 2, 1, W8)):
        mm = m[mask]; t = ((mul * mm - r) // 3) % Mc
        j = np.stack([(t - 2) // 3, ((t + Mc) - 2) // 3, ((t + 2 * Mc) - 2) // 3])
        vals = C[j]; vs = np.sort(vals, axis=0)
        g2[mask] = (vs[1] - vs[0]) / vs[0]; g3[mask] = (vs[2] - vs[0]) / vs[0]
        o = np.argmin(vals, axis=0)
        tri[mask] = j[o, np.arange(j.shape[1])]
        Wb[mask] = W
        del mm, t, j, vals, vs, o
    hm = tri >= 0
    u = np.full(N, 1.0 / N)
    for _ in range(150):
        nu = np.bincount(i4, weights=W0 * u, minlength=N)
        nu += np.bincount(tri[hm], weights=Wb[hm] * u[hm], minlength=N)
        nu /= nu.sum()
        if np.abs(nu - u).sum() < 1e-12:
            u = nu; break
        u = nu
    return dict(N=N, C=C, u=u, hm=hm, g2=g2, g3=g3, W0=W0, W2=W2, W8=W8,
                lam=lam, K=K, Mc=Mc, M3=M3)

def part_ivb(S, nsamp=2_000_000, seed=23):
    rng = np.random.default_rng(seed)
    hm, u, g2, g3 = S['hm'], S['u'], S['g2'], S['g3']
    w = u[hm] / u[hm].sum()
    rows = rng.choice(np.where(hm)[0], size=nsamp, p=w)
    G2, G3 = g2[rows], g3[rows]
    pool = np.concatenate([g2[hm], g3[hm] - g2[hm]])
    pool = pool - pool.mean()
    xi = pool[rng.integers(0, len(pool), size=(3, nsamp))]
    c2 = xi[0] - xi[1] - G2
    c3 = xi[0] - xi[2] - G3
    X = np.maximum(0, np.maximum(c2, c3))
    sw = X > 0
    j3 = (c3 >= c2)          # switching partner (3 if its clip is larger)
    xij = np.where(j3, xi[2], xi[1])
    gj = np.where(j3, G3, G2)
    L = (X * (2 * xi[0] - X)).mean()
    main = (X * gj)[sw].sum() / nsamp
    corr_term = (X * (xi[0] + xij))[sw].sum() / nsamp
    E1 = (xi[0] ** 2).mean()
    print(f"k={S['K']}: energy loss L/E[xi^2] = {L / E1:.4f}")
    print(f"        main (gap) term  = {main / E1:+.4f}   correction (sum) term = {corr_term / E1:+.4f}")
    print(f"        identity check: main+corr = {(main + corr_term) / E1:.4f} vs L = {L / E1:.4f}")
    print(f"        correction/main ratio = {corr_term / main:+.3f}")

def part_iva(S, nsamp=400_000, seed=29):
    K, C, u, Mc, M3f = S['K'], S['C'], S['u'], S['Mc'], S['M3']
    N = S['N']; lam = S['lam']; W0 = S['W0']
    A = math.log2(3); W2, W8 = S['W2'], S['W8']
    M = 3 ** K; M3s = 3 ** (K - 2)
    rng = np.random.default_rng(seed)
    def cls(mm): return ((mm % M) - 2) // 3
    def tri_min(b):
        return np.minimum(np.minimum(C[b % M3s], C[b % M3s + M3s]), C[b % M3s + 2 * M3s])
    probs = u / u.sum()
    Vs = {}; corrs = {}
    for P in range(3, 8):
        d = 2 * 3 ** P
        js = rng.choice(N, size=nsamp, p=probs)
        m = 3 * js + 2; mp = m + d
        LHS = C[cls(mp)] - C[cls(m)]
        y = 4 * m; yp = 4 * mp
        T = W0 * ((C[cls(y + d)] - C[cls(y)]) + (C[cls(yp)] - C[cls(y + d)]))
        mod9 = m % 9
        B = np.zeros_like(LHS)
        for mv, mul, sub, W in ((2, 4, 2, W2), (8, 2, 1, W8)):
            sel = mod9 == mv
            t = ((mul * m[sel] - sub) // 3) % Mc
            tp = ((mul * mp[sel] - sub) // 3) % Mc
            B[sel] = W * (tri_min((tp - 2) // 3) - tri_min((t - 2) // 3))
        resid = LHS - (T + B)
        Vs[P] = np.var(LHS)
        corrs[P] = np.corrcoef(T, B)[0, 1]
        print(f"  P={P}: V = {Vs[P]:.3e}  corr(T,B) = {corrs[P]:+.4f}  "
              f"identity resid/std = {resid.std() / LHS.std():.4f}")
    # two-term closure fit on interior scales
    Ps = [4, 5, 6]
    Amat = np.array([[Vs[P - 1], Vs[P + 1]] for P in Ps])
    b = np.array([Vs[P] for P in Ps])
    (a, c), *_ = np.linalg.lstsq(Amat, b, rcond=None)
    print(f"  closure fit: a = {a:.4f}, c = {c:.4f}  (c - a = {c - a:+.4f})")
    for P in Ps:
        e = Vs[P] - a * Vs[P - 1] - c * Vs[P + 1]
        print(f"    P={P}: residual e/V = {e / Vs[P]:+.4f}   margin (c-a) = {c - a:+.4f}  "
              f"{'OK margin > |e/V|' if abs(e / Vs[P]) < abs(c - a) else 'MARGIN FAILS'}")

if __name__ == "__main__":
    for K, lam, path in ((13, 1.818, r"E:\projects\collatz\research\certificates\cert_k13.npy"),
                         (15, 1.841, r"E:\projects\collatz\research\certificates\cert_k15.npy")):
        S = load(K, lam, path)
        part_ivb(S)
        if K == 13:
            part_iva(S)
