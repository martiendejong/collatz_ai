"""R1466-1490 TRACK A (WALL 1): two-point Hilbert contraction of the TRUE
nonlinear K-L operator vs its FROZEN linearization, across k.
Wall 1 = prove strict Hilbert-metric contraction uniform in k, supplied by
switching. Measure: (1) asymptotic two-point rate mu_true vs mu_frozen,
(2) switching fraction along extremal trajectories, (3) near-tie density
(the competition condition), (4) amplitude dependence of strictness."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
rng = np.random.default_rng(7)
ALPHA = math.log2(3)
lam = 2.0
W0 = lam ** -2; W2 = lam ** (ALPHA - 2); W8 = lam ** (ALPHA - 1)

def build(k):
    N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
    i = np.arange(N, dtype=np.int64); m = 3 * i + 2
    i4 = (((4 * m) % M) - 2) // 3
    r9 = m % 9
    m2 = np.where(r9 == 2)[0]; m8 = np.where(r9 == 8)[0]
    j2 = ((((4 * m[m2] - 2) // 3) % Mc) - 2) // 3
    j8 = ((((2 * m[m8] - 1) // 3) % Mc) - 2) // 3
    return N, M3, i4, m2, m8, j2, j8

def make_apply(N, M3, i4, m2, m8, j2, j8):
    def apply_T(c, frozen_arg=None, want_arg=False):
        out = W0 * c[i4]
        args = None
        if want_arg:
            args = np.empty(len(m2) + len(m8), dtype=np.int8)
        for off, (mi, ji, W) in enumerate(((m2, j2, W2), (m8, j8, W8))):
            tri = np.stack((c[ji], c[ji + M3], c[ji + 2 * M3]))
            if frozen_arg is not None:
                a = frozen_arg[off]
                b = tri[a, np.arange(len(ji))]
            else:
                b = tri.min(axis=0)
                if want_arg:
                    sl = slice(0, len(m2)) if off == 0 else slice(len(m2), None)
                    args[sl] = tri.argmin(axis=0)
            out[mi] += W * b
        return (out, args) if want_arg else out
    return apply_T

def hilbert(x, y):
    r = np.log(x) - np.log(y)
    return r.max() - r.min()

print(f"{'k':>3s} {'mu_true':>9s} {'mu_frozen':>9s} {'gap':>8s} {'switch%':>8s} {'neartie%':>9s}")
for k in (9, 11, 13):
    N, M3, i4, m2, m8, j2, j8 = build(k)
    ap = make_apply(N, M3, i4, m2, m8, j2, j8)
    # eigenvector by power iteration (Hilbert-normalized)
    v = np.ones(N)
    for _ in range(300):
        v = ap(v); v /= v.max()
    v = np.maximum(v, 1e-300)
    # argmin pattern of v (for frozen operator) + near-tie density
    tri2 = np.stack((v[j2], v[j2 + M3], v[j2 + 2 * M3]))
    tri8 = np.stack((v[j8], v[j8 + M3], v[j8 + 2 * M3]))
    arg_v = (tri2.argmin(axis=0).astype(np.int8), tri8.argmin(axis=0).astype(np.int8))
    def neartie(tri, tol):
        s = np.sort(tri, axis=0)
        return float(((s[1] - s[0]) / s[0] < tol).mean())
    nt = 0.5 * (neartie(tri2, 0.05) + neartie(tri8, 0.05))
    # two-point trajectories: x,y = v*exp(+-g), g random, amplitude eps
    eps = 0.5
    g = rng.standard_normal(N)
    x = v * np.exp(eps * g); y = v * np.exp(-eps * g)
    xf, yf = x.copy(), y.copy()
    d_prev = hilbert(x, y); df_prev = d_prev
    mus, muf, sw = [], [], []
    for t in range(60):
        x, ax = ap(x, want_arg=True); y, ay = ap(y, want_arg=True)
        x /= x.max(); y /= y.max()
        d = hilbert(x, y); mus.append(d / d_prev); d_prev = d
        sw.append(float((ax != ay).mean()))
        xf = ap(xf, frozen_arg=arg_v); yf = ap(yf, frozen_arg=arg_v)
        xf /= xf.max(); yf /= yf.max()
        df = hilbert(xf, yf); muf.append(df / df_prev); df_prev = df
    mu_t = float(np.mean(mus[-20:])); mu_f = float(np.mean(muf[-20:]))
    sw_t = float(np.mean(sw[-20:]))
    print(f"{k:3d} {mu_t:9.4f} {mu_f:9.4f} {mu_f - mu_t:8.4f} {100*sw_t:7.2f}% {100*nt:8.2f}%")

# amplitude dependence at k=11: strictness vs eps
print("\namplitude dependence (k=11): eps -> mu_true, switch%")
N, M3, i4, m2, m8, j2, j8 = build(11)
ap = make_apply(N, M3, i4, m2, m8, j2, j8)
v = np.ones(N)
for _ in range(300):
    v = ap(v); v /= v.max()
for eps in (0.01, 0.05, 0.2, 0.5, 1.0, 2.0):
    g = rng.standard_normal(N)
    x = v * np.exp(eps * g); y = v * np.exp(-eps * g)
    d_prev = hilbert(x, y); mus, sws = [], []
    for t in range(50):
        x, ax = ap(x, want_arg=True); y, ay = ap(y, want_arg=True)
        x /= x.max(); y /= y.max()
        d = hilbert(x, y); mus.append(d / d_prev); d_prev = d
        sws.append(float((ax != ay).mean()))
    print(f"  eps={eps:4.2f}: mu_true={np.mean(mus[-15:]):.4f}  switch={100*np.mean(sws[-15:]):.2f}%")
