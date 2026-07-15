"""R586-594: CV_1 SATURATION TEST. Compute the K-L eigenvector at k=8..12
(own critical lambda via bisection), measure fine-end CV_1; combine with
k=13..20 measured values; fit CV_1(k) = CVinf - b*rho^k."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)

def build(k, lam, iters=400):
    N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
    idx = np.arange(N, dtype=np.int64); m = 3 * idx + 2
    i4 = (((4 * m) % M) - 2) // 3
    mod9 = m % 9
    sel2, sel8 = mod9 == 2, mod9 == 8
    t2 = (((4 * m[sel2] - 2) // 3) % Mc); b2 = (t2 - 2) // 3
    t8 = (((2 * m[sel8] - 1) // 3) % Mc); b8 = (t8 - 2) // 3
    w0 = lam ** -2.0; w2 = lam ** (ALPHA - 2.0); w8 = lam ** (ALPHA - 1.0)
    c = np.ones(N)
    for _ in range(iters):
        f = w0 * c[i4]
        f[sel2] += w2 * np.minimum(np.minimum(c[b2], c[b2 + M3]), c[b2 + 2 * M3])
        f[sel8] += w8 * np.minimum(np.minimum(c[b8], c[b8 + M3]), c[b8 + 2 * M3])
        c = f / f.max()
    f = w0 * c[i4]
    f[sel2] += w2 * np.minimum(np.minimum(c[b2], c[b2 + M3]), c[b2 + 2 * M3])
    f[sel8] += w8 * np.minimum(np.minimum(c[b8], c[b8 + M3]), c[b8 + 2 * M3])
    return c, float((f / c).min())

def cv1(c, k):
    N = 3 ** (k - 1); B = 3
    idx = np.arange(N, dtype=np.int64)
    sel = idx[(idx // B) % 3 == 0]; sel = sel[sel + 2 * B < N]
    t = np.stack([c[sel], c[sel + B], c[sel + 2 * B]])
    return float((t.std(0) / t.mean(0)).mean())

pts = {13: 0.4144, 15: 0.4313, 17: 0.4453, 20: 0.4629}
for k in (8, 9, 10, 11, 12):
    lo, hi = 1.5, 2.0
    for _ in range(12):
        mid = (lo + hi) / 2
        c, mr = build(k, mid, 300)
        if mr >= 1.0: lo = mid
        else: hi = mid
    c, _ = build(k, lo, 400)
    pts[k] = cv1(c, k)
    print(f"k={k}: lam*={lo:.4f} CV_1={pts[k]:.4f}", flush=True)

ks = sorted(pts)
print("\nall points:", " ".join(f"k{k}:{pts[k]:.4f}" for k in ks))
# saturation fit CV1(k) = A - b*rho^k  (grid over rho)
best = None
for rho in np.arange(0.70, 0.99, 0.005):
    X = np.column_stack([np.ones(len(ks)), -np.array([rho**k for k in ks])])
    y = np.array([pts[k] for k in ks])
    coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = ((X @ coef - y)**2).sum()
    if best is None or r < best[0]: best = (r, rho, coef)
r, rho, (A, b) = best
print(f"saturation fit: CV_1(k) = {A:.4f} - {b:.3f}*{rho:.3f}^k   (residual {r:.2e})")
print(f"CV_1(infinity) = {A:.4f}  -> BOUNDED: {'YES' if A < 1 else 'suspect'}")
