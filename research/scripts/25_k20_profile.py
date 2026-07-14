"""R198: K-L eigenvector at k=20 (1,162,261,467 classes), chunked float32, single lambda.
Outputs: q, intra-triple CV, profile CV_p, lattice fit (a,c,theta)."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)
k = 20
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
LAM = 1.885
W0 = np.float32(LAM ** -2.0); W2 = np.float32(LAM ** (ALPHA - 2.0)); W8 = np.float32(LAM ** (ALPHA - 1.0))
CH = 3 ** 14
ITERS = 150

c = np.ones(N, dtype=np.float32)
print(f"k={k}: N={N:,} classes, lambda={LAM}, {ITERS} iterations", flush=True)
for it in range(ITERS):
    f = np.empty(N, dtype=np.float32)
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        j = np.arange(lo, hi, dtype=np.int64)
        m = 3 * j + 2
        i4 = (((4 * m) % M) - 2) // 3
        out = W0 * c[i4]
        mod9 = m % 9
        for mask_val, mul, sub, W in ((2, 4, 2, W2), (8, 2, 1, W8)):
            sel = mod9 == mask_val
            if sel.any():
                mm = m[sel]
                t = ((mul * mm - sub) // 3) % Mc
                b = (t - 2) // 3
                cb = np.minimum(np.minimum(c[b], c[b + M3]), c[b + 2 * M3])
                out[sel] += W * cb
        f[lo:hi] = out
    del c
    mx = f.max()
    f /= mx
    c = f
    if it % 25 == 0: print(f"  iter {it}: max-ratio norm {mx:.6f}", flush=True)

# measurements
tri = np.stack([c[:M3], c[M3:2*M3], c[2*M3:3*M3]])
q = float(3 * tri.min(axis=0).sum(dtype=np.float64) / c.sum(dtype=np.float64))
cv_top = float((tri.std(axis=0) / tri.mean(axis=0)).mean(dtype=np.float64))
del tri
prof = []
for p in range(1, k - 1):
    B = 3 ** p
    acc_n = 0.0; acc_c = 0.0
    for lo in range(0, N, 3 * CH):
        hi = min(lo + 3 * CH, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        d = (idx // B) % 3
        sel = idx[d == 0]
        sel = sel[sel + 2 * B < N]
        if sel.size == 0: continue
        t = np.stack([c[sel], c[sel + B], c[sel + 2 * B]]).astype(np.float64)
        cv = t.std(axis=0) / t.mean(axis=0)
        acc_c += cv.sum(); acc_n += cv.size
    prof.append(acc_c / acc_n)
    print(f"  profile p={p}: CV={prof[-1]:.4f}", flush=True)
prof = np.array(prof)
lo_i, hi_i = 2, len(prof) - 2
A = np.stack([prof[lo_i-1:hi_i-1], prof[lo_i+1:hi_i+1]], axis=1)
y = prof[lo_i:hi_i]
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
a, cc = coef
disc = 1 - 4 * a * cc
th = (1 - math.sqrt(disc)) / (2 * cc) if disc > 0 else float('nan')
print(f"\nFINAL k=20 (lam={LAM}): q = {q:.5f}   CV_top = {cv_top:.5f}")
print(f"lattice fit: a = {a:.4f}  c = {cc:.4f}  a+c = {a+cc:.4f}  theta = {th:.4f}")
print(f"PREDICTIONS were: (a,c) trend -> k20 ~ (0.472, 0.521), theta ~ 0.849")
