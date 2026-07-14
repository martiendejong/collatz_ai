"""E33: K-L system at k=17 (43,046,721 classes), float32 memory-lean."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)
k = 19
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1)

idx = np.arange(N, dtype=np.int64)
m = 3 * idx + 2
i4 = ((((4 * m) % M) - 2) // 3).astype(np.int32)
mod9 = (m % 9).astype(np.int8)
del idx
w2_rows = np.where(mod9 == 2)[0].astype(np.int32)
w8_rows = np.where(mod9 == 8)[0].astype(np.int32)
del mod9
m2 = 3 * w2_rows.astype(np.int64) + 2
t = ((4 * m2 - 2) // 3) % Mc
r2 = np.stack([((t - 2) // 3), (((t + Mc) - 2) // 3), (((t + 2 * Mc) - 2) // 3)]).astype(np.int32)
del m2, t
m8 = 3 * w8_rows.astype(np.int64) + 2
t = ((2 * m8 - 1) // 3) % Mc
r8 = np.stack([((t - 2) // 3), (((t + Mc) - 2) // 3), (((t + 2 * Mc) - 2) // 3)]).astype(np.int32)
del m8, t, m
print("setup done", flush=True)

def min_ratio(lam, iters=150):
    w0 = np.float32(lam ** -2.0)
    w2 = np.float32(lam ** (ALPHA - 2.0))
    w8 = np.float32(lam ** (ALPHA - 1.0))
    c = np.ones(N, dtype=np.float32)
    for it in range(iters):
        f = w0 * c[i4]
        cb = np.minimum(np.minimum(c[r2[0]], c[r2[1]]), c[r2[2]])
        f[w2_rows] += w2 * cb
        cb = np.minimum(np.minimum(c[r8[0]], c[r8[1]]), c[r8[2]])
        f[w8_rows] += w8 * cb
        c = f / f.max()
    f = w0 * c[i4]
    cb = np.minimum(np.minimum(c[r2[0]], c[r2[1]]), c[r2[2]])
    f[w2_rows] += w2 * cb
    cb = np.minimum(np.minimum(c[r8[0]], c[r8[1]]), c[r8[2]])
    f[w8_rows] += w8 * cb
    return float((f / c).min())

lo, hi = 1.858, 1.892
for step in range(8):
    mid = (lo + hi) / 2
    mr = min_ratio(mid)
    print(f"lam={mid:.5f} min_ratio={mr:.6f}", flush=True)
    if mr >= 1.0: lo = mid
    else: hi = mid
print(f"FINAL k=17: lambda={lo:.5f} gamma={math.log2(lo):.5f}")
