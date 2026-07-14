"""R367-373: EXCLUDED-SET ANATOMY. At time t, which numbers in a window are
still outside the wavefront? Predictors: trailing-ones k, n mod small powers
of 2, family membership. Measures P(excluded|k) and the mod-64 profile, plus
overlap with the coefficient-survivor classes of R284 (mod 2^12)."""
import sys
import numpy as np
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

# window [2^20, 2^21), odd numbers; t = median d in window
lo, hi = 1 << 20, 1 << 21
nn = np.arange(lo | 1, hi, 2)
dw = d[nn]
t_med = int(np.median(dw))
for t, label in ((t_med, "median"), (int(np.percentile(dw, 90)), "p90"), (int(np.percentile(dw, 99)), "p99")):
    exc = dw > t
    # trailing ones of n
    k = np.zeros(len(nn), dtype=int)
    x = nn.copy()
    while (x & 1).any():
        odd = (x & 1).astype(bool)
        k[odd] += 1
        x = np.where(odd, x >> 1, x)
        x[~odd] = 0
    print(f"\nt = {t} ({label}): excluded fraction {exc.mean():.4f}")
    print("  P(excluded | trailing-ones k):")
    for kk in range(1, 9):
        sel = k == kk
        if sel.sum() > 300:
            print(f"    k={kk}: {exc[sel].mean():.4f}  (n={sel.sum()})")
    # mod 64 concentration
    m64 = nn % 64
    ps = [exc[m64 == r].mean() for r in range(1, 64, 2)]
    ps = np.array(ps)
    print(f"  P(excluded | n mod 64): min {ps.min():.4f}, max {ps.max():.4f}, "
          f"ratio {ps.max()/max(ps.min(),1e-9):.1f}x; top classes:",
          [int(r) for r in np.array(range(1, 64, 2))[np.argsort(-ps)[:4]]])
