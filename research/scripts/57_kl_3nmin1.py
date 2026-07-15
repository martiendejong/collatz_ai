"""R523-530: K-L ANALOG FOR 3n-1 — does the spectral machinery SEE the cycles?
Same functional system, recursion via (3m-1). For 3n+1 lambda* -> 2 (gamma->1).
For 3n-1 the tree of 1 has density ~1/3, so lambda* should stay AWAY from 2.
Bisect lambda_max at k=11 for both maps and compare. Also alpha (tempering) for both."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ALPHA = math.log2(3)
k = 11
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)

def min_ratio(lam, sign, iters=250):
    # classes m == 2 mod 3 (3n+1) or m == 1 mod 3 (3n-1): preimage of odd step:
    # 3n+1: m odd, 3n+1=m*2^w -> n=(m*2^w-1)/3 needs m*2^w==1 mod 3
    # 3n-1: n=(m*2^w+1)/3 needs m*2^w==-1==2 mod 3
    idx = np.arange(N, dtype=np.int64)
    m = 3 * idx + 2
    i4 = (((4 * m) % M) - 2) // 3
    mod9 = m % 9
    w0 = np.float64(lam ** -2.0); w2 = np.float64(lam ** (ALPHA - 2.0)); w8 = np.float64(lam ** (ALPHA - 1.0))
    # branch classes: for 3n+1 (sign=+1): m==2 mod 9 -> w even branch (mul 4), m==8 -> mul 2
    # for 3n-1 (sign=-1): preimage n=(m*2^w+1)/3: m*2^w==2 mod 3: m==2: w even; m==1: w odd.
    # we keep the same class lattice (m==2 mod 3) and adjust the t-map:
    c = np.ones(N, dtype=np.float64)
    for it in range(iters):
        f = w0 * c[i4]
        for mv, mul in ((2, 4), (8, 2)):
            sel = mod9 == mv
            mm = m[sel]
            t = ((mul * mm - sign * (2 if mul == 4 else 1) * 1) // 3) % Mc if sign > 0 else \
                ((mul * mm + (1 if (mul * mv + 1) % 3 == 0 else 2)) // 3) % Mc
            if sign < 0:
                # n = (m*2^w + 1)/3 integer when m*2^w == 2 mod 3
                x = mul * mm + 1
                good = (x % 3 == 0)
                if not good.all():
                    x = mul * mm + 2  # fallback shouldn't happen
                t = (x // 3) % Mc
            b = (t - 2) // 3
            W = w2 if mul == 4 else w8
            cb = np.minimum(np.minimum(c[b], c[b + M3]), c[b + 2 * M3])
            f[sel] += W * cb
        r = f / c
        c = f / f.max()
    return float(r.min() / r.max()), c  # feasibility proxy: uniformity of ratio

def lam_star(sign):
    lo, hi = 1.5, 2.0
    for _ in range(14):
        mid = (lo + hi) / 2
        # feasibility: does iteration keep min ratio >= 1 after normalization?
        # use direct min-ratio of converged vector:
        idxr, c = min_ratio(mid, sign)
        # apply once more and take true min ratio
        lam = mid
        w0 = np.float64(lam ** -2.0); w2 = np.float64(lam ** (ALPHA - 2.0)); w8 = np.float64(lam ** (ALPHA - 1.0))
        idx = np.arange(N, dtype=np.int64); m = 3 * idx + 2
        i4 = (((4 * m) % M) - 2) // 3; mod9 = m % 9
        f = w0 * c[i4]
        for mv, mul in ((2, 4), (8, 2)):
            sel = mod9 == mv; mm = m[sel]
            if sign > 0:
                t = ((mul * mm - (2 if mul == 4 else 1)) // 3) % Mc
            else:
                t = ((mul * mm + 1) // 3) % Mc
            b = (t - 2) // 3
            W = w2 if mul == 4 else w8
            f[sel] += W * np.minimum(np.minimum(c[b], c[b + M3]), c[b + 2 * M3])
        mr = float((f / c).min())
        if mr >= 1.0: lo = mid
        else: hi = mid
    return lo

lp = lam_star(+1)
lm = lam_star(-1)
print(f"k={k}: lambda*(3n+1) = {lp:.4f} -> gamma = {math.log2(lp):.4f}")
print(f"k={k}: lambda*(3n-1) = {lm:.4f} -> gamma = {math.log2(lm):.4f}")
print(f"difference: {lp - lm:+.4f}  (positive = machinery SEES the 3n-1 cycles)")
