"""R522: k=20 certification, attempt 2. Resume from checkpoint: 60 more power
iterations at lam=1.885, then monotone polish c <- min(c, f(c)) until pointwise
feasible (or 40 polish rounds), then integer-certify + exact verify."""
import sys, math, decimal
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
decimal.getcontext().prec = 80
ALPHA = math.log2(3)
k = 20
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
CH = 3 ** 13
LAM = 1.885

c = np.load("certificates/k20_eig.npy")
print(f"resumed from checkpoint: {N:,} values", flush=True)

def apply_f(c, lam):
    w0 = np.float32(lam ** -2.0); w2 = np.float32(lam ** (ALPHA - 2.0)); w8 = np.float32(lam ** (ALPHA - 1.0))
    f = np.empty(N, dtype=np.float32)
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        j = np.arange(lo, hi, dtype=np.int64)
        m = 3 * j + 2
        i4 = (((4 * m) % M) - 2) // 3
        out = w0 * c[i4]
        mod9 = m % 9
        for mv, mul, sub, W in ((2, 4, 2, w2), (8, 2, 1, w8)):
            sel = mod9 == mv
            if sel.any():
                mm = m[sel]
                t = ((mul * mm - sub) // 3) % Mc
                b = (t - 2) // 3
                out[sel] += W * np.minimum(np.minimum(c[b], c[b + M3]), c[b + 2 * M3])
        f[lo:hi] = out
    return f

for it in range(60):
    f = apply_f(c, LAM)
    mx = f.max(); f /= mx
    c = f
    if it % 15 == 0:
        print(f"  extra iter {it}: norm {mx:.6f}", flush=True)
        np.save("certificates/k20_eig.npy", c)
np.save("certificates/k20_eig.npy", c)

for lam in (1.885, 1.883, 1.881, 1.879):
    cp = c.copy()
    ok = False
    for p in range(40):
        f = apply_f(cp, lam)
        mr = float((f / cp).min())
        if mr >= 1.0005:
            ok = True
            print(f"lam={lam}: feasible after {p} polish rounds, min_ratio={mr:.6f}, "
                  f"collapse check min={cp.min():.3e}", flush=True)
            break
        np.minimum(cp, f, out=cp)
        cp /= cp.max()
    if ok and cp.min() > 1e-7:
        S = 10 ** 9
        Ci = np.empty(N, dtype=np.int64)
        for lo in range(0, N, CH):
            hi = min(lo + CH, N)
            Ci[lo:hi] = np.floor(cp[lo:hi].astype(np.float64) * S).astype(np.int64)
        np.save("certificates/cert_k20.npy", Ci)
        del cp
        lam_num, lam_den = int(round(lam * 1000)), 1000
        Q = 10 ** 18
        dl = decimal.Decimal(lam_num) / decimal.Decimal(lam_den)
        ln_l = dl.ln()
        W0 = (lam_den ** 2 * Q) // (lam_num ** 2)
        W2 = int(((decimal.Decimal('-0.41503749927884390') * ln_l).exp() * Q).to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
        W8 = int(((decimal.Decimal('0.58496250072115610') * ln_l).exp() * Q).to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
        viol = 0
        for lo in range(0, N, CH):
            hi = min(lo + CH, N)
            idx = np.arange(lo, hi, dtype=np.int64)
            m = 3 * idx + 2
            i4 = (((4 * m) % M) - 2) // 3
            mod9 = m % 9
            rhs = W0 * Ci[i4].astype(object)
            for mask, mul in ((mod9 == 2, 4), (mod9 == 8, 2)):
                if mask.any():
                    mm = m[mask]
                    t = (((mul * mm - (2 if mul == 4 else 1)) // 3) % Mc)
                    j0 = (t - 2) // 3
                    cb = np.minimum(np.minimum(Ci[j0], Ci[j0 + M3]), Ci[j0 + 2 * M3]).astype(object)
                    rhs[mask] = rhs[mask] + (W2 if mul == 4 else W8) * cb
            viol += int((Ci[lo:hi].astype(object) * Q > rhs).sum())
        g = math.log2(lam_num / lam_den)
        print(f"\nk=20 lam={lam_num}/{lam_den}: violations={viol}: "
              f"{'VERIFIED pi(x) >= x^%.4f *** NEW RECORD ***' % g if viol == 0 else 'FAILED'}", flush=True)
        break
    elif ok:
        print(f"lam={lam}: polish collapsed (min {cp.min():.3e}), trying lower lam", flush=True)
