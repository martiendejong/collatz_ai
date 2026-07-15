"""R509: EXACT CERTIFICATION OF k=20. Steps:
(1) load saved eigenvector (float32, 1.16B), chunked min-ratio test at candidate
    rational lambdas; pick the largest with margin >= 5e-4;
(2) scale by S=10^9, floor to int64 -> cert_k20.npy;
(3) exact-integer verification (object dtype, chunked) a la verify_certificates.py.
Result: verified pi(x) >= x^gamma with gamma = log2(lambda)."""
import sys, math, decimal
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
decimal.getcontext().prec = 80
ALPHA = math.log2(3)
k = 20
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
CH = 3 ** 13

c = np.load("certificates/k20_eig.npy")
print(f"loaded eigenvector: {c.shape[0]:,} float32; min {c.min():.3e} max {c.max():.3e}", flush=True)

def min_ratio(lam):
    w0 = np.float64(lam ** -2.0); w2 = np.float64(lam ** (ALPHA - 2.0)); w8 = np.float64(lam ** (ALPHA - 1.0))
    worst = np.inf
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        j = np.arange(lo, hi, dtype=np.int64)
        m = 3 * j + 2
        i4 = (((4 * m) % M) - 2) // 3
        out = w0 * c[i4].astype(np.float64)
        mod9 = m % 9
        for mv, mul, sub, W in ((2, 4, 2, w2), (8, 2, 1, w8)):
            sel = mod9 == mv
            if sel.any():
                mm = m[sel]
                t = ((mul * mm - sub) // 3) % Mc
                b = (t - 2) // 3
                cb = np.minimum(np.minimum(c[b], c[b + M3]), c[b + 2 * M3]).astype(np.float64)
                out[sel] += W * cb
        worst = min(worst, float((out / c[lo:hi]).min()))
    return worst

best = None
for num in (1885, 1883, 1881, 1879, 1877, 1875):
    lam = num / 1000
    mr = min_ratio(lam)
    print(f"lam={lam}: min_ratio={mr:.6f}", flush=True)
    if mr >= 1.0005: best = num; break
if best is None:
    print("no lambda with sufficient margin; aborting"); sys.exit(1)
lam_num, lam_den = best, 1000
print(f"certifying at lambda={lam_num}/{lam_den}, gamma={math.log2(lam_num/lam_den):.4f}", flush=True)

S = 10 ** 9
Ci = np.empty(N, dtype=np.int64)
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    Ci[lo:hi] = np.floor(c[lo:hi].astype(np.float64) * S).astype(np.int64)
np.save("certificates/cert_k20.npy", Ci)
print("cert_k20.npy saved", flush=True)
del c

# exact verification (mirrors verify_certificates.py)
Q = 10 ** 18
dl = decimal.Decimal(lam_num) / decimal.Decimal(lam_den)
ln_l = dl.ln()
W0 = (lam_den ** 2 * Q) // (lam_num ** 2)
e2 = decimal.Decimal('-0.41503749927884390')
e8 = decimal.Decimal('0.58496250072115610')
W2 = int(((e2 * ln_l).exp() * Q).to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
W8 = int(((e8 * ln_l).exp() * Q).to_integral_value(rounding=decimal.ROUND_FLOOR)) - 1
viol = 0
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    idx = np.arange(lo, hi, dtype=np.int64)
    m = 3 * idx + 2
    i4 = (((4 * m) % M) - 2) // 3
    mod9 = m % 9
    Cc = Ci[lo:hi].astype(object)
    rhs = W0 * Ci[i4].astype(object)
    for mask, mul, Wb in ((mod9 == 2, 4, W2), (mod9 == 8, 2, W8)):
        if mask.any():
            mm = m[mask]
            t = (((mul * mm - (2 if mul == 4 else 1)) // 3) % Mc)
            j0 = (t - 2) // 3
            cb = np.minimum(np.minimum(Ci[j0], Ci[j0 + M3]), Ci[j0 + 2 * M3]).astype(object)
            rhs[mask] = rhs[mask] + Wb * cb
    viol += int((Cc * Q > rhs).sum())
    if (lo // CH) % 50 == 0: print(f"  verified {hi:,}/{N:,} violations so far {viol}", flush=True)
g = math.log2(lam_num / lam_den)
print(f"\nk=20, lambda={lam_num}/{lam_den}: violations={viol}: "
      f"{'VERIFIED => pi(x) >= x^%.4f  *** NEW RECORD ***' % g if viol == 0 else 'FAILED'}", flush=True)
