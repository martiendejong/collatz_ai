"""R576b: k=20 certification, attempt 4. Polish to fixed point at lam=1.885,
save; margin comes from evaluating at lower lam'; S=1e10 flooring so margin
1.00005 suffices. Certify at the largest feasible lam'."""
import sys, math, decimal
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
decimal.getcontext().prec = 80
ALPHA = math.log2(3)
k = 20
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
CH = 3 ** 13

c = np.load("certificates/k20_eig.npy")
print(f"loaded {N:,}", flush=True)

def f_chunk(lo, hi, lam):
    w0 = np.float32(lam ** -2.0); w2 = np.float32(lam ** (ALPHA - 2.0)); w8 = np.float32(lam ** (ALPHA - 1.0))
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
    return out

def min_ratio(lam):
    worst = np.inf
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        worst = min(worst, float((f_chunk(lo, hi, lam) / c[lo:hi]).min()))
    return worst

LAM = 1.885
prev = 0.0
for rnd in range(30):
    mr = min_ratio(LAM)
    print(f"polish {rnd}: min_ratio {mr:.7f} min_c {c.min():.3e}", flush=True)
    if mr >= 1.0 - 1e-6: break
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        np.minimum(c[lo:hi], f_chunk(lo, hi, LAM), out=c[lo:hi])
    prev = mr
np.save("certificates/k20_polished.npy", c)
print("polished vector saved", flush=True)

target = None
for lamp in (1.884, 1.883, 1.882, 1.881, 1.879, 1.877):
    mr = min_ratio(lamp)
    print(f"margin at lam'={lamp}: {mr:.7f}", flush=True)
    if mr >= 1.00005: target = lamp; break
if target is None:
    print("no margin found", flush=True); sys.exit(1)

S = 10 ** 10
Ci = np.lib.format.open_memmap("certificates/cert_k20.npy", mode="w+", dtype=np.int64, shape=(N,))
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    Ci[lo:hi] = np.floor(c[lo:hi].astype(np.float64) * S).astype(np.int64)
Ci.flush(); del c
lam_num, lam_den = int(round(target * 1000)), 1000
print(f"int certificate written; verifying at {lam_num}/{lam_den}", flush=True)

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
    rhs = W0 * np.asarray(Ci[i4]).astype(object)
    for mask, mul in ((mod9 == 2, 4), (mod9 == 8, 2)):
        if mask.any():
            mm = m[mask]
            t = (((mul * mm - (2 if mul == 4 else 1)) // 3) % Mc)
            j0 = (t - 2) // 3
            cb = np.minimum(np.minimum(np.asarray(Ci[j0]), np.asarray(Ci[j0 + M3])), np.asarray(Ci[j0 + 2 * M3])).astype(object)
            rhs[mask] = rhs[mask] + (W2 if mul == 4 else W8) * cb
    viol += int((np.asarray(Ci[lo:hi]).astype(object) * Q > rhs).sum())
    if (lo // CH) % 40 == 0: print(f"  verify {hi:,}/{N:,} viol {viol}", flush=True)
g = math.log2(lam_num / lam_den)
print(f"\nk=20 lam={lam_num}/{lam_den}: viol={viol}: "
      f"{'VERIFIED pi(x) >= x^%.4f *** NEW RECORD ***' % g if viol == 0 else 'FAILED'}", flush=True)
