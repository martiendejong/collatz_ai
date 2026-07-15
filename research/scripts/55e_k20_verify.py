"""R576c: k=20 FINAL — the vector is feasible AT lam=1.885 with margin 1.94e-4
(>> flooring error ~1e-6 at S=1e10). Floor and exactly verify at 1885/1000."""
import sys, math, decimal
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
decimal.getcontext().prec = 80
k = 20
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
CH = 3 ** 13
c = np.load("certificates/k20_polished.npy", mmap_mode="r")
S = 10 ** 10
Ci = np.lib.format.open_memmap("certificates/cert_k20.npy", mode="w+", dtype=np.int64, shape=(N,))
for lo in range(0, N, CH):
    hi = min(lo + CH, N)
    Ci[lo:hi] = np.floor(np.asarray(c[lo:hi], dtype=np.float64) * S).astype(np.int64)
Ci.flush(); del c
print("int certificate written (S=1e10); exact verify at 1885/1000", flush=True)
lam_num, lam_den = 1885, 1000
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
    if (lo // CH) % 40 == 0: print(f"  {hi:,}/{N:,} viol {viol}", flush=True)
g = math.log2(lam_num / lam_den)
print(f"\nk=20 lam={lam_num}/{lam_den}: violations={viol}: "
      f"{'VERIFIED pi(x) >= x^%.4f *** NEW RECORD ***' % g if viol == 0 else 'FAILED'}", flush=True)
