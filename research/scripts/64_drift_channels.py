"""R601-610: CHANNEL-RESOLVED LATTICE COEFFICIENTS. The naive L1 mass split
(transport 1/4 -> {P,P+1}, branch 3/4 -> {P-1,P}) predicts a > c, but measured
a < c. Measure the actual per-channel transfer on cert_k13's difference field:
for offsets delta = 2*3^P, decompose Delta_P(j) into its transport and branch
source amplitudes and regress. Where does the inversion come from?"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
k = 13
C = np.load("certificates/cert_k13.npy").astype(np.float64)
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
lam = 1.818
A = math.log2(3)
W0 = lam**-2; W2 = lam**(A-2); W8 = lam**(A-1)
rng = np.random.default_rng(601)

print(f"{'P':>2} {'|D_P|':>9} {'transport |D_P|,|D_P+1|':>24} {'branch |D_P-1|,|D_P|':>21}")
for P in (3, 5, 7):
    d = 3 ** P
    js = rng.integers(0, N - 2*d, 40000)
    m = 3 * js + 2
    # own difference at scale P
    D = np.abs(C[js + d] - C[js])
    # transport source: i4(j) offset 4d = d + 3d -> measure both parts
    i4 = (((4 * m) % M) - 2) // 3
    ok = (i4 + 4*d < N)
    Dt_P  = np.abs(C[(i4 + d) % N] - C[i4])
    Dt_P1 = np.abs(C[(i4 + 3*d) % N] - C[i4])
    # branch source (only for m==2 or 8 mod 9): t(j), offset 4d/3 = d/3 + d
    mod9 = m % 9
    out = []
    for mv, mul, sub in ((2, 4, 2), (8, 2, 1)):
        sel = mod9 == mv
        mm = m[sel]
        t = ((mul * mm - sub) // 3) % Mc
        b = (t - 2) // 3
        db = d // 3
        Db_Pm1 = np.abs(C[(b + db) % M3 + 0] - C[b % M3])
        Db_P   = np.abs(C[(b + d) % M3] - C[b % M3])
        out.append((Db_Pm1.mean(), Db_P.mean()))
    bm1 = np.mean([o[0] for o in out]); bP = np.mean([o[1] for o in out])
    print(f"{P:>2} {D.mean():>9.5f}  tp: {W0*Dt_P.mean():>8.5f},{W0*Dt_P1.mean():>8.5f}"
          f"   br: {bm1:>8.5f},{bP:>8.5f} (x w~)")
    # effective shares: which source scale dominates D_P?
    r_up = np.corrcoef(D[ok], np.abs(C[(i4 + 3*d) % N] - C[i4])[ok])[0,1]
    r_self_t = np.corrcoef(D[ok], np.abs(C[(i4 + d) % N] - C[i4])[ok])[0,1]
    print(f"    corr(D_P, transport-3delta-part) = {r_up:.3f} | corr(D_P, transport-delta-part) = {r_self_t:.3f}")
