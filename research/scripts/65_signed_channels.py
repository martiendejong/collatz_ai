"""R611-620: SIGNED CHANNEL COVARIANCE. Exact difference identity (Thm 16):
Delta_P(m) = W0*[D_delta(y+3d) + D_3delta(y)]  (transport telescoped, y=i4(m))
           + w(m)*[min-triple difference at offset (4/3 or 2/3)*delta].
Branch offset telescopes to scales P-1 and P. Measure Var/Cov of all signed
parts on cert_k13 -> effective (a,c) reconstruction -> where the inversion
a < c comes from."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
k = 13
C = np.load("certificates/cert_k13.npy").astype(np.float64)
C /= C.mean()
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
lam = 1.818; A = math.log2(3)
W0 = lam**-2; W2 = lam**(A-2); W8 = lam**(A-1)
rng = np.random.default_rng(611)

def cls(m): return ((m % M) - 2) // 3
def tri_min(b):
    return np.minimum(np.minimum(C[b % M3], C[b % M3 + M3]), C[b % M3 + 2*M3])

P = 5
d = 2 * 3 ** P     # class-space offset, divisible by 9 -> same mod-9 type
js = rng.integers(0, N, 60000)
m = 3 * js + 2
mp = m + d
LHS = C[cls(mp)] - C[cls(m)]
y = 4 * m; yp = 4 * mp        # transport images (offset 4d)
T_low  = C[cls(y + d)] - C[cls(y)]            # D_delta at y   (scale P)
T_high = C[cls(yp)] - C[cls(y + d)]           # D_3delta       (scale P+1)
T = W0 * (T_low + T_high)
mod9 = m % 9
B = np.zeros_like(LHS); Bl = np.zeros_like(LHS); Bh = np.zeros_like(LHS)
for mv, mul, sub, W in ((2, 4, 2, W2), (8, 2, 1, W8)):
    sel = mod9 == mv
    mm = m[sel]; mmp = mp[sel]
    t  = ((mul * mm - sub) // 3) % Mc
    tp = ((mul * mmp - sub) // 3) % Mc
    b, bp = (t - 2) // 3, (tp - 2) // 3
    full = tri_min(bp) - tri_min(b)
    off = (mul * d) // 3   # branch offset: (4/3)d or (2/3)d
    dl = off - (off // (3**P)) * 0  # telescope: off = low(P-1 scale part) + rest
    # telescope at scale boundary 3^P: off = q*3^P + r
    q, r = divmod(off, 3 ** P)
    mid = (t + r) % Mc
    bm = (mid - 2) // 3
    low  = tri_min(bm) - tri_min(b)     # scale P-1 part (offset r < 3^P)
    high = tri_min(bp) - tri_min(bm)    # scale >= P part
    B[sel] = W * full; Bl[sel] = W * low; Bh[sel] = W * high
resid = LHS - (T + B)
print(f"identity check: corr(LHS, T+B) = {np.corrcoef(LHS, T+B)[0,1]:.6f}, "
      f"residual std/LHS std = {resid.std()/LHS.std():.4f}")
V = np.var
print(f"\nvariance shares of LHS (Var {V(LHS):.3e}):")
for name, x in (("transport-low (P)", W0*T_low), ("transport-high (P+1)", W0*T_high),
                ("branch-low (P-1)", Bl), ("branch-high (P/P+)", Bh)):
    print(f"  {name:>22}: Var {V(x):.3e} share {V(x)/V(LHS):.3f} corr(LHS) {np.corrcoef(LHS,x)[0,1]:+.3f}")
print(f"\ncross-covariances (normalized by Var LHS):")
parts = {"Tl": W0*T_low, "Th": W0*T_high, "Bl": Bl, "Bh": Bh}
ks_ = list(parts)
for i in range(4):
    for j2 in range(i+1, 4):
        cv = np.cov(parts[ks_[i]], parts[ks_[j2]])[0,1]/V(LHS)
        print(f"  cov({ks_[i]},{ks_[j2]}) = {cv:+.3f}")
