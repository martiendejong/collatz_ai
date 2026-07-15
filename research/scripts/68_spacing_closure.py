"""R641-650: THE SPACING CLOSURE. Fine triples (m, m+9, m+18) map through the
K-L equation to image triples with spacing 36 (transport) and 12 (branch).
Measure CV(s) = mean triple-CV at class-spacing s; test the closure
   CV(s)^2 = A*CV(4s)^2 + B*CV(4s/3)^2
with fitted (A,B) vs theory (W0^2, w~^2) x incoherence. If it closes, the
fine-end is a self-consistent bounded map -> Prop 23 semi-analytic."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
k = 13
C = np.load("certificates/cert_k13.npy").astype(np.float64)
C /= C.mean()
N = 3 ** (k - 1); M = 3 ** k
rng = np.random.default_rng(641)

def cv_spacing(s_m):
    """mean CV of triples (m, m+s, m+2s) in class space (s_m divisible by 3)."""
    si = s_m // 3
    js = rng.integers(0, N - 2*si, 150000)
    t = np.stack([C[js], C[js + si], C[js + 2*si]])
    return float((t.std(0) / t.mean(0)).mean())

spac = {}
for s in (9, 12, 16*3, 18, 24, 36, 48, 72, 108, 144, 216):
    if s % 3 == 0:
        spac[s] = cv_spacing(s)
print("CV by class-spacing s:")
for s in sorted(spac): print(f"  s={s:>4}: CV = {spac[s]:.4f}")

# closure test on pairs (s, 4s, 4s/3) where all measured
lam = 1.818; A_ = math.log2(3)
W0 = lam**-2; wbar = (lam**(A_-2) + lam**(A_-1)) / 3 * 1.5  # rough mean branch weight
rows = []
for s in (9, 12, 18, 36):
    if 4*s in spac and (4*s) % 3 == 0 and (4*s//3) in spac:
        rows.append((s, spac[s], spac[4*s], spac[4*s//3]))
X = np.array([[r[2]**2, r[3]**2] for r in rows])
y = np.array([r[1]**2 for r in rows])
coef, *_ = np.linalg.lstsq(X, y, rcond=None)
A, B = coef
pred = X @ coef
print(f"\nclosure fit CV(s)^2 = A*CV(4s)^2 + B*CV(4s/3)^2:")
print(f"  A = {A:.4f} (theory W0^2 = {W0**2:.4f})   B = {B:.4f} (theory ~w2 incoh)")
for (s, cs, c4, c43), p in zip(rows, pred):
    print(f"  s={s:>3}: measured {cs**2:.5f} predicted {p:.5f} ratio {cs**2/p:.3f}")
