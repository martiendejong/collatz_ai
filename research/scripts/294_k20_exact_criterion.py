"""
294_k20_exact_criterion.py
==========================
Post-hoc analysis of the SAVED k=20 eigenvector (Script 289 memmap, final in w.dat
after 35 iterations = odd number of ping-pong swaps).

Verifies the EXACT criterion of step (3b): c2/c0 < R (column-min means), equivalent
g2 > R*g0 (Obs 490) — not the Gaussian sigma-proxy used in Scripts 285/289.
Cross-checks the sigma-ratio 1.07769 reported by Script 289.
"""
import numpy as np
from math import log2, sqrt

ALPHA = log2(3.0)
lam = 1.05
k = 20
N = 3**(k-1)
Nl = N // 3
Nl3 = Nl // 3
A = lam**-2
rho = 1.576710          # from Script 289 (converged wmax)
t = A / rho
R = (t**2 + lam) / (1 + t*lam)

v = np.memmap("E:/temp/collatz_k20/w.dat", dtype=np.float32, mode="r", shape=(N,))

CH = 10_000_000
# accumulators: class means, column-min means, sigma-stats (cross-check)
sum_mu0 = sum_mu2 = 0.0
sum_c0 = sum_c2 = 0.0
sum_s0 = sum_s2 = 0.0
for a in range(0, Nl3, CH):
    b = min(a + CH, Nl3)
    # class-0 column elements: v[3j], v[3j+Nl], v[3j+2Nl] for j in [a,b)
    x1 = np.asarray(v[3*a       : 3*b       : 3], dtype=np.float64)
    x2 = np.asarray(v[3*a+Nl    : 3*b+Nl    : 3], dtype=np.float64)
    x3 = np.asarray(v[3*a+2*Nl  : 3*b+2*Nl  : 3], dtype=np.float64)
    mu = (x1+x2+x3)/3.0
    mn = np.minimum(np.minimum(x1, x2), x3)
    var = ((x1-mu)**2 + (x2-mu)**2 + (x3-mu)**2)/3.0
    sum_mu0 += mu.sum(); sum_c0 += mn.sum(); sum_s0 += np.sqrt(var).sum()
    # class-2
    y1 = np.asarray(v[3*a+2      : 3*b+2      : 3], dtype=np.float64)
    y2 = np.asarray(v[3*a+2+Nl   : 3*b+2+Nl   : 3], dtype=np.float64)
    y3 = np.asarray(v[3*a+2+2*Nl : 3*b+2+2*Nl : 3], dtype=np.float64)
    mu2c = (y1+y2+y3)/3.0
    mn2 = np.minimum(np.minimum(y1, y2), y3)
    var2 = ((y1-mu2c)**2 + (y2-mu2c)**2 + (y3-mu2c)**2)/3.0
    sum_mu2 += mu2c.sum(); sum_c2 += mn2.sum(); sum_s2 += np.sqrt(var2).sum()

mu0 = sum_mu0/Nl3; mu2 = sum_mu2/Nl3
c0 = sum_c0/Nl3;  c2 = sum_c2/Nl3
g0 = mu0 - c0;    g2 = mu2 - c2

print(f"k={k} lam={lam}: rho={rho} t={t:.6f} R={R:.6f}")
print(f"  mu2/mu0 = {mu2/mu0:.8f}  (identity: should equal R={R:.8f}; diff={mu2/mu0-R:+.2e})")
print(f"  c2/c0   = {c2/c0:.8f}  (< R? {'YES' if c2/c0 < R else 'NO — FAIL'})  margin R - c2/c0 = {R-c2/c0:+.3e}")
print(f"  g2/(R*g0) = {g2/(R*g0):.6f}  (> 1? {'YES' if g2 > R*g0 else 'NO — FAIL'})")
sig_ratio = (sum_s2/sum_mu2) / (sum_s0/sum_mu0)
print(f"  sigma-ratio cross-check: {sig_ratio:.5f} (Script 289 reported 1.07769)")
print("DONE")
