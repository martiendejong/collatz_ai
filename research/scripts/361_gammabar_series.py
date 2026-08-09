# 361: C6 from the 2026-08-09 analysis — extend the gammabar(2,k) series to
# k=17 using cached endpoint vectors (lam=2.00) plus cold runs for the 12-14 gap.
# Mass identity (Obs 511/514, exact): 1 - rho(2,k) = (3/4)*gammabar(2,k).
# Root law (Obs 513): amplitude rates = sqrt(c); c(2) via rate(gammabar)^2.
# This is the third instrument on c(2) next to TR and the rho-rate, extended
# two depths beyond the Obs 513 data.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
LAM = 2.0
A = LAM**-2.0; B1 = LAM**(ALPHA-2.0); B3 = LAM**(ALPHA-1.0)
CACHE = r"E:\projects\collatz\research\cache"

def rho_k(k, iters_cold=2000, tol=1e-13):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s, r = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4*s) % Nl; R3 = (2*s+1) % Nl
    fn = os.path.join(CACHE, f"v_lam2.00_k{k}.npy")
    if os.path.exists(fn):
        v = np.load(fn).astype(np.float64)
        iters, min_it = 600, 200
    else:
        v = np.ones(N, dtype=np.float64)
        iters, min_it = iters_cold, 1200
    rho_prev = 0.0
    stable = 0
    for it in range(iters):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max())
        v = w/rho
        stable = stable + 1 if abs(rho - rho_prev) < tol else 0
        if it > min_it and stable >= 30:
            break
        rho_prev = rho
    # Rayleigh-type growth on converged vector (max-norm growth = Perron rho)
    return rho, it

print("361: gammabar(2,k)-reeks via massa-identiteit 1-rho = (3/4)*gammabar")
print(f"{'k':>3} {'rho(2,k)':>12} {'gammabar':>12} {'ratio':>8} {'iters':>6}")
gb_prev = None
series = []
for k in range(9, 18):
    rho, it = rho_k(k)
    gb = 4.0*(1.0-rho)/3.0
    ratio = gb/gb_prev if gb_prev else float('nan')
    series.append((k, rho, gb, ratio))
    print(f"{k:>3} {rho:>12.8f} {gb:>12.8f} {ratio:>8.4f} {it:>6}", flush=True)
    gb_prev = gb
print()
rates = [r for _, _, _, r in series[1:]]
print(f"rate-reeks: {[f'{r:.4f}' for r in rates]}")
print(f"c(2)-schattingen (rate^2): {[f'{r*r:.4f}' for r in rates]}")
print("vergelijk: TR-instrument c(2)-reeks eindigde op 0.8357 en steeg nog (Obs 526)")
