# 360: G3 from the 2026-08-09 analysis — resolve the ~0.018 Fourier-Hoelder
# calibration discrepancy between Scripts 219 and 227 (Obs 433).
# Suspect: 219 fits log(MEAN |vhat|) per v3-level; 227 fits MEAN(log |vhat|)
# (geometric mean). Jensen gap varies per level -> slope shifts.
# Method: ONE vector per k (identical iterations, float64, full fft/N),
# then three aggregations on the same data:
#   A = log of arithmetic mean per level   (219 convention)
#   G = mean of logs per level             (227 convention)
#   M = log of max per level               (true Hoelder envelope)
# Report alpha for each and the per-k discrepancy A-G.
import numpy as np
from math import log2, log

ALPHA_MATH = log2(3.0)
LAM = 1.70
A_ = LAM**-2.0; B1 = LAM**(ALPHA_MATH-2.0); B3 = LAM**(ALPHA_MATH-1.0)

def perron(k, n_iter=400):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s, r = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4*s) % Nl; R3 = (2*s+1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A_*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    return v/v.mean()

def levels_v3(n_count):
    ns = np.arange(n_count, dtype=np.int64)
    v3 = np.zeros(n_count, dtype=np.int32)
    tmp = ns.copy(); act = tmp > 0
    while True:
        na = act & (tmp % 3 == 0)
        if not na.any():
            break
        tmp[na] //= 3; v3[na] += 1; act = na
    return v3

def fit(xs, ys):
    xs = np.array(xs); ys = np.array(ys)
    ok = np.isfinite(ys)
    return np.polyfit(xs[ok], ys[ok], 1)[0] if ok.sum() >= 3 else float('nan')

print(f"360: Hoelder-kalibratie, lam={LAM}, identieke vectoren, full fft/N")
print(f"{'k':>3} {'alpha_A(219)':>12} {'alpha_G(227)':>12} {'alpha_M(max)':>12} {'A-G':>8}")
serA, serG, serM, ks = [], [], [], []
for k in range(10, 17):
    v = perron(k)
    N = v.size
    vhat = np.abs(np.fft.fft(v))/N
    v3 = levels_v3(N)
    xs, ya, yg, ym = [], [], [], []
    for lev in range(k-1):
        m = (np.arange(N) > 0) & (v3 == lev)
        if m.sum() == 0:
            continue
        mags = vhat[m]
        xs.append(-lev*log(3.0))
        ya.append(log(mags.mean()))
        yg.append(float(np.log(mags + 1e-300).mean()))
        ym.append(log(mags.max()))
    aA, aG, aM = fit(xs, ya), fit(xs, yg), fit(xs, ym)
    serA.append(aA); serG.append(aG); serM.append(aM); ks.append(k)
    print(f"{k:>3} {aA:>12.4f} {aG:>12.4f} {aM:>12.4f} {aA-aG:>8.4f}", flush=True)

print()
for name, ser in [("A(219)", serA), ("G(227)", serG), ("M(max)", serM)]:
    d = np.diff(ser)
    ratios = d[1:]/(np.abs(d[:-1])+1e-12)
    avg_r = float(np.mean(ratios[-3:]))
    if 0 < avg_r < 1:
        inf_est = ser[-1] - abs(d[-1])*avg_r/(1-avg_r)
    else:
        inf_est = float('nan')
    print(f"conventie {name}: reeks {[f'{a:.4f}' for a in ser]}")
    print(f"   decrements {[f'{x:+.4f}' for x in d]}  geo-extrapolatie alpha_inf ~ {inf_est:.4f}")
