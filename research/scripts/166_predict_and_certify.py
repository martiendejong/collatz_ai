"""
166_predict_and_certify.py
===========================
EXTENDING THE RECORD: k = 17, 18 via predict-and-certify.

The gamma(k) model is now predictive (k=16 predicted 0.889, observed 0.8893),
so bisection is unnecessary: aim directly at the predicted exponent minus a
safety margin and run ONE certificate attempt per level.

Memory-lean implementation for N up to 3^17 = 129M states:
  - index maps int32 (N < 2^31)
  - Perron iteration in float32 (margin targets ~4e-4 >> float32 noise per
    step; the iteration only needs to CONVERGE, accuracy is irrelevant)
  - final verification pass in float64 with directed-rounding (4-ulp lowered)
    coefficients: the certificate check itself is rigorous.

Targets: k=17: gamma_cert = 0.8950 (predicted lambda* ~ 0.8968)
         k=18: gamma_cert = 0.9020 (predicted ~ 0.9043) -- the x^0.90 line.
"""
import numpy as np
import mpmath as mp
from math import log2
import sys

mp.mp.dps = 60
ALPHA_MP = mp.log(3)/mp.log(2)

def lowered_coeffs(lam_mp):
    vals = [lam_mp**-2, lam_mp**(ALPHA_MP-2), lam_mp**(ALPHA_MP-1)]
    out = []
    for x in vals:
        f = float(x)
        for _ in range(4):
            f = np.nextafter(f, 0.0)
        assert mp.mpf(f) <= x
        out.append(f)
    return out

def make_maps(k):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = ((4*i + 2) % N).astype(np.int32)
    s = i // 3
    r = (i % 3).astype(np.int8)
    Nl = N // 3
    m1 = (r == 0); m3 = (r == 2)
    R1 = ((4*s) % Nl).astype(np.int32)[m1]      # compact: only where needed
    R3 = ((2*s + 1) % Nl).astype(np.int32)[m3]
    idx1 = np.nonzero(m1)[0].astype(np.int32)
    idx3 = np.nonzero(m3)[0].astype(np.int32)
    del i, s, r, m1, m3
    return N, Nl, T4, idx1, R1, idx3, R3

def F(v, co, maps, dtype):
    A, B1, B3 = (dtype(c) for c in co)
    N, Nl, T4, idx1, R1, idx3, R3 = maps
    cbar = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:3*Nl])
    w = v[T4] * A
    w[idx1] += cbar[R1] * B1
    w[idx3] += cbar[R3] * B3
    return w

def predict_and_certify(k, gamma_target, iters):
    lam_mp = mp.mpf(2)**mp.mpf(str(gamma_target))
    co = lowered_coeffs(lam_mp)
    maps = make_maps(k)
    N = maps[0]
    v = np.ones(N, dtype=np.float32)
    for it in range(iters):
        w = F(v, co, maps, np.float32)
        v = w / w.max()
        if it % 200 == 0:
            print(f"    iter {it}", flush=True)
    # rigorous final check in float64
    v64 = v.astype(np.float64)
    del v
    w64 = F(v64, co, maps, np.float64)
    margin = float((w64 / v64).min())
    ok = margin >= 1.0 + 1e-6
    if ok:
        np.save(f"certificate_k{k}.npy", v64.astype(np.float32))
    return float(lam_mp), margin, ok

if __name__ == "__main__":
    k = int(sys.argv[1]); gt = float(sys.argv[2]); iters = int(sys.argv[3])
    print(f"k={k}: certifying gamma = {gt} (N = {3**(k-1):,}) ...", flush=True)
    lam, margin, ok = predict_and_certify(k, gt, iters)
    print(f"k={k}  lambda={lam:.6f}  gamma={gt}  min F(v)/v = {margin:.6f}  "
          f"{'CERTIFIED' if ok else 'FAILED'}")
