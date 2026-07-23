"""
166b_certify_k19.py — memory-lean k=19 certificate (chunked map build).
Same rigor as 166: float32 iteration, float64 + 4-ulp-lowered final check.
Steady-state memory ~7.1 GB (maps 3.5 + vectors 3.1 + cbar 0.5), no int64
full-size intermediates.
"""
import numpy as np
import mpmath as mp
from math import log2
import sys

mp.mp.dps = 60
ALPHA_MP = mp.log(3)/mp.log(2)
CHUNK = 1 << 23          # 8M

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

def make_maps_chunked(k):
    N = 3**(k-1)
    Nl = N // 3
    T4 = np.empty(N, dtype=np.int32)
    n1 = n3 = 0
    # first pass: count D1/D3 sizes = N/3 each (exact), allocate directly
    idx1 = np.empty(Nl, dtype=np.int32)
    R1 = np.empty(Nl, dtype=np.int32)
    idx3 = np.empty(Nl, dtype=np.int32)
    R3 = np.empty(Nl, dtype=np.int32)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        i = np.arange(lo, hi, dtype=np.int64)
        T4[lo:hi] = ((4*i + 2) % N).astype(np.int32)
        s, r = np.divmod(i, 3)
        sel1 = (r == 0)
        sel3 = (r == 2)
        c1 = int(sel1.sum()); c3 = int(sel3.sum())
        idx1[n1:n1+c1] = i[sel1].astype(np.int32)
        R1[n1:n1+c1] = ((4*s[sel1]) % Nl).astype(np.int32)
        idx3[n3:n3+c3] = i[sel3].astype(np.int32)
        R3[n3:n3+c3] = ((2*s[sel3] + 1) % Nl).astype(np.int32)
        n1 += c1; n3 += c3
    assert n1 == Nl and n3 == Nl
    return N, Nl, T4, idx1, R1, idx3, R3

def F(v, co, maps, dtype):
    A, B1, B3 = (dtype(c) for c in co)
    N, Nl, T4, idx1, R1, idx3, R3 = maps
    cbar = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:3*Nl])
    w = v[T4]
    w *= A
    w[idx1] += cbar[R1] * B1
    w[idx3] += cbar[R3] * B3
    return w

if __name__ == "__main__":
    k, gt, iters = 19, 0.9070, int(sys.argv[1]) if len(sys.argv) > 1 else 1800
    lam_mp = mp.mpf(2)**mp.mpf(str(gt))
    co = lowered_coeffs(lam_mp)
    print(f"k={k}: building maps (N = {3**(k-1):,}) ...", flush=True)
    maps = make_maps_chunked(k)
    N = maps[0]
    print("maps built; iterating ...", flush=True)
    v = np.ones(N, dtype=np.float32)
    for it in range(iters):
        w = F(v, co, maps, np.float32)
        mx = w.max()
        np.divide(w, mx, out=w)
        v, w = w, v
        del w
        if it % 100 == 0:
            print(f"    iter {it}", flush=True)
    v64 = v.astype(np.float64)
    del v                      # free 1.5GB before the float64 gather peak
    w64 = F(v64, co, maps, np.float64)
    np.divide(w64, v64, out=w64)
    margin = float(w64.min())
    del w64
    ok = margin >= 1.0 + 1e-6
    if ok:
        np.save(f"certificate_k{k}.npy", v64.astype(np.float32))
    print(f"k={k}  lambda={float(lam_mp):.6f}  gamma={gt}  "
          f"min F(v)/v = {margin:.6f}  {'CERTIFIED' if ok else 'FAILED'}")
