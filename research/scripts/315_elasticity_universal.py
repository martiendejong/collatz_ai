"""A: is e -> 2 universal across lambda? B: mechanism test — multifractal coupling:
slope of log(local roughness) vs log(local mean) across cells should be e-1 -> 1."""
import numpy as np
from math import log2
import os

ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

def get_v(lam, k=13, iters=400):
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        return np.load(fn)
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N)
    for _ in range(iters):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    np.save(fn, v)
    return v

for lam in [1.05, 1.30, 1.70, 1.90]:
    k = 13
    v = get_v(lam)
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar3 = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    G = vbar3 - cb
    es = []; slopes = []
    for q in [2, 3, 4, 5, 6]:
        M = 3**(q+1); Mp = 3**q
        cmV = vbar3.reshape(Nl//M, M).mean(axis=0)
        cmG = G.reshape(Nl//M, M).mean(axis=0)
        U = cmV.reshape(3, Mp).T; W = cmG.reshape(3, Mp).T
        Uc = U - U.mean(axis=1, keepdims=True); Wc = W - W.mean(axis=1, keepdims=True)
        vU = (Uc**2).mean(axis=1); ok = vU > 1e-300
        bC = (Uc*Wc).mean(axis=1)[ok]/vU[ok]
        gamC = (W.mean(axis=1)/U.mean(axis=1))[ok]
        e = np.median(bC/(gamC+1e-300))
        es.append(e)
        # B: multifractal coupling across cells at this level:
        # slope of log(gamma_C) on log(mean_C) should be ~ e-1
        mC = U.mean(axis=1)[ok]
        x = np.log(mC); y = np.log(np.abs(gamC)+1e-300)
        sl = np.polyfit(x, y, 1)[0]
        slopes.append(sl)
    print(f"lam={lam}: e(q=2..6) = {[f'{x:.2f}' for x in es]} | dlog(gam)/dlog(m) = {[f'{s:+.2f}' for s in slopes]}")
