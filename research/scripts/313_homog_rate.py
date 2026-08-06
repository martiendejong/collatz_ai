"""
Homogenization rate: distribution of local proportionality b_C = cov(u_C,w_C)/var(u_C)
per cell, across depth q. If sd(b_C)/|mean| shrinks geometrically, the field
homogenizes with that rate (prediction: the cascade ratio itself).
Also: local correlation corr_C(u,w) distribution.
"""
import numpy as np
from math import log2

CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar3 = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    G = vbar3 - cb
    print(f"\n=== lam={lam} k={k} ===")
    print("q | mean(b_C) | sd(b_C) | sd/mean | frac corr_C>0 | median corr_C")
    prev_sd = None
    for q in range(0, 9):
        M = 3**(q+1)
        Mp = 3**q
        cmV = vbar3.reshape(Nl//M, M).mean(axis=0)
        cmG = G.reshape(Nl//M, M).mean(axis=0)
        U = cmV.reshape(3, Mp).T
        W = cmG.reshape(3, Mp).T
        Uc = U - U.mean(axis=1, keepdims=True)
        Wc = W - W.mean(axis=1, keepdims=True)
        vU = (Uc**2).mean(axis=1)
        cUW = (Uc*Wc).mean(axis=1)
        vW = (Wc**2).mean(axis=1)
        ok = vU > 0
        bC = cUW[ok]/vU[ok]
        # weight by var (aggregate b = weighted mean)
        bw = cUW[ok].sum()/vU[ok].sum()
        corrC = cUW[ok]/np.sqrt(vU[ok]*vW[ok] + 1e-300)
        sd = bC.std()
        ratio = f" sd-ratio {sd/prev_sd:.3f}" if prev_sd else ""
        print(f"{q} | {bw:+.5f} | {sd:.5f} | {sd/abs(bw):8.2f} | {(corrC>0).mean():.3f} | {np.median(corrC):+.3f}{ratio}")
        prev_sd = sd
