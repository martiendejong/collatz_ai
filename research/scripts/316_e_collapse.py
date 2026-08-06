"""Does e(q) collapse as a function of q (intrinsic) or k-q (finite-size)?"""
import numpy as np
from math import log2
ALPHA = log2(3.0)
CACHE = "E:/projects/collatz/research/cache"

def eprofile(lam, k):
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar3 = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    G = vbar3 - cb
    out = {}
    for q in range(2, k-4):
        M = 3**(q+1); Mp = 3**q
        cmV = vbar3.reshape(Nl//M, M).mean(axis=0)
        cmG = G.reshape(Nl//M, M).mean(axis=0)
        U = cmV.reshape(3, Mp).T; W = cmG.reshape(3, Mp).T
        Uc = U - U.mean(axis=1, keepdims=True); Wc = W - W.mean(axis=1, keepdims=True)
        vU = (Uc**2).mean(axis=1); ok = vU > 1e-300
        bC = (Uc*Wc).mean(axis=1)[ok]/vU[ok]
        gamC = (W.mean(axis=1)/U.mean(axis=1))[ok]
        out[q] = float(np.median(bC/(gamC+1e-300)))
    return out

lam = 1.70
for k in [10, 12, 13, 14]:
    prof = eprofile(lam, k)
    line_q = " ".join(f"q{q}:{prof[q]:.2f}" for q in sorted(prof))
    print(f"k={k}: {line_q}")
print("\nby distance from top (k-2-q = digits above):")
for k in [10, 12, 13, 14]:
    prof = eprofile(lam, k)
    line_d = " ".join(f"d{k-2-q}:{prof[q]:.2f}" for q in sorted(prof, reverse=True))
    print(f"k={k}: {line_d}")
