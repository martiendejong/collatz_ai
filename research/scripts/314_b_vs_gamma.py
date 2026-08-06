"""Is b_C = local relative gap gamma_C = E[G|C]/E[vbar|C]? Then the c<1 machinery
reduces entirely to the distribution of cell relative gaps."""
import numpy as np

CACHE = "E:/projects/collatz/research/cache"
for lam in [1.05, 1.70]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar3 = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    G = vbar3 - cb
    print(f"\n=== lam={lam} ===")
    print("q | corr(b_C, gamma_C) | mean b/gamma | median b/gamma")
    for q in [1, 2, 3, 4, 5, 6]:
        M = 3**(q+1); Mp = 3**q
        cmV = vbar3.reshape(Nl//M, M).mean(axis=0)
        cmG = G.reshape(Nl//M, M).mean(axis=0)
        U = cmV.reshape(3, Mp).T; W = cmG.reshape(3, Mp).T
        Uc = U - U.mean(axis=1, keepdims=True)
        Wc = W - W.mean(axis=1, keepdims=True)
        vU = (Uc**2).mean(axis=1)
        ok = vU > 1e-300
        bC = (Uc*Wc).mean(axis=1)[ok]/vU[ok]
        gamC = (W.mean(axis=1)/U.mean(axis=1))[ok]   # cell relative gap
        cc = np.corrcoef(bC, gamC)[0, 1]
        rat = bC/(gamC + 1e-300)
        print(f"{q} | {cc:+.3f} | {rat.mean():+.3f} | {np.median(rat):+.3f}")
