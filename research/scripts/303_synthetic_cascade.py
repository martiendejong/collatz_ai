"""
Corrected mean-field cascade: class-pure cb columns + uniform sub-class mixing.
cb-col of class q = elementwise min of 3 draws from class-q population.
Receiving columns get sub-class q uniform (matching m mod 3 equidistribution).
Fluctuation scale renormalized each level (exact in linearized-min regime);
decay ratio tracked -> compare with plateau sqrt(0.41)=0.64 (variance ratio 0.41).
"""
import numpy as np
from math import log2

rng = np.random.default_rng(7)
ALPHA = log2(3.0)

for lam, rho_meas, meas in [(1.05, 1.5765, (1.093, 0.73, 0.41)), (1.70, 1.0469, (1.15, 0.31, 0.70))]:
    A = lam**-2
    t = A / rho_meas
    w = t * lam**ALPHA
    R = (t**2+lam)/(1+t*lam)
    P = 150_000
    v0 = 1 + 0.001*rng.standard_normal((P,3))
    v2 = R + 0.001*rng.standard_normal((P,3))
    ratios = []
    for lev in range(50):
        v1 = t * v0
        pops = [v0, v1, v2]
        def cbdraw():
            # uniform sub-class q; 3 same-class draws, elementwise min
            q = rng.integers(0, 3, size=P)
            out = np.empty((P,3))
            for qq in range(3):
                sel = q == qq
                n = sel.sum()
                idx = rng.integers(0, P, size=(n,3))
                pop = pops[qq]
                out[sel] = np.minimum(np.minimum(pop[idx[:,0]], pop[idx[:,1]]), pop[idx[:,2]])
            return out
        nv0 = t*v2[rng.integers(0,P,P)] + w*cbdraw()
        nv2 = (t**2)*v0[rng.integers(0,P,P)] + lam*w*cbdraw()
        sc = nv0.mean()
        v0 = nv0/sc; v2 = nv2/sc
        # fluctuation scale (within-column relative std, class 0)
        eps = (v0.std(axis=1)/v0.mean(axis=1)).mean()
        ratios.append(eps)
        # renormalize fluctuations to keep numerics healthy (scale-equivariant regime)
        if eps < 3e-4:
            for arr in (v0, v2):
                mcol = arr.mean(axis=1, keepdims=True)
                arr *= 1  # in-place via slicing:
            v0 = v0.mean(axis=1, keepdims=True) + (v0 - v0.mean(axis=1, keepdims=True))*(1e-3/eps)
            v2 = v2.mean(axis=1, keepdims=True) + (v2 - v2.mean(axis=1, keepdims=True))*(1e-3/eps)
    mu0 = v0.mean(); mu2 = v2.mean()
    def gap(c): return c.mean(axis=1)-c.min(axis=1)
    def sig(c): return c.std(axis=1)
    s_ratio = (sig(v2).mean()/mu2)/(sig(v0).mean()/mu0)
    v1 = t*v0; pops = [v0, v1, v2]
    cb1 = cbdraw(); cb2 = cbdraw()
    a = v2[rng.integers(0,P,P)]; b = v1[rng.integers(0,P,P)]
    S0 = t*gap(a) + w*gap(cb1) - gap(t*a+w*cb1)
    S2 = t*gap(b) + lam*w*gap(cb2) - gap(t*b+lam*w*cb2)
    dec = np.array(ratios)
    # decay ratio in the clean mid-range (levels 10-25, between transient and renorm kicks)
    rr = dec[11:26]/dec[10:25]
    print(f"lam={lam}: mu2/mu0={mu2/mu0:.5f} (R={R:.5f}) sigma-ratio={s_ratio:.4f} (meas {meas[0]})")
    print(f"   s2/s0={S2.mean()/S0.mean():.4f} (meas {meas[1]}) | eps-decay/level={np.median(rr):.4f} "
          f"-> var-ratio={np.median(rr)**2:.4f} (meas plateau {meas[2]})", flush=True)
