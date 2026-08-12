# 379: Route A step 15 — quantify the excess contraction of L for the Doeblin
# lemma: (a) k-flatness of the factor; (b) linearity in perturbation amplitude
# (small-perturbation limit = the linearised smoothing the lemma must bound).
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def setup(lam, k):
    N = 3**(k-1); Nl = N//3; Nll = Nl//3
    s = np.arange(Nl, dtype=np.int64)
    sig4 = (4*s+2) % Nl
    sl = s//3
    R1p = (4*sl) % Nll; R3p = (2*sl+1) % Nll
    r = s % 3
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    v = np.load(os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")).astype(np.float64)
    rho = 1.0
    for _ in range(200):
        i = np.arange(N); T4 = (4*i+2) % N
        cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        ri = i % 3; si = i//3
        w[ri == 2] += B3*cbv[((2*si+1) % Nl)[ri == 2]]
        w[ri == 0] += B1*cbv[((4*si) % Nl)[ri == 0]]
        rho = float(w.max()); v = w/rho
    t = A/rho
    def L(h):
        f = np.zeros(Nl)
        f[r == 0] = (B1/rho)*h[R1p[r == 0]]
        f[r == 2] = (B3/rho)*h[R3p[r == 2]]
        g = np.zeros(Nl)
        for _ in range(int(np.log(1e-14)/np.log(t))+1):
            g = t*g[sig4] + f
        return g
    cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    cb3 = (cbv[:Nll] + cbv[Nll:2*Nll] + cbv[2*Nll:])/3.0
    return L, cb3

def osc(x, jmax=4):
    lx = np.log(x)
    return np.array([float(np.sqrt((lx.reshape(len(lx)//3**(j+1), 3**(j+1)).std(axis=1)**2).mean()))
                     for j in range(jmax)])

print(f"{'lam':>5} {'k':>3} {'eps':>5}  excess-contractie per schaal (uit-echt)/(in-echt)")
for lam in [1.70, 2.00]:
    for k in [10, 12, 14]:
        L, cb3 = setup(lam, k)
        base_in = osc(cb3); base_out = osc(L(cb3))
        rng = np.random.default_rng(379)
        for eps in ([0.5, 0.1] if k == 12 else [0.5]):
            h = cb3*np.exp(eps*rng.standard_normal(len(cb3)))
            e_in = osc(h) - base_in
            e_out = osc(L(h)) - base_out
            ratio = e_out/e_in
            print(f"{lam:>5} {k:>3} {eps:>5}  {np.array2string(ratio, precision=3)}", flush=True)
